// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include "cmpmc/core/core.hpp"
#include "cmpmc/transport/engine.hpp"
#include "cmpmc/superconductor/kaplan_film.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define CMPMC_HOST_DEVICE __host__ __device__
#else
#define CMPMC_HOST_DEVICE
#endif

namespace cmpmc::experimental {

// Experimental: queue-based ("wavefront") phonon transport pipeline, structured
// for a later GPU port. Guarded by require_validated(); see docs/experimental.md.
//
// Rather than one loop that handles every interaction type per track, all active
// tracks pass through four stages connected by flat work queues, so each stage
// runs a uniform workload:
//   1. stage1_free_flight_and_raytrace: advance each track by
//      min(d_phys, d_geom) and sort it into q_isotope, q_decay or q_boundary;
//   2. stage2_isotope_scatter: new DOS-weighted mode, isotropic k, v_g lookup;
//   3. stage3_anharmonic_decay: inverse-CDF energy split, two daughters pushed
//      back to the active queue;
//   4. stage4_boundary_and_kaplan: sensor absorption via KaplanSurrogateLUT, or
//      surface reflection.
// Known simplification relative to SimulationEngine::transport: decay daughters
// are emitted isotropically.
struct WavefrontQueues {
    std::vector<PhononTrack> q_active;
    std::vector<PhononTrack> q_next_active;
    std::vector<PhononTrack> q_isotope;
    std::vector<PhononTrack> q_decay;

    struct BoundaryEvent {
        PhononTrack trk;
        SurfaceId surface_id;
        Vec3 inward_normal;
    };
    std::vector<BoundaryEvent> q_boundary;

    void clear_all() noexcept {
        q_active.clear();
        q_next_active.clear();
        q_isotope.clear();
        q_decay.clear();
        q_boundary.clear();
    }

    void reserve_all(size_t cap) {
        q_active.reserve(cap);
        q_next_active.reserve(cap);
        q_isotope.reserve(cap);
        q_decay.reserve(cap / 4 + 64);
        q_boundary.reserve(cap / 2 + 64);
    }
};

class WavefrontPipelineEngine {
public:
    WavefrontPipelineEngine(const CrystalLattice& lattice, const DetectorModel& detector)
        : lattice_(lattice),
          detector_(detector),
          kaplan_surrogate_((require_validated(detector.allow_unvalidated, "experimental::WavefrontPipelineEngine",
                                               "isotropic decay daughters"),
                             KaplanFilm::aluminum_sensor())) {}

    void transport_wavefront(std::vector<PhononTrack>& initial_tracks,
                             std::vector<SensorHit>& hits,
                             SimStats& stats,
                             Rng& rng) const {
        WavefrontQueues queues;
        queues.reserve_all(std::max<size_t>(4096, initial_tracks.size() * 8));
        queues.q_active = std::move(initial_tracks);
        initial_tracks.clear();

        while (!queues.q_active.empty()) {
            queues.q_next_active.clear();
            queues.q_isotope.clear();
            queues.q_decay.clear();
            queues.q_boundary.clear();

            // Stage 1: free flight and boundary intersection (sorts into 3 queues)
            stage1_free_flight_and_raytrace(queues, stats, rng);

            // Stage 2: isotope and polycrystal elastic scattering
            stage2_isotope_scatter(queues, stats, rng);

            // Stage 3: anharmonic down-conversion (L -> TT / LT)
            stage3_anharmonic_decay(queues, stats, rng);

            // Stage 4: boundary reflection and sensor absorption (Kaplan surrogate table)
            stage4_boundary_and_kaplan(queues, hits, stats, rng);

            std::swap(queues.q_active, queues.q_next_active);
        }
    }

private:
    const CrystalLattice& lattice_;
    const DetectorModel& detector_;
    KaplanSurrogateLUT kaplan_surrogate_;

    void stage1_free_flight_and_raytrace(WavefrontQueues& q, SimStats& stats, Rng& rng) const {
        const double min_e = detector_.min_phonon_eV;
        const uint32_t max_b = detector_.max_bounces;

        for (PhononTrack& trk : q.q_active) {
            if (trk.energy < min_e || trk.bounces >= max_b) {
                stats.killed_tracks++;
                stats.escaped_or_killed_eV += trk.energy * trk.weight;
                continue;
            }

            stats.total_steps++;

            const double r_iso  = lattice_.isotope_rate(trk.energy);
            const double r_poly = lattice_.polycrystal_rate(trk.vg_mag);
            const double r_down = lattice_.downconversion_rate(trk.mode, trk.energy);
            const double r_tot  = r_iso + r_poly + r_down;

            const double d_phys = (r_tot > 1e-12)
                ? (-std::log(rng.uniform()) * (trk.vg_mag / r_tot))
                : std::numeric_limits<double>::infinity();

            const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, trk.vg_dir);

            if (d_phys < bhit.distance) {
                trk.pos  += trk.vg_dir * d_phys;
                trk.time += d_phys / trk.vg_mag;
                if (r_down > 0.0 && rng.uniform() * r_tot < r_down) {
                    q.q_decay.push_back(trk);
                } else {
                    q.q_isotope.push_back(trk);
                }
            } else {
                trk.pos  += trk.vg_dir * bhit.distance;
                trk.time += bhit.distance / trk.vg_mag;
                trk.bounces++;
                q.q_boundary.push_back({trk, bhit.surface_id, bhit.inward_normal});
            }
        }
    }

    void stage2_isotope_scatter(WavefrontQueues& q, SimStats& stats, Rng& rng) const {
        stats.isotope_scatters += q.q_isotope.size();
        for (PhononTrack& trk : q.q_isotope) {
            trk.mode  = lattice_.sample_dos_mode(rng);
            trk.k_dir = rng.isotropic_dir();
            const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
            trk.vg_dir = kin.vg_dir;
            trk.vg_mag = kin.vg_mag;
            q.q_next_active.push_back(trk);
        }
    }

    void stage3_anharmonic_decay(WavefrontQueues& q, SimStats& stats, Rng& rng) const {
        const double fTT = lattice_.config().fTT_frac;
        for (const PhononTrack& trk : q.q_decay) {
            const bool is_tt = (rng.uniform() < fTT);
            if (is_tt) stats.decays_L_to_TT++;
            else       stats.decays_L_to_LT++;

            const double x = is_tt
                ? lattice_.sample_tt_energy_fraction(rng)
                : lattice_.sample_lt_energy_fraction(rng);

            const PhononMode m1 = is_tt ? lattice_.sample_transverse_mode(rng) : PhononMode::L;
            const PhononMode m2 = lattice_.sample_transverse_mode(rng);

            const Vec3 kd1 = rng.isotropic_dir();
            const Vec3 kd2 = rng.isotropic_dir();
            const auto kin1 = lattice_.map_k_to_vg_oct(m1, kd1);
            const auto kin2 = lattice_.map_k_to_vg_oct(m2, kd2);

            q.q_next_active.push_back(PhononTrack{
                .pos = trk.pos, .k_dir = kd1, .vg_dir = kin1.vg_dir,
                .energy = x * trk.energy, .vg_mag = kin1.vg_mag,
                .time = trk.time, .weight = trk.weight, .bounces = trk.bounces, .mode = m1
            });
            q.q_next_active.push_back(PhononTrack{
                .pos = trk.pos, .k_dir = kd2, .vg_dir = kin2.vg_dir,
                .energy = (1.0 - x) * trk.energy, .vg_mag = kin2.vg_mag,
                .time = trk.time, .weight = trk.weight, .bounces = trk.bounces, .mode = m2
            });
        }
    }

    void stage4_boundary_and_kaplan(WavefrontQueues& q,
                                    std::vector<SensorHit>& hits,
                                    SimStats& stats,
                                    Rng& rng) const {
        KaplanBuffer kaplan_refl;
        for (auto& ev : q.q_boundary) {
            PhononTrack& trk = ev.trk;
            const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(ev.surface_id)];

            int32_t ch_id = -1;
            bool hit_active_pad = true;
            if (surf.electrode_mask.has_value()) {
                ch_id = surf.electrode_mask->find_channel(trk.pos, rng);
                hit_active_pad = (ch_id >= 0);
            }

            if (hit_active_pad && surf.kaplan_sensor.has_value() && rng.uniform() < surf.prob_abs) {
                if (rng.uniform() < surf.kaplan_sensor->film_absorption) {
                    const double edep = kaplan_surrogate_.absorb_phonon_fast(trk.energy, kaplan_refl, rng);
                    if (edep > 0.0) {
                        const double w_edep = edep * trk.weight;
                        stats.sensor_hits++;
                        stats.total_deposited_eV += w_edep;
                        hits.push_back(SensorHit{
                            .pos = trk.pos, .time = trk.time, .edep = w_edep,
                            .incident_energy = trk.energy, .surface_id = ev.surface_id,
                            .excitation = ExcitationType::Phonon, .mode = trk.mode, .channel_id = ch_id
                        });
                    }
                    double sum_refl = 0.0;
                    for (double re : kaplan_refl) sum_refl += re;
                    stats.escaped_or_killed_eV += std::max(0.0, trk.energy - edep - sum_refl) * trk.weight;
                    if (!kaplan_refl.empty()) {
                        trk.pos    += ev.inward_normal * 1e-11;
                        trk.energy  = kaplan_refl[0];
                        trk.mode    = lattice_.sample_dos_mode(rng);
                        trk.k_dir   = rng.lambertian_dir(ev.inward_normal);
                        const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
                        trk.vg_dir  = kin.vg_dir;
                        trk.vg_mag  = kin.vg_mag;
                        q.q_next_active.push_back(trk);
                    }
                    continue;
                }
            }

            if (rng.uniform() > surf.prob_refl) {
                stats.killed_tracks++;
                stats.escaped_or_killed_eV += trk.energy * trk.weight;
                continue;
            }

            stats.surface_refls++;
            trk.pos   += ev.inward_normal * 1e-11;
            trk.k_dir  = rng.lambertian_dir(ev.inward_normal);
            const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
            trk.vg_dir = kin.vg_dir;
            trk.vg_mag = kin.vg_mag;
            q.q_next_active.push_back(trk);
        }
    }
};

} // namespace cmpmc::experimental
