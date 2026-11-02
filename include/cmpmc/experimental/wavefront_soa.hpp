// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include "cmpmc/core/core.hpp"
#include "cmpmc/transport/engine.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cmpmc::experimental {

// Experimental: structure-of-arrays (SoA) phonon transport.
// Guarded by require_validated(); see docs/experimental.md.
//
// Phonon tracks are transposed into 16-wide, cache-line-aligned batches
// (PhononBatchSoA16) so that the cylinder intersection and free-flight sampling
// can be auto-vectorised across 16 tracks. Known simplifications relative to
// SimulationEngine::transport: only the L -> TT decay channel, isotropic decay
// daughters, diffuse walls only, and only the first re-emitted film phonon.
struct alignas(64) PhononBatchSoA16 {
    static constexpr size_t W = 16;

    alignas(64) double px[W]{};
    alignas(64) double py[W]{};
    alignas(64) double pz[W]{};

    alignas(64) double vx[W]{};
    alignas(64) double vy[W]{};
    alignas(64) double vz[W]{};
    alignas(64) double vmag[W]{};

    alignas(64) double kx[W]{};
    alignas(64) double ky[W]{};
    alignas(64) double kz[W]{};

    alignas(64) double energy[W]{};
    alignas(64) double time[W]{};
    alignas(64) double weight[W]{};
    alignas(64) uint32_t bounces[W]{};
    alignas(64) uint8_t mode[W]{};

    // Outputs computed in lockstep across all 16 lanes
    alignas(64) double d_geom[W]{};
    alignas(64) double d_phys[W]{};
    alignas(64) double rate_iso[W]{};
    alignas(64) double rate_down[W]{};
    alignas(64) uint8_t hit_surf[W]{};

    // Vectorized 16-wide ray-cylinder exit intersection (branch-free arithmetic)
    inline void intersect_cylinder_16(double radius, double half_z) noexcept {
        const double r2 = radius * radius;
        for (size_t i = 0; i < W; ++i) {
            // 1. Axial endcap distance
            const double vz_i = vz[i];
            const double target_z = (vz_i >= 0.0) ? half_z : -half_z;
            const double denom_z  = ( std::abs(vz_i) > 1e-15 ) ? vz_i : ((vz_i >= 0.0) ? 1e-15 : -1e-15);
            const double tz       = std::max(0.0, (target_z - pz[i]) / denom_z);

            // 2. Radial cylinder sidewall distance
            const double a    = std::max(1e-18, vx[i] * vx[i] + vy[i] * vy[i]);
            const double b    = px[i] * vx[i] + py[i] * vy[i];
            const double c    = px[i] * px[i] + py[i] * py[i] - r2;
            const double disc = std::max(0.0, b * b - a * c);
            const double tr   = std::max(0.0, (-b + std::sqrt(disc)) / a);

            const bool hit_endcap = (tz <= tr);
            d_geom[i]   = hit_endcap ? tz : tr;
            hit_surf[i] = hit_endcap
                ? static_cast<uint8_t>((vz_i >= 0.0) ? SurfaceId::TopFace : SurfaceId::BottomFace)
                : static_cast<uint8_t>(SurfaceId::SideWall);
        }
    }

    // Vectorized 16-wide free-flight sampling and rate evaluation
    inline void compute_free_flight_16(double B_iso, double A_anh, const double u_rand[W]) noexcept {
        constexpr double inv_h = 1.0 / units::h_Planck;
        for (size_t i = 0; i < W; ++i) {
            const double nu  = energy[i] * inv_h;
            const double nu2 = nu * nu;
            const double nu4 = nu2 * nu2;
            const double r_i = B_iso * nu4;
            const double r_d = (mode[i] == static_cast<uint8_t>(PhononMode::L)) ? (A_anh * nu4 * nu) : 0.0;
            const double r_t = r_i + r_d;

            rate_iso[i]  = r_i;
            rate_down[i] = r_d;
            d_phys[i]    = (r_t > 1e-12) ? (-std::log(u_rand[i]) * (vmag[i] / r_t)) : 1e30;
        }
    }
};

class WavefrontSoAEngine {
public:
    WavefrontSoAEngine(const CrystalLattice& lattice, const DetectorModel& detector)
        : lattice_(lattice),
          detector_(detector),
          scalar_engine_(lattice, detector),
          kaplan_surrogate_((require_validated(detector.allow_unvalidated, "experimental::WavefrontSoAEngine",
                                               "L->TT only, isotropic daughters, diffuse-only walls, first re-emitted phonon only"),
                             KaplanFilm::aluminum_sensor())) {}

    // Transport phonon stack using 16-wide SoA SIMD batches when queue >= 16
    void transport_soa(std::vector<PhononTrack>& stack,
                       std::vector<SensorHit>& hits,
                       SimStats& stats,
                       Rng& rng) const {
        PhononBatchSoA16 batch{};
        alignas(64) double u_rand[PhononBatchSoA16::W];
        KaplanBuffer kaplan_refl;

        const double R      = detector_.solid.radius;
        const double H      = detector_.solid.half_z;
        const double B_iso  = lattice_.config().B_iso;
        const double A_anh  = lattice_.config().A_anh;

        while (stack.size() >= PhononBatchSoA16::W) {
            const size_t base = stack.size() - PhononBatchSoA16::W;
            // Pack 16 tracks into SoA registers
            for (size_t i = 0; i < PhononBatchSoA16::W; ++i) {
                const PhononTrack& t = stack[base + i];
                batch.px[i]      = t.pos.x;
                batch.py[i]      = t.pos.y;
                batch.pz[i]      = t.pos.z;
                batch.vx[i]      = t.vg_dir.x;
                batch.vy[i]      = t.vg_dir.y;
                batch.vz[i]      = t.vg_dir.z;
                batch.vmag[i]    = t.vg_mag;
                batch.kx[i]      = t.k_dir.x;
                batch.ky[i]      = t.k_dir.y;
                batch.kz[i]      = t.k_dir.z;
                batch.energy[i]  = t.energy;
                batch.time[i]    = t.time;
                batch.weight[i]  = t.weight;
                batch.bounces[i] = t.bounces;
                batch.mode[i]    = static_cast<uint8_t>(t.mode);
                u_rand[i]        = rng.uniform();
            }
            stack.resize(base);

            // Execute 16-wide vectorized geometry & physics step calculation
            batch.intersect_cylinder_16(R, H);
            batch.compute_free_flight_16(B_iso, A_anh, u_rand);
            stats.total_steps += PhononBatchSoA16::W;

            // Advance all 16 lanes in lockstep
            for (size_t i = 0; i < PhononBatchSoA16::W; ++i) {
                const double d_step = std::min(batch.d_phys[i], batch.d_geom[i]);
                batch.px[i]   += batch.vx[i] * d_step;
                batch.py[i]   += batch.vy[i] * d_step;
                batch.pz[i]   += batch.vz[i] * d_step;
                batch.time[i] += d_step / batch.vmag[i];
            }

            // Resolve interaction outcome for each lane
            for (size_t i = 0; i < PhononBatchSoA16::W; ++i) {
                if (batch.energy[i] < detector_.min_phonon_eV || batch.bounces[i] >= detector_.max_bounces) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += batch.energy[i] * batch.weight[i];
                    continue;
                }

                PhononTrack trk{
                    .pos     = {batch.px[i], batch.py[i], batch.pz[i]},
                    .k_dir   = {batch.kx[i], batch.ky[i], batch.kz[i]},
                    .vg_dir  = {batch.vx[i], batch.vy[i], batch.vz[i]},
                    .energy  = batch.energy[i],
                    .vg_mag  = batch.vmag[i],
                    .time    = batch.time[i],
                    .weight  = batch.weight[i],
                    .bounces = batch.bounces[i],
                    .mode    = static_cast<PhononMode>(batch.mode[i])
                };

                if (batch.d_phys[i] < batch.d_geom[i]) {
                    const double r_down = batch.rate_down[i];
                    const double r_tot  = batch.rate_iso[i] + r_down;
                    if (r_down > 0.0 && rng.uniform() * r_tot < r_down) {
                        // Anharmonic decay via the inverse-CDF tables and the octahedral v_g table
                        stats.decays_L_to_TT++;
                        const double x = lattice_.sample_tt_energy_fraction(rng);
                        const PhononMode m1 = lattice_.sample_transverse_mode(rng);
                        const PhononMode m2 = lattice_.sample_transverse_mode(rng);
                        const Vec3 kd1 = rng.isotropic_dir();
                        const Vec3 kd2 = rng.isotropic_dir();
                        const auto kin1 = lattice_.map_k_to_vg_oct(m1, kd1);
                        const auto kin2 = lattice_.map_k_to_vg_oct(m2, kd2);

                        stack.push_back(PhononTrack{
                            .pos = trk.pos, .k_dir = kd2, .vg_dir = kin2.vg_dir,
                            .energy = (1.0 - x) * trk.energy, .vg_mag = kin2.vg_mag,
                            .time = trk.time, .weight = trk.weight, .bounces = trk.bounces, .mode = m2
                        });
                        trk.energy = x * trk.energy;
                        trk.mode   = m1;
                        trk.k_dir  = kd1;
                        trk.vg_dir = kin1.vg_dir;
                        trk.vg_mag = kin1.vg_mag;
                        stack.push_back(trk);
                    } else {
                        stats.isotope_scatters++;
                        trk.mode  = lattice_.sample_dos_mode(rng);
                        trk.k_dir = rng.isotropic_dir();
                        const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
                        trk.vg_dir = kin.vg_dir;
                        trk.vg_mag = kin.vg_mag;
                        stack.push_back(trk);
                    }
                } else {
                    // Boundary hit
                    trk.bounces++;
                    const SurfaceId sid = static_cast<SurfaceId>(batch.hit_surf[i]);
                    const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(sid)];

                    Vec3 inward_normal{0.0, 0.0, -1.0};
                    if (sid == SurfaceId::BottomFace) {
                        inward_normal = {0.0, 0.0, 1.0};
                    } else if (sid == SurfaceId::SideWall) {
                        const double r_xy = std::hypot(trk.pos.x, trk.pos.y);
                        inward_normal = (r_xy > 1e-15) ? Vec3{-trk.pos.x / r_xy, -trk.pos.y / r_xy, 0.0} : Vec3{-1.0, 0.0, 0.0};
                    }

                    if (surf.kaplan_sensor.has_value() && rng.uniform() < surf.prob_abs) {
                        if (rng.uniform() < surf.kaplan_sensor->film_absorption) {
                            const double edep = kaplan_surrogate_.absorb_phonon_fast(trk.energy, kaplan_refl, rng);
                            if (edep > 0.0) {
                                stats.sensor_hits++;
                                stats.total_deposited_eV += edep * trk.weight;
                                hits.push_back(SensorHit{
                                    .pos = trk.pos, .time = trk.time, .edep = edep * trk.weight,
                                    .incident_energy = trk.energy, .surface_id = sid,
                                    .excitation = ExcitationType::Phonon, .mode = trk.mode
                                });
                            }
                            double sum_refl = 0.0;
                            for (double re : kaplan_refl) sum_refl += re;
                            stats.escaped_or_killed_eV += std::max(0.0, trk.energy - edep - sum_refl) * trk.weight;
                            if (!kaplan_refl.empty()) {
                                trk.pos    += inward_normal * 1e-11;
                                trk.energy  = kaplan_refl[0];
                                trk.mode    = lattice_.sample_dos_mode(rng);
                                trk.k_dir   = rng.lambertian_dir(inward_normal);
                                const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
                                trk.vg_dir  = kin.vg_dir;
                                trk.vg_mag  = kin.vg_mag;
                                stack.push_back(trk);
                            }
                            continue;
                        }
                    }

                    stats.surface_refls++;
                    trk.pos   += inward_normal * 1e-11;
                    trk.k_dir  = rng.lambertian_dir(inward_normal);
                    const auto kin = lattice_.map_k_to_vg_oct(trk.mode, trk.k_dir);
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                    stack.push_back(trk);
                }
            }
        }

        // Drain remaining (<16) tracks with scalar engine
        if (!stack.empty()) {
            scalar_engine_.transport(stack, hits, stats, rng);
        }
    }

private:
    const CrystalLattice& lattice_;
    const DetectorModel& detector_;
    SimulationEngine scalar_engine_;
    KaplanSurrogateLUT kaplan_surrogate_;
};

} // namespace cmpmc::experimental
