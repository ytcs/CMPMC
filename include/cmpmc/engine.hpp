#pragma once

#include "cmpmc/core.hpp"
#include "cmpmc/geometry.hpp"
#include "cmpmc/kaplan.hpp"
#include "cmpmc/lattice.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace cmpmc {

// ============================================================================
// PhononTrack: Flat 64-byte Cache-Aligned State Struct
// Replaces: G4Track + G4DynamicParticle + G4CMPPhononTrackInfo + G4Step
// ============================================================================
struct alignas(64) PhononTrack {
    Vec3 pos;            // 24B: Position in detector frame (m)
    Vec3 k_dir;          // 24B: Wavevector unit direction
    Vec3 vg_dir;         // 24B: Group velocity unit direction
    double energy;       //  8B: Phonon energy (eV)
    double vg_mag;       //  8B: Group velocity magnitude (m/s)
    double time;         //  8B: Global track time (s)
    uint32_t bounces;    //  4B: Boundary reflection counter
    PhononMode mode;     //  1B: Polarization mode (L, FT, ST)
};

struct SensorHit {
    Vec3 pos;
    double time;
    double edep;
    double incident_energy;
    SurfaceId surface_id;
    PhononMode mode;
};

struct SimStats {
    uint64_t total_steps        = 0;
    uint64_t isotope_scatters   = 0;
    uint64_t decays_L_to_TT     = 0;
    uint64_t decays_L_to_LT     = 0;
    uint64_t surface_refls      = 0;
    uint64_t sensor_hits        = 0;
    uint64_t killed_tracks      = 0;
    double total_injected_eV    = 0.0;
    double total_deposited_eV   = 0.0;
    double escaped_or_killed_eV = 0.0;

    SimStats& operator+=(const SimStats& r) noexcept {
        total_steps          += r.total_steps;
        isotope_scatters     += r.isotope_scatters;
        decays_L_to_TT       += r.decays_L_to_TT;
        decays_L_to_LT       += r.decays_L_to_LT;
        surface_refls        += r.surface_refls;
        sensor_hits          += r.sensor_hits;
        killed_tracks        += r.killed_tracks;
        total_injected_eV    += r.total_injected_eV;
        total_deposited_eV   += r.total_deposited_eV;
        escaped_or_killed_eV += r.escaped_or_killed_eV;
        return *this;
    }
};

// ============================================================================
// DetectorModel & SimulationEngine
// ============================================================================
struct DetectorModel {
    CylinderSolid solid{};
    std::array<SurfaceProperty, static_cast<size_t>(SurfaceId::Count)> surfaces{
        SurfaceProperty::g4cmp_al_sensor_surface(), // TopFace
        SurfaceProperty::g4cmp_al_sensor_surface(), // BottomFace
        SurfaceProperty::g4cmp_bare_wall_surface()  // SideWall
    };
    uint32_t max_bounces    = 5000;
    double min_phonon_eV    = 10.0 * units::ueV; // Cutoff below 2*Delta_W
};

class SimulationEngine {
public:
    SimulationEngine(const CrystalLattice& lattice, const DetectorModel& detector)
        : lattice_(lattice), detector_(detector) {}

    // Create a primary phonon track with group velocity initialized from the LUT
    inline PhononTrack make_track(const Vec3& pos, const Vec3& k_dir,
                                  PhononMode mode, double energy_eV, double time_s = 0.0) const noexcept {
        const Vec3 k_unit = k_dir.unit();
        const auto kin = lattice_.map_k_to_vg(mode, k_unit);
        return PhononTrack{
            .pos     = pos,
            .k_dir   = k_unit,
            .vg_dir  = kin.vg_dir,
            .energy  = energy_eV,
            .vg_mag  = kin.vg_mag,
            .time    = time_s,
            .bounces = 0,
            .mode    = mode
        };
    }

    // Transport a batch of primary phonons using a zero-allocation per-thread stack
    void transport(std::vector<PhononTrack>& stack,
                   std::vector<SensorHit>& hits,
                   SimStats& stats,
                   Rng& rng) const {
        SmallVec<double, 64> kaplan_reflected;

        while (!stack.empty()) {
            PhononTrack trk = stack.back();
            stack.pop_back();

            // Follow active track in registers until absorbed or killed
            for (;;) {
                if (trk.energy < detector_.min_phonon_eV || trk.bounces >= detector_.max_bounces) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += trk.energy;
                    break;
                }

                stats.total_steps++;

                // 1. Combined bulk interaction rate (1/s) -> 1 logarithm per step!
                const double rate_iso  = lattice_.isotope_rate(trk.energy);
                const double rate_down = lattice_.downconversion_rate(trk.mode, trk.energy);
                const double rate_tot  = rate_iso + rate_down;

                const double d_phys = (rate_tot > 1e-12)
                    ? (-std::log(rng.uniform()) * (trk.vg_mag / rate_tot))
                    : std::numeric_limits<double>::infinity();

                // 2. Analytic ray-cylinder exit distance
                const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, trk.vg_dir);

                // 3A. Bulk Interaction wins (Isotope Scatter or Anharmonic Decay)
                if (d_phys < bhit.distance) {
                    trk.pos  += trk.vg_dir * d_phys;
                    trk.time += d_phys / trk.vg_mag;

                    if (rate_down > 0.0 && rng.uniform() * rate_tot < rate_down) {
                        // Anharmonic downconversion of L-phonon
                        perform_anharmonic_decay(trk, stack, stats, rng);
                    } else {
                        // Elastic isotope scattering (mode conversion + isotropic k)
                        stats.isotope_scatters++;
                        trk.mode  = lattice_.sample_dos_mode(rng);
                        trk.k_dir = rng.isotropic_dir();
                        const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                        trk.vg_dir = kin.vg_dir;
                        trk.vg_mag = kin.vg_mag;
                    }
                    continue;
                }

                // 3B. Boundary Interaction wins
                trk.pos  += trk.vg_dir * bhit.distance;
                trk.time += bhit.distance / trk.vg_mag;
                trk.bounces++;

                const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];

                // Check for thin-film sensor interaction (G4CMPPhononElectrode + KaplanQP)
                if (surf.kaplan_sensor.has_value() && rng.uniform() < surf.prob_abs) {
                    const KaplanFilm& film = *surf.kaplan_sensor;
                    if (rng.uniform() < film.film_absorption) {
                        const double edep = film.absorb_phonon(trk.energy, kaplan_reflected, rng);
                        if (edep > 0.0) {
                            stats.sensor_hits++;
                            stats.total_deposited_eV += edep;
                            hits.push_back(SensorHit{
                                .pos             = trk.pos,
                                .time            = trk.time,
                                .edep            = edep,
                                .incident_energy = trk.energy,
                                .surface_id      = bhit.surface_id,
                                .mode            = trk.mode
                            });
                        }
                        // Energy lost to sub-threshold thermalization inside film
                        double sum_refl = 0.0;
                        for (double re : kaplan_reflected) sum_refl += re;
                        stats.escaped_or_killed_eV += std::max(0.0, trk.energy - edep - sum_refl);

                        // Re-emit any escaped Kaplan phonons back into the crystal
                        if (kaplan_reflected.empty()) {
                            break; // Track fully absorbed
                        }
                        // Nudge inward to prevent boundary re-trigger
                        const Vec3 reemit_pos = trk.pos + bhit.inward_normal * 1e-11;
                        for (size_t i = 1; i < kaplan_reflected.size(); ++i) {
                            const PhononMode m = lattice_.sample_dos_mode(rng);
                            const Vec3 kd = sample_inward_wavevector(m, bhit.inward_normal, rng);
                            stack.push_back(make_track(reemit_pos, kd, m, kaplan_reflected[i], trk.time));
                        }
                        // Continue first reflected phonon in-place
                        trk.pos    = reemit_pos;
                        trk.energy = kaplan_reflected[0];
                        trk.mode   = lattice_.sample_dos_mode(rng);
                        trk.k_dir  = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                        const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                        trk.vg_dir = kin.vg_dir;
                        trk.vg_mag = kin.vg_mag;
                        continue;
                    }
                }

                // Surface Reflection (Specular vs. Lambertian Diffuse)
                if (rng.uniform() > surf.prob_refl) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += trk.energy;
                    break;
                }

                stats.surface_refls++;
                trk.pos += bhit.inward_normal * 1e-11; // Nudge inside crystal

                const double p_spec = surf.specular_prob(trk.energy);
                if (rng.uniform() < p_spec) {
                    // Specular reflection of wavevector around surface normal
                    Vec3 k_refl = trk.k_dir - bhit.inward_normal * (2.0 * trk.k_dir.dot(bhit.inward_normal));
                    auto kin = lattice_.map_k_to_vg(trk.mode, k_refl);
                    // In strongly anisotropic crystals, specular k can occasionally have vg pointing outward;
                    // fall back to inward Lambertian wavevector if needed (replaces G4CMP 100-trial loop)
                    if (kin.vg_dir.dot(bhit.inward_normal) <= 1e-6) {
                        k_refl = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                        kin = lattice_.map_k_to_vg(trk.mode, k_refl);
                    }
                    trk.k_dir  = k_refl;
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                } else {
                    // Diffuse Lambertian reflection
                    trk.k_dir  = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                    const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                }
            }
        }
    }

private:
    const CrystalLattice& lattice_;
    const DetectorModel& detector_;

    // Sample a Lambertian wavevector guaranteed to have an inward-pointing group velocity
    inline Vec3 sample_inward_wavevector(PhononMode mode, const Vec3& inward_normal, Rng& rng) const noexcept {
        for (int trial = 0; trial < 16; ++trial) {
            const Vec3 kd = rng.lambertian_dir(inward_normal);
            const auto kin = lattice_.map_k_to_vg(mode, kd);
            if (kin.vg_dir.dot(inward_normal) > 1e-6) {
                return kd;
            }
        }
        return inward_normal;
    }

    // Anharmonic decay: L -> T + T (74%) or L -> L' + T (26%)
    // Updates `trk` in-place with Daughter 1 and pushes Daughter 2 onto `stack`.
    inline void perform_anharmonic_decay(PhononTrack& trk,
                                         std::vector<PhononTrack>& stack,
                                         SimStats& stats,
                                         Rng& rng) const {
        const Vec3 dir1 = trk.k_dir;
        const Vec3 perp = dir1.orthogonal_unit();
        const double phi = units::two_pi * rng.uniform();
        const double c_phi = std::cos(phi);
        const double s_phi = std::sin(phi);

        if (rng.uniform() < lattice_.config().fTT_frac) {
            // Channel 1: L -> T1 + T2
            stats.decays_L_to_TT++;
            const double x   = lattice_.sample_tt_energy_fraction(rng);
            const double th1 = lattice_.theta_TT(x);
            const double th2 = lattice_.theta_TT(1.0 - x);

            const Vec3 k1 = dir1.rotate_around(perp, std::cos(th1), std::sin(th1)).rotate_around(dir1, c_phi, s_phi).unit();
            const Vec3 k2 = dir1.rotate_around(perp, std::cos(-th2), std::sin(-th2)).rotate_around(dir1, c_phi, s_phi).unit();

            const PhononMode m1 = lattice_.sample_transverse_mode(rng);
            const PhononMode m2 = lattice_.sample_transverse_mode(rng);

            // Push second transverse phonon onto stack
            stack.push_back(make_track(trk.pos, k2, m2, (1.0 - x) * trk.energy, trk.time));

            // Continue first transverse phonon in-place
            trk.energy = x * trk.energy;
            trk.mode   = m1;
            trk.k_dir  = k1;
            const auto kin1 = lattice_.map_k_to_vg(m1, k1);
            trk.vg_dir = kin1.vg_dir;
            trk.vg_mag = kin1.vg_mag;
        } else {
            // Channel 2: L -> L' + T
            stats.decays_L_to_LT++;
            const double x   = lattice_.sample_lt_energy_fraction(rng);
            const double thL = lattice_.theta_L_to_LT(x);
            const double thT = lattice_.theta_T_from_LT(x);

            const Vec3 kL = dir1.rotate_around(perp, std::cos(thL), std::sin(thL)).rotate_around(dir1, c_phi, s_phi).unit();
            const Vec3 kT = dir1.rotate_around(perp, std::cos(-thT), std::sin(-thT)).rotate_around(dir1, c_phi, s_phi).unit();

            const PhononMode mT = lattice_.sample_transverse_mode(rng);

            // Push transverse phonon onto stack
            stack.push_back(make_track(trk.pos, kT, mT, (1.0 - x) * trk.energy, trk.time));

            // Continue L' phonon in-place
            trk.energy = x * trk.energy;
            trk.mode   = PhononMode::L;
            trk.k_dir  = kL;
            const auto kinL = lattice_.map_k_to_vg(PhononMode::L, kL);
            trk.vg_dir = kinL.vg_dir;
            trk.vg_mag = kinL.vg_mag;
        }
    }
};

} // namespace cmpmc
