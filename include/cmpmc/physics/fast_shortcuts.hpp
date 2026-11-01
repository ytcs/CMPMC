#pragma once

#include "cmpmc/core.hpp"
#include "cmpmc/lattice.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cmpmc::physics {

// ============================================================================
// Tier-4 Algorithmic Physics Shortcuts (10x - 100x Step Reduction)
//   1. JumpToBallisticKernel: Asymptotic Tamura anharmonic downconversion
//      shortcut that jumps high-frequency diffusive phonons (where
//      lambda_MFP << clearance to boundary) directly to their quasi-ballistic
//      daughter phonon spectrum + 3D diffusive displacement (dr, dt).
//   2. ChordLukeStepper: Macroscopic chord stepping for high-voltage (HV /
//      CDMSlite) Neganov-Trofimov-Luke phonon emission along E-field lines.
//   3. RussianRoulette: Unbiased statistical survival weighting for high-bounce
//      phonons in sparse sensor geometries (qubits / KIDs).
// ============================================================================

struct JumpToBallisticConfig {
    bool enabled              = false;
    double safety_margin_mfp  = 8.0;   // Require clearance to wall >= 8 * lambda_tot(E)
    double target_freq_GHz    = 600.0; // Quasi-ballistic transition frequency (~0.6 THz in Ge/Si)
    uint32_t max_cascade_gen  = 16;    // Max in-register cascade generations per jump
};

// Perform an in-register diffusive + anharmonic downconversion cascade for a
// high-frequency phonon located at least `safety_margin_mfp * lambda_tot` inside
// the crystal bulk, without touching the geometry navigator on intermediate steps.
// Preserves exact energy conservation (sum of daughter energies == trk.energy).
inline uint32_t jump_to_ballistic_cascade(PhononTrack& trk,
                                          double wall_clearance_m,
                                          const CrystalLattice& lattice,
                                          const JumpToBallisticConfig& cfg,
                                          std::vector<PhononTrack>& out_stack,
                                          uint64_t& stat_decays_TT,
                                          uint64_t& stat_decays_LT,
                                          uint64_t& stat_iso_scatters,
                                          Rng& rng) noexcept {
    const double e_cut = cfg.target_freq_GHz * units::GHz * units::h_Planck;
    if (!cfg.enabled || trk.energy <= e_cut) return 0;

    const double r_iso0  = lattice.isotope_rate(trk.energy);
    const double r_down0 = lattice.downconversion_rate(trk.mode, trk.energy);
    const double r_tot0  = r_iso0 + r_down0;
    if (r_tot0 <= 1e-6) return 0;

    const double mfp0 = trk.vg_mag / r_tot0;
    if (wall_clearance_m < cfg.safety_margin_mfp * mfp0) {
        return 0; // Too close to crystal boundary; use exact ray-traced stepping
    }

    struct MiniPhonon {
        Vec3 pos;
        double energy;
        double time;
        PhononMode mode;
    };

    SmallVec<MiniPhonon, 64> local_q;
    SmallVec<MiniPhonon, 64> final_q;
    local_q.push_back({trk.pos, trk.energy, trk.time, trk.mode});

    uint32_t fast_steps = 0;
    const double fTT   = lattice.config().fTT_frac;
    const double dos_L = std::max(0.05, lattice.config().fLDOS);
    const double v_avg = 0.3333333333333333 * (lattice.config().v_sound + 2.0 * lattice.config().v_trans);

    while (!local_q.empty() && (final_q.size() + local_q.size() < 48)) {
        MiniPhonon p = local_q.back();
        local_q.pop_back();

        if (p.energy <= e_cut) {
            final_q.push_back(p);
            continue;
        }

        const double r_iso  = lattice.isotope_rate(p.energy);
        // Above e_cut, rapid isotope scattering mixes modes (T <-> L with fraction dos_L),
        // so the effective downconversion rate of the mode-mixed state is dos_L * r_down(L)
        const double r_down_L = lattice.downconversion_rate(PhononMode::L, p.energy);
        const double r_down   = (p.mode == PhononMode::L) ? r_down_L : (dos_L * r_down_L);
        if (r_down <= 1e-6) {
            final_q.push_back(p);
            continue;
        }

        const double dt_decay = -std::log(rng.uniform()) / r_down;
        const double n_iso_expected = r_iso * dt_decay;
        const uint64_t n_iso = static_cast<uint64_t>(n_iso_expected + 0.5);
        stat_iso_scatters += n_iso;
        fast_steps        += static_cast<uint32_t>(1 + n_iso);

        const double mfp_eff = v_avg / std::max(r_iso + r_down, 1.0);
        const double rms_dr  = std::min(0.25 * wall_clearance_m,
                                        std::sqrt(std::max(0.0, n_iso_expected)) * mfp_eff);
        p.pos  += rng.isotropic_dir() * rms_dr;
        p.time += dt_decay;

        if (rng.uniform() < fTT) {
            stat_decays_TT++;
            const double x  = lattice.sample_tt_energy_fraction(rng);
            const double e1 = x * p.energy;
            const double e2 = (1.0 - x) * p.energy;
            const PhononMode m1 = lattice.sample_transverse_mode(rng);
            const PhononMode m2 = lattice.sample_transverse_mode(rng);
            if (e1 > e_cut) local_q.push_back({p.pos, e1, p.time, m1});
            else            final_q.push_back({p.pos, e1, p.time, m1});
            if (e2 > e_cut) local_q.push_back({p.pos, e2, p.time, m2});
            else            final_q.push_back({p.pos, e2, p.time, m2});
        } else {
            stat_decays_LT++;
            const double x  = lattice.sample_lt_energy_fraction(rng);
            const double e1 = x * p.energy;
            const double e2 = (1.0 - x) * p.energy;
            const PhononMode m2 = lattice.sample_transverse_mode(rng);
            if (e1 > e_cut) local_q.push_back({p.pos, e1, p.time, PhononMode::L});
            else            final_q.push_back({p.pos, e1, p.time, PhononMode::L});
            if (e2 > e_cut) local_q.push_back({p.pos, e2, p.time, m2});
            else            final_q.push_back({p.pos, e2, p.time, m2});
        }
    }

    while (!local_q.empty()) {
        final_q.push_back(local_q.back());
        local_q.pop_back();
    }

    if (final_q.empty()) return 0;

    // Emit daughters 1..N-1 onto stack, and overwrite `trk` in-place with daughter 0
    for (size_t i = 1; i < final_q.size(); ++i) {
        const Vec3 kd = rng.isotropic_dir();
        const auto kin = lattice.map_k_to_vg_oct(final_q[i].mode, kd);
        out_stack.push_back(PhononTrack{
            .pos     = final_q[i].pos,
            .k_dir   = kd,
            .vg_dir  = kin.vg_dir,
            .energy  = final_q[i].energy,
            .vg_mag  = kin.vg_mag,
            .time    = final_q[i].time,
            .weight  = trk.weight,
            .bounces = trk.bounces,
            .mode    = final_q[i].mode
        });
    }

    const Vec3 kd0 = rng.isotropic_dir();
    const auto kin0 = lattice.map_k_to_vg_oct(final_q[0].mode, kd0);
    trk.pos    = final_q[0].pos;
    trk.k_dir  = kd0;
    trk.vg_dir = kin0.vg_dir;
    trk.energy = final_q[0].energy;
    trk.vg_mag = kin0.vg_mag;
    trk.time   = final_q[0].time;
    trk.mode   = final_q[0].mode;

    return fast_steps;
}

// ============================================================================
// 2. Russian Roulette Variance Reduction for High-Bounce Tracks
//    Returns true if track survives (with weight multiplied by 1/p_surv),
//    or false if track is terminated.
// ============================================================================
inline bool apply_russian_roulette(double& weight,
                                   uint32_t bounces,
                                   uint32_t threshold_bounces,
                                   double survival_prob,
                                   Rng& rng) noexcept {
    if (threshold_bounces == 0 || bounces <= threshold_bounces) return true;
    if ((bounces - threshold_bounces) % 32 != 0) return true;
    if (rng.uniform() < survival_prob) {
        weight /= survival_prob;
        return true;
    }
    return false;
}

} // namespace cmpmc::physics
