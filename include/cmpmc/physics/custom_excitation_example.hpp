#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/core/units.hpp"
#include <cmath>
#include <concepts>
#include <vector>

namespace cmpmc::physics {

// ============================================================================
// How to Add a New Excitation Species in CMPMC
//
// Any struct that satisfies `ExcitationProcessConcept` can be plugged directly
// into `SimulationEngine::transport_custom_excitation<Model>(...)` with ZERO
// virtual-dispatch overhead!
//
// Below are two ready-to-use examples:
//   1. `BogoliubovQPModel` : 3D Bogoliubov quasiparticle diffusion & trapping
//   2. `MagnonModel`       : Ferromagnetic spin-wave dispersion & Gilbert damping
// ============================================================================

template <typename Model>
concept ExcitationProcessConcept = requires(const Model& m,
                                            GenericExcitationTrack& trk,
                                            std::vector<GenericExcitationTrack>& secondaries,
                                            Rng& rng) {
    // 1. Bulk interaction rate Gamma(E) in [1/s]
    { m.bulk_rate(trk) } -> std::convertible_to<double>;

    // 2. Execute bulk scattering / decay / pair-breaking in-place
    { m.perform_bulk_interaction(trk, secondaries, rng) } -> std::same_as<void>;

    // 3. Minimum survival energy cutoff (eV)
    { m.min_energy_eV() } -> std::convertible_to<double>;
};

// ----------------------------------------------------------------------------
// Example 1: 3D Superconducting Bogoliubov Quasiparticle Model
// Reference: G4CMPBogoliubovQP / BCS Quasiparticle Group Velocity
//   v_g(E) = v_F * sqrt(1 - (Delta / E)^2)
// ----------------------------------------------------------------------------
struct BogoliubovQPModel {
    double gap_energy_eV     = 173.715 * units::ueV; // Aluminum Delta_0
    double fermi_velocity    = 2.03e6;               // Al Fermi velocity (m/s)
    double elastic_mfp_m     = 60.0 * units::nm;     // Impurity/boundary elastic MFP
    double trap_rate_Hz      = 1.0e5;                // Quasiparticle trapping rate (1/s)

    constexpr double min_energy_eV() const noexcept { return gap_energy_eV; }

    inline double group_velocity(double energy_eV) const noexcept {
        if (energy_eV <= gap_energy_eV) return 0.0;
        const double ratio = gap_energy_eV / energy_eV;
        return fermi_velocity * std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
    }

    inline double bulk_rate(const GenericExcitationTrack& trk) const noexcept {
        const double vg = group_velocity(trk.energy);
        const double r_elastic = (elastic_mfp_m > 0.0) ? (vg / elastic_mfp_m) : 0.0;
        return r_elastic + trap_rate_Hz;
    }

    inline void perform_bulk_interaction(GenericExcitationTrack& trk,
                                         std::vector<GenericExcitationTrack>& /*secondaries*/,
                                         Rng& rng) const noexcept {
        const double tot = bulk_rate(trk);
        if (tot > 0.0 && rng.uniform() * tot < trap_rate_Hz) {
            // Quasiparticle trapped below gap
            trk.energy = 0.0;
            return;
        }
        // Diffusive elastic impurity scatter
        trk.mom_dir  = rng.isotropic_dir();
        trk.velocity = group_velocity(trk.energy);
    }
};

// ----------------------------------------------------------------------------
// Example 2: Ferromagnetic Magnon (Spin-Wave) Model (e.g. YIG Calorimeters)
// Quadratic exchange dispersion: E(k) = D_ex * k^2 + E_Zeeman
// ----------------------------------------------------------------------------
struct MagnonModel {
    double zeeman_gap_eV      = 50.0 * units::ueV;
    double exchange_stiff_eVm2= 5.0e-18; // D_ex in eV * m^2
    double magnon_mfp_m       = 1.0 * units::mm;

    constexpr double min_energy_eV() const noexcept { return zeeman_gap_eV; }

    inline double group_velocity(double energy_eV) const noexcept {
        if (energy_eV <= zeeman_gap_eV) return 0.0;
        return (2.0 / units::hbar_Planck) * std::sqrt(exchange_stiff_eVm2 * (energy_eV - zeeman_gap_eV));
    }

    inline double bulk_rate(const GenericExcitationTrack& trk) const noexcept {
        return group_velocity(trk.energy) / magnon_mfp_m;
    }

    inline void perform_bulk_interaction(GenericExcitationTrack& trk,
                                         std::vector<GenericExcitationTrack>& /*secondaries*/,
                                         Rng& rng) const noexcept {
        trk.mom_dir  = rng.isotropic_dir();
        trk.velocity = group_velocity(trk.energy);
    }
};

static_assert(ExcitationProcessConcept<BogoliubovQPModel>);
static_assert(ExcitationProcessConcept<MagnonModel>);

} // namespace cmpmc::physics
