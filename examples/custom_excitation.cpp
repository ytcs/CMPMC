// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Defining a new excitation species for SimulationEngine::transport_custom_excitation.
//
// A model is any type satisfying cmpmc::ExcitationModel: it provides a bulk
// interaction rate, the state update at a bulk interaction, and an energy cutoff.
// The engine handles free flight, boundary intersection and sensor absorption.
//
// The two models below are illustrative, not validated physics.

#include "cmpmc/cmpmc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace cmpmc;

namespace {

// ----------------------------------------------------------------------------
// Example 1: 3D Superconducting Bogoliubov Quasiparticle Model
// BCS quasiparticle group velocity:
//   v_g(E) = v_F * sqrt(1 - (Delta / E)^2)
// ----------------------------------------------------------------------------
struct BogoliubovQPModel {
    double gap_energy_eV     = 173.715 * units::ueV; // Aluminum Delta_0
    double fermi_velocity    = 2.03e6;               // Al Fermi velocity (m/s)
    double elastic_mfp_m     = 60.0 * units::nm;     // Impurity/boundary elastic MFP
    double trap_rate_Hz      = 1.0e5;                // Quasiparticle trapping rate (1/s)

    constexpr double min_energy_eV() const noexcept { return gap_energy_eV; }

    double group_velocity(double energy_eV) const noexcept {
        if (energy_eV <= gap_energy_eV) return 0.0;
        const double ratio = gap_energy_eV / energy_eV;
        return fermi_velocity * std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
    }

    double bulk_rate(const GenericExcitationTrack& trk) const noexcept {
        const double vg = group_velocity(trk.energy);
        const double r_elastic = (elastic_mfp_m > 0.0) ? (vg / elastic_mfp_m) : 0.0;
        return r_elastic + trap_rate_Hz;
    }

    void perform_bulk_interaction(GenericExcitationTrack& trk,
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

    double group_velocity(double energy_eV) const noexcept {
        if (energy_eV <= zeeman_gap_eV) return 0.0;
        return (2.0 / units::hbar_Planck) * std::sqrt(exchange_stiff_eVm2 * (energy_eV - zeeman_gap_eV));
    }

    double bulk_rate(const GenericExcitationTrack& trk) const noexcept {
        return group_velocity(trk.energy) / magnon_mfp_m;
    }

    void perform_bulk_interaction(GenericExcitationTrack& trk,
                                         std::vector<GenericExcitationTrack>& /*secondaries*/,
                                         Rng& rng) const noexcept {
        trk.mom_dir  = rng.isotropic_dir();
        trk.velocity = group_velocity(trk.energy);
    }
};

static_assert(cmpmc::ExcitationModel<BogoliubovQPModel>);
static_assert(cmpmc::ExcitationModel<MagnonModel>);

} // namespace

template <ExcitationModel Model>
void run(const char* label, const Model& model, const SimulationEngine& engine, double energy_eV) {
    Rng rng(42);
    std::vector<GenericExcitationTrack> stack;
    for (int i = 0; i < 200; ++i) {
        GenericExcitationTrack trk{};
        trk.pos      = {0.0, 0.0, 0.0};
        trk.mom_dir  = rng.isotropic_dir();
        trk.energy   = energy_eV;
        trk.velocity = model.group_velocity(energy_eV);
        trk.weight   = 1.0;
        stack.push_back(trk);
    }
    std::vector<SensorHit> hits;
    SimStats stats;
    engine.transport_custom_excitation(model, stack, hits, stats, rng);
    std::printf("%-14s steps %10llu  sensor hits %6zu  absorbed in bulk %6llu\n", label,
                static_cast<unsigned long long>(stats.total_steps), hits.size(),
                static_cast<unsigned long long>(stats.killed_tracks));
}

int main() {
    const CrystalLattice lattice(materials::Germanium());
    DetectorModel detector;
    const SimulationEngine engine(lattice, detector);

    // Mean free paths are chosen so that a track crosses the 2.54 cm crystal in a
    // modest number of steps; a realistic 60 nm QP mean free path would need ~1e12.
    BogoliubovQPModel qp;
    qp.elastic_mfp_m = 2.0 * units::mm;
    qp.trap_rate_Hz  = 2.0e4;
    run("Bogoliubov QP", qp, engine, 300.0 * units::ueV);
    run("Magnon", MagnonModel{}, engine, 120.0 * units::ueV);
    return 0;
}
