// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Energy closure of the phonon transport loop: for a short Ge run with the default
// detector, deposited + escaped/killed energy equals the injected energy.
#include "test_support.hpp"
#include "cmpmc/transport/engine.hpp"

using namespace cmpmc;

int main() {
    std::printf("# test_engine_energy_closure\n");
    const CrystalLattice lat(materials::Germanium());
    DetectorModel det;
    const SimulationEngine eng(lat, det);
    Rng rng(1);
    std::vector<PhononTrack> st; std::vector<SensorHit> hits; SimStats S;
    double Ein = 0;
    for (int i = 0; i < 2000; ++i) { st.push_back(eng.make_track({0, 0, 0}, rng.isotropic_dir(), PhononMode::L, 0.004)); Ein += 0.004; }
    eng.transport(st, hits, S, rng);
    cmpmc_test::check("Ge phonon run: deposited + escaped/killed = injected", S.total_deposited_eV + S.escaped_or_killed_eV, Ein, 1e-9);
    return cmpmc_test::result();
}
