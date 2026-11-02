// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors
//
// CMPMC side of the CMPMC-vs-G4CMP speed comparison (see benchmarks/README.md).
//
// The detector reproduces G4CMP's examples/phonon configuration, which is also
// what benchmarks/g4cmp/g4cmp_bench.cc builds: Ge cylinder R = 3.81 cm,
// half-height 1.27 cm, Al Kaplan films on both faces entered with probability 0.3,
// the example's polynomial reflection model on every surface, and G4CMP's
// default limit of 100 surface bounces per phonon. Single thread; the timed
// region is transport only.
//
// Usage: cmpmc_bench <mode> <events> [seed] [voltage_V] [kinematics]
//   mode       phonon | charge | chargefull (as in g4cmp_bench)
//   kinematics octahedral (default) | bilinear | exact

#include "cmpmc/cmpmc.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace cmpmc;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s phonon|charge|chargefull <events> [seed] [voltage_V] [kinematics]\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    const int events = std::atoi(argv[2]);
    const unsigned long long seed = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 12345ULL;
    const double voltage = argc > 4 ? std::atof(argv[4]) : 4.0;
    const std::string kin = argc > 5 ? argv[5] : "octahedral";

    const auto t0 = std::chrono::steady_clock::now();
    CrystalLattice lattice(materials::Germanium());
    if (kin == "bilinear") lattice.set_kinematics_mode(CrystalLattice::KinematicsMode::Bilinear);
    else if (kin == "exact") lattice.set_kinematics_mode(CrystalLattice::KinematicsMode::Exact);

    DetectorModel det;
    det.allow_unvalidated = true;   // the example's polynomial surface coefficients are guarded
    SurfaceProperty sensor = SurfaceProperty::demo_polynomial_sensor_surface();
    sensor.prob_abs = 0.3;
    // Charges are absorbed at every surface, the default of both codes (G4CMP's
    // charge-bounce limit of 0 turns any surface contact into absorption).
    const SurfaceProperty wall = SurfaceProperty::demo_polynomial_bare_surface();
    det.surfaces = {sensor, sensor, wall};
    det.max_bounces = 100;
    det.e_field = UniformElectricField::from_bias(mode == "phonon" ? 0.0 : voltage, 2.54 * units::cm);
    det.luke_sample_fraction = (mode == "charge") ? 0.0 : 1.0;

    const SimulationEngine engine(lattice, det);
    const ChargePhysics charges(lattice);
    Rng rng(seed);
    const auto t1 = std::chrono::steady_clock::now();

    SimStats stats;
    std::vector<SensorHit> hits;
    std::vector<PhononTrack> phonons;
    std::vector<ChargeTrack> carriers;
    unsigned long long tracks = 0;
    const size_t n_valleys = lattice.config().valley_dirs.size();

    for (int ev = 0; ev < events; ++ev) {
        if (mode == "phonon") {
            const double u = rng.uniform();
            const PhononMode m = u < 0.53539 ? PhononMode::ST : u < 0.90217 ? PhononMode::FT : PhononMode::L;
            phonons.push_back(engine.make_track({0, 0, 0}, rng.isotropic_dir(), m, 7.5e-3));
            stats.total_injected_eV += 7.5e-3;
            engine.transport(phonons, hits, stats, rng);
            ++tracks;
        } else {
            const auto valley = static_cast<uint16_t>(rng.uniform() * n_valleys) % n_valleys;
            carriers.push_back(charges.make_charge(ExcitationType::Electron, {0, 0, 0}, rng.isotropic_dir(), 1e-6, valley));
            carriers.push_back(charges.make_charge(ExcitationType::Hole, {0, 0, 0}, rng.isotropic_dir(), 1e-6, 0));
            engine.transport_charges(carriers, phonons, hits, stats, rng);
            tracks += 2 + phonons.size();
            if (!phonons.empty()) engine.transport(phonons, hits, stats, rng);
        }
    }
    const auto t2 = std::chrono::steady_clock::now();

    const double init_s = std::chrono::duration<double>(t1 - t0).count();
    const double run_s = std::chrono::duration<double>(t2 - t1).count();
    std::printf("code=CMPMC mode=%s kinematics=%s events=%d init_s=%.3f run_s=%.3f us_per_event=%.1f "
                "steps=%llu phonon_steps=%llu charge_steps=%llu tracks=%llu edep_eV=%.6f luke_phonons=%llu\n",
                mode.c_str(), kin.c_str(), events, init_s, run_s, 1e6 * run_s / events,
                static_cast<unsigned long long>(stats.total_steps + stats.charge_steps),
                static_cast<unsigned long long>(stats.total_steps),
                static_cast<unsigned long long>(stats.charge_steps), tracks, stats.total_deposited_eV,
                static_cast<unsigned long long>(stats.luke_phonons_emitted));
    return 0;
}
