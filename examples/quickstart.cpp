// Minimal CMPMC program: propagate phonons in a germanium detector with
// aluminium sensor films on both faces and report where the energy went.
//
//   cmake -B build && cmake --build build && ./build/examples/quickstart

#include "cmpmc/cmpmc.hpp"

#include <cstdio>
#include <vector>

int main() {
    using namespace cmpmc;

    // 1. Material: built-in germanium preset. Building the lattice tabulates the
    //    phonon kinematics (group velocity vs. wavevector) and decay samplers.
    const CrystalLattice lattice(materials::Germanium());

    // 2. Detector: the default is a 76.2 mm x 25.4 mm cylinder, diffuse walls, and
    //    an Al sensor film (Kaplan quasiparticle cascade) on the top and bottom faces.
    DetectorModel detector;

    // 3. Engine and initial state: 200 longitudinal phonons at the Debye energy,
    //    emitted isotropically from the centre.
    const SimulationEngine engine(lattice, detector);
    Rng rng(2026);

    std::vector<PhononTrack> phonons;
    SimStats stats;
    const double energy = lattice.debye_energy();
    for (int i = 0; i < 200; ++i) {
        phonons.push_back(engine.make_track({0.0, 0.0, 0.0}, rng.isotropic_dir(), PhononMode::L, energy));
        stats.total_injected_eV += energy;
    }

    // 4. Transport until every track is absorbed or terminated.
    std::vector<SensorHit> hits;
    engine.transport(phonons, hits, stats, rng);

    std::printf("CMPMC %s\n", version_string);
    std::printf("injected            %10.4f eV\n", stats.total_injected_eV);
    std::printf("deposited in films  %10.4f eV (%zu sensor hits)\n", stats.total_deposited_eV, hits.size());
    std::printf("lost / below cutoff %10.4f eV\n", stats.escaped_or_killed_eV);
    std::printf("anharmonic decays   %10llu\n",
                static_cast<unsigned long long>(stats.decays_L_to_TT + stats.decays_L_to_LT));
    return 0;
}
