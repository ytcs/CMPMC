// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Fano binomial sampler and the energy partition.
//  (a) sample_fano_binomial: <N> = n_bar and Var N / n_bar = F for n_bar in [2, 1e4],
//      F in {0.119, 0.15, 0.30} (Fano's definition, arXiv:2004.10709, Eq. (3));
//  (b) partition_deposit energy bookkeeping at zero field:
//      omega = N_eh E_g + heat, heat = prompt phonons + carrier kinetic energy).
#include "test_support.hpp"
#include "cmpmc/charge/energy_partition.hpp"
#include <vector>

using namespace cmpmc;

int main() {
    std::printf("# test_fano_partition\n");
    Rng rng(77);
    for (double F : {0.119, 0.15, 0.30}) {
        for (double nb : {2.0, 3.5, 7.3, 20.0, 100.0, 1.0e4}) {
            const int NS = 400000; double s = 0, s2 = 0;
            for (int i = 0; i < NS; ++i) { const double n = double(physics::sample_fano_binomial(nb, F, rng)); s += n; s2 += n * n; }
            const double m = s / NS, v = s2 / NS - m * m;
            const std::string tag = "F=" + std::to_string(F).substr(0, 5) + " nbar=" + std::to_string(nb).substr(0, 6);
            cmpmc_test::check(tag + " <N>", m, nb, 0.01);
            cmpmc_test::check(tag + " Var/<N> = F", v / m, F, 0.05);
        }
    }
    const CrystalLattice lat(materials::Silicon());
    const EnergyPartition part(lat);
    PartitionOptions opt; opt.max_sampled_pairs = 1u << 30; opt.max_sampled_phonons = 1u << 30;
    for (double E : {50.0, 1000.0}) {
        std::vector<ChargeTrack> q; std::vector<PhononTrack> ph;
        const auto r = part.partition_deposit(Vec3{0, 0, 0}, 0.0, E, false, q, ph, rng, opt);
        double eq = 0, ep = 0; for (auto& c : q) eq += c.ekin * c.weight; for (auto& p : ph) ep += p.energy * p.weight;
        cmpmc_test::check("Si ER " + std::to_string(int(E)) + " eV: N E_g + E_kin + E_ph = omega",
                     r.n_pairs_true * lat.config().bandgap + eq + ep, E, 1e-9);
    }
    return cmpmc_test::result();
}
