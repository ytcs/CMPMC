// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Frequency- and film-loss-dependent pair-breaking efficiency eta_pb(Omega, tau_l)
// of the Kaplan cascade with rate-competition film loss (FilmLossModel::KaplanRateCompetition).
//  (a) exact T = 0 Kaplan pair-breaking rate 1/tau_pb(Omega) = I(Omega)/(pi D tau0):
//      1/tau_pb(2D) = 1/tau0 (closed form) and I(Omega) against an independent brute-force
//      quadrature (sin^2 substitution, 2e5 nodes); large-Omega limit I -> Omega - 2D + O(D ln).
//  (b) Omega < 4D: N_qp = 2 exactly, eta = 2D/Omega for every tau_l (only one pair can be
//      broken; arXiv:1401.1937, Fig. 4: all curves coincide below 4D).
//  (c) eta_s(h nu/D; tau_l/tau_0^phi) against arXiv:1401.1937, Fig. 4 (Chang-Scalapino kinetic
//      equations, Al, T/Tc = 0.1, probe off), digitised from the published figure (pixel
//      read-out, +-0.005). Photon seed (case II coherence factor, arXiv:1401.1937, Eq. (1)).
//      Tolerance 0.015 absolute for tau_l/tau_0^phi = 0.5-4; the thick-film curve (8.0) is
//      checked at 0.045: the kinetic-equation solution sits at its own no-loss value there
//      (text: "tau_l/tau_0^phi = 8 ... eta_s ~ 0.59"), the event cascade still loses 6-11%
//      of 2D-4D phonons (documented model difference: E = D assumption for the excess, +4% per
//      arXiv:1401.1937, Sec. V).
//  (d) no-loss asymptote 0.59 (arXiv:1401.1937; literature range 0.57-0.6) and the F_qp report.
//  (e) eta decreases monotonically as the film gets thinner (tau_l down) at fixed Omega >= 5D;
//  (f) cascade energy closure at 50-350 D (checks the cascade buffer capacity).
#include "test_support.hpp"
#include "cmpmc/superconductor/kaplan_film.hpp"
#include <string>
#include <vector>

using namespace cmpmc;

static double brute_I(double om, double d) {
    const int N = 200000; double s = 0.0;
    for (int i = 0; i < N; ++i) {
        const double u = (i + 0.5) / N, sn = std::sin(units::half_pi * u);
        const double E = d + (om - 2 * d) * sn * sn;
        const double jac = (om - 2 * d) * units::pi * sn * std::cos(units::half_pi * u) / N;
        s += jac * (E * (om - E) + d * d) / std::sqrt((E * E - d * d) * ((om - E) * (om - E) - d * d));
    }
    return s;
}

int main() {
    std::printf("# test_pair_breaking\n");
    KaplanFilm f = KaplanFilm::aluminum_sensor();
    f.loss_model = FilmLossModel::KaplanRateCompetition;
    const double D = f.gap_energy, t0 = f.phonon_lifetime;

    // (a)
    cmpmc_test::check("tau0/tau_pb(2D) (closed form 1)", f.kaplan_pair_breaking_rate(2.0 * D * (1 + 1e-12)) * t0, 1.0, 1e-9);
    for (double W : {2.5, 3.0, 5.0, 10.0, 40.0})
        cmpmc_test::check("I(Omega)/D vs brute quadrature, Omega=" + std::to_string(W).substr(0, 4) + "D",
                     f.kaplan_pb_integral(W * D) / D, brute_I(W * D, D) / D, 1e-5);
    cmpmc_test::info("tau0/tau_pb(100D) * pi / (Omega/D) (-> 1)", f.kaplan_pair_breaking_rate(100 * D) * t0 * units::pi / 100.0);

    Rng rng(20261101);
    const int N = 40000;
    // (b)
    for (double tl : {0.5, 8.0}) {
        f.phonon_loss_time = tl * t0;
        for (double W : {2.2, 3.0, 3.9}) {
            const auto r = f.pair_breaking_efficiency(W * D, N, rng, KaplanFilm::Seed::Photon);
            cmpmc_test::check("eta photon, tau_l/tau0=" + std::to_string(tl).substr(0, 3) + ", W=" + std::to_string(W).substr(0, 3) + " vs 2D/Omega",
                         r.eta, 2.0 / W, 1e-12, 1e-12);
            cmpmc_test::check("  F_qp (N_qp = 2 exactly)", r.F_qp, 0.0, 0.0, 1e-12);
        }
    }
    // (c) arXiv:1401.1937, Fig. 4, digitised: rows tau_l/tau_0^phi, columns h nu/D = 5, 6, 7, 8, 10
    const double hv[5] = {5, 6, 7, 8, 10};
    const double tls[5] = {0.5, 1.0, 2.0, 4.0, 8.0};
    const double fig4[5][5] = {{0.475, 0.437, 0.409, 0.385, 0.354},
                               {0.504, 0.479, 0.457, 0.442, 0.416},
                               {0.534, 0.520, 0.506, 0.495, 0.479},
                               {0.564, 0.558, 0.545, 0.540, 0.533},
                               {0.605, 0.605, 0.597, 0.594, 0.594}};
    double eta_tab[5][5];
    for (int a = 0; a < 5; ++a) {
        f.phonon_loss_time = tls[a] * t0;
        for (int b = 0; b < 5; ++b) {
            const auto r = f.pair_breaking_efficiency(hv[b] * D, N, rng, KaplanFilm::Seed::Photon);
            eta_tab[a][b] = r.eta;
            cmpmc_test::check("eta_s tau_l/tau0=" + std::to_string(tls[a]).substr(0, 3) + " hv=" + std::to_string((int)hv[b]) + "D vs arXiv:1401.1937 Fig.4",
                         r.eta, fig4[a][b], 0.0, (tls[a] < 5.0) ? 0.015 : 0.045);
        }
    }
    // (d) no loss
    f.phonon_loss_time = 0.0;
    for (double W : {10.0, 20.0}) {
        const auto r = f.pair_breaking_efficiency(W * D, N, rng, KaplanFilm::Seed::Photon);
        cmpmc_test::check("eta_s no loss, hv=" + std::to_string((int)W) + "D vs 0.59 (arXiv:1401.1937)", r.eta, 0.59, 0.0, 0.01);
        cmpmc_test::info("  F_qp no loss", r.F_qp);
        const auto p = f.pair_breaking_efficiency(W * D, N, rng, KaplanFilm::Seed::Phonon);
        cmpmc_test::info("  eta (phonon seed, case I) no loss", p.eta);
    }
    f.phonon_loss_time = 1.0 * t0;
    {
        const auto p = f.pair_breaking_efficiency(10 * D, N, rng, KaplanFilm::Seed::Phonon);
        cmpmc_test::info("eta phonon seed incl. entry loss, tau_l/tau0=1, 10D", p.eta);
        cmpmc_test::info("  F_qp", p.F_qp);
    }
    // (e) monotonic in tau_l
    for (int b = 0; b < 5; ++b) {
        bool mono = true;
        for (int a = 1; a < 5; ++a) mono = mono && (eta_tab[a][b] > eta_tab[a - 1][b]);
        cmpmc_test::check("eta increases with tau_l at hv=" + std::to_string((int)hv[b]) + "D", mono ? 1.0 : 0.0, 1.0, 0.0);
    }
    // (f) energy closure of the cascade at Debye-scale energies (buffer capacity): absorbed +
    //     re-emitted = incident for every event (mean-free-path loss law, thin film, sub-gap re-emission)
    {
        KaplanFilm g = KaplanFilm::aluminum_sensor();
        g.subgap_absorption = 0.0; g.temperature = 0.0;
        KaplanBuffer refl;
        for (double W : {50.0, 150.0, 210.0, 350.0}) {
            double worst = 0.0;
            for (int i = 0; i < 2000; ++i) {
                double e = g.absorb_phonon(W * D, refl, rng);
                for (double r : refl) e += r;
                worst = std::max(worst, std::abs(e / (W * D) - 1.0));
            }
            cmpmc_test::check("cascade energy closure, Omega=" + std::to_string((int)W) + "D (max over 2000)", worst, 0.0, 0.0, 1e-12);
        }
    }
    return cmpmc_test::result();
}
