// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Isotropic AMM / DMM phonon interface transmission (phonon/interface_transmission.hpp)
//  (a) P-SV oblique AMM: energy-flux conservation R_P + R_S + T_P + T_S = 1 at every angle and
//      both incident modes (not imposed by the solver: an independent check of the 4x4 boundary
//      problem), for Si/Al, Ge/Al, Si/W, Al2O3/Al;
//  (b) normal-incidence limit of the 4x4 solve = 4 Z1 Z2/(Z1 + Z2)^2 (closed form) and no mode
//      conversion at theta = 0;
//  (c) identical media: T = 1, R = 0 at all angles;
//  (d) reciprocity of SH transmission T_12(theta1) = T_21(theta2);
//  (e) reference table (AMM L-normal, T-normal, SH flux average, DMM) for six substrate/film
//      pairs, computed independently from the same inputs and quoted to 3 digits (tolerance
//      0.0006). Inputs: Si, Ge from Leman, arXiv:1109.1193, Table I; Al2O3 from arXiv:2408.04732,
//      Table 3, with the AFLOW density (database value); W, Ti: AFLOW (database values);
//      Al, Ir: moduli not independently verified, used as test inputs only;
//  (f) DMM detailed balance t_12 + t_21 = 1.
#include "test_support.hpp"
#include "cmpmc/phonon/interface_transmission.hpp"
#include <string>
#include <vector>

using namespace cmpmc::physics;

int main() {
    std::printf("# test_interface_amm\n");
    const IsotropicMedium Si{2330., 9000., 5400.}, Ge{5320., 5310., 3250.};
    const IsotropicMedium Al2O3 = IsotropicMedium::from_lame(3866., 139.8e9, 166.24e9);
    const IsotropicMedium Al = IsotropicMedium::from_bulk_shear(2720., 76e9, 26e9);
    const IsotropicMedium W = IsotropicMedium::from_bulk_shear(18560., 297.672e9, 150.536e9);
    const IsotropicMedium Ti = IsotropicMedium::from_bulk_shear(4265., 113.529e9, 63.5361e9);
    const IsotropicMedium Ir = IsotropicMedium::from_bulk_shear(21980., 320e9, 210e9);

    struct Pair { const char* n; IsotropicMedium s, f; double L, T, SH, D; };
    const std::vector<Pair> pairs = {
        {"Si/Al", Si, Al, 0.991, 0.960, 0.960, 0.743}, {"Si/W", Si, W, 0.588, 0.621, 0.482, 0.778},
        {"Si/Ti", Si, Ti, 0.974, 0.982, 0.896, 0.658}, {"Si/Ir", Si, Ir, 0.522, 0.527, 0.406, 0.752},
        {"Ge/Al", Ge, Al, 0.943, 0.881, 0.905, 0.510}, {"Al2O3/Al", Al2O3, Al, 0.822, 0.748, 0.858, 0.810}};

    for (const auto& pr : pairs) {
        const std::string n = pr.n;
        // (a)
        double worst = 0.0;
        for (int i = 0; i < 400; ++i) {
            const double th = (i + 0.5) / 400 * 1.5707963267948966;
            worst = std::max(worst, std::abs(amm_psv(pr.s, pr.f, true, th).sum() - 1.0));
            worst = std::max(worst, std::abs(amm_psv(pr.s, pr.f, false, th).sum() - 1.0));
            worst = std::max(worst, std::abs(amm_psv(pr.f, pr.s, true, th).sum() - 1.0));
            worst = std::max(worst, std::abs(amm_psv(pr.f, pr.s, false, th).sum() - 1.0));
        }
        cmpmc_test::check(n + " max |R+T-1| over angles, P & SV, both sides", worst, 0.0, 0.0, 1e-10);
        // (b)
        const auto P0 = amm_psv(pr.s, pr.f, true, 0.0), S0 = amm_psv(pr.s, pr.f, false, 0.0);
        cmpmc_test::check(n + " 4x4 P normal vs 4Z1Z2/(Z1+Z2)^2", P0.T_P, amm_normal(pr.s.rho * pr.s.vL, pr.f.rho * pr.f.vL), 1e-12);
        cmpmc_test::check(n + " 4x4 SV normal vs 4Z1Z2/(Z1+Z2)^2", S0.T_S, amm_normal(pr.s.rho * pr.s.vT, pr.f.rho * pr.f.vT), 1e-12);
        cmpmc_test::check(n + " no conversion at normal incidence", P0.T_S + P0.R_S + S0.T_P + S0.R_P, 0.0, 0.0, 1e-14);
        // (e)
        double sh = 0.0; const int NA = 20000;
        for (int i = 0; i < NA; ++i) { const double th = (i + 0.5) / NA * 1.5707963267948966;
            sh += amm_sh(pr.s, pr.f, th) * 2 * std::sin(th) * std::cos(th); }
        sh *= 1.5707963267948966 / NA;
        cmpmc_test::check(n + " AMM L normal vs reference table", P0.T_P, pr.L, 0.0, 6e-4);
        cmpmc_test::check(n + " AMM T normal vs reference table", S0.T_S, pr.T, 0.0, 6e-4);
        cmpmc_test::check(n + " AMM SH flux avg vs reference table", sh, pr.SH, 0.0, 6e-4);
        cmpmc_test::check(n + " DMM vs reference table", dmm(pr.s, pr.f), pr.D, 0.0, 6e-4);
        // (f)
        cmpmc_test::check(n + " DMM t12 + t21", dmm(pr.s, pr.f) + dmm(pr.f, pr.s), 1.0, 1e-14);
        const IsotropicAMMInterface I{pr.s, pr.f};
        cmpmc_test::info(n + " engine AMM flux-avg t_L (P-SV, oblique)", I.flux_average(true));
        cmpmc_test::info(n + " engine AMM flux-avg t_T ((SV+SH)/2)", I.flux_average(false));
    }
    // (c)
    double worst = 0.0;
    for (int i = 0; i < 100; ++i) {
        const double th = (i + 0.5) / 100 * 1.5;
        worst = std::max(worst, std::abs(amm_psv(Si, Si, true, th).T_P - 1.0));
        worst = std::max(worst, std::abs(amm_psv(Si, Si, false, th).T_S - 1.0));
        worst = std::max(worst, std::abs(amm_sh(Si, Si, th) - 1.0));
    }
    cmpmc_test::check("identical media: max |T - 1|", worst, 0.0, 0.0, 1e-12);
    // (d)
    const double th1 = 0.6, th2 = std::asin(Al.vT / Si.vT * std::sin(th1));
    cmpmc_test::check("SH reciprocity T_12(th1) = T_21(th2), Si/Al", amm_sh(Si, Al, th1), amm_sh(Al, Si, th2), 1e-12);
    // P-SV reciprocity: T_P(1->2, P at th1) equals T_P(2->1, P at th2) (Snell-related angles)
    const double thp2 = std::asin(Al.vL / Si.vL * std::sin(th1));
    cmpmc_test::check("P-SV reciprocity T_PP 12 vs 21, Si/Al", amm_psv(Si, Al, true, th1).T_P, amm_psv(Al, Si, true, thp2).T_P, 1e-10);
    return cmpmc_test::result();
}
