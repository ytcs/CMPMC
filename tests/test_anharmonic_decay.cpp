// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Anharmonic L-phonon decay in the isotropic approximation.
//  (a) daughter angles conserve crystal momentum (k_L = k_1 + k_2) for every allowed x;
//  (b) daughter energy-fraction distributions sampled from CMPMC's inverse-CDF tables
//      match Leman, arXiv:1109.1193, Eqs. (9) and (12), integrated independently here;
//  (c) energy conservation per decay (exact by construction, checked on samples).
#include "test_support.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include <vector>

using namespace cmpmc;

// Leman (9): L -> L' + T, x = E_L'/E_L, d = v_L/v_T, (d-1)/(d+1) < x < 1
static double leman_lt(double d, double x) {
    const double a = 1.0 - x * x;
    const double b = (1.0 + x) * (1.0 + x) - d * d * (1.0 - x) * (1.0 - x);
    const double c = 1.0 + x * x - d * d * (1.0 - x) * (1.0 - x);
    return (a * a / (x * x)) * b * c * c;
}
// Leman (12): L -> T + T in Leman's variable X = d * E_T1/E_L, (d-1)/2 < X < (d+1)/2;
// beta, gamma, lambda, mu in units of 1e11 Pa (the convention of the crystal configuration).
static double leman_tt(double d, double X, const CrystalConfig& c) {
    const double be = c.beta / 1e11, ga = c.gamma / 1e11, la = c.lambda / 1e11, mu = c.mu / 1e11;
    const double A = 0.5 * (1 - d * d) * (be + la + (1 + d * d) * (ga + mu));
    const double B = be + la + 2 * d * d * (ga + mu);
    const double C = be + la + 2 * (ga + mu);
    const double D = (1 - d * d) * (2 * be + 4 * ga + la + 3 * mu);
    const double t1 = A + B * d * X - B * X * X;
    const double t2 = C * X * (d - X) - D / (d - X) * (X - d - (1 - d * d) / (4 * X));
    return t1 * t1 + t2 * t2;
}

struct Moments { double m1 = 0, m2 = 0; };
template <class F>
static Moments pdf_moments(F&& p, double lo, double hi) {
    const int N = 200000; double s0 = 0, s1 = 0, s2 = 0; const double dx = (hi - lo) / N;
    for (int i = 0; i < N; ++i) { const double x = lo + (i + 0.5) * dx; const double w = std::max(0.0, p(x));
        s0 += w; s1 += w * x; s2 += w * x * x; }
    return {s1 / s0, s2 / s0 - (s1 / s0) * (s1 / s0)};
}

static void run(const CrystalConfig& cfg) {
    const CrystalLattice lat(cfg);
    const double vL = cfg.v_sound, vT = cfg.v_trans, d = vL / vT;
    const std::string n = cfg.name;

    // (a) momentum closure. Units: E_L = 1, k = E / v.
    double worst_lt = 0, worst_tt = 0;
    for (int i = 1; i < 400; ++i) {
        const double xl = (d - 1) / (d + 1) + ((1.0 - (d - 1) / (d + 1)) * i) / 400.0;
        const double tL = lat.theta_L_to_LT(xl), tT = lat.theta_T_from_LT(xl);
        const double kL = xl / vL, kT = (1 - xl) / vT, k0 = 1 / vL;
        const double par = kL * std::cos(tL) + kT * std::cos(tT) - k0;
        const double per = kL * std::sin(tL) - kT * std::sin(tT);
        worst_lt = std::max(worst_lt, std::hypot(par, per) / k0);
        const double xt = 0.5 * (1 - 1 / d) + (1 / d) * i / 400.0;
        const double t1 = lat.theta_TT(xt), t2 = lat.theta_TT(1 - xt);
        const double k1 = xt / vT, k2 = (1 - xt) / vT;
        const double par2 = k1 * std::cos(t1) + k2 * std::cos(t2) - k0;
        const double per2 = k1 * std::sin(t1) - k2 * std::sin(t2);
        worst_tt = std::max(worst_tt, std::hypot(par2, per2) / k0);
    }
    cmpmc_test::check(n + " L->L'T max |k1+k2-k0|/k0", worst_lt, 0.0, 0.0, 1e-9);
    cmpmc_test::check(n + " L->TT  max |k1+k2-k0|/k0", worst_tt, 0.0, 0.0, 1e-9);
    {
        const double cth = std::cos(lat.theta_TT(0.5));
        cmpmc_test::check(n + " L->TT cos(theta) at x=0.5 (Leman (13): 1/d)", cth, 1.0 / d, 1e-9);
    }

    // (b) energy-fraction distributions
    Rng rng(12345);
    const int NS = 2000000;
    double s1 = 0, s2 = 0, t1 = 0, t2 = 0, econs = 0;
    for (int i = 0; i < NS; ++i) {
        const double x = lat.sample_lt_energy_fraction(rng); s1 += x; s2 += x * x;
        const double y = lat.sample_tt_energy_fraction(rng); t1 += y; t2 += y * y;
        const double E = 1.0 + rng.uniform();
        econs = std::max(econs, std::abs((x * E + (1 - x) * E) - E) / E);
    }
    const double mlt = s1 / NS, vlt = s2 / NS - mlt * mlt, mtt = t1 / NS, vtt = t2 / NS - mtt * mtt;
    const Moments lt = pdf_moments([&](double x) { return leman_lt(d, x); }, (d - 1) / (d + 1), 1.0);
    const Moments tt = pdf_moments([&](double f) { return leman_tt(d, d * f, cfg); }, 0.5 * (1 - 1 / d), 0.5 * (1 + 1 / d));
    cmpmc_test::check(n + " L->L'T <x> vs Leman (9)", mlt, lt.m1, 2e-3);
    cmpmc_test::check(n + " L->L'T Var(x) vs Leman (9)", vlt, lt.m2, 2e-2);
    cmpmc_test::check(n + " L->TT  <x> vs Leman (12)", mtt, tt.m1, 2e-3);
    cmpmc_test::check(n + " L->TT  Var(x) vs Leman (12)", vtt, tt.m2, 2e-2);
    cmpmc_test::check(n + " energy conservation per decay", econs, 0.0, 0.0, 1e-14);
}

int main() {
    std::printf("# test_anharmonic_decay\n");
    run(materials::Germanium());
    run(materials::Silicon());
    return cmpmc_test::result();
}
