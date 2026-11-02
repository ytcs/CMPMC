// Material presets against independent references.
//  (a) isotope constant B rebuilt independently, B = 4 pi^3 Gamma_md Omega <v^-3>,
//      <v^-3> = (1/3) sum_modes <v_p^-3>_{4pi} from CMPMC's own Christoffel solver and
//      Gamma_md from IUPAC isotope abundances (2408.04732 Eqs. (4)-(9) [FT]);
//  (b) published values: Si eps_eh, F_eh (WP5-v4 Sec. 6.4: 3.749 eV, 0.119 at 5.9 keV),
//      Debye (max optical) frequencies (2408.04732 Table 5 [FT]: Si 15.28, Ge 9.02,
//      GaAs 8.58 THz), GaAs A, F_TT (2408.04732 Table 3).
#include "test_support.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include <vector>

using namespace cmpmc;

struct Iso { double m, f; };
static double gamma_md(const std::vector<std::vector<Iso>>& sites) {
    // 2408.04732 Eqs. (5)-(9) with equal stoichiometry per listed site
    double Mbar = 0, dM2 = 0;
    for (const auto& s : sites) { double Mn = 0; for (auto i : s) Mn += i.f * i.m;
        double d = 0; for (auto i : s) d += i.f * (i.m - Mn) * (i.m - Mn); Mbar += Mn; dM2 += d; }
    Mbar /= sites.size(); dM2 /= sites.size();
    return dM2 / (Mbar * Mbar);
}

static double v3inv(const CrystalConfig& cfg) {
    const auto V = cfg.build_voigt_matrix(); double C[3][3][3][3];
    auto vi = [](int i, int j) { if (i == j) return i; if (i + j == 3) return 3; if (i + j == 2) return 4; return 5; };
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m)
        C[i][j][l][m] = V[vi(i, j)][vi(l, m)];
    Rng rng(5); double s = 0; const int N = 400000;
    for (int t = 0; t < N; ++t) { Vec3 vg[3]; double vp[3];
        physics::solve_christoffel_exact(C, cfg.density, rng.isotropic_dir(), vg, vp);
        for (int m = 0; m < 3; ++m) s += 1.0 / (vp[m] * vp[m] * vp[m]); }
    return s / (3.0 * N);
}

int main() {
    std::printf("# test_presets\n");
    const double pi3 = units::pi * units::pi * units::pi;
    // IUPAC 2013 representative abundances and atomic masses (u)
    const std::vector<Iso> Si = {{27.97693, 0.92223}, {28.97649, 0.04685}, {29.97377, 0.03092}};
    const std::vector<Iso> Ge = {{69.92425, 0.2057}, {71.92208, 0.2745}, {72.92346, 0.0775},
                                 {73.92118, 0.3650}, {75.92140, 0.0773}};
    const std::vector<Iso> Ga = {{68.92557, 0.60108}, {70.92470, 0.39892}};
    const std::vector<Iso> As = {{74.92160, 1.0}};
    struct Case { CrystalConfig c; double a; std::vector<std::vector<Iso>> s; };
    for (auto& k : {Case{materials::Silicon(), 5.431e-10, {Si}}, Case{materials::Germanium(), 5.658e-10, {Ge}},
                    Case{materials::GalliumArsenide(), 5.653e-10, {Ga, As}}}) {
        const double g = gamma_md(k.s), Om = k.a * k.a * k.a / 8.0, B = 4 * pi3 * g * Om * v3inv(k.c);
        cmpmc_test::info(k.c.name + " Gamma_md (independent)", g);
        cmpmc_test::check(k.c.name + " B_iso preset vs 4pi^3 Gamma Omega <v^-3>", k.c.B_iso, B, 0.10);
    }

    const auto si = materials::Silicon(), ge = materials::Germanium(), ga = materials::GalliumArsenide();
    cmpmc_test::check("Si pair energy vs campaign 3.749 eV (WP5-v4 6.4)", si.pair_energy, 3.749, 0.005);
    cmpmc_test::check("Si Fano vs campaign 0.119 (WP5-v4 6.4)", si.fano_factor, 0.119, 0.05);
    cmpmc_test::check("Si Debye (THz) vs 2408.04732 T5 15.28", si.debye_freq / units::THz, 15.28, 0.03);
    cmpmc_test::check("Ge Debye (THz) vs campaign 8.64 (WP8-R T1.1)", ge.debye_freq / units::THz, 8.64, 0.01);
    cmpmc_test::check("GaAs Debye (THz) vs 2408.04732 T5 8.58", ga.debye_freq / units::THz, 8.58, 0.05);
    cmpmc_test::check("GaAs A_anh vs 2408.04732 T3 7.77e-55", ga.A_anh, 7.77e-55, 0.02);
    cmpmc_test::check("GaAs F_TT vs 2408.04732 T3 0.77", ga.fTT_frac, 0.77, 0.02);
    return cmpmc_test::result();
}
