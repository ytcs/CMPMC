// Audit test: material presets.
//  (a) every phonon/charge constant of the presets vs G4CMP CrystalMaps/<mat>/config.txt
//      (the port claim), read with CMPMC's own load_g4cmp_config_file; skipped if the
//      G4CMP tree is absent (set G4CMP_CRYSTALMAPS);
//  (b) isotope constant B rebuilt independently, B = 4 pi^3 Gamma_md Omega <v^-3>,
//      <v^-3> = (1/3) sum_modes <v_p^-3>_{4pi} from CMPMC's own Christoffel solver and
//      Gamma_md from IUPAC isotope abundances (2408.04732 Eqs. (4)-(9) [FT]);
//  (c) campaign values: Si eps_eh, F_eh (WP5-v4 Sec. 6.4: 3.749 eV, 0.119 at 5.9 keV),
//      Debye (max optical) frequencies (2408.04732 Table 5 [FT]: Si 15.28, Ge 9.02,
//      GaAs 8.58 THz), GaAs A, F_TT (2408.04732 Table 3).
#include "audit_common.hpp"
#include "cmpmc/lattice.hpp"
#include "cmpmc/materials/config_reader.hpp"
#include <cstdlib>
#include <filesystem>
#include <vector>

using namespace cmpmc;

static void cmp_g4cmp(const CrystalConfig& p, const std::string& dir) {
    const std::string f = dir + "/" + p.name + "/config.txt";
    if (!std::filesystem::exists(f)) { std::printf("# skip %s (no %s)\n", p.name.c_str(), f.c_str()); return; }
    const CrystalConfig g = materials::load_g4cmp_config_file(f, p.density, p.name);
    const std::string n = p.name + " preset vs G4CMP config: ";
    const auto Vp = p.build_voigt_matrix(), Vg = g.build_voigt_matrix();
    double dv = 0, vmax = 0;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { dv = std::max(dv, std::abs(Vp[i][j] - Vg[i][j])); vmax = std::max(vmax, std::abs(Vg[i][j])); }
    audit::check(n + "max |C_ij - C_ij(G4CMP)| / max C", dv / vmax, 0.0, 0.0, 1e-6);
    audit::check(n + "B_iso", p.B_iso, g.B_iso, 1e-6, 1e-60);
    audit::check(n + "A_anh", p.A_anh, g.A_anh, 1e-6);
    audit::check(n + "fTT", p.fTT_frac, g.fTT_frac, 1e-6);
    audit::check(n + "LDOS", p.fLDOS, g.fLDOS, 1e-6);
    audit::info(n + "Debye (preset; G4CMP value next line)", p.debye_freq);
    audit::info(n + "Debye (G4CMP; deliberate deviation for Ge, GaAs)", g.debye_freq);
    audit::check(n + "v_L", p.v_sound, g.v_sound, 1e-6);
    audit::check(n + "v_T", p.v_trans, g.v_trans, 1e-6);
}

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
    const char* env = std::getenv("G4CMP_CRYSTALMAPS");
    const std::string dir = env ? env : "/home/user/ytcs/G4CMP-ref/CrystalMaps";
    for (const auto& c : {materials::Germanium(), materials::Silicon(), materials::Sapphire(),
                          materials::CalciumTungstate(), materials::GalliumArsenide(),
                          materials::LithiumFluoride(), materials::CalciumFluoride(),
                          materials::AluminumCrystal(), materials::NiobiumCrystal()})
        cmp_g4cmp(c, dir);

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
        audit::info(k.c.name + " Gamma_md (independent)", g);
        audit::check(k.c.name + " B_iso preset vs 4pi^3 Gamma Omega <v^-3>", k.c.B_iso, B, 0.10);
    }

    const auto si = materials::Silicon(), ge = materials::Germanium(), ga = materials::GalliumArsenide();
    audit::check("Si pair energy vs campaign 3.749 eV (WP5-v4 6.4)", si.pair_energy, 3.749, 0.005);
    audit::check("Si Fano vs campaign 0.119 (WP5-v4 6.4)", si.fano_factor, 0.119, 0.05);
    audit::check("Si Debye (THz) vs 2408.04732 T5 15.28", si.debye_freq / units::THz, 15.28, 0.03);
    audit::check("Ge Debye (THz) vs campaign 8.64 (WP8-R T1.1)", ge.debye_freq / units::THz, 8.64, 0.01);
    audit::check("GaAs Debye (THz) vs 2408.04732 T5 8.58", ga.debye_freq / units::THz, 8.58, 0.05);
    audit::check("GaAs A_anh vs 2408.04732 T3 7.77e-55", ga.A_anh, 7.77e-55, 0.02);
    audit::check("GaAs F_TT vs 2408.04732 T3 0.77", ga.fTT_frac, 0.77, 0.02);
    return audit::result();
}
