// G4CMP config.txt import, and agreement of the built-in presets with the
// corresponding G4CMP CrystalMaps entries.
//
// Requires a CrystalMaps directory, passed through the CMake cache variable
// CMPMC_G4CMP_CRYSTALMAPS_DIR (or the environment variable of the same name).
// Without it the test reports SKIPPED.
#include "test_support.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include "cmpmc/io/g4cmp_config_import.hpp"
#include <cstdlib>
#include <filesystem>
#include <string>

using namespace cmpmc;

static void compare_with_imported_config(const CrystalConfig& p, const std::string& dir) {
    const std::string f = dir + "/" + p.name + "/config.txt";
    if (!std::filesystem::exists(f)) { std::printf("# no %s, not compared\n", f.c_str()); return; }
    const CrystalConfig g = io::import_g4cmp_config(f, p.density, p.name);
    const std::string n = p.name + " preset vs imported config: ";
    const auto Vp = p.build_voigt_matrix(), Vg = g.build_voigt_matrix();
    double dv = 0, vmax = 0;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { dv = std::max(dv, std::abs(Vp[i][j] - Vg[i][j])); vmax = std::max(vmax, std::abs(Vg[i][j])); }
    cmpmc_test::check(n + "max |C_ij - C_ij(imported)| / max C", dv / vmax, 0.0, 0.0, 1e-6);
    cmpmc_test::check(n + "B_iso", p.B_iso, g.B_iso, 1e-6, 1e-60);
    cmpmc_test::check(n + "A_anh", p.A_anh, g.A_anh, 1e-6);
    cmpmc_test::check(n + "fTT", p.fTT_frac, g.fTT_frac, 1e-6);
    cmpmc_test::check(n + "LDOS", p.fLDOS, g.fLDOS, 1e-6);
    cmpmc_test::info(n + "Debye (preset; imported value next line)", p.debye_freq);
    cmpmc_test::info(n + "Debye (imported; presets use the max optical frequency)", g.debye_freq);
    cmpmc_test::check(n + "v_L", p.v_sound, g.v_sound, 1e-6);
    cmpmc_test::check(n + "v_T", p.v_trans, g.v_trans, 1e-6);
}

int main() {
    std::printf("# test_g4cmp_import\n");
    const char* env = std::getenv("CMPMC_G4CMP_CRYSTALMAPS_DIR");
    if (env == nullptr || *env == '\0' || !std::filesystem::is_directory(env))
        return cmpmc_test::skip("CMPMC_G4CMP_CRYSTALMAPS_DIR is not set to a CrystalMaps directory");
    const std::string dir = env;
    for (const auto& c : {materials::Germanium(), materials::Silicon(), materials::Sapphire(),
                          materials::CalciumTungstate(), materials::GalliumArsenide(),
                          materials::LithiumFluoride(), materials::CalciumFluoride(),
                          materials::AluminumCrystal(), materials::NiobiumCrystal()})
        compare_with_imported_config(c, dir);
    return cmpmc_test::result();
}
