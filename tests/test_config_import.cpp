// io::import_g4cmp_config on a synthetic file: keywords, units, comments, unknown keys
// and the defaults of fields the file does not set.
#include "test_support.hpp"
#include "cmpmc/io/g4cmp_config_import.hpp"
#include <filesystem>
#include <fstream>
#include <stdexcept>

using namespace cmpmc;

int main() {
    std::printf("# test_config_import\n");
    const auto path = std::filesystem::temp_directory_path() / "cmpmc_test_config_import.txt";
    {
        std::ofstream out(path);
        out << "# synthetic lattice file\n"
               "tetragonal\n"
               "stiffness 1 1 150.0 GPa   # C11\n"
               "stiffness 3 3 120.0 GPa\n"
               "stiffness 6 6 40.0 GPa\n"
               "scat 1.5e-42\n"
               "decay 2.5e-55\n"
               "decayTT 0.6\n"
               "LDOS 0.1\nSTDOS 0.5\nFTDOS 0.4\n"
               "Debye 5.0 THz\n"
               "vsound 6.0 km/s\n"
               "vtrans 3500 m/s\n"
               "bandgap 1.1 eV\n"
               "pairEnergy 3.8 eV\n"
               "fanoFactor 0.12\n"
               "emass 0.9 0.2 0.2\n"
               "valleyDir 1 0 0\n"
               "valleyDir 0 1 0\n"
               "someFutureKey 42\n";
    }
    const CrystalConfig c = io::import_g4cmp_config(path.string(), 4000.0, "Synthetic");
    std::filesystem::remove(path);

    const auto V = c.build_voigt_matrix();
    cmpmc_test::check("crystal system is tetragonal", c.group == BravaisGroup::Tetragonal ? 1.0 : 0.0, 1.0, 0.0);
    cmpmc_test::check("name", c.name == "Synthetic" ? 1.0 : 0.0, 1.0, 0.0);
    cmpmc_test::check("density (argument)", c.density, 4000.0, 0.0);
    cmpmc_test::check("C11 (GPa unit)", V[0][0], 150.0 * units::GPa, 1e-12);
    cmpmc_test::check("C33", V[2][2], 120.0 * units::GPa, 1e-12);
    cmpmc_test::check("C66", V[5][5], 40.0 * units::GPa, 1e-12);
    cmpmc_test::check("B_iso", c.B_iso, 1.5e-42, 1e-12);
    cmpmc_test::check("A_anh", c.A_anh, 2.5e-55, 1e-12);
    cmpmc_test::check("fTT", c.fTT_frac, 0.6, 1e-12);
    cmpmc_test::check("LDOS + STDOS + FTDOS", c.fLDOS + c.fSTDOS + c.fFTDOS, 1.0, 1e-12);
    cmpmc_test::check("Debye (THz unit)", c.debye_freq, 5.0 * units::THz, 1e-12);
    cmpmc_test::check("v_L (km/s unit)", c.v_sound, 6000.0, 1e-12);
    cmpmc_test::check("v_T (m/s unit)", c.v_trans, 3500.0, 1e-12);
    cmpmc_test::check("band gap (eV)", c.bandgap, 1.1, 1e-12);
    cmpmc_test::check("pair energy (eV)", c.pair_energy, 3.8, 1e-12);
    cmpmc_test::check("Fano factor", c.fano_factor, 0.12, 1e-12);
    cmpmc_test::check("electron mass tensor m_L", c.e_mass_tensor.x, 0.9, 1e-12);
    cmpmc_test::check("valley count", static_cast<double>(c.valley_dirs.size()), 2.0, 0.0);

    const CrystalConfig defaults{};
    cmpmc_test::check("unset field keeps default (hole mass)", c.h_mass_ratio, defaults.h_mass_ratio, 0.0);

    bool threw = false;
    try { (void)io::import_g4cmp_config("/nonexistent/cmpmc/config.txt"); } catch (const std::runtime_error&) { threw = true; }
    cmpmc_test::check("missing file throws std::runtime_error", threw ? 1.0 : 0.0, 1.0, 0.0);
    return cmpmc_test::result();
}
