// Audit test: bulk phonon rate laws (source form, nu in Hz):
//   Gamma_iso = B nu^4 for every mode (Tamura; 2408.04732 Eq. (3)), Gamma_anh = A nu^5 for L
//   only (ICA; 2408.04732 Eq. (1); G4CMP 2302.04315), campaign WP8-R Sec. 2.1.
// Preset A values are reported against the literature spread (2408.04732 Table 3 [FT]).
#include "test_support.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"

using namespace cmpmc;

int main() {
    std::printf("# test_rates\n");
    for (const auto& c : {materials::Silicon(), materials::Germanium()}) {
        const CrystalLattice lat(c);
        for (double nu : {0.3e12, 1e12, 3e12, 10e12}) {
            const double E = nu * units::h_Planck;
            const std::string t = c.name + " nu=" + std::to_string(nu / 1e12).substr(0, 4) + " THz ";
            cmpmc_test::check(t + "Gamma_iso = B nu^4", lat.isotope_rate(E), c.B_iso * std::pow(nu, 4), 1e-12);
            cmpmc_test::check(t + "Gamma_anh(L) = A nu^5", lat.downconversion_rate(PhononMode::L, E), c.A_anh * std::pow(nu, 5), 1e-12);
            cmpmc_test::check(t + "Gamma_anh(FT) = 0", lat.downconversion_rate(PhononMode::FT, E), 0.0, 0.0, 0.0);
        }
    }
    cmpmc_test::info("Si A preset (1e-55 s^4); 2408.04732 T3: calc 1.15, Tamura 0.741", materials::Silicon().A_anh / 1e-55);
    cmpmc_test::info("Ge A preset (1e-55 s^4); 2408.04732 T3: calc 6.8, Tamura 16.5; 1109.4343: 1.61", materials::Germanium().A_anh / 1e-55);
    return cmpmc_test::result();
}
