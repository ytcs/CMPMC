// Kaplan thin-film QP cascade (KaplanFilm), its surrogate table, and the
// 3D Bogoliubov-QP rates (BulkSuperconductor).
//  (a) pair-breaking QP-energy sampler and QP phonon-emission sampler vs the Kaplan
//      densities [E1 E2 + D^2]/sqrt((E1^2 - D^2)(E2^2 - D^2)) and
//      (E - E')^2 (E' - D^2/E)/sqrt(E'^2 - D^2) (Kaplan et al.,
//      Phys. Rev. B 14, 4854 (1976)), by moments;
//  (b) pair-breaking efficiency eta_pb = D <N_qp>/Omega vs Omega/D in the no-loss
//      limit (thick film, no sub-gap absorption): exact value 2D/Omega for Omega < 4D
//      (only one pair can be broken); structure above 4D and the asymptote (0.57-0.6,
//      arXiv:1401.1937) are reported;
//  (c) the film escape model's linear lifetime law tau(Omega) = tau0/(1 + 0.29 (Omega/D - 2))
//      vs the exact Kaplan pair-breaking rate (BulkSuperconductor LUT);
//  (d) KaplanSurrogateLUT vs the full cascade in a re-entry cavity (reflected phonons
//      return to the film with probability p_ret): collected-energy fraction.
#include "test_support.hpp"
#include "cmpmc/superconductor/kaplan_film.hpp"
#include "cmpmc/superconductor/bulk_superconductor.hpp"
#include <vector>

using namespace cmpmc;

template <class F> static double mean_pdf(F&& p, double lo, double hi) {
    // integrate in u with x = lo + (hi - lo) sin^2(pi u / 2) to tame endpoint singularities
    const int N = 400000; double s0 = 0, s1 = 0;
    for (int i = 0; i < N; ++i) { const double u = (i + 0.5) / N, sn = std::sin(units::half_pi * u);
        const double x = lo + (hi - lo) * sn * sn, jac = (hi - lo) * units::pi * sn * std::cos(units::half_pi * u);
        const double w = p(x) * jac; s0 += w; s1 += w * x; }
    return s1 / s0;
}

static double cascade_with_reentry(const KaplanFilm& film, const KaplanSurrogateLUT* lut, double E, double p_ret, Rng& rng) {
    std::vector<double> todo{E}; double edep = 0; KaplanBuffer refl;
    while (!todo.empty()) {
        const double e = todo.back(); todo.pop_back();
        edep += lut ? lut->absorb_phonon_fast(e, refl, rng) : film.absorb_phonon(e, refl, rng);
        for (double r : refl) if (rng.uniform() < p_ret) todo.push_back(r);
    }
    return edep / E;
}

int main() {
    std::printf("# test_kaplan\n");
    KaplanFilm film = KaplanFilm::aluminum_sensor();
    const double D = film.gap_energy;
    Rng rng(31);
    // (a)
    for (double W : {2.5, 5.0, 20.0}) {
        const double Om = W * D; double s = 0; const int NS = 400000;
        for (int i = 0; i < NS; ++i) s += film.sample_qp_energy(Om, rng);
        const double ref = mean_pdf([&](double x) { return (x * (Om - x) + D * D) / std::sqrt((x * x - D * D) * ((Om - x) * (Om - x) - D * D)); }, D, Om - D);
        cmpmc_test::check("pair-breaking <E1>/D, Omega=" + std::to_string(W).substr(0, 4) + "D", s / NS / D, ref / D, 2e-3);
        double s2 = 0; for (int i = 0; i < NS; ++i) s2 += film.sample_phonon_energy(Om, rng);
        const double ref2 = mean_pdf([&](double x) { return (Om - x) * (Om - x) * (x - D * D / Om) / std::sqrt(x * x - D * D); }, D, Om);
        cmpmc_test::check("QP emission <Omega_ph>/D, E_qp=" + std::to_string(W).substr(0, 4) + "D", s2 / NS / D, (Om - ref2) / D, 2e-3);
    }
    // (b) eta_pb in the no-loss limit
    KaplanFilm thick = film;
    thick.film_thickness = 1.0; thick.subgap_absorption = 0.0; thick.keep_all_phonons = false;
    thick.temperature = 0.0; thick.absorber_eff = 1.0;
    for (double W : {2.0001, 2.5, 3.0, 3.9, 4.5, 5.0, 6.0, 8.0, 10.0, 20.0, 40.0, 100.0}) {
        double s = 0; const int NS = (W > 50) ? 20000 : 100000; KaplanBuffer refl;
        for (int i = 0; i < NS; ++i) s += thick.absorb_phonon(W * D, refl, rng);
        const double eta = s / NS / (W * D);
        const std::string tag = "eta_pb no-loss, Omega=" + std::to_string(W).substr(0, 5) + "D";
        if (W < 4.0) cmpmc_test::check(tag + " vs 2D/Omega (W13 limit)", eta, 2.0 / W, 2e-3);
        else cmpmc_test::info(tag, eta);
    }
    // with the default 600 nm Al film (escape on): eta including loss
    {
        KaplanFilm f = film; f.subgap_absorption = 0.0; f.keep_all_phonons = false;
        for (double W : {3.0, 5.0, 10.0, 20.0, 40.0}) {
            double s = 0; const int NS = 100000; KaplanBuffer refl;
            for (int i = 0; i < NS; ++i) s += f.absorb_phonon(W * D, refl, rng);
            cmpmc_test::info("eta_pb 600 nm Al film (escape on), Omega=" + std::to_string(W).substr(0, 4) + "D", s / NS / (W * D));
        }
    }
    // (c) linear lifetime law vs exact Kaplan pair-breaking rate
    {
        CrystalConfig al = materials::AluminumCrystal();
        const physics::BulkSuperconductor sc(al);
        double worst = 0;
        for (double W : {2.5, 4.0, 10.0, 20.0, 40.0}) {
            const double exact = sc.phonon_pair_breaking_rate(W * sc.delta_T_eV) * al.sc_tau0_ph_s;
            const double lin = 1.0 + film.phonon_lifetime_slope * (W - 2.0);
            cmpmc_test::info("tau0/tau_pb exact Kaplan, Omega=" + std::to_string(W).substr(0, 4) + "D", exact);
            cmpmc_test::info("tau0/tau_pb film linear law, same Omega", lin);
            worst = std::max(worst, std::abs(lin / exact - 1.0));
        }
        cmpmc_test::info("max |linear law / exact Kaplan - 1| on [2.5, 40] D", worst);
    }
    // (d) surrogate vs full cascade with re-entry
    {
        KaplanFilm f = film;
        const KaplanSurrogateLUT lut(f);
        for (double p_ret : {0.0, 0.9}) {
            for (double W : {3.0, 6.0, 15.0, 35.0}) {
                double a = 0, b = 0; const int NS = 50000;
                for (int i = 0; i < NS; ++i) { a += cascade_with_reentry(f, nullptr, W * D, p_ret, rng);
                                               b += cascade_with_reentry(f, &lut, W * D, p_ret, rng); }
                cmpmc_test::check("surrogate/full collected fraction, p_ret=" + std::to_string(p_ret).substr(0, 3) +
                             " Omega=" + std::to_string(W).substr(0, 4) + "D", b / NS, a / NS, 0.03);
            }
        }
    }
    return cmpmc_test::result();
}
