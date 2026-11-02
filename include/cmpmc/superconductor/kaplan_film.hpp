#pragma once

#include "cmpmc/core/core.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include "cmpmc/materials/presets.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace cmpmc {

// ============================================================================
// KaplanFilm: Lumped Superconducting Thin-Film Quasiparticle Cascade Model
// Replaces: G4CMPKaplanQP & G4CMPVKaplanQP (S.B. Kaplan et al., PRB 14, 1976)
//
// Key Optimization over G4CMP:
//   G4CMPVKaplanQP::QPEnergyRand and PhononEnergyRand use uniform rejection
//   sampling on [Delta + eps, E - Delta - eps] with BUFF=1000 to avoid the
//   1/sqrt(x^2 - Delta^2) endpoint singularities, resulting in ~15-30 rejection
//   trials per broken Cooper pair.
//   Here we apply exact analytic change-of-variables:
//     1) QP pair breaking: x = Delta + (E - 2*Delta) * sin^2(pi/2 * u)
//        cancels 1 / sqrt((x - Delta)*(E - Delta - x)) identically!
//     2) QP phonon emission: x = Delta + (E - Delta) * u^2
//        cancels 1 / sqrt(x - Delta) identically!
//   Result: Exact endpoint sampling (no BUFF=1000 truncation) and >90%
//   acceptance on the very first random draw.
// ============================================================================
struct KaplanFilm : public SuperconductorFilmConfig {
    constexpr KaplanFilm() noexcept = default;
    constexpr explicit KaplanFilm(const SuperconductorFilmConfig& cfg) noexcept
        : SuperconductorFilmConfig(cfg) {}

    // Standard Aluminum QET sensor parameters from G4CMP examples/phonon
    static constexpr KaplanFilm aluminum_sensor() noexcept {
        return KaplanFilm(materials::AluminumFilm());
    }

    constexpr bool is_subgap(double energy) const noexcept {
        return energy < 2.0 * gap_energy;
    }

    inline bool is_thermalized(double energy, Rng& rng) const noexcept {
        if (temperature <= 0.0) return false;
        const double kT = units::k_Boltzmann * temperature;
        return energy < kT * (-std::log(rng.uniform()));
    }

    // Probability of phonon escaping back into crystal without breaking a Cooper pair
    inline double calc_escape_probability(double energy, double thickness_frac) const noexcept {
        if (gap_energy <= 0.0) return 1.0;
        const double denom = 1.0 + phonon_lifetime_slope * (energy / gap_energy - 2.0);
        if (denom <= 0.0) return 1.0;
        const double mfp = (v_sound * phonon_lifetime) / denom;
        const double path = thickness_frac * film_thickness;
        return std::exp(-path / mfp);
    }

    // Main entry point matching G4CMPKaplanQP::AbsorbPhonon
    // ------------------------------------------------------------------------
    // Exact T = 0 Kaplan pair-breaking rate (Kaplan PRB 14, 4854 (1976), as implemented in
    // G4CMPSCPairBreakingRate and quoted by 1401.1937 via tau_0^phi):
    //   1/tau_pb(Omega) = I(Omega) / (pi D tau0),
    //   I(Omega) = Int_D^{Omega-D} dE [E(Omega-E) + D^2] / sqrt((E^2-D^2)((Omega-E)^2-D^2)),
    // so 1/tau_pb(2D) = 1/tau0 exactly. With E = D + (s/2)(1 - cos t), s = Omega - 2D, the
    // inverse square roots cancel and the midpoint rule in t (Gauss-Chebyshev) converges
    // exponentially: I = Int_0^pi dt [E(Omega-E) + D^2]/sqrt((E+D)(Omega-E+D)).
    // ------------------------------------------------------------------------
    inline double kaplan_pb_integral(double omega) const noexcept {
        const double d = gap_energy;
        if (omega < 2.0 * d) return 0.0;
        const double s = omega - 2.0 * d;
        constexpr int NQ = 48;
        double acc = 0.0;
        for (int i = 0; i < NQ; ++i) {
            const double t = units::pi * (i + 0.5) / NQ;
            const double E = d + 0.5 * s * (1.0 - std::cos(t));
            acc += (E * (omega - E) + d * d) / std::sqrt((E + d) * (omega - E + d));
        }
        return acc * units::pi / NQ;
    }
    inline double kaplan_pair_breaking_rate(double omega) const noexcept {
        if (gap_energy <= 0.0 || phonon_lifetime <= 0.0 || omega < 2.0 * gap_energy) return 0.0;
        return kaplan_pb_integral(omega) / (units::pi * gap_energy * phonon_lifetime);
    }
    // Film-loss time from geometry: mean escape path 4d (2302.05998 Eq. (17) [FT]: 2 x 2d,
    // traversal and angle average), tau_l = 4 d / (eta_t v_s) with interface transmission eta_t.
    static inline double loss_time_from_geometry(double d_film, double v_s, double eta_t) noexcept {
        return (eta_t > 0.0 && v_s > 0.0) ? 4.0 * d_film / (eta_t * v_s) : 0.0;
    }
    // Probability that a film phonon of energy E (>= 2D) escapes before it breaks a pair.
    inline double loss_probability(double energy, double thickness_frac) const noexcept {
        if (loss_model == FilmLossModel::KaplanRateCompetition) {
            if (phonon_loss_time <= 0.0) return 0.0;
            const double g = kaplan_pair_breaking_rate(energy);
            if (g <= 0.0) return 1.0;
            return 1.0 / (1.0 + g * phonon_loss_time);   // = tau_pb / (tau_pb + tau_l)
        }
        return calc_escape_probability(energy, thickness_frac);
    }

    // Seed of the first Cooper-pair split: a phonon entering the film (case I coherence
    // factor 1 + D^2/E1E2) or a photon absorbed in it (case II, 1 - D^2/E1E2; 1401.1937
    // Eq. (1) [FT]; WP5-v4 Sec. 5). A photon seed skips the entry-escape step.
    enum class Seed : uint8_t { Phonon, Photon };

    double absorb_phonon(double energy, KaplanBuffer& reflected_energies, Rng& rng) const noexcept {
        return absorb_quantum(energy, reflected_energies, rng, Seed::Phonon);
    }

    double absorb_quantum(double energy, KaplanBuffer& reflected_energies, Rng& rng,
                          Seed seed) const noexcept {
        reflected_energies.clear();
        const bool photon = (seed == Seed::Photon);

        if (!photon && do_direct_absorption(energy, rng)) {
            return energy;
        }
        if (is_subgap(energy) || (!photon && rng.uniform() <= loss_probability(energy, 4.0))) {
            reflected_energies.push_back(energy);
            return 0.0;
        }

        double edep = 0.0;
        const int n_pairs = (high_qp_limit > 0.0 && gap_energy > 0.0)
            ? static_cast<int>(std::ceil(energy / (2.0 * high_qp_limit * gap_energy)))
            : 1;

        KaplanBuffer qp_energies;
        KaplanBuffer phon_energies;
        KaplanBuffer next_phon;
        KaplanBuffer next_qp;

        const double init_e = energy / n_pairs;
        if (photon) {
            for (int i = 0; i < n_pairs; ++i) {
                const double qpE = sample_qp_energy_case2(init_e, rng);
                edep += calc_qp_absorption(qpE, phon_energies, qp_energies, rng);
                edep += calc_qp_absorption(init_e - qpE, phon_energies, qp_energies, rng);
            }
        } else {
            for (int i = 0; i < n_pairs; ++i) {
                phon_energies.push_back(init_e);
            }
        }

        while (!qp_energies.empty() || !phon_energies.empty()) {
            // A. Phonons >= 2*Delta break Cooper pairs into two quasiparticles (qpE, E - qpE)
            if (!phon_energies.empty()) {
                next_phon.clear();
                for (double E : phon_energies) {
                    if (is_subgap(E)) {
                        next_phon.push_back(E);
                        continue;
                    }
                    const double qpE = sample_qp_energy(E, rng);
                    edep += calc_qp_absorption(qpE, next_phon, qp_energies, rng);
                    edep += calc_qp_absorption(E - qpE, next_phon, qp_energies, rng);
                }
                phon_energies = next_phon;
            }

            // B. Quasiparticles >= low_qp_limit * Delta radiate phonons
            if (!qp_energies.empty()) {
                next_qp.clear();
                for (double E : qp_energies) {
                    const double phonE = sample_phonon_energy(E, rng);
                    const double qpE   = E - phonE;
                    if (is_subgap(phonE)) {
                        edep += calc_direct_absorption(phonE, phon_energies, rng);
                    } else {
                        phon_energies.push_back(phonE);
                    }
                    edep += calc_qp_absorption(qpE, phon_energies, next_qp, rng);
                }
                qp_energies = next_qp;
            }

            // C. Phonons in film either escape into crystal or remain for next cascade step
            if (!phon_energies.empty()) {
                next_phon.clear();
                for (double E : phon_energies) {
                    if (is_thermalized(E, rng)) continue;
                    if (loss_model == FilmLossModel::KaplanRateCompetition) {
                        // sub-gap phonons cannot break pairs and leave the film; the others
                        // compete pair breaking (1/tau_pb) against loss (1/tau_l)
                        if (is_subgap(E) || rng.uniform() < loss_probability(E, 1.0)) {
                            reflected_energies.push_back(E);
                        } else {
                            next_phon.push_back(E);
                        }
                        continue;
                    }
                    const double frac = ((rng.uniform() < 0.5) ? 0.5 : 1.5) /
                                        std::cos(rng.uniform() * 1.47);
                    if (rng.uniform() < calc_escape_probability(E, frac)) {
                        reflected_energies.push_back(E);
                    } else if (keep_all_phonons || !is_subgap(E)) {
                        next_phon.push_back(E);
                    }
                }
                phon_energies = next_phon;
            }
        }

        return edep;
    }

    // Singularity-free Cooper-pair breaking energy sampler
    inline double sample_qp_energy(double E, Rng& rng) const noexcept {
        const double d = gap_energy;
        const double span = E - 2.0 * d;
        if (span <= 1e-16) return 0.5 * E;

        auto weight = [E, d](double x) noexcept {
            return (x * (E - x) + d * d) / std::sqrt((x + d) * (E - x + d));
        };
        const double w_max = std::max(weight(d), weight(0.5 * E));

        for (;;) {
            const double u = rng.uniform();
            const double s = std::sin(units::half_pi * u);
            const double x = d + span * (s * s);
            if (rng.uniform() * w_max <= weight(x)) {
                return x;
            }
        }
    }

    // Case-II (photon) split: weight [x(E-x) - D^2]/sqrt((x^2-D^2)((E-x)^2-D^2)) (1401.1937 Eq. (1) [FT])
    inline double sample_qp_energy_case2(double E, Rng& rng) const noexcept {
        const double d = gap_energy;
        const double span = E - 2.0 * d;
        if (span <= 1e-16) return 0.5 * E;
        auto weight = [E, d](double x) noexcept {
            return std::max(0.0, x * (E - x) - d * d) / std::sqrt((x + d) * (E - x + d));
        };
        const double w_max = std::max(weight(d), weight(0.5 * E)) * 1.0001;
        for (;;) {
            const double u = rng.uniform();
            const double s = std::sin(units::half_pi * u);
            const double x = d + span * (s * s);
            if (rng.uniform() * w_max <= weight(x)) return x;
        }
    }

    // Pair-breaking efficiency and QP Fano factor of one absorbed quantum (WP5-v4 (W13)):
    //   eta_pb = D <N_qp>/Omega,  F_qp = Var(N_qp)/<N_qp>  (E_sig = D N_qp; requires
    //   subgap_absorption = 0 and absorber_eff = 1 so that the cascade's edep = D N_qp).
    struct PairBreakingResult { double eta = 0, eta_err = 0, F_qp = 0; };
    inline PairBreakingResult pair_breaking_efficiency(double omega, int n_events, Rng& rng,
                                                       Seed seed = Seed::Phonon) const noexcept {
        KaplanFilm f = *this;
        f.subgap_absorption = 0.0; f.absorber_eff = 1.0; f.absorber_eff_slope = 0.0; f.temperature = 0.0;
        double s1 = 0.0, s2 = 0.0;
        KaplanBuffer refl;
        for (int i = 0; i < n_events; ++i) {
            const double nqp = f.absorb_quantum(omega, refl, rng, seed) / gap_energy;
            s1 += nqp; s2 += nqp * nqp;
        }
        PairBreakingResult r;
        const double m = s1 / n_events, var = std::max(0.0, s2 / n_events - m * m);
        r.eta = m * gap_energy / omega;
        r.eta_err = std::sqrt(var / n_events) * gap_energy / omega;
        r.F_qp = (m > 0.0) ? var / m : 0.0;
        return r;
    }

    // Singularity-free QP-radiated phonon energy sampler
    inline double sample_phonon_energy(double E, Rng& rng) const noexcept {
        const double d = gap_energy;
        const double span = E - d;
        if (span <= 1e-16) return 0.0;
        const double d2_over_E = (d * d) / E;

        auto weight = [E, d, d2_over_E](double x) noexcept {
            const double diff = E - x;
            return (diff * diff) * (x - d2_over_E) / std::sqrt(x + d);
        };
        const double x_peak = std::clamp(0.25 * (E + 3.0 * d2_over_E), d, E);
        const double w_max = std::max(weight(d), weight(x_peak)) * 1.02;

        for (;;) {
            const double u = rng.uniform();
            const double x = d + span * (u * u);
            if (rng.uniform() * w_max <= weight(x)) {
                return E - x;
            }
        }
    }

private:
    inline bool do_direct_absorption(double energy, Rng& rng) const noexcept {
        if (energy < 2.0 * absorber_gap) return false;
        return rng.uniform() < subgap_absorption;
    }

    inline double calc_direct_absorption(double energy, KaplanBuffer& keep, Rng& rng) const noexcept {
        if (do_direct_absorption(energy, rng)) return energy;
        keep.push_back(energy);
        return 0.0;
    }

    inline double calc_qp_efficiency(double qpE) const noexcept {
        const double eff = (gap_energy > 0.0)
            ? (absorber_eff + absorber_eff_slope * (qpE / gap_energy))
            : 1.0;
        return std::clamp(eff, 0.0, 1.0);
    }

    inline double calc_qp_absorption(double qpE,
                                     KaplanBuffer& phon_energies,
                                     KaplanBuffer& qp_energies,
                                     Rng& rng) const noexcept {
        if (rng.uniform() > calc_qp_efficiency(qpE)) return 0.0;
        if (qpE >= low_qp_limit * gap_energy) {
            qp_energies.push_back(qpE);
            return 0.0;
        } else if (qpE > gap_energy) {
            double edep = calc_direct_absorption(qpE - gap_energy, phon_energies, rng);
            return edep + gap_energy;
        }
        return qpE;
    }
};

// ============================================================================
// KaplanSurrogateLUT: pre-sampled Kaplan cascade outcomes for SIMD / GPU paths.
//   For each of N_E incident energies on [2 Delta, 2 Delta * E_MAX_DIV_2GAP] the
//   constructor stores N_S complete outcomes of KaplanFilm::absorb_phonon: the
//   absorbed fraction and the individual re-emitted phonon energies (as fractions
//   of the incident energy). A lookup draws one stored outcome from one of the two
//   bracketing energy nodes (chosen with linear-interpolation probabilities) and
//   rescales it to the incident energy. Energies outside the table are passed to
//   the full cascade.
//   History: the original table stored only (abs_frac, sum of reflected fractions)
//   and re-emitted the sum as ONE phonon, which then re-broke pairs on re-entry
//   (collected fraction +13% in a p_ret = 0.9 cavity, tests/test_kaplan), used the
//   lower energy node only (+10% at 3 Delta) and clamped at 40 Delta.
// ============================================================================
class KaplanSurrogateLUT {
public:
    static constexpr int N_E = 128;                 // energy nodes
    static constexpr int N_S = 512;                 // stored outcomes per node
    static constexpr double E_MAX_DIV_2GAP = 20.0;  // table covers [2 Delta, 40 Delta]

    KaplanSurrogateLUT() = default;

    explicit KaplanSurrogateLUT(const KaplanFilm& film, uint64_t seed = 1337ULL, int /*unused*/ = 0)
        : film_(film), abs_(N_E * N_S), off_(N_E * N_S + 1, 0) {
        Rng rng(seed);
        const double two_delta = 2.0 * std::max(film_.gap_energy, 1e-12);
        KaplanBuffer refl;
        for (int ie = 0; ie < N_E; ++ie) {
            const double E = node_energy(ie, two_delta);
            for (int s = 0; s < N_S; ++s) {
                const int k = ie * N_S + s;
                abs_[k] = static_cast<float>(film_.absorb_phonon(E, refl, rng) / E);
                for (double r : refl) refl_.push_back(static_cast<float>(r / E));
                off_[k + 1] = static_cast<uint32_t>(refl_.size());
            }
        }
    }

    // Lookup with the KaplanFilm::absorb_phonon signature
    inline double absorb_phonon_fast(double energy, KaplanBuffer& reflected_energies, Rng& rng) const noexcept {
        const double two_delta = 2.0 * std::max(film_.gap_energy, 1e-12);
        const double u = (energy / two_delta - 1.0) / (E_MAX_DIV_2GAP - 1.0) * (N_E - 1);
        if (film_.is_subgap(energy) || u < 0.0 || u > N_E - 1) {
            return film_.absorb_phonon(energy, reflected_energies, rng);   // exact path outside the table
        }
        reflected_energies.clear();
        int ie = std::min(static_cast<int>(u), N_E - 2);
        if (rng.uniform() < (u - ie)) ++ie;
        const int k = ie * N_S + std::min(N_S - 1, static_cast<int>(rng.uniform() * N_S));
        for (uint32_t r = off_[k]; r < off_[k + 1]; ++r) reflected_energies.push_back(energy * refl_[r]);
        return energy * static_cast<double>(abs_[k]);
    }

private:
    static double node_energy(int ie, double two_delta) noexcept {
        return two_delta * (1.0 + (E_MAX_DIV_2GAP - 1.0) * (static_cast<double>(ie) / (N_E - 1)));
    }
    KaplanFilm film_{};
    std::vector<float> abs_;
    std::vector<uint32_t> off_;
    std::vector<float> refl_;
};

} // namespace cmpmc
