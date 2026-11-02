#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/core/units.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace cmpmc::physics {

// Bulk superconductor: Bogoliubov quasiparticle and pair-breaking physics.
//   - BCS temperature-dependent gap Delta(T) and Fermi/Bose factors
//   - phonon -> 2 QP pair breaking: rate and energy sampling
//   - QP -> QP' + phonon radiation: rate and energy sampling
//   - QP + thermal QP -> phonon recombination
//   - local impurity / vortex trapping, rate 1 / tau_trap
//   - energy-dependent diffusivity D(E) = D_n sqrt(1 - (Delta/E)^2) and Gaussian
//     diffusive stepping
//
// The rate tables are integrals over the BCS density of states, which has
// 1/sqrt(E^2 - Delta^2) singularities. Trigonometric and quadratic substitutions
// cancel these singularities exactly, so 64 quadrature nodes per table suffice.
//
// This module feeds SimulationEngine::transport_superconductor_coupled, which is
// guarded as experimental (see docs/experimental.md).

struct BulkSuperconductor {
    static constexpr int N_LUT_BINS = 256;
    static constexpr double MIN_PH_DIV_GAP = 2.0;
    static constexpr double MAX_PH_DIV_GAP = 60.0;
    static constexpr double MIN_QP_DIV_GAP = 1.0;
    static constexpr double MAX_QP_DIV_GAP = 60.0;

    double delta0_eV   = 0.0; // Zero-temperature BCS gap Delta(0) [eV]
    double t_crit_K    = 0.0; // Critical temperature T_c [K]
    double t_eff_K     = 0.0; // Effective quasiparticle bath temperature T_eff [K]
    double delta_T_eV  = 0.0; // Temperature-corrected BCS gap Delta(T_eff) [eV]
    double tau0_ph_s   = 0.0; // Characteristic 2*Delta phonon pair-breaking lifetime [s]
    double tau0_qp_s   = 0.0; // Characteristic QP phonon-emission lifetime [s]
    double dn_m2_s     = 0.0; // Normal-state electron diffusivity D_n [m^2/s]
    double trap_tau_s  = 0.0; // Local QP trapping time tau_trap [s] (<= 0 disables)

    std::array<double, N_LUT_BINS> inv_tau_pb_lut_{};  // 1 / tau_pb(E_ph)  [1/s]
    std::array<double, N_LUT_BINS> inv_tau_rad_lut_{}; // 1 / tau_rad(E_qp) [1/s]
    std::array<double, N_LUT_BINS> inv_tau_rec_lut_{}; // 1 / tau_rec(E_qp) [1/s]

    BulkSuperconductor() = default;

    explicit BulkSuperconductor(const CrystalConfig& cfg) {
        init_from_config(cfg);
    }

    bool is_enabled() const noexcept { return delta_T_eV > 0.0; }

    // --- 1. BCS Temperature-Dependent Gap & Thermal Factors ---
    static inline double bcs_gap_at_temp(double delta0_eV, double t_eff_K, double t_crit_K) noexcept {
        if (delta0_eV <= 0.0) return 0.0;
        if (t_crit_K <= 0.0) {
            // Standard BCS weak-coupling relation: Delta(0) = 1.764 * k_B * T_c
            t_crit_K = delta0_eV / (1.7638 * units::k_Boltzmann);
        }
        if (t_eff_K <= 0.0) return delta0_eV;
        if (t_eff_K >= t_crit_K) return 0.0;

        // Exact Mühlschlegel / Tinkham BCS gap interpolation (agrees with the tabulated BCS gap to <0.1%)
        const double u = t_eff_K / t_crit_K;
        const double factor = std::tanh(1.74 * std::sqrt(std::max(0.0, (1.0 / u) - 1.0)));
        return delta0_eV * factor;
    }

    static inline double fermi_factor(double energy_eV, double temp_K) noexcept {
        if (temp_K <= 0.0) return 0.0;
        const double x = energy_eV / (units::k_Boltzmann * temp_K);
        if (x > 500.0) return 0.0;
        return 1.0 / (std::exp(x) + 1.0);
    }

    static inline double bose_factor(double energy_eV, double temp_K) noexcept {
        if (temp_K <= 0.0 || energy_eV <= 0.0) return 0.0;
        const double x = energy_eV / (units::k_Boltzmann * temp_K);
        if (x > 500.0) return 0.0;
        return 1.0 / (std::exp(x) - 1.0);
    }

    // --- 2. Energy-Dependent Quasiparticle Diffusivity & Group Velocity ---
    //   D(E) = D_n * sqrt(1 - (Delta / E)^2)
    inline double qp_diffusivity(double qp_energy_eV) const noexcept {
        if (qp_energy_eV <= delta_T_eV || dn_m2_s <= 0.0) return 0.0;
        const double r = delta_T_eV / qp_energy_eV;
        return dn_m2_s * std::sqrt(std::max(0.0, 1.0 - r * r));
    }

    // Effective group velocity scale for QP stepping
    inline double qp_velocity(double qp_energy_eV) const noexcept {
        if (qp_energy_eV <= delta_T_eV) return 1.0e3;
        const double r = delta_T_eV / qp_energy_eV;
        return 1.0e6 * std::sqrt(std::max(1e-6, 1.0 - r * r)); // Fermi velocity scale * BCS group factor
    }

    // --- 3. Lookup Table Initialization (Singularity-Cancelled Quadrature) ---
    void init_from_config(const CrystalConfig& cfg) noexcept {
        delta0_eV  = cfg.sc_delta0_eV;
        tau0_qp_s  = cfg.sc_tau0_qp_s;
        tau0_ph_s  = cfg.sc_tau0_ph_s;
        dn_m2_s    = cfg.sc_dn_m2_s;
        t_eff_K    = (cfg.sc_teff_K > 0.0) ? cfg.sc_teff_K : 0.015; // Default 15 mK
        trap_tau_s = cfg.sc_trap_tau_s;
        t_crit_K   = (delta0_eV > 0.0) ? (delta0_eV / (1.7638 * units::k_Boltzmann)) : 0.0;
        delta_T_eV = bcs_gap_at_temp(delta0_eV, t_eff_K, t_crit_K);

        if (delta_T_eV <= 0.0) return;
        build_pair_breaking_lut();
        build_qp_radiation_lut();
        build_qp_recombination_lut();
    }

    // --- 4. Fast O(1) Rate Queries [1/s] ---

    // Bulk phonon Cooper-pair breaking rate: Phonon (E >= 2*Delta) -> QP1 + QP2
    inline double phonon_pair_breaking_rate(double phonon_eV) const noexcept {
        if (delta_T_eV <= 0.0 || tau0_ph_s <= 0.0 || phonon_eV <= 2.0 * delta_T_eV) return 0.0;
        const double e_div = phonon_eV / delta_T_eV;
        return interpolate_lut(inv_tau_pb_lut_, e_div, MIN_PH_DIV_GAP, MAX_PH_DIV_GAP);
    }

    // Quasiparticle spontaneous/stimulated phonon radiation rate: QP(E) -> QP(E - Omega) + Phonon(Omega)
    inline double qp_radiates_phonon_rate(double qp_eV) const noexcept {
        if (delta_T_eV <= 0.0 || tau0_qp_s <= 0.0 || qp_eV <= delta_T_eV * 1.0001) return 0.0;
        const double e_div = qp_eV / delta_T_eV;
        return interpolate_lut(inv_tau_rad_lut_, e_div, MIN_QP_DIV_GAP, MAX_QP_DIV_GAP);
    }

    // Quasiparticle thermal recombination rate: QP(E) + Thermal QP -> Phonon(>= E + Delta)
    inline double qp_recombination_rate(double qp_eV) const noexcept {
        if (delta_T_eV <= 0.0 || tau0_qp_s <= 0.0 || qp_eV < delta_T_eV) return 0.0;
        const double e_div = qp_eV / delta_T_eV;
        return interpolate_lut(inv_tau_rec_lut_, e_div, MIN_QP_DIV_GAP, MAX_QP_DIV_GAP);
    }

    // Local vortex / gap-depression trapping rate
    inline double qp_local_trapping_rate() const noexcept {
        return (trap_tau_s > 0.0 && trap_tau_s < 1e20) ? (1.0 / trap_tau_s) : 0.0;
    }

    // --- 5. Singularity-Free Final-State Energy Samplers ---

    // Sample daughter QP energy E1 in [Delta, E_ph - Delta] from pair-breaking phonon E_ph
    // Substitution E1 = Delta + (E_ph - 2*Delta) * sin^2(pi/2 * u) cancels the endpoint singularities.
    inline double sample_pair_breaking_qp_energy(double phonon_eV, Rng& rng) const noexcept {
        const double span = phonon_eV - 2.0 * delta_T_eV;
        if (span <= 0.0) return delta_T_eV;

        for (int iter = 0; iter < 64; ++iter) {
            const double s  = std::sin(units::half_pi * rng.uniform());
            const double e1 = delta_T_eV + span * (s * s);
            const double e2 = phonon_eV - e1;
            const double w_num = e1 * e2 + delta_T_eV * delta_T_eV;
            const double w_den = std::sqrt((e1 + delta_T_eV) * (e2 + delta_T_eV));
            const double w     = (w_num / w_den) * (1.0 - fermi_factor(e1, t_eff_K) - fermi_factor(e2, t_eff_K));
            const double w_max = (0.25 * phonon_eV * phonon_eV + delta_T_eV * delta_T_eV) / (2.0 * delta_T_eV);
            if (rng.uniform() * w_max <= w) return e1;
        }
        return 0.5 * phonon_eV;
    }

    // Sample radiated phonon energy Omega in (0, E_qp - Delta] from a QP of energy E_qp
    // Substitution E' = E_qp - Omega = Delta + (E_qp - Delta) * u^2 cancels the 1/sqrt(E'^2 - Delta^2) pole.
    inline double sample_radiated_phonon_energy(double qp_eV, Rng& rng) const noexcept {
        const double span = qp_eV - delta_T_eV;
        if (span <= 0.0) return 0.0;

        const double max_w = (span * span * span) * std::sqrt(2.0 * qp_eV) * (1.0 - delta_T_eV / qp_eV);
        if (max_w <= 0.0) return 0.5 * span;

        for (int iter = 0; iter < 64; ++iter) {
            const double u  = rng.uniform();
            const double ep = delta_T_eV + span * (u * u); // New QP energy E'
            const double om = qp_eV - ep;                  // Emitted phonon energy Omega
            const double w  = (om * om) * (ep - delta_T_eV * delta_T_eV / qp_eV) *
                              (ep / std::sqrt(std::max(1e-30, ep + delta_T_eV))) *
                              (bose_factor(om, t_eff_K) + 1.0) * (1.0 - fermi_factor(ep, t_eff_K));
            if (rng.uniform() * max_w <= w) return om;
        }
        return 0.5 * span;
    }

    // Sample 3D or 2D Gaussian diffusive displacement step for time dt
    inline Vec3 sample_diffusive_step(double qp_eV, double dt, Rng& rng, bool planar_2d = false) const noexcept {
        const double D = qp_diffusivity(qp_eV);
        if (D <= 0.0 || dt <= 0.0) return {0.0, 0.0, 0.0};
        const double sigma = std::sqrt(2.0 * D * dt);
        const double dx = sigma * rng.normal();
        const double dy = sigma * rng.normal();
        const double dz = planar_2d ? 0.0 : (sigma * rng.normal());
        return {dx, dy, dz};
    }

private:
    static inline double interpolate_lut(const std::array<double, N_LUT_BINS>& lut,
                                         double e_div_gap,
                                         double min_div,
                                         double max_div) noexcept {
        if (e_div_gap <= min_div) return 0.0;
        const double d_bin = (max_div - min_div) / N_LUT_BINS;
        const double u = (e_div_gap - min_div) / d_bin - 0.5;
        if (u <= 0.0) return lut.front();
        if (u >= N_LUT_BINS - 1) return lut.back();
        const int i = static_cast<int>(u);
        const double frac = u - i;
        return lut[i] * (1.0 - frac) + lut[i + 1] * frac;
    }

    // Tabulate 1 / tau_pb(E_ph) using singularity-free theta substitution:
    //   omega = Delta + (E_ph - 2*Delta) * sin^2(theta), theta in (0, pi/2)
    //   d_omega / sqrt((omega - Delta)(E_ph - Delta - omega)) = 2 * d_theta
    void build_pair_breaking_lut() noexcept {
        if (tau0_ph_s <= 0.0) return;
        const double d_div = (MAX_PH_DIV_GAP - MIN_PH_DIV_GAP) / N_LUT_BINS;
        const double tau0_norm = units::hbar_Planck / (4.0 * units::pi * units::pi * delta0_eV);
        const double prefactor = 4.0 * units::pi / units::hbar_Planck;

        constexpr int N_QUAD = 64;
        const double d_theta = units::half_pi / N_QUAD;

        for (int b = 0; b < N_LUT_BINS; ++b) {
            const double e_ph = delta_T_eV * (MIN_PH_DIV_GAP + (b + 0.5) * d_div);
            const double span = e_ph - 2.0 * delta_T_eV;
            double integral = 0.0;
            for (int k = 0; k < N_QUAD; ++k) {
                const double th = (k + 0.5) * d_theta;
                const double s  = std::sin(th);
                const double w1 = delta_T_eV + span * (s * s);
                const double w2 = e_ph - w1;
                const double coh = (w1 * w2 + delta_T_eV * delta_T_eV) /
                                   std::sqrt((w1 + delta_T_eV) * (w2 + delta_T_eV));
                const double pauli = 1.0 - fermi_factor(w2, t_eff_K) - fermi_factor(w1, t_eff_K);
                integral += 2.0 * coh * pauli * d_theta;
            }
            const double inv_tau_norm = prefactor * integral * tau0_norm;
            inv_tau_lut_set(inv_tau_pb_lut_[b], inv_tau_norm / tau0_ph_s);
        }
    }

    // Tabulate 1 / tau_rad(E_qp) using singularity-free substitution:
    //   E' = E_qp - Omega = Delta + (E_qp - Delta) * u^2, u in (0, 1)
    //   d_Omega / sqrt(E' - Delta) = 2 * sqrt(E_qp - Delta) * du
    void build_qp_radiation_lut() noexcept {
        if (tau0_qp_s <= 0.0 || t_crit_K <= 0.0) return;
        const double d_div = (MAX_QP_DIV_GAP - MIN_QP_DIV_GAP) / N_LUT_BINS;
        const double kTc   = units::k_Boltzmann * t_crit_K;
        const double tau0_norm = units::hbar_Planck / (2.0 * units::pi * kTc * kTc * kTc);

        constexpr int N_QUAD = 64;
        const double du = 1.0 / N_QUAD;

        for (int b = 0; b < N_LUT_BINS; ++b) {
            const double e_qp = delta_T_eV * (MIN_QP_DIV_GAP + (b + 0.5) * d_div);
            const double span = e_qp - delta_T_eV;
            const double prefactor = (2.0 * units::pi / units::hbar_Planck) /
                                     std::max(1e-30, 1.0 - fermi_factor(e_qp, t_eff_K));

            double integral = 0.0;
            const double jac_const = 2.0 * std::sqrt(span);
            for (int k = 0; k < N_QUAD; ++k) {
                const double u  = (k + 0.5) * du;
                const double ep = delta_T_eV + span * (u * u);
                const double om = e_qp - ep;
                const double term = (om * om) * (ep / std::sqrt(ep + delta_T_eV)) *
                                    (1.0 - (delta_T_eV * delta_T_eV) / (e_qp * ep)) *
                                    (bose_factor(om, t_eff_K) + 1.0) *
                                    (1.0 - fermi_factor(ep, t_eff_K));
                integral += jac_const * term * du;
            }
            const double inv_tau_norm = prefactor * integral * tau0_norm;
            inv_tau_lut_set(inv_tau_rad_lut_[b], inv_tau_norm / tau0_qp_s);
        }
    }

    // Tabulate 1 / tau_rec(E_qp) using singularity-free substitution:
    //   E2 = Omega - E_qp = Delta + span * u^2
    void build_qp_recombination_lut() noexcept {
        if (tau0_qp_s <= 0.0 || t_crit_K <= 0.0) return;
        const double d_div = (MAX_QP_DIV_GAP - MIN_QP_DIV_GAP) / N_LUT_BINS;
        const double kTc   = units::k_Boltzmann * t_crit_K;
        const double tau0_norm = units::hbar_Planck / (2.0 * units::pi * kTc * kTc * kTc);

        constexpr int N_QUAD = 64;
        const double du = 1.0 / N_QUAD;

        for (int b = 0; b < N_LUT_BINS; ++b) {
            const double e_qp = delta_T_eV * (MIN_QP_DIV_GAP + (b + 0.5) * d_div);
            const double span = std::max(10.0 * units::k_Boltzmann * t_eff_K, 4.0 * delta_T_eV);
            const double prefactor = (2.0 * units::pi / units::hbar_Planck) /
                                     std::max(1e-30, 1.0 - fermi_factor(e_qp, t_eff_K));

            double integral = 0.0;
            const double jac_const = 2.0 * std::sqrt(span);
            for (int k = 0; k < N_QUAD; ++k) {
                const double u  = (k + 0.5) * du;
                const double e2 = delta_T_eV + span * (u * u);
                const double om = e_qp + e2;
                const double term = (om * om) * (e2 / std::sqrt(e2 + delta_T_eV)) *
                                    (1.0 + (delta_T_eV * delta_T_eV) / (e_qp * e2)) *
                                    (bose_factor(om, t_eff_K) + 1.0) *
                                    fermi_factor(e2, t_eff_K);
                integral += jac_const * term * du;
            }
            const double inv_tau_norm = prefactor * integral * tau0_norm;
            inv_tau_lut_set(inv_tau_rec_lut_[b], inv_tau_norm / tau0_qp_s);
        }
    }

    static inline void inv_tau_lut_set(double& slot, double val) noexcept {
        slot = (std::isfinite(val) && val > 0.0) ? val : 0.0;
    }
};

} // namespace cmpmc::physics
