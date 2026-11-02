#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/core/units.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cmpmc::physics {

// Charge-carrier physics kernels: stateless functions for
//   1. the non-parabolic band energy function,
//   2. electron intervalley scattering rates:
//      - quadratic field parameterisation (arXiv:1807.07986, Eq. 1),
//      - linear field parameterisation    (arXiv:1807.07986, Eq. 2),
//      - deformation-potential matrix elements plus Erginsoy neutral-impurity
//        scattering (C. Jacoboni and L. Reggiani, Rev. Mod. Phys. 55, 645 (1983));
//   3. bulk trapping and trap impact-ionisation rates;
//   4. the spatial distribution of the initial charge cloud.

// --- 1. Non-Parabolic Band Energy Function (Jacoboni 1983) ---
//   gamma(E) = sqrt(E * (1 + alpha * E)) * (1 + 2 * alpha * E)
inline double nonparabolic_energy_func_SI(double e_eV, double alpha_inv_eV) noexcept {
    if (e_eV <= 0.0) return 0.0;
    const double e_J = e_eV * units::e_SI;
    const double aE  = alpha_inv_eV * e_eV;
    return std::sqrt(e_J * (1.0 + aE)) * (1.0 + 2.0 * aE);
}

// --- 2. Matrix-Element Optical + Neutral Impurity Intervalley Rate ---
// Rate in SI units [1/s].
inline double intervalley_matrix_rate(double ekin_eV, const CrystalConfig& cfg) noexcept {
    const size_t n_valleys = cfg.valley_dirs.size();
    if (n_valleys < 2 || ekin_eV <= 0.0) return 0.0;

    const double hbar_SI  = units::hbar_Planck * units::e_SI;
    const double hbar2_SI = hbar_SI * hbar_SI;

    // Density-of-states effective mass: m_DOS = (m_xx * m_yy * m_zz)^(1/3) * m_e
    const double m_prod_ratio = cfg.e_mass_tensor.x * cfg.e_mass_tensor.y * cfg.e_mass_tensor.z;
    const double m_dos_ratio  = std::cbrt(m_prod_ratio);
    const double m_dos32_SI   = std::sqrt(m_prod_ratio) * std::pow(units::m_e_kg, 1.5);

    // (a) Optical deformation-potential intervalley phonon scattering
    const double n_final_valleys = static_cast<double>(2 * n_valleys - 1);
    const double opt_scale = (n_final_valleys * m_dos32_SI) /
                             (std::sqrt(2.0) * units::pi * hbar2_SI * cfg.density);

    double opt_rate = 0.0;
    const size_t n_op = std::min(cfg.iv_energy_eV.size(), cfg.iv_deform_eV_m.size());
    for (size_t i = 0; i < n_op; ++i) {
        const double e_min_eV = cfg.iv_energy_eV[i];
        if (ekin_eV <= e_min_eV || e_min_eV <= 0.0) continue;

        const double d_op_SI  = cfg.iv_deform_eV_m[i] * units::e_SI; // [J/m]
        const double e_min_SI = e_min_eV * units::e_SI;              // [J]
        const double e_func   = nonparabolic_energy_func_SI(ekin_eV - e_min_eV, cfg.alpha_nonparabolic);
        opt_rate += opt_scale * ((d_op_SI * d_op_SI) / e_min_SI) * e_func;
    }

    // (b) Erginsoy neutral impurity scattering rate (n_I in m^-3 = neutral_impurity_cm3 * 1e6)
    double imp_rate = 0.0;
    if (cfg.neutral_impurity_cm3 > 0.0 && cfg.epsilon_rel > 0.0) {
        const double n_I_m3  = cfg.neutral_impurity_cm3 * 1.0e6;
        const double e_T_eV  = 0.75 * m_dos_ratio / cfg.epsilon_rel;
        const double e_trk_J = ekin_eV * units::e_SI;
        const double e_T_J   = e_T_eV * units::e_SI;
        imp_rate = (4.0 * std::sqrt(2.0) * n_I_m3 * hbar2_SI * std::sqrt(e_trk_J)) /
                   (m_dos32_SI * (e_trk_J + e_T_J));
    }

    return opt_rate + imp_rate;
}

// --- 3. Unified Intervalley Scattering Rate Dispatcher ---
inline double intervalley_scattering_rate(double ekin_eV,
                                          const Vec3& E_field_Vm,
                                          const CrystalConfig& cfg) noexcept {
    if (cfg.valley_dirs.size() < 2) return 0.0;

    if (cfg.iv_model == IVRateModel::Matrix) {
        return intervalley_matrix_rate(ekin_eV, cfg);
    }
    if (ekin_eV < cfg.iv_min_energy) return 0.0;

    if (cfg.iv_model == IVRateModel::Quadratic) {
        const double e0 = cfg.iv_quad_field;
        return cfg.iv_quad_rate * std::pow(e0 * e0 + E_field_Vm.norm2(), 0.5 * cfg.iv_quad_power);
    } else {
        const double E_Vcm = E_field_Vm.norm() * 0.01;
        return cfg.iv_lin_rate0 + cfg.iv_lin_rate1 * std::pow(E_Vcm, cfg.iv_lin_power);
    }
}

// --- 4. Bulk Charge Trapping & Trap Impact Ionization Rates [1/s] ---
// Constant mean-free-path model: rate = v / lambda for each enabled channel.
inline constexpr double charge_trapping_rate(ExcitationType type,
                                             double v_mag,
                                             const CrystalConfig& cfg) noexcept {
    const double mfp = (type == ExcitationType::Electron) ? cfg.e_trapping_mfp_m : cfg.h_trapping_mfp_m;
    return (mfp > 0.0 && mfp < 1e20) ? (v_mag / mfp) : 0.0;
}

inline constexpr double trap_ionization_e_rate(ExcitationType projectile_type,
                                               double v_mag,
                                               const CrystalConfig& cfg) noexcept {
    const double mfp = (projectile_type == ExcitationType::Electron)
        ? cfg.e_trap_ion_d0_mfp_m
        : cfg.h_trap_ion_d_mfp_m;
    return (mfp > 0.0 && mfp < 1e20) ? (v_mag / mfp) : 0.0;
}

inline constexpr double trap_ionization_h_rate(ExcitationType projectile_type,
                                               double v_mag,
                                               const CrystalConfig& cfg) noexcept {
    const double mfp = (projectile_type == ExcitationType::Electron)
        ? cfg.e_trap_ion_a_mfp_m
        : cfg.h_trap_ion_a0_mfp_m;
    return (mfp > 0.0 && mfp < 1e20) ? (v_mag / mfp) : 0.0;
}

// --- 5. 3D Initial Charge Cloud Spatial Generator ---
// Distributes `n_pairs` e-/h+ pairs in a sphere around `center` with linear
// radial falloff and radius scaled to match the crystal unit-cell density.
struct ChargeCloudGenerator {
    static inline double compute_cloud_radius(uint64_t n_pairs, const CrystalConfig& cfg) noexcept {
        // Estimate unit-cell lattice constant a0 from atomic mass and crystal density
        // (8 atoms per diamond unit cell -> V_cell = 8 * A / (N_A * rho))
        constexpr double N_Avogadro = 6.02214076e23;
        const double vol_cell = (cfg.density > 0.0 && cfg.atomic_A > 0.0)
            ? (8.0 * (cfg.atomic_A * 1e-3) / (N_Avogadro * cfg.density))
            : std::pow(0.55 * units::nm, 3.0);
        const double a0      = std::cbrt(vol_cell);
        const double r_scale = 0.25 * a0; // 8 atoms per diamond unit cell -> r = a0/4 * N^(1/3)
        return std::max(a0, r_scale * std::cbrt(static_cast<double>(n_pairs)));
    }

    static inline Vec3 sample_offset(double r_max, Rng& rng) noexcept {
        if (r_max <= 0.0) return {0.0, 0.0, 0.0};
        const double u = rng.uniform();
        const double r = r_max * (1.0 - std::sqrt(std::max(0.0, 1.0 - u * u)));
        return rng.isotropic_dir() * r;
    }
};

} // namespace cmpmc::physics
