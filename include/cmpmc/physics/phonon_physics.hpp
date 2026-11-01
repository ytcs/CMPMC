#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/core/units.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include <algorithm>
#include <cmath>

namespace cmpmc::physics {

// ============================================================================
// Pure Phonon Physics Equations
//
// This file contains the human-readable analytical equations for:
//   1. Anisotropic Christoffel Wave Equation & Group Velocity
//   2. Tamura Isotope Elastic Scattering Rate (PRB 31, 1985)
//   3. Tamura Three-Phonon Anharmonic Downconversion Rates & PDFs (PRB 31, 1985)
//
// Note for Physicists:
//   If you modify any rate or PDF formula in this file, the engine's startup
//   LUT builder (`include/cmpmc/engine/lut.hpp`) will automatically tabulate
//   your new equations into fast O(1) lookup tables at initialization!
// ============================================================================

// --- 1. Bulk Scattering Rates ---

// Rayleigh-like elastic scattering off isotopic mass fluctuations:
//   Gamma_iso(nu) = B_iso * nu^4   [1/s]
inline constexpr double isotope_scattering_rate(double energy_eV, double B_iso) noexcept {
    const double nu  = energy_eV / units::h_Planck;
    const double nu2 = nu * nu;
    return B_iso * (nu2 * nu2);
}

// Three-phonon anharmonic decay of Longitudinal (L) phonons (L -> T+T, L -> L'+T):
//   Gamma_anh(nu) = A_anh * nu^5   [1/s] (for L-mode; 0 for FT/ST modes)
inline constexpr double anharmonic_decay_rate(PhononMode mode, double energy_eV, double A_anh) noexcept {
    if (mode != PhononMode::L) return 0.0;
    const double nu  = energy_eV / units::h_Planck;
    const double nu2 = nu * nu;
    return A_anh * (nu2 * nu2 * nu);
}

// Polycrystalline grain-boundary elastic scattering rate (G4CMPPhononPolycrystalElasticScattering):
//   Gamma_poly = v_group / lambda_polycrystal   [1/s]
inline constexpr double polycrystal_scattering_rate(double vg_mag, double polycrystal_mfp_m) noexcept {
    if (polycrystal_mfp_m <= 0.0) return 0.0;
    return vg_mag / polycrystal_mfp_m;
}

// --- 2. Density-of-States Polarization Mode Sampling ---

inline PhononMode sample_dos_mode(const CrystalConfig& cfg, Rng& rng) noexcept {
    const double u = rng.uniform();
    if (u < cfg.fLDOS) return PhononMode::L;
    if (u < cfg.fLDOS + cfg.fSTDOS) return PhononMode::ST;
    return PhononMode::FT;
}

inline PhononMode sample_transverse_mode(const CrystalConfig& cfg, Rng& rng) noexcept {
    const double p_st = cfg.fSTDOS / (cfg.fSTDOS + cfg.fFTDOS);
    return (rng.uniform() < p_st) ? PhononMode::ST : PhononMode::FT;
}

// --- 3. Two-Body Anharmonic Decay Kinematics & Energy Partition PDFs ---
// Reference: S. Tamura, Phys. Rev. B 31, 2574 (1985)
// Let x = E_daughter1 / E_parent_L, and m = v_L / v_T.

// Opening angle theta_L' of outgoing L' phonon in L -> L' + T decay
inline double decay_angle_L_to_LT(double x, double v_L, double v_T) noexcept {
    const double m = v_L / v_T;
    const double cos_th = (1.0 + (x * x) * (1.0 - m * m) + 2.0 * x * m * m - m * m) / (2.0 * x);
    return std::acos(std::clamp(cos_th, -1.0, 1.0));
}

// Opening angle theta_T of outgoing T phonon in L -> L' + T decay
inline double decay_angle_T_from_LT(double x, double v_L, double v_T) noexcept {
    const double m = v_L / v_T;
    const double cos_th = (1.0 - x * x + m * m * (1.0 - x) * (1.0 - x)) / (2.0 * m * (1.0 - x));
    return std::acos(std::clamp(cos_th, -1.0, 1.0));
}

// Opening angle theta_T of outgoing T phonon with energy fraction x in L -> T + T decay
inline double decay_angle_TT(double x, double v_L, double v_T) noexcept {
    const double m = v_L / v_T;
    const double cos_th = (1.0 + x * x - (1.0 - x) * (1.0 - x) - 1.0 / (m * m)) / (2.0 * x / m);
    return std::acos(std::clamp(cos_th, -1.0, 1.0));
}

// Unnormalized energy partition PDF for L -> L' + T (G4CMP FV4)
// Valid on x in [(m - 1)/(m + 1), 1]
inline double tamura_lt_pdf(double x, const CrystalConfig& cfg) noexcept {
    const double m = cfg.v_sound / cfg.v_trans;
    const double xmin = (m - 1.0) / (m + 1.0);
    if (x <= xmin || x >= 1.0) return 0.0;

    const double s = (1.0 + (x * x) * (1.0 - m * m) + 2.0 * x * m * m - m * m) / (2.0 * x);
    if (std::abs(s) >= 1.0) return 0.0;
    const double thL = std::acos(s);
    const double thT = decay_angle_T_from_LT(x, cfg.v_sound, cfg.v_trans);

    const double ca = cfg.beta + cfg.lambda + (1.0 - m * m) * (cfg.gamma + cfg.mu);
    const double cb = cfg.beta + cfg.lambda + 2.0 * m * m * (cfg.gamma + cfg.mu);
    const double cc = cfg.beta + cfg.lambda;
    const double cd = cfg.gamma + cfg.lambda - m * m * (cfg.gamma + cfg.mu);

    const double u = (ca * x * x + cb * x * (1.0 - x) + cc * (1.0 - x) * (1.0 - x) +
                      cd * (1.0 - x * x - m * m * (1.0 - x) * (1.0 - x)) * std::cos(thL) * std::cos(thL - thT))
                     * std::sin(thL - thT);
    return u * u * (x / (1.0 - x));
}

// Unnormalized energy partition PDF for L -> T + T (G4CMP FTV5)
// Valid on x in [(1 - 1/m)/2, (1 + 1/m)/2]
inline double tamura_tt_pdf(double x, const CrystalConfig& cfg) noexcept {
    const double m = cfg.v_sound / cfg.v_trans;
    const double xmin = 0.5 * (1.0 - 1.0 / m);
    const double xmax = 0.5 * (1.0 + 1.0 / m);
    if (x <= xmin || x >= xmax) return 0.0;

    const double d = (1.0 + x * x - (1.0 - x) * (1.0 - x) - 1.0 / (m * m)) / (2.0 * x / m);
    if (std::abs(d) >= 1.0) return 0.0;

    const double g    = cfg.gamma;
    const double l    = cfg.lambda;
    const double mu   = cfg.mu;
    const double d2   = d * d;
    const double s2   = 1.0 - d2;
    const double s    = std::sqrt(s2);
    const double s2th = 2.0 * s * d;
    const double c2th = d2 - s2;

    const double a = 0.5 * (1.0 - m * m) * (g + mu) * d * (s2th - s);
    const double b = ((g + l) * d + m * (g + mu) * (c2th - d2));
    const double c = m * m * (g + mu) * d * s * c2th + 0.5 * m * (g + l) * (2.0 * d2 * s - s2th);

    const double term1 = (a + b * x + c * x * x);
    const double term2 = 0.25 * (g + l - m * m * (g + mu)) * s2th +
                         m * (g + l - (g + mu)) * s * x -
                         m * m * (g + mu) * s * c2th * x * x;

    return (term1 * term1 + term2 * term2) * x * (1.0 - x);
}

// --- 4. Exact Christoffel Acoustic Wave Solver ---
// Solves the 3x3 real symmetric eigenvalue problem:
//   D_il(k) * e_l = v_phase^2 * e_i,   where D_il(k) = (1/rho) * C_ijlm * k_j * k_m
// and computes the anisotropic group velocity:
//   v_group,d = (1 / (rho * v_phase)) * e_i * C_ijld * k_j * e_l
inline void solve_christoffel_exact(const double C[3][3][3][3],
                                    double density,
                                    const Vec3& k_unit,
                                    Vec3 vg_out[3],
                                    double vphase_out[3]) noexcept {
    double D[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const double inv_rho = 1.0 / density;
    for (int i = 0; i < 3; ++i) {
        for (int l = 0; l < 3; ++l) {
            double sum = 0.0;
            for (int j = 0; j < 3; ++j) {
                for (int m = 0; m < 3; ++m) {
                    sum += C[i][j][l][m] * k_unit[j] * k_unit[m];
                }
            }
            D[i][l] = sum * inv_rho;
        }
    }

    // Cyclic Jacobi eigensolver for 3x3 symmetric matrix D
    double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 12; ++sweep) {
        const double off = std::abs(D[0][1]) + std::abs(D[0][2]) + std::abs(D[1][2]);
        if (off < 1e-15 * (std::abs(D[0][0]) + std::abs(D[1][1]) + std::abs(D[2][2]))) break;

        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(D[p][q]) < 1e-18) continue;
                const double tau = (D[q][q] - D[p][p]) / (2.0 * D[p][q]);
                const double t = (tau >= 0.0)
                    ? (1.0 / (tau + std::sqrt(1.0 + tau * tau)))
                    : (-1.0 / (-tau + std::sqrt(1.0 + tau * tau)));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = t * c;

                const double d_pp = D[p][p];
                const double d_qq = D[q][q];
                const double d_pq = D[p][q];
                D[p][p] = c * c * d_pp - 2.0 * s * c * d_pq + s * s * d_qq;
                D[q][q] = s * s * d_pp + 2.0 * s * c * d_pq + c * c * d_qq;
                D[p][q] = D[q][p] = 0.0;

                for (int r = 0; r < 3; ++r) {
                    if (r != p && r != q) {
                        const double d_rp = D[r][p];
                        const double d_rq = D[r][q];
                        D[r][p] = D[p][r] = c * d_rp - s * d_rq;
                        D[r][q] = D[q][r] = s * d_rp + c * d_rq;
                    }
                    const double v_rp = V[r][p];
                    const double v_rq = V[r][q];
                    V[r][p] = c * v_rp - s * v_rq;
                    V[r][q] = s * v_rp + c * v_rq;
                }
            }
        }
    }

    const double evals[3] = {D[0][0], D[1][1], D[2][2]};
    const Vec3 evecs[3] = {
        {V[0][0], V[1][0], V[2][0]},
        {V[0][1], V[1][1], V[2][1]},
        {V[0][2], V[1][2], V[2][2]}
    };

    // Mode classification: L = polarization most parallel to k_unit; FT/ST sorted by phase speed
    int long_idx = 0;
    double max_dot = -1.0;
    for (int i = 0; i < 3; ++i) {
        const double d = std::abs(evecs[i].dot(k_unit));
        if (d > max_dot) {
            max_dot = d;
            long_idx = i;
        }
    }
    const int rem1 = (long_idx + 1) % 3;
    const int rem2 = (long_idx + 2) % 3;
    const int ft_idx = (evals[rem1] >= evals[rem2]) ? rem1 : rem2;
    const int st_idx = (evals[rem1] >= evals[rem2]) ? rem2 : rem1;
    const int mode_indices[3] = {long_idx, ft_idx, st_idx};

    for (int m = 0; m < 3; ++m) {
        const int idx = mode_indices[m];
        const double vp = std::sqrt(std::max(evals[idx], 1.0));
        vphase_out[m] = vp;
        const Vec3& e = evecs[idx];
        const Vec3 slow = k_unit / vp;

        Vec3 vg{0.0, 0.0, 0.0};
        for (int dim = 0; dim < 3; ++dim) {
            double sum = 0.0;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int l = 0; l < 3; ++l) {
                        sum += e[i] * C[i][j][l][dim] * slow[j] * e[l];
                    }
                }
            }
            vg[dim] = sum * inv_rho;
        }
        vg_out[m] = vg;
    }
}

} // namespace cmpmc::physics
