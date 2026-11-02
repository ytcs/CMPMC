#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/core/units.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include <algorithm>
#include <cmath>

namespace cmpmc::physics {

// Phonon physics kernels: stateless functions for
//   1. the anisotropic Christoffel wave equation and group velocity,
//   2. isotope (mass-defect) elastic scattering (S. Tamura, Phys. Rev. B 31, 2574 (1985)),
//   3. three-phonon anharmonic down-conversion rates and energy-partition PDFs.
//
// CrystalLattice tabulates the rate and PDF functions defined here into O(1)
// lookup tables at construction, so a change to a formula in this file takes
// effect without touching the sampling code.

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

// Polycrystalline grain-boundary elastic scattering rate:
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

// Opening angle theta_T of outgoing T phonon with energy fraction x in L -> T + T decay.
// Momentum closure k_L = k_T1 + k_T2 with k_L = E/v_L, k_T1 = xE/v_T, k_T2 = (1-x)E/v_T gives
//   cos(theta_T1) = (1 - m^2 (1-x)^2 + m^2 x^2) / (2 m x)
// (Leman, arXiv:1109.1193, Eq. (13)). The second
// daughter uses x -> 1 - x (Leman Eq. (14)).
inline double decay_angle_TT(double x, double v_L, double v_T) noexcept {
    const double m = v_L / v_T;
    const double cos_th = (1.0 - m * m * (1.0 - x) * (1.0 - x) + m * m * x * x) / (2.0 * m * x);
    return std::acos(std::clamp(cos_th, -1.0, 1.0));
}

// Unnormalized energy partition PDF for L -> L' + T, x = E_L'/E_L, m = v_L/v_T:
//   P(x) ~ (1/x^2) (1 - x^2)^2 [(1 + x)^2 - m^2 (1 - x)^2] [1 + x^2 - m^2 (1 - x)^2]^2
// on (m - 1)/(m + 1) < x < 1 (Leman, arXiv:1109.1193 Eq. (9), quoting Tamura PRB 31,
// 2574 (1985)).
inline double tamura_lt_pdf(double x, const CrystalConfig& cfg) noexcept {
    const double m = cfg.v_sound / cfg.v_trans;
    const double xmin = (m - 1.0) / (m + 1.0);
    if (x <= xmin || x >= 1.0) return 0.0;
    const double a = 1.0 - x * x;
    const double b = (1.0 + x) * (1.0 + x) - m * m * (1.0 - x) * (1.0 - x);
    const double c = 1.0 + x * x - m * m * (1.0 - x) * (1.0 - x);
    return (a * a / (x * x)) * b * c * c;
}

// Unnormalized energy partition PDF for L -> T + T in the energy fraction x = E_T1/E_L,
// on (1 - 1/m)/2 < x < (1 + 1/m)/2. Leman, arXiv:1109.1193 Eq. (12) is written in
// X = m x, (m - 1)/2 < X < (m + 1)/2 (the Jacobian dX/dx = m is constant):
//   P ~ [A + B m X - B X^2]^2 + [C X (m - X) - D/(m - X) (X - m - (1 - m^2)/(4X))]^2,
//   A = (1 - m^2)[beta + lambda + (1 + m^2)(gamma + mu)]/2, B = beta + lambda + 2 m^2 (gamma + mu),
//   C = beta + lambda + 2 (gamma + mu), D = (1 - m^2)(2 beta + 4 gamma + lambda + 3 mu)
// Symmetric under x -> 1 - x.
inline double tamura_tt_pdf(double x, const CrystalConfig& cfg) noexcept {
    const double m = cfg.v_sound / cfg.v_trans;
    const double xmin = 0.5 * (1.0 - 1.0 / m);
    const double xmax = 0.5 * (1.0 + 1.0 / m);
    if (x <= xmin || x >= xmax) return 0.0;
    const double be = cfg.beta / 1e11, ga = cfg.gamma / 1e11, la = cfg.lambda / 1e11, mu = cfg.mu / 1e11;
    const double A = 0.5 * (1.0 - m * m) * (be + la + (1.0 + m * m) * (ga + mu));
    const double B = be + la + 2.0 * m * m * (ga + mu);
    const double C = be + la + 2.0 * (ga + mu);
    const double D = (1.0 - m * m) * (2.0 * be + 4.0 * ga + la + 3.0 * mu);
    const double X = m * x;
    const double t1 = A + B * m * X - B * X * X;
    const double t2 = C * X * (m - X) - D / (m - X) * (X - m - (1.0 - m * m) / (4.0 * X));
    return t1 * t1 + t2 * t2;
}

// Mode labelling and group velocities from the Christoffel eigen-decomposition:
//   L = polarization most parallel to k_unit; FT/ST ordered by phase speed;
//   v_group,d = (1/(rho v_p)) e_i C_ijld k_j e_l.
inline void christoffel_modes_from_eigen(const double C[3][3][3][3], double density, const Vec3& k_unit,
                                         const double evals[3], const Vec3 evecs[3],
                                         Vec3 vg_out[3], double vphase_out[3]) noexcept {
    const double inv_rho = 1.0 / density;
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
    christoffel_modes_from_eigen(C, density, k_unit, evals, evecs, vg_out, vphase_out);
}

// --- 5. Closed-form Christoffel solve (trigonometric Cardano roots + cross-product eigenvectors) ---
// Eigenvalues of the symmetric 3x3 D: with q = tr D/3, p = sqrt(tr((D - qI)^2)/6), B = (D - qI)/p,
// r = det(B)/2: lambda_k = q + 2p cos(acos(r)/3 + 2 pi k/3) (Smith, CACM 4, 168 (1961); standard
// algebra). Eigenvector of lambda: the largest cross product of two rows of D - lambda I.
// Near a degeneracy (relative eigenvalue gap < gap_tol, e.g. FT/ST along [100], [111]) the
// eigenvectors are ill-conditioned: fall back to the Jacobi solve.
inline void solve_christoffel_cardano(const double C[3][3][3][3], double density, const Vec3& k_unit,
                                      Vec3 vg_out[3], double vphase_out[3], double gap_tol = 1e-4) noexcept {
    double D[3][3];
    const double inv_rho = 1.0 / density;
    for (int i = 0; i < 3; ++i)
        for (int l = i; l < 3; ++l) {
            double sum = 0.0;
            for (int j = 0; j < 3; ++j)
                for (int m = 0; m < 3; ++m) sum += C[i][j][l][m] * k_unit[j] * k_unit[m];
            D[i][l] = D[l][i] = sum * inv_rho;
        }
    const double q = (D[0][0] + D[1][1] + D[2][2]) / 3.0;
    const double p1 = D[0][1] * D[0][1] + D[0][2] * D[0][2] + D[1][2] * D[1][2];
    const double p2 = (D[0][0] - q) * (D[0][0] - q) + (D[1][1] - q) * (D[1][1] - q) + (D[2][2] - q) * (D[2][2] - q) + 2.0 * p1;
    const double p = std::sqrt(p2 / 6.0);
    if (!(p > 1e-12 * std::abs(q))) { solve_christoffel_exact(C, density, k_unit, vg_out, vphase_out); return; }
    double B[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = (D[i][j] - (i == j ? q : 0.0)) / p;
    const double detB = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
                      - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
                      + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
    const double r = std::clamp(0.5 * detB, -1.0, 1.0);
    const double phi = std::acos(r) / 3.0;
    double ev[3];
    ev[0] = q + 2.0 * p * std::cos(phi);
    ev[2] = q + 2.0 * p * std::cos(phi + 2.0943951023931957);
    ev[1] = 3.0 * q - ev[0] - ev[2];
    const double scale = std::abs(ev[0]);
    if (std::abs(ev[0] - ev[1]) < gap_tol * scale || std::abs(ev[1] - ev[2]) < gap_tol * scale) {
        solve_christoffel_exact(C, density, k_unit, vg_out, vphase_out);
        return;
    }
    Vec3 evecs[3];
    for (int k = 0; k < 3; ++k) {
        const Vec3 r0{D[0][0] - ev[k], D[0][1], D[0][2]};
        const Vec3 r1{D[1][0], D[1][1] - ev[k], D[1][2]};
        const Vec3 r2{D[2][0], D[2][1], D[2][2] - ev[k]};
        const Vec3 c01 = r0.cross(r1), c02 = r0.cross(r2), c12 = r1.cross(r2);
        const double n01 = c01.norm2(), n02 = c02.norm2(), n12 = c12.norm2();
        const Vec3& best = (n01 >= n02 && n01 >= n12) ? c01 : ((n02 >= n12) ? c02 : c12);
        evecs[k] = best.unit();
    }
    christoffel_modes_from_eigen(C, density, k_unit, ev, evecs, vg_out, vphase_out);
}

} // namespace cmpmc::physics
