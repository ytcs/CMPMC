#pragma once

// Phonon interface transmission models (AMM / DMM); used by SurfaceProperty for
// substrate-to-film entry.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

namespace cmpmc::physics {

// Phonon transmission across a welded interface between two isotropic solids
// (substrate 1 -> film 2): acoustic mismatch (AMM) and diffuse mismatch (DMM) models.
//
//  AMM (elastodynamic boundary conditions: continuity of the
//  displacement u_x, u_z and of the tractions tau_xz, tau_zz across z = 0):
//    * normal incidence, mode nu:   t = 4 Z1 Z2 / (Z1 + Z2)^2,  Z = rho v_nu;
//    * SH (transverse, polarised perpendicular to the plane of incidence), Snell refraction:
//      t(theta1) = 4 Z1 cos(theta1) Z2 cos(theta2) / (Z1 cos(theta1) + Z2 cos(theta2))^2,
//      0 beyond the critical angle;
//    * P-SV (incident P = L or SV), oblique, with mode conversion: the 4x4 complex linear
//      system for (R_P, R_S, T_P, T_S); energy-flux coefficients
//      E_j = rho_j v_j^2 Re(s_zj)|A_j|^2 / (rho_1 v_inc^2 s_z,inc), s_z the vertical slowness
//      (evanescent waves carry no flux). sum_j E_j = 1 is not imposed, so it serves as
//      an independent check of the solution.
//  DMM (Debye density of states):  t_12 = S2 / (S1 + S2),  S = v_L^-2 + 2 v_T^-2.
//
//  Scope: isotropic media (the substrate is replaced by its isotropic (v_L, v_T)); half-space
//  limit (no film interference; valid for wavelength << d_film).
//  The AMM of an ideal interface is an upper bound: measured Si->Al transmission is 0.3-0.55
//  (arXiv:1805.02495, p. 1).
struct IsotropicMedium {
    double rho = 0.0;   // kg/m^3
    double vL  = 0.0;   // m/s
    double vT  = 0.0;   // m/s

    static IsotropicMedium from_bulk_shear(double rho, double K_Pa, double G_Pa) noexcept {
        return {rho, std::sqrt((K_Pa + 4.0 * G_Pa / 3.0) / rho), std::sqrt(G_Pa / rho)};
    }
    static IsotropicMedium from_lame(double rho, double lam_Pa, double mu_Pa) noexcept {
        return {rho, std::sqrt((lam_Pa + 2.0 * mu_Pa) / rho), std::sqrt(mu_Pa / rho)};
    }
    double mu() const noexcept { return rho * vT * vT; }
    double lam() const noexcept { return rho * (vL * vL - 2.0 * vT * vT); }
};

inline double amm_normal(double Z1, double Z2) noexcept {
    return 4.0 * Z1 * Z2 / ((Z1 + Z2) * (Z1 + Z2));
}

inline double amm_sh(const IsotropicMedium& m1, const IsotropicMedium& m2, double theta1) noexcept {
    const double s2 = m2.vT / m1.vT * std::sin(theta1);
    if (s2 >= 1.0) return 0.0;
    const double c1 = std::cos(theta1), c2 = std::sqrt(1.0 - s2 * s2);
    const double a = m1.rho * m1.vT * c1, b = m2.rho * m2.vT * c2;
    return (a + b > 0.0) ? 4.0 * a * b / ((a + b) * (a + b)) : 0.0;
}

inline double dmm(const IsotropicMedium& m1, const IsotropicMedium& m2) noexcept {
    const double s1 = 1.0 / (m1.vL * m1.vL) + 2.0 / (m1.vT * m1.vT);
    const double s2 = 1.0 / (m2.vL * m2.vL) + 2.0 / (m2.vT * m2.vT);
    return s2 / (s1 + s2);
}

// Energy-flux reflection / transmission coefficients for P or SV incidence.
struct PSVCoefficients { double R_P = 0, R_S = 0, T_P = 0, T_S = 0;
                         double transmitted() const noexcept { return T_P + T_S; }
                         double sum() const noexcept { return R_P + R_S + T_P + T_S; } };

inline PSVCoefficients amm_psv(const IsotropicMedium& m1, const IsotropicMedium& m2,
                               bool incident_P, double theta1) noexcept {
    using C = std::complex<double>;
    const double v_inc = incident_P ? m1.vL : m1.vT;
    const double p = std::sin(theta1) / v_inc;                          // horizontal slowness
    auto sz = [p](double v) { return std::sqrt(C(1.0 / (v * v) - p * p, 0.0)); };   // Im >= 0 branch
    const C a1 = sz(m1.vL), b1 = sz(m1.vT), a2 = sz(m2.vL), b2 = sz(m2.vT);
    // displacement polarisation and tractions (tau_xz, tau_zz), common factor i*omega dropped
    struct Wave { C dx, dz, txz, tzz; };
    auto wave = [p](const IsotropicMedium& m, bool P, C s) -> Wave {
        const double v = P ? m.vL : m.vT;
        const C dx = P ? C(v * p) : v * s;
        const C dz = P ? v * s : C(-v * p);
        const double mu = m.mu(), la = m.lam();
        return {dx, dz, mu * (s * dx + p * dz), la * (p * dx + s * dz) + 2.0 * mu * s * dz};
    };
    const Wave inc = wave(m1, incident_P, incident_P ? a1 : b1);
    const Wave rP = wave(m1, true, -a1), rS = wave(m1, false, -b1);
    const Wave tP = wave(m2, true, a2), tS = wave(m2, false, b2);
    // rows: u_x, u_z, tau_xz, tau_zz continuity; unknowns (R_P, R_S, T_P, T_S)
    std::array<std::array<C, 5>, 4> M{};
    auto col = [&](int j, const Wave& w, double sgn) {
        M[0][j] = sgn * w.dx; M[1][j] = sgn * w.dz; M[2][j] = sgn * w.txz; M[3][j] = sgn * w.tzz; };
    col(0, rP, 1.0); col(1, rS, 1.0); col(2, tP, -1.0); col(3, tS, -1.0);
    M[0][4] = -inc.dx; M[1][4] = -inc.dz; M[2][4] = -inc.txz; M[3][4] = -inc.tzz;
    for (int c = 0; c < 4; ++c) {                                         // Gauss-Jordan, partial pivot
        int piv = c;
        for (int r = c + 1; r < 4; ++r) if (std::abs(M[r][c]) > std::abs(M[piv][c])) piv = r;
        std::swap(M[c], M[piv]);
        const C d = M[c][c];
        if (std::abs(d) == 0.0) return {};
        for (int k = c; k < 5; ++k) M[c][k] /= d;
        for (int r = 0; r < 4; ++r) if (r != c) {
            const C f = M[r][c];
            for (int k = c; k < 5; ++k) M[r][k] -= f * M[c][k];
        }
    }
    const double flux_inc = m1.rho * v_inc * v_inc * std::real(incident_P ? a1 : b1);
    auto E = [&](int j, const IsotropicMedium& m, double v, C s) {
        return m.rho * v * v * std::abs(std::real(s)) * std::norm(M[j][4]) / flux_inc; };
    PSVCoefficients out;
    out.R_P = E(0, m1, m1.vL, a1); out.R_S = E(1, m1, m1.vT, b1);
    out.T_P = E(2, m2, m2.vL, a2); out.T_S = E(3, m2, m2.vT, b2);
    return out;
}

// Mode- and angle-resolved AMM substrate -> film transmission for the phonon engine.
// L phonons are P waves; a transverse phonon's polarisation relative to the plane of incidence
// is not tracked, so t_T = (t_SV + t_SH)/2 (isotropic-substrate approximation).
struct IsotropicAMMInterface {
    IsotropicMedium substrate{}, film{};
    double transmission(bool longitudinal, double cos_theta) const noexcept {
        const double th = std::acos(std::clamp(std::abs(cos_theta), 0.0, 1.0));
        if (longitudinal) return amm_psv(substrate, film, true, th).transmitted();
        return 0.5 * (amm_psv(substrate, film, false, th).transmitted() + amm_sh(substrate, film, th));
    }
    // flux (cos theta sin theta) average over the incident hemisphere
    double flux_average(bool longitudinal, int n = 4000) const noexcept {
        double s = 0.0;
        for (int i = 0; i < n; ++i) {
            const double th = (i + 0.5) / n * 0.5 * 3.14159265358979323846;
            s += transmission(longitudinal, std::cos(th)) * 2.0 * std::sin(th) * std::cos(th);
        }
        return s * 0.5 * 3.14159265358979323846 / n;
    }
};

} // namespace cmpmc::physics
