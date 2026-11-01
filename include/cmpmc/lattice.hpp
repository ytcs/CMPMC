#pragma once

#include "cmpmc/core.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace cmpmc {

// ============================================================================
// CrystalLattice: Anisotropic Elasticity, Christoffel Solver, Group Velocity
// LUT, and Precomputed Anharmonic Decay Inverse-CDF Tables
// Replaces: G4LatticeLogical, G4LatticePhysical, G4CMPPhononKinematics,
//           G4CMPPhononKinTable, and G4CMPPhononDownconversionProcess PDFs
// ============================================================================
class CrystalLattice {
public:
    static constexpr int N_CTHETA = 181; // cos(theta) in [-1, 1]
    static constexpr int N_PHI    = 361; // phi in [0, 2*pi]
    static constexpr int N_CDF    = 1024; // Inverse-CDF bins for anharmonic decay

    struct Config {
        std::string name = "Ge";
        double density   = 5323.0;              // kg / m^3
        // Voigt 6x6 stiffness matrix in Pa (cubic uses C11, C12, C44)
        double C11       = 126.0 * units::GPa;
        double C12       =  44.0 * units::GPa;
        double C44       =  67.0 * units::GPa;

        // Tamura anharmonic & isotope scattering parameters (S. Tamura, PRB 31, 1985)
        double beta      = -73.2 * units::GPa;
        double gamma     = -70.8 * units::GPa;
        double lambda    =  37.6 * units::GPa;
        double mu        =  56.1 * units::GPa;

        double B_iso     = 3.67e-41;            // s^3 (isotope scattering rate coeff)
        double A_anh     = 6.43e-55;            // s^4 (anharmonic decay rate coeff)
        double fTT_frac  = 0.74;                // Branching fraction for L -> T + T

        // Phonon density of states (S. Tamura et al., PRB 44, 1991)
        double fLDOS     = 0.097834;
        double fSTDOS    = 0.535390;
        double fFTDOS    = 0.366770;

        double debye_freq = 2.0 * units::THz;
        double v_sound    = 5324.2077;          // m/s (isotropic longitudinal sound speed)
        double v_trans    = 3258.7879;          // m/s (isotropic transverse sound speed)

        // Crystal frame to Detector frame rotation (Miller orientation)
        Mat3 orient       = Mat3::identity();
    };

    struct Kinematics {
        Vec3 vg_dir;   // Unit group velocity direction in detector frame
        double vg_mag; // Group velocity magnitude (m/s)
    };

    explicit CrystalLattice(const Config& cfg) : cfg_(cfg) {
        build_stiffness_tensor();
        build_group_velocity_lut();
        build_decay_cdf_tables();
    }

    // Standard G4CMP CrystalMaps/Ge/config.txt preset
    static Config germanium_config() noexcept {
        return Config{};
    }

    // Standard G4CMP CrystalMaps/Si/config.txt preset
    static Config silicon_config() noexcept {
        Config c;
        c.name       = "Si";
        c.density    = 2329.0; // kg/m^3
        c.C11        = 165.6 * units::GPa;
        c.C12        =  63.9 * units::GPa;
        c.C44        =  79.5 * units::GPa;
        c.beta       = -42.9 * units::GPa;
        c.gamma      = -94.5 * units::GPa;
        c.lambda     =  52.4 * units::GPa;
        c.mu         =  68.0 * units::GPa;
        c.B_iso      = 2.43e-42;
        c.A_anh      = 7.41e-56;
        c.fTT_frac   = 0.74;
        c.fLDOS      = 0.093;
        c.fSTDOS     = 0.531;
        c.fFTDOS     = 0.376;
        c.debye_freq = 15.0 * units::THz;
        c.v_sound    = 9000.0;
        c.v_trans    = 5400.0;
        return c;
    }

    const Config& config() const noexcept { return cfg_; }
    double debye_energy() const noexcept { return cfg_.debye_freq * units::h_Planck; }

    // ========================================================================
    // Fast Bilinear LUT Query: wavevector k_dir -> (vg_dir, vg_mag)
    // Replaces runtime Christoffel matrix solve on every step (~5 ns vs ~280 ns)
    // ========================================================================
    inline Kinematics map_k_to_vg(PhononMode mode, const Vec3& k_det) const noexcept {
        // Rotate wavevector from detector frame into crystal lattice frame
        const Vec3 k = cfg_.orient.mul_transpose(k_det).unit();

        const double ctheta = std::clamp(k.z, -1.0, 1.0);
        double phi = std::atan2(k.y, k.x);
        if (phi < 0.0) phi += units::two_pi;

        const double u = (ctheta + 1.0) * 0.5 * (N_CTHETA - 1);
        const double v = (phi / units::two_pi) * (N_PHI - 1);

        const int i0 = std::clamp(static_cast<int>(u), 0, N_CTHETA - 2);
        const int j0 = std::clamp(static_cast<int>(v), 0, N_PHI - 2);
        const double du = u - i0;
        const double dv = v - j0;

        const auto& grid = vg_lut_[static_cast<int>(mode)];
        const Vec3 vg_cryst =
            grid[i0][j0]         * ((1.0 - du) * (1.0 - dv)) +
            grid[i0 + 1][j0]     * (du * (1.0 - dv)) +
            grid[i0][j0 + 1]     * ((1.0 - du) * dv) +
            grid[i0 + 1][j0 + 1] * (du * dv);

        const double mag = vg_cryst.norm();
        // Rotate group velocity back into detector frame
        const Vec3 vg_det = cfg_.orient * (vg_cryst / std::max(mag, 1e-12));
        return {vg_det.unit(), mag};
    }

    // Exact Christoffel matrix solver (used to populate LUT or for reference checks)
    void solve_christoffel_exact(const Vec3& k_unit, Vec3 vg_out[3], double vphase_out[3]) const noexcept {
        // 1. Build symmetric 3x3 Christoffel matrix D_il = (1/rho) * C_ijlm * k_j * k_m
        double D[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        const double inv_rho = 1.0 / cfg_.density;
        for (int i = 0; i < 3; ++i) {
            for (int l = 0; l < 3; ++l) {
                double sum = 0.0;
                for (int j = 0; j < 3; ++j) {
                    for (int m = 0; m < 3; ++m) {
                        sum += C_[i][j][l][m] * k_unit[j] * k_unit[m];
                    }
                }
                D[i][l] = sum * inv_rho;
            }
        }

        // 2. Solve 3x3 real symmetric eigensystem via exact cyclic Jacobi rotations
        double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int sweep = 0; sweep < 12; ++sweep) {
            const double off = std::abs(D[0][1]) + std::abs(D[0][2]) + std::abs(D[1][2]);
            if (off < 1e-15 * ( std::abs(D[0][0]) + std::abs(D[1][1]) + std::abs(D[2][2]) )) break;

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

        double evals[3] = {D[0][0], D[1][1], D[2][2]};
        Vec3 evecs[3] = {
            {V[0][0], V[1][0], V[2][0]},
            {V[0][1], V[1][1], V[2][1]},
            {V[0][2], V[1][2], V[2][2]}
        };

        // 3. Classify modes matching G4CMPPhononKinematics:
        // Longitudinal (L) = eigenvector most parallel to k_unit (largest |e . k|)
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

        // 4. Compute group velocity vg[dim] = (1 / (rho * v_phase)) * e_i * C_ijld * k_j * e_l
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
                            sum += e[i] * C_[i][j][l][dim] * slow[j] * e[l];
                        }
                    }
                }
                vg[dim] = sum * inv_rho;
            }
            vg_out[m] = vg;
        }
    }

    // ========================================================================
    // Bulk Scattering Rates & Mode Sampling (Tamura PRB 31, 1985)
    // ========================================================================
    // Isotope scattering rate (1/s): Gamma_iso = B_iso * nu^4
    inline double isotope_rate(double energy_eV) const noexcept {
        const double nu = energy_eV / units::h_Planck;
        const double nu2 = nu * nu;
        return cfg_.B_iso * (nu2 * nu2);
    }

    // Anharmonic downconversion rate (1/s): Gamma_anh = A_anh * nu^5 (only for L-mode)
    inline double downconversion_rate(PhononMode mode, double energy_eV) const noexcept {
        if (mode != PhononMode::L) return 0.0;
        const double nu = energy_eV / units::h_Planck;
        const double nu2 = nu * nu;
        return cfg_.A_anh * (nu2 * nu2 * nu);
    }

    // Sample new polarization mode after elastic isotope scattering according to DOS
    inline PhononMode sample_dos_mode(Rng& rng) const noexcept {
        const double u = rng.uniform();
        if (u < cfg_.fLDOS) return PhononMode::L;
        if (u < cfg_.fLDOS + cfg_.fSTDOS) return PhononMode::ST;
        return PhononMode::FT;
    }

    // Sample transverse polarization mode (ST vs FT) for anharmonic decay products
    inline PhononMode sample_transverse_mode(Rng& rng) const noexcept {
        const double p_st = cfg_.fSTDOS / (cfg_.fSTDOS + cfg_.fFTDOS);
        return (rng.uniform() < p_st) ? PhononMode::ST : PhononMode::FT;
    }

    // ========================================================================
    // O(1) Inverse-CDF Sampling for Anharmonic Downconversion
    // Replaces 30-iteration rejection loops in G4CMPPhononDownconversionProcess
    // ========================================================================
    // Sample x = E_L' / E_L for L -> L' + T decay (FV4 distribution)
    inline double sample_lt_energy_fraction(Rng& rng) const noexcept {
        const double u = rng.uniform() * (N_CDF - 1);
        const int i = std::min(static_cast<int>(u), N_CDF - 2);
        const double frac = u - i;
        return cdf_lt_[i] * (1.0 - frac) + cdf_lt_[i + 1] * frac;
    }

    // Sample x = E_T1 / E_L for L -> T + T decay (FTV5 distribution)
    inline double sample_tt_energy_fraction(Rng& rng) const noexcept {
        const double u = rng.uniform() * (N_CDF - 1);
        const int i = std::min(static_cast<int>(u), N_CDF - 2);
        const double frac = u - i;
        return cdf_tt_[i] * (1.0 - frac) + cdf_tt_[i + 1] * frac;
    }

    // Opening angles for decay secondaries (exact energy-momentum conservation)
    inline double theta_L_to_LT(double x) const noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        const double cos_th = (1.0 + (x * x) * (1.0 - m * m) + 2.0 * x * m * m - m * m) / (2.0 * x);
        return std::acos(std::clamp(cos_th, -1.0, 1.0));
    }

    inline double theta_T_from_LT(double x) const noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        const double cos_th = (1.0 - x * x + m * m * (1.0 - x) * (1.0 - x)) / (2.0 * m * (1.0 - x));
        return std::acos(std::clamp(cos_th, -1.0, 1.0));
    }

    inline double theta_TT(double x) const noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        const double cos_th = (1.0 + x * x - (1.0 - x) * (1.0 - x) - 1.0 / (m * m)) / (2.0 * x / m);
        return std::acos(std::clamp(cos_th, -1.0, 1.0));
    }

private:
    Config cfg_;
    double C_[3][3][3][3]{};
    // vg_lut_[mode][ctheta_idx][phi_idx] stores velocity vector in crystal frame
    std::array<std::array<std::array<Vec3, N_PHI>, N_CTHETA>, 3> vg_lut_{};
    std::array<double, N_CDF> cdf_lt_{};
    std::array<double, N_CDF> cdf_tt_{};

    static constexpr int voigt_index(int i, int j) noexcept {
        if (i == j) return i;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
        return 5;
    }

    void build_stiffness_tensor() noexcept {
        double V[6][6]{};
        V[0][0] = V[1][1] = V[2][2] = cfg_.C11;
        V[0][1] = V[1][0] = V[0][2] = V[2][0] = V[1][2] = V[2][1] = cfg_.C12;
        V[3][3] = V[4][4] = V[5][5] = cfg_.C44;

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int l = 0; l < 3; ++l)
                    for (int m = 0; m < 3; ++m)
                        C_[i][j][l][m] = V[voigt_index(i, j)][voigt_index(l, m)];
    }

    void build_group_velocity_lut() noexcept {
        for (int i = 0; i < N_CTHETA; ++i) {
            const double ctheta = -1.0 + 2.0 * (static_cast<double>(i) / (N_CTHETA - 1));
            const double stheta = std::sqrt(std::max(0.0, 1.0 - ctheta * ctheta));
            for (int j = 0; j < N_PHI; ++j) {
                const double phi = units::two_pi * (static_cast<double>(j) / (N_PHI - 1));
                Vec3 k_dir{stheta * std::cos(phi), stheta * std::sin(phi), ctheta};
                if (stheta < 1e-7) {
                    // Slightly perturb away from exact pole to avoid degenerate transverse eigenvectors
                    k_dir = Vec3{1e-6, 0.0, ctheta}.unit();
                }
                Vec3 vg[3];
                double vp[3];
                solve_christoffel_exact(k_dir, vg, vp);
                for (int m = 0; m < 3; ++m) {
                    vg_lut_[m][i][j] = vg[m];
                }
            }
        }
    }

    // Exact Tamura FV4 PDF (G4CMPPhononDownconversionProcess::FV4)
    double pdf_fv4(double x) const noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        const double xmin = (m - 1.0) / (m + 1.0);
        if (x <= xmin || x >= 1.0) return 0.0;

        const double s = (1.0 + (x * x) * (1.0 - m * m) + 2.0 * x * m * m - m * m) / (2.0 * x);
        if (std::abs(s) >= 1.0) return 0.0;
        const double thL = std::acos(s);
        const double thT = theta_T_from_LT(x);

        const double ca = cfg_.beta + cfg_.lambda + (1.0 - m * m) * (cfg_.gamma + cfg_.mu);
        const double cb = cfg_.beta + cfg_.lambda + 2.0 * m * m * (cfg_.gamma + cfg_.mu);
        const double cc = cfg_.beta + cfg_.lambda;
        const double cd = cfg_.gamma + cfg_.lambda - m * m * (cfg_.gamma + cfg_.mu);

        const double u = (ca * x * x + cb * x * (1.0 - x) + cc * (1.0 - x) * (1.0 - x) +
                          cd * (1.0 - x * x - m * m * (1.0 - x) * (1.0 - x)) * std::cos(thL) * std::cos(thL - thT))
                         * std::sin(thL - thT);
        return u * u * (x / (1.0 - x));
    }

    // Exact Tamura FTV5 PDF (G4CMPPhononDownconversionProcess::FTV5)
    double pdf_ftv5(double x) const noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        const double xmin = (1.0 - 1.0 / m) * 0.5;
        const double xmax = (1.0 + 1.0 / m) * 0.5;
        if (x <= xmin || x >= xmax) return 0.0;

        const double d = (1.0 + x * x - (1.0 - x) * (1.0 - x) - 1.0 / (m * m)) / (2.0 * x / m);
        if (std::abs(d) >= 1.0) return 0.0;

        const double g = cfg_.gamma;
        const double l = cfg_.lambda;
        const double mu = cfg_.mu;
        const double d2 = d * d;
        const double s2 = 1.0 - d2;
        const double s = std::sqrt(s2);
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

    void build_decay_cdf_tables() noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        build_single_cdf(cdf_lt_, (m - 1.0) / (m + 1.0), 1.0,
                         [this](double x) { return pdf_fv4(x); });
        build_single_cdf(cdf_tt_, 0.5 * (1.0 - 1.0 / m), 0.5 * (1.0 + 1.0 / m),
                         [this](double x) { return pdf_ftv5(x); });
    }

    template <typename Func>
    static void build_single_cdf(std::array<double, N_CDF>& inv_cdf,
                                 double xmin, double xmax, Func&& pdf_fn) noexcept {
        constexpr int N_INT = 4096;
        std::array<double, N_INT + 1> cdf{};
        const double dx = (xmax - xmin) / N_INT;
        double sum = 0.0;
        cdf[0] = 0.0;
        for (int i = 0; i < N_INT; ++i) {
            const double x_mid = xmin + (i + 0.5) * dx;
            sum += std::max(0.0, pdf_fn(x_mid)) * dx;
            cdf[i + 1] = sum;
        }
        if (sum <= 0.0) {
            for (int k = 0; k < N_CDF; ++k) inv_cdf[k] = 0.5 * (xmin + xmax);
            return;
        }
        for (double& v : cdf) v /= sum;

        inv_cdf[0] = xmin;
        inv_cdf[N_CDF - 1] = xmax;
        int j = 0;
        for (int k = 1; k < N_CDF - 1; ++k) {
            const double target = static_cast<double>(k) / (N_CDF - 1);
            while (j < N_INT && cdf[j + 1] < target) ++j;
            const double denom = cdf[j + 1] - cdf[j];
            const double t = (denom > 1e-16) ? ((target - cdf[j]) / denom) : 0.0;
            inv_cdf[k] = xmin + (j + t) * dx;
        }
    }
};

} // namespace cmpmc
