#pragma once

#include "cmpmc/core.hpp"
#include "cmpmc/materials/config_reader.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include "cmpmc/materials/presets.hpp"
#include "cmpmc/physics/charge_physics.hpp"
#include "cmpmc/physics/phonon_physics.hpp"
#include "cmpmc/physics/superconductor_3d.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace cmpmc {

// ============================================================================
// CrystalLattice: Precomputes 2D Group-Velocity LUT, 1D Inverse-CDF Decay
// Tables, and Herring-Vogt Valley Tensors from pure physics equations in
// `include/cmpmc/physics/phonon_physics.hpp` and `materials/crystal_config.hpp`.
// ============================================================================
class CrystalLattice {
public:
    using Config = CrystalConfig;

    static constexpr int N_CTHETA = 181; // cos(theta) in [-1, 1]
    static constexpr int N_PHI    = 361; // phi in [0, 2*pi]
    static constexpr int N_OCT    = 256; // 256x256 zero-trig Octahedral FP32 LUT
    static constexpr int N_CDF    = 1024; // Inverse-CDF bins for anharmonic decay

    struct alignas(16) OctLutEntry {
        float vx;
        float vy;
        float vz;
        float vmag;
    };

    struct Kinematics {
        Vec3 vg_dir;   // Unit group velocity direction in detector frame
        double vg_mag; // Group velocity magnitude (m/s)
    };

    explicit CrystalLattice(const Config& cfg)
        : cfg_(cfg),
          vg_lut_(3 * N_CTHETA * N_PHI),
          vg_oct_lut_(3 * N_OCT * N_OCT),
          sc_3d_(cfg) {
        build_stiffness_tensor();
        build_group_velocity_lut();
        build_octahedral_lut();
        build_decay_cdf_tables();
        build_charge_tensors();
    }

    static Config germanium_config() { return materials::Germanium(); }
    static Config silicon_config()   { return materials::Silicon(); }

    const Config& config() const noexcept { return cfg_; }
    const physics::Superconductor3D& superconductor() const noexcept { return sc_3d_; }
    double debye_energy() const noexcept { return cfg_.debye_freq * units::h_Planck; }

    // Kinematics query: wavevector k_dir -> (vg_dir, vg_mag).
    // The legacy 181 x 361 bilinear (cos theta, phi) table is uniform in cos theta, so
    // near the poles (the crystal z axis, e.g. [001] towards the sensor faces) its polar
    // step is ~8 degrees: tests/test_christoffel measured a 99th-percentile direction
    // error of 14-16 mrad for L phonons in Si/Ge (max 84 mrad; 0.3-0.4 rad for FT/ST).
    // Default is now the 256 x 256 octahedral table (q99 0.2 mrad for L, 2.5 mrad for
    // FT/ST; errors > 1 degree only on FT/ST degeneracy lines, 0.1% of directions; and
    // faster). Exact per-call Christoffel solve (~1.7 us) and the legacy table remain
    // selectable.
    enum class KinematicsMode : uint8_t { Octahedral = 0, Exact = 1, Bilinear = 2 };
    void set_kinematics_mode(KinematicsMode m) noexcept { kin_mode_ = m; }
    KinematicsMode kinematics_mode() const noexcept { return kin_mode_; }

    inline Kinematics map_k_to_vg_exact(PhononMode mode, const Vec3& k_det) const noexcept {
        const Vec3 k = cfg_.orient.mul_transpose(k_det).unit();
        Vec3 vg[3];
        double vp[3];
        physics::solve_christoffel_cardano(C_, cfg_.density, k, vg, vp);
        const Vec3& v = vg[static_cast<int>(mode)];
        const double mag = v.norm();
        return {(cfg_.orient * (v / std::max(mag, 1e-12))).unit(), mag};
    }

    inline Kinematics map_k_to_vg(PhononMode mode, const Vec3& k_det) const noexcept {
        if (kin_mode_ == KinematicsMode::Octahedral) return map_k_to_vg_oct(mode, k_det);
        if (kin_mode_ == KinematicsMode::Exact) return map_k_to_vg_exact(mode, k_det);
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

        const int m_idx = static_cast<int>(mode);
        const Vec3 vg_cryst =
             lut_at(m_idx, i0,     j0)     * ((1.0 - du) * (1.0 - dv)) +
            lut_at(m_idx, i0 + 1, j0)     * (du * (1.0 - dv)) +
            lut_at(m_idx, i0,     j0 + 1) * ((1.0 - du) * dv) +
            lut_at(m_idx, i0 + 1, j0 + 1) * (du * dv);

        const double mag = vg_cryst.norm();
        const Vec3 vg_det = cfg_.orient * (vg_cryst / std::max(mag, 1e-12));
        return {vg_det.unit(), mag};
    }

    // Zero-Trigonometry Octahedral Equal-Area Mapping: unit 3D vector -> (u, v) in [-1, 1]^2
    static inline void oct_encode(const Vec3& k, double& u_out, double& v_out) noexcept {
        const double l1 = std::abs(k.x) + std::abs(k.y) + std::abs(k.z);
        const double inv_l1 = (l1 > 1e-18) ? (1.0 / l1) : 1.0;
        double px = k.x * inv_l1;
        double py = k.y * inv_l1;
        if (k.z < 0.0) {
            const double ox = (1.0 - std::abs(py)) * (px >= 0.0 ? 1.0 : -1.0);
            const double oy = (1.0 - std::abs(px)) * (py >= 0.0 ? 1.0 : -1.0);
            px = ox;
            py = oy;
        }
        u_out = std::clamp(px, -1.0, 1.0);
        v_out = std::clamp(py, -1.0, 1.0);
    }

    static inline Vec3 oct_decode(double u, double v) noexcept {
        double kz = 1.0 - std::abs(u) - std::abs(v);
        double kx = u;
        double ky = v;
        if (kz < 0.0) {
            kx = (1.0 - std::abs(v)) * (u >= 0.0 ? 1.0 : -1.0);
            ky = (1.0 - std::abs(u)) * (v >= 0.0 ? 1.0 : -1.0);
        }
        return Vec3{kx, ky, kz}.unit();
    }

    // Tier-2 Zero-Trig FP32 Octahedral LUT Query (eliminates std::atan2 from inner loop)
    inline Kinematics map_k_to_vg_oct(PhononMode mode, const Vec3& k_det) const noexcept {
        const Vec3 k = cfg_.orient.mul_transpose(k_det);
        double u_oct = 0.0, v_oct = 0.0;
        oct_encode(k, u_oct, v_oct);

        const double u = (u_oct + 1.0) * 0.5 * (N_OCT - 1);
        const double v = (v_oct + 1.0) * 0.5 * (N_OCT - 1);

        const int i0 = std::clamp(static_cast<int>(u), 0, N_OCT - 2);
        const int j0 = std::clamp(static_cast<int>(v), 0, N_OCT - 2);
        const float du = static_cast<float>(u - i0);
        const float dv = static_cast<float>(v - j0);

        const float w00 = (1.0f - du) * (1.0f - dv);
        const float w10 = du * (1.0f - dv);
        const float w01 = (1.0f - du) * dv;
        const float w11 = du * dv;

        const int m_idx = static_cast<int>(mode);
        const OctLutEntry& e00 = oct_lut_at(m_idx, i0,     j0);
        const OctLutEntry& e10 = oct_lut_at(m_idx, i0 + 1, j0);
        const OctLutEntry& e01 = oct_lut_at(m_idx, i0,     j0 + 1);
        const OctLutEntry& e11 = oct_lut_at(m_idx, i0 + 1, j0 + 1);

        const Vec3 vg_cryst{
            static_cast<double>(e00.vx * w00 + e10.vx * w10 + e01.vx * w01 + e11.vx * w11),
            static_cast<double>(e00.vy * w00 + e10.vy * w10 + e01.vy * w01 + e11.vy * w11),
            static_cast<double>(e00.vz * w00 + e10.vz * w10 + e01.vz * w01 + e11.vz * w11)
        };
        const double mag = vg_cryst.norm();
        const Vec3 vg_det = cfg_.orient * (vg_cryst / std::max(mag, 1e-12));
        return {vg_det.unit(), mag};
    }

    // Forward to readable equations in `include/cmpmc/physics/phonon_physics.hpp`
    inline double isotope_rate(double energy_eV) const noexcept {
        return physics::isotope_scattering_rate(energy_eV, cfg_.B_iso);
    }

    inline double downconversion_rate(PhononMode mode, double energy_eV) const noexcept {
        return physics::anharmonic_decay_rate(mode, energy_eV, cfg_.A_anh);
    }

    inline double polycrystal_rate(double vg_mag) const noexcept {
        return physics::polycrystal_scattering_rate(vg_mag, cfg_.polycrystal_mfp_m);
    }

    inline PhononMode sample_dos_mode(Rng& rng) const noexcept {
        return physics::sample_dos_mode(cfg_, rng);
    }

    inline PhononMode sample_transverse_mode(Rng& rng) const noexcept {
        return physics::sample_transverse_mode(cfg_, rng);
    }

    // O(1) Inverse-CDF Sampling for Anharmonic Downconversion
    inline double sample_lt_energy_fraction(Rng& rng) const noexcept {
        const double u = rng.uniform() * (N_CDF - 1);
        const int i = std::min(static_cast<int>(u), N_CDF - 2);
        const double frac = u - i;
        return cdf_lt_[i] * (1.0 - frac) + cdf_lt_[i + 1] * frac;
    }

    inline double sample_tt_energy_fraction(Rng& rng) const noexcept {
        const double u = rng.uniform() * (N_CDF - 1);
        const int i = std::min(static_cast<int>(u), N_CDF - 2);
        const double frac = u - i;
        return cdf_tt_[i] * (1.0 - frac) + cdf_tt_[i + 1] * frac;
    }

    inline double theta_L_to_LT(double x) const noexcept {
        return physics::decay_angle_L_to_LT(x, cfg_.v_sound, cfg_.v_trans);
    }

    inline double theta_T_from_LT(double x) const noexcept {
        return physics::decay_angle_T_from_LT(x, cfg_.v_sound, cfg_.v_trans);
    }

    inline double theta_TT(double x) const noexcept {
        return physics::decay_angle_TT(x, cfg_.v_sound, cfg_.v_trans);
    }

    // ========================================================================
    // Charge Carrier Herring-Vogt Kinematics
    // ========================================================================
    size_t num_valleys() const noexcept { return valley_rot_.size(); }
    double electron_mass_eV() const noexcept { return e_cond_mass_eV_; }
    double hole_mass_eV() const noexcept { return h_mass_eV_; }
    double l0_electron() const noexcept { return l0_e_; }
    double l0_hole() const noexcept { return l0_h_; }
    double k_sound_electron() const noexcept { return k_sound_e_; }
    double k_sound_hole() const noexcept { return k_sound_h_; }
    // Largest stretch factor of the Herring-Vogt map k -> k* (sqrt(m_c / m_min))
    double hv_max_scale() const noexcept { return std::max({hv_sqrt_inv_.x, hv_sqrt_inv_.y, hv_sqrt_inv_.z}); }

    inline Vec3 rotate_to_valley(size_t iv, const Vec3& v) const noexcept {
        return valley_rot_[iv] * v;
    }

    inline Vec3 rotate_from_valley(size_t iv, const Vec3& v) const noexcept {
        return valley_rot_[iv].mul_transpose(v);
    }

    inline Vec3 ellipsoidal_to_spherical(size_t iv, const Vec3& v) const noexcept {
        const Vec3 v_val = rotate_to_valley(iv, v);
        return {v_val.x * hv_sqrt_inv_.x, v_val.y * hv_sqrt_inv_.y, v_val.z * hv_sqrt_inv_.z};
    }

    inline Vec3 spherical_to_ellipsoidal(size_t iv, const Vec3& v) const noexcept {
        const Vec3 scaled{v.x * hv_sqrt_.x, v.y * hv_sqrt_.y, v.z * hv_sqrt_.z};
        return rotate_from_valley(iv, scaled);
    }

    inline double map_p_to_ekin_electron(size_t iv, const Vec3& p_eV) const noexcept {
        const Vec3 p_hv = ellipsoidal_to_spherical(iv, p_eV);
        return std::sqrt(p_hv.norm2() + e_cond_mass_eV_ * e_cond_mass_eV_) - e_cond_mass_eV_;
    }

    inline double map_p_to_ekin_hole(const Vec3& p_eV) const noexcept {
        return std::sqrt(p_eV.norm2() + h_mass_eV_ * h_mass_eV_) - h_mass_eV_;
    }

    inline Vec3 map_p_to_vel_electron(size_t iv, const Vec3& p_eV) const noexcept {
        const double ekin = map_p_to_ekin_electron(iv, p_eV);
        const Vec3 pv = rotate_to_valley(iv, p_eV);
        const Vec3 v_valley{pv.x / mass_ratio_.x, pv.y / mass_ratio_.y, pv.z / mass_ratio_.z};
        return rotate_from_valley(iv, v_valley) * (units::c_light / (ekin + e_cond_mass_eV_));
    }

    inline Vec3 map_p_to_vel_hole(const Vec3& p_eV) const noexcept {
        const double ekin = map_p_to_ekin_hole(p_eV);
        return p_eV * (units::c_light / (ekin + h_mass_eV_));
    }

    inline Vec3 map_p_to_k_electron(size_t /*iv*/, const Vec3& p_eV) const noexcept {
        return p_eV / units::hbarc;
    }

    inline Vec3 map_k_to_p_electron(size_t /*iv*/, const Vec3& k) const noexcept {
        return k * units::hbarc;
    }

    inline Vec3 map_ekin_to_p_electron(size_t iv, const Vec3& p_dir, double ekin) const noexcept {
        const Vec3 u = p_dir.unit();
        const Vec3 u_hv = ellipsoidal_to_spherical(iv, u);
        const double band_u = u_hv.norm2();
        if (band_u <= 0.0 || ekin <= 0.0) return {0.0, 0.0, 0.0};
        const double p_mag = std::sqrt((ekin * ekin + 2.0 * ekin * e_cond_mass_eV_) / band_u);
        return u * p_mag;
    }

    inline double intervalley_rate(double ekin_eV, const Vec3& E_field_Vm) const noexcept {
        return physics::intervalley_scattering_rate(ekin_eV, E_field_Vm, cfg_);
    }

private:
    Config cfg_;
    double C_[3][3][3][3]{};
    std::vector<Vec3> vg_lut_;
    std::vector<OctLutEntry> vg_oct_lut_;
    std::array<double, N_CDF> cdf_lt_{};
    std::array<double, N_CDF> cdf_tt_{};
    physics::Superconductor3D sc_3d_;
    KinematicsMode kin_mode_ = KinematicsMode::Octahedral;

    inline const Vec3& lut_at(int m, int i, int j) const noexcept {
        return vg_lut_[(m * N_CTHETA + i) * N_PHI + j];
    }
    inline Vec3& lut_at(int m, int i, int j) noexcept {
        return vg_lut_[(m * N_CTHETA + i) * N_PHI + j];
    }
    inline const OctLutEntry& oct_lut_at(int m, int i, int j) const noexcept {
        return vg_oct_lut_[(m * N_OCT + i) * N_OCT + j];
    }
    inline OctLutEntry& oct_lut_at(int m, int i, int j) noexcept {
        return vg_oct_lut_[(m * N_OCT + i) * N_OCT + j];
    }

    double e_cond_mass_eV_ = 0.12 * units::m_e_eV;
    double h_mass_eV_      = 0.35 * units::m_e_eV;
    Vec3 mass_ratio_{1.0, 1.0, 1.0};
    Vec3 hv_sqrt_{1.0, 1.0, 1.0};
    Vec3 hv_sqrt_inv_{1.0, 1.0, 1.0};
    std::vector<Mat3> valley_rot_;
    double l0_e_           = 1e-7;
    double l0_h_           = 1e-7;
    double k_sound_e_      = 1e7;
    double k_sound_h_      = 1e7;

    static constexpr int voigt_index(int i, int j) noexcept {
        if (i == j) return i;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
        return 5;
    }

    void build_stiffness_tensor() noexcept {
        const auto V = cfg_.build_voigt_matrix();

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
                    k_dir = Vec3{1e-6, 0.0, ctheta}.unit();
                }
                Vec3 vg[3];
                double vp[3];
                physics::solve_christoffel_exact(C_, cfg_.density, k_dir, vg, vp);
                for (int m = 0; m < 3; ++m) {
                    lut_at(m, i, j) = vg[m];
                }
            }
        }
    }

    void build_octahedral_lut() noexcept {
        for (int i = 0; i < N_OCT; ++i) {
            const double u = -1.0 + 2.0 * (static_cast<double>(i) / (N_OCT - 1));
            for (int j = 0; j < N_OCT; ++j) {
                const double v = -1.0 + 2.0 * (static_cast<double>(j) / (N_OCT - 1));
                Vec3 k_dir = oct_decode(u, v);
                if (std::abs(k_dir.x) < 1e-6 && std::abs(k_dir.y) < 1e-6) {
                    k_dir = Vec3{1e-6, 0.0, k_dir.z}.unit();
                }
                Vec3 vg[3];
                double vp[3];
                physics::solve_christoffel_exact(C_, cfg_.density, k_dir, vg, vp);
                for (int m = 0; m < 3; ++m) {
                    oct_lut_at(m, i, j) = OctLutEntry{
                        static_cast<float>(vg[m].x),
                        static_cast<float>(vg[m].y),
                        static_cast<float>(vg[m].z),
                        static_cast<float>(vg[m].norm())
                    };
                }
            }
        }
    }

    void build_charge_tensors() noexcept {
        const double m_xx = cfg_.e_mass_tensor.x;
        const double m_yy = cfg_.e_mass_tensor.y;
        const double m_zz = cfg_.e_mass_tensor.z;
        const double m_cond_ratio = 3.0 / (1.0 / m_xx + 1.0 / m_yy + 1.0 / m_zz);

        e_cond_mass_eV_ = m_cond_ratio * units::m_e_eV;
        h_mass_eV_      = cfg_.h_mass_ratio * units::m_e_eV;

        mass_ratio_ = {m_xx / m_cond_ratio, m_yy / m_cond_ratio, m_zz / m_cond_ratio};
        hv_sqrt_    = {std::sqrt(mass_ratio_.x), std::sqrt(mass_ratio_.y), std::sqrt(mass_ratio_.z)};
        hv_sqrt_inv_= {1.0 / hv_sqrt_.x, 1.0 / hv_sqrt_.y, 1.0 / hv_sqrt_.z};

        valley_rot_.clear();
        for (const Vec3& raw_dir : cfg_.valley_dirs) {
            // valley_dirs are crystal-frame axes; all charge vectors (p, E, v) live in
            // the detector frame, so rotate the axis with the crystal orientation.
            const Vec3 vdir = (cfg_.orient * raw_dir.unit()).unit();
            const double vx = vdir.x;
            const double vy = vdir.y;
            const double vz = vdir.z;
            const double f  = std::sqrt(std::max(0.0, 1.0 - vz * vz));
            const double a  = (std::abs(vy) < 1e-12) ? 0.0 : ((std::abs(vx) < 1e-12) ? -1.0 : -vy / f);
            const double b  = (std::abs(vy) < 1e-12) ? 1.0 : ((std::abs(vx) < 1e-12) ?  0.0 :  vx / f);
            const double d  = -vz * b;
            const double e  =  vz * a;
            valley_rot_.push_back(Mat3::from_rows(vdir, {a, b, 0.0}, {d, e, f}));
        }

        auto compute_l0 = [&](double mass_ratio_me, double ac_deform_eV) noexcept -> double {
            if (ac_deform_eV <= 0.0) return 1e-6;
            const double hbar_SI = units::hbar_Planck * units::e_SI;
            const double m_SI    = mass_ratio_me * units::m_e_kg;
            const double ac_SI   = ac_deform_eV * units::e_SI;
            const double hbar4   = (hbar_SI * hbar_SI) * (hbar_SI * hbar_SI);
            const double m3      = m_SI * m_SI * m_SI;
            return 0.5 * units::pi * hbar4 * cfg_.density / (m3 * ac_SI * ac_SI);
        };

        l0_e_ = compute_l0(m_cond_ratio, cfg_.ac_deform_e);
        l0_h_ = compute_l0(cfg_.h_mass_ratio, cfg_.ac_deform_h);

        k_sound_e_ = (e_cond_mass_eV_ * cfg_.v_sound) / (units::c_light * units::hbarc);
        k_sound_h_ = (h_mass_eV_      * cfg_.v_sound) / (units::c_light * units::hbarc);
    }

    void build_decay_cdf_tables() noexcept {
        const double m = cfg_.v_sound / cfg_.v_trans;
        build_single_cdf(cdf_lt_, (m - 1.0) / (m + 1.0), 1.0,
                         [this](double x) { return physics::tamura_lt_pdf(x, cfg_); });
        build_single_cdf(cdf_tt_, 0.5 * (1.0 - 1.0 / m), 0.5 * (1.0 + 1.0 / m),
                         [this](double x) { return physics::tamura_tt_pdf(x, cfg_); });
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
