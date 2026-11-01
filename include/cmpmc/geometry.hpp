#pragma once

#include "cmpmc/core.hpp"
#include "cmpmc/geometry/mesh_field.hpp"
#include "cmpmc/geometry/surface_mask.hpp"
#include "cmpmc/kaplan.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace cmpmc {

struct BoundaryHit {
    double distance;
    Vec3 inward_normal;
    SurfaceId surface_id;
};

// ============================================================================
// SurfaceProperty: Boundary Reflection, Absorption & Kaplan Sensor Model
// Replaces: G4CMPSurfaceProperty & G4CMPPhononElectrode
// ============================================================================
struct SurfaceProperty {
    double prob_abs      = 0.0;   // Probability of phonon entering sensor / boundary absorption
    double prob_refl     = 1.0;   // Probability of phonon reflection if not absorbed
    double refl_cutoff   = 350.0; // GHz cutoff for specular/diffuse polynomials
    double anh_cutoff    = 520.0; // GHz cutoff for surface anharmonicity

    // Charge boundary properties (G4CMPDriftBoundaryProcess)
    double prob_charge_abs  = 1.0; // Probability of charge absorption/recombination at boundary
    double prob_charge_spec = 0.0; // Specular reflection probability for reflected charges
    double min_k_elec       = 0.0; // Minimum perpendicular wavevector k_perp [1/m] for e- absorption
    double min_k_hole       = 0.0; // Minimum perpendicular wavevector k_perp [1/m] for h+ absorption

    // Polynomial coefficients vs frequency in GHz (evaluated via Horner's method)
    std::vector<double> spec_coeffs = {1.0};
    std::vector<double> diff_coeffs = {0.0};
    std::vector<double> anh_coeffs  = {0.0};

    std::optional<KaplanFilm> kaplan_sensor = std::nullopt;
    std::optional<geometry::ElectrodePatternMask> electrode_mask = std::nullopt;

    // Evaluate polynomial c[0] + c[1]*f + c[2]*f^2 + ... via Horner's rule
    static inline double eval_poly(const std::vector<double>& coeffs, double f_ghz) noexcept {
        if (coeffs.empty()) return 0.0;
        double acc = 0.0;
        for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
            acc = acc * f_ghz + (*it);
        }
        return acc;
    }

    inline double specular_prob(double energy_eV) const noexcept {
        const double f_ghz = std::min((energy_eV / units::h_Planck) / units::GHz, refl_cutoff);
        return std::clamp(eval_poly(spec_coeffs, f_ghz), 0.0, 1.0);
    }

    // Standard G4CMP examples/phonon Top/Bottom Aluminum sensor surface
    static SurfaceProperty g4cmp_al_sensor_surface() {
        SurfaceProperty s;
        s.prob_abs      = 1.0;
        s.prob_refl     = 1.0;
        s.anh_cutoff    = 520.0;
        s.refl_cutoff   = 350.0;
        s.spec_coeffs   = {0.928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};
        s.diff_coeffs   = {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
        s.anh_coeffs    = {0.0, 0.0, 0.0, 0.0, 0.0, 1.51e-14};
        s.kaplan_sensor = KaplanFilm::aluminum_sensor();
        return s;
    }

    // Standard G4CMP examples/phonon bare crystal sidewall surface
    static SurfaceProperty g4cmp_bare_wall_surface() {
        SurfaceProperty s;
        s.prob_abs      = 0.0;
        s.prob_refl     = 1.0;
        s.anh_cutoff    = 520.0;
        s.refl_cutoff   = 350.0;
        s.spec_coeffs   = {0.928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};
        s.diff_coeffs   = {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
        s.anh_coeffs    = {0.0, 0.0, 0.0, 0.0, 0.0, 1.51e-14};
        s.kaplan_sensor = std::nullopt;
        return s;
    }
};

// ============================================================================
// Analytic Crystal Primitives (Replaces G4Tubs / G4Box + G4Navigator)
// ============================================================================
struct CylinderSolid {
    double radius = 3.81 * units::cm; // Default G4CMP examples/phonon radius
    double half_z = 1.27 * units::cm; // Default G4CMP examples/phonon half-height

    constexpr bool contains(const Vec3& p) const noexcept {
        return (p.x * p.x + p.y * p.y <= radius * radius + 1e-12) &&
               (std::abs(p.z) <= half_z + 1e-12);
    }

    // Reflect exterior points back into the cylinder volume (G4CMPChargeCloud::AdjustToVolume)
    inline void fold_inside(Vec3& p) const noexcept {
        if (p.z > half_z)  p.z =  2.0 * half_z - p.z;
        if (p.z < -half_z) p.z = -2.0 * half_z - p.z;
        const double r_xy = std::hypot(p.x, p.y);
        if (r_xy > radius && r_xy > 0.0) {
            const double r_new = std::max(0.0, 2.0 * radius - r_xy);
            const double scale = r_new / r_xy;
            p.x *= scale;
            p.y *= scale;
        }
    }

    // Compute exact distance to exit cylinder from interior point `p` along unit vector `v`
    inline BoundaryHit intersect_exit(const Vec3& p, const Vec3& v) const noexcept {
        double t_min = std::numeric_limits<double>::infinity();
        SurfaceId surf = SurfaceId::SideWall;

        // 1. Planar endcaps z = +half_z (TopFace) and z = -half_z (BottomFace)
        if (v.z > 1e-15) {
            const double tz = (half_z - p.z) / v.z;
            if (tz >= 0.0 && tz < t_min) {
                t_min = tz;
                surf = SurfaceId::TopFace;
            }
        } else if (v.z < -1e-15) {
            const double tz = (-half_z - p.z) / v.z;
            if (tz >= 0.0 && tz < t_min) {
                t_min = tz;
                surf = SurfaceId::BottomFace;
            }
        }

        // 2. Cylindrical sidewall x^2 + y^2 = R^2
        const double a = v.x * v.x + v.y * v.y;
        if (a > 1e-18) {
            const double b = p.x * v.x + p.y * v.y;
            const double c = p.x * p.x + p.y * p.y - radius * radius;
            const double disc = std::max(0.0, b * b - a * c);
            const double tr = (-b + std::sqrt(disc)) / a;
            if (tr >= 0.0 && tr < t_min) {
                t_min = tr;
                surf = SurfaceId::SideWall;
            }
        }

        t_min = std::max(0.0, t_min);
        Vec3 inward_normal{0.0, 0.0, -1.0};
        if (surf == SurfaceId::TopFace) {
            inward_normal = {0.0, 0.0, -1.0};
        } else if (surf == SurfaceId::BottomFace) {
            inward_normal = {0.0, 0.0, 1.0};
        } else {
            const Vec3 hit_p = p + v * t_min;
            const double r_xy = std::hypot(hit_p.x, hit_p.y);
            inward_normal = (r_xy > 1e-15)
                ? Vec3{-hit_p.x / r_xy, -hit_p.y / r_xy, 0.0}
                : Vec3{-1.0, 0.0, 0.0};
        }

        return {t_min, inward_normal, surf};
    }
};

struct BoxSolid {
    double half_x = 0.5 * units::cm;
    double half_y = 0.5 * units::cm;
    double half_z = 0.5 * units::cm;

    inline BoundaryHit intersect_exit(const Vec3& p, const Vec3& v) const noexcept {
        double t_min = std::numeric_limits<double>::infinity();
        Vec3 inward_normal{0.0, 0.0, -1.0};
        SurfaceId surf = SurfaceId::SideWall;

        auto check_axis = [&](double pos_i, double vel_i, double half_i,
                              const Vec3& norm_pos, const Vec3& norm_neg,
                              SurfaceId s_pos, SurfaceId s_neg) noexcept {
            if (vel_i > 1e-15) {
                const double t = (half_i - pos_i) / vel_i;
                if (t >= 0.0 && t < t_min) {
                    t_min = t;
                    inward_normal = norm_pos;
                    surf = s_pos;
                }
            } else if (vel_i < -1e-15) {
                const double t = (-half_i - pos_i) / vel_i;
                if (t >= 0.0 && t < t_min) {
                    t_min = t;
                    inward_normal = norm_neg;
                    surf = s_neg;
                }
            }
        };

        check_axis(p.x, v.x, half_x, {-1, 0, 0}, {1, 0, 0}, SurfaceId::SideWall, SurfaceId::SideWall);
        check_axis(p.y, v.y, half_y, {0, -1, 0}, {0, 1, 0}, SurfaceId::SideWall, SurfaceId::SideWall);
        check_axis(p.z, v.z, half_z, {0, 0, -1}, {0, 0, 1}, SurfaceId::TopFace, SurfaceId::BottomFace);

        return {std::max(0.0, t_min), inward_normal, surf};
    }
};

} // namespace cmpmc
