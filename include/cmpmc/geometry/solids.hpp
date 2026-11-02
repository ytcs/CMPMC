// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include "cmpmc/core/core.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace cmpmc {

// Result of a boundary query from inside a solid.
struct BoundaryHit {
    double distance;
    Vec3 inward_normal;
    SurfaceId surface_id;
};

// Analytic crystal solids. Each solid returns, for an interior point and a unit
// direction, the exit distance, the inward surface normal and the surface id, so
// reflections are handled in place without a volume-navigation step.
struct CylinderSolid {
    double radius = 3.81 * units::cm; // 3-inch diameter
    double half_z = 1.27 * units::cm; // 1-inch thickness

    constexpr bool contains(const Vec3& p) const noexcept {
        return (p.x * p.x + p.y * p.y <= radius * radius + 1e-12) &&
               (std::abs(p.z) <= half_z + 1e-12);
    }

    // Reflect exterior points back into the cylinder volume (mirror fold)
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
