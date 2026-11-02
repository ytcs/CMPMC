// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

// Experimental: next-event estimation onto sparse sensor pads.
// Guarded by require_validated(); see docs/experimental.md.

#include "cmpmc/core/core.hpp"
#include <algorithm>
#include <cmath>

namespace cmpmc::experimental {

// Next-event estimation (shadow rays) onto a sparse sensor pad.
//    At a diffuse boundary reflection at `refl_pos` with `inward_normal`,
//    projects a deterministic weighted contribution directly onto a small
//    circular sensor pad (center `pad_center`, radius `pad_radius`, normal
//    `pad_inward_normal`) using the Lambertian solid-angle view factor:
//      dP_NEE = (cos(theta_emit) * cos(theta_inc) * Area_pad) / (pi * R^2)
//    Returns the weighted energy contribution (eV) and sets out_flight_time_s.
//    Known limitation: a real-space Lambertian estimator that ignores phonon
//    focusing (biased whenever v_g is not parallel to k), and it is not called by
//    SimulationEngine.
inline double evaluate_sparse_sensor_nee(const Vec3& refl_pos,
                                         const Vec3& refl_inward_normal,
                                         const Vec3& pad_center,
                                         const Vec3& pad_inward_normal,
                                         double pad_radius_m,
                                         double phonon_energy_eV,
                                         double track_weight,
                                         double& out_flight_time_s,
                                         double v_group_ms,
                                         bool allow_unvalidated = false) {
    require_validated(allow_unvalidated, "evaluate_sparse_sensor_nee",
                      "Lambertian real-space estimator ignores phonon focusing; not wired into the engine");
    const Vec3 disp = pad_center - refl_pos;
    const double d2 = disp.norm2();
    const double area = units::pi * pad_radius_m * pad_radius_m;
    if (d2 <= area * 4.0) return 0.0; // Avoid near-field singularity

    const double dist = std::sqrt(d2);
    const Vec3 dir    = disp / dist;

    const double cos_emit = dir.dot(refl_inward_normal);
    const double cos_inc  = (-dir).dot(pad_inward_normal);
    if (cos_emit <= 0.0 || cos_inc <= 0.0) return 0.0;

    // Exact Lambertian view factor dF = (cos_emit * cos_inc * pi * r^2) / (pi * d^2)
    const double view_prob = std::min(0.25, (cos_emit * cos_inc * pad_radius_m * pad_radius_m) / d2);
    out_flight_time_s = dist / std::max(v_group_ms, 1.0);
    return phonon_energy_eV * track_weight * view_prob;
}

} // namespace cmpmc::experimental
