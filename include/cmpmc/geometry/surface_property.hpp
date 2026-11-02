#pragma once

#include "cmpmc/core/core.hpp"
#include "cmpmc/geometry/solids.hpp"
#include "cmpmc/geometry/electrode_pattern.hpp"
#include "cmpmc/superconductor/kaplan_film.hpp"
#include "cmpmc/phonon/interface_transmission.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace cmpmc {

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
    // Mode- and angle-resolved substrate -> film transmission (isotropic AMM). When set, a phonon
    // reaching an active pad enters the film with probability t_AMM(mode, cos theta_k), which
    // replaces KaplanFilm::film_absorption (campaign WP8-R Sec. 3.2; IF-AUD-4).
    std::optional<physics::IsotropicAMMInterface> amm_interface = std::nullopt;
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

    // Frequency-independent specular probability used where the polynomials do not
    // apply (G4CMP MaterialPropertiesTable "specProb"; G4CMPSurfaceProperty ctor pSpecProb)
    double spec_prob_const = 0.0;
    // true for the G4CMP example coefficients ("not well-motivated" demonstration values,
    // R15 CMPMC-AUD IF-AUD-4); the engine refuses them unless allow_unvalidated is set.
    bool demo_values = false;

    struct ReflectionSplit { double p_spec, p_diff, p_anh; };

    // Phonon reflection branching exactly as G4CMPSurfaceProperty::{Specular,Diffuse,
    // Anharmonic}ReflProb + G4CMPPhononBoundaryProcess::DoReflection (normalised):
    //   f > anh_cutoff:                p_anh = 0, p_diff = 1 - specProb, p_spec = specProb
    //   refl_cutoff < f <= anh_cutoff: p_anh = P_anh(f), p_diff = P_diff(refl_cutoff),
    //                                  p_spec = 1 - p_diff - p_anh
    //   f <= refl_cutoff:              polynomials P_spec, P_diff, P_anh at f
    inline ReflectionSplit reflection_split(double energy_eV) const noexcept {
        const double f = (energy_eV / units::h_Planck) / units::GHz;
        const double anh = (anh_coeffs.empty() || f > anh_cutoff) ? 0.0 : eval_poly(anh_coeffs, f);
        double diff;
        if (diff_coeffs.empty() || f > anh_cutoff) diff = 1.0 - spec_prob_const;
        else diff = eval_poly(diff_coeffs, std::min(f, refl_cutoff));
        double spec;
        if (spec_coeffs.empty() || f > refl_cutoff) spec = 1.0 - diff - anh;
        else spec = eval_poly(spec_coeffs, f);
        spec = std::max(0.0, spec); diff = std::max(0.0, diff);
        const double a = std::max(0.0, anh);
        const double norm = spec + diff + a;
        if (norm <= 0.0) return {0.0, 1.0, 0.0};
        return {spec / norm, diff / norm, a / norm};
    }

    inline double specular_prob(double energy_eV) const noexcept {
        return reflection_split(energy_eV).p_spec;
    }

    // ---- Campaign defaults (no demo constants) ----------------------------------------------
    // Fully diffuse (Lambertian) reflection and no surface down-conversion: "Generally diffusive
    // scattering has been found to be a good model for phonon-surface reflections" (Leman
    // 1109.1193 Sec. IV.G [FT]). Surface down-conversion has no sourced coefficient (R15
    // CMPMC-AUD IF-AUD-4), so p_anh = 0; a campaign run that needs it, or a specular fraction,
    // must supply sourced spec/diff/anh coefficients (and amm_interface for film entry).
    static SurfaceProperty diffuse_al_sensor_surface() {
        SurfaceProperty s = diffuse_bare_surface();
        s.prob_abs      = 1.0;
        s.kaplan_sensor = KaplanFilm::aluminum_sensor();
        return s;
    }
    static SurfaceProperty diffuse_bare_surface() {
        SurfaceProperty s;
        s.prob_abs    = 0.0;
        s.prob_refl   = 1.0;
        s.spec_coeffs = {0.0};
        s.diff_coeffs = {1.0};
        s.anh_coeffs  = {};
        s.spec_prob_const = 0.0;
        s.kaplan_sensor = std::nullopt;
        return s;
    }

    // ---- G4CMP examples/phonon DEMONSTRATION values (guarded: demo_values = true) ------------
    // G4CMP's own documentation calls these polynomials "not well-motivated"; they reproduce the
    // G4CMP example (R15 CMPMC-AUD T9) and must not feed a campaign number.
    static SurfaceProperty demo_polynomial_sensor_surface() {
        SurfaceProperty s = demo_polynomial_bare_surface();
        s.prob_abs      = 1.0;
        s.kaplan_sensor = KaplanFilm::aluminum_sensor();
        return s;
    }
    static SurfaceProperty demo_polynomial_bare_surface() {
        SurfaceProperty s;
        s.prob_abs      = 0.0;
        s.prob_refl     = 1.0;
        s.anh_cutoff    = 520.0;
        s.refl_cutoff   = 350.0;
        s.spec_coeffs   = {0.928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};
        s.diff_coeffs   = {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
        s.anh_coeffs    = {0.0, 0.0, 0.0, 0.0, 0.0, 1.51e-14};
        s.kaplan_sensor = std::nullopt;
        s.demo_values   = true;
        return s;
    }
};

} // namespace cmpmc
