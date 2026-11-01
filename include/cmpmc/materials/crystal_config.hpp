#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/units.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace cmpmc {

enum class IVRateModel : uint8_t {
    Quadratic = 0, // Edelweiss parametrization: rate * (E0^2 + E^2)^(power/2)
    Linear    = 1  // Power-law parametrization: rate0 + rate1 * E^power
};

// ============================================================================
// CrystalConfig: Pure Physical Datasheet for a Semiconductor / Insulator Crystal
//
// Physicists can define a new crystal simply by filling this struct (see
// `include/cmpmc/materials/presets.hpp` for Germanium and Silicon examples).
// ============================================================================
struct CrystalConfig {
    std::string name     = "Ge";
    double atomic_Z      = 32.0;                 // Atomic number Z (for Lindhard NIEL)
    double atomic_A      = 72.63;                // Atomic mass A in g/mol (for Lindhard NIEL)
    double density       = 5323.0;               // Mass density rho (kg / m^3)

    // --- 2nd-Order Elastic Stiffness Tensor (Voigt C11, C12, C44 in Pa) ---
    double C11           = 126.0 * units::GPa;
    double C12           =  44.0 * units::GPa;
    double C44           =  67.0 * units::GPa;

    // --- 3rd-Order Anharmonic & Isotope Parameters (S. Tamura, PRB 31, 1985) ---
    double beta          = -73.2 * units::GPa;
    double gamma         = -70.8 * units::GPa;
    double lambda        =  37.6 * units::GPa;
    double mu            =  56.1 * units::GPa;

    double B_iso         = 3.67e-41;             // s^3 (Isotope scattering: Gamma = B_iso * nu^4)
    double A_anh         = 6.43e-55;             // s^4 (L-phonon decay: Gamma = A_anh * nu^5)
    double fTT_frac      = 0.74;                 // Branching fraction for L -> T + T (vs L -> L + T)

    // --- Phonon Density of States Weights (S. Tamura et al., PRB 44, 1991) ---
    double fLDOS         = 0.097834;
    double fSTDOS        = 0.535390;
    double fFTDOS        = 0.366770;

    double debye_freq    = 2.0 * units::THz;     // Characteristic primary phonon frequency
    double v_sound       = 5324.2077;            // Isotropic longitudinal sound speed (m/s)
    double v_trans       = 3258.7879;            // Isotropic transverse sound speed (m/s)

    // --- Charge Carrier Band Structure & Ionization Parameters ---
    double bandgap       = 0.74 * units::eV;     // Bandgap E_g (eV)
    double pair_energy   = 2.96 * units::eV;     // Mean energy per e-/h+ pair epsilon_eh (eV)
    double fano_factor   = 0.1057;               // Fano factor F for ionization fluctuations
    double h_mass_ratio  = 0.350;                // Effective hole mass m_h / m_e
    Vec3 e_mass_tensor   = {1.588, 0.081, 0.081};// Valley electron mass tensor (m_L, m_T, m_T) / m_e
    std::vector<Vec3> valley_dirs = {            // Conduction band valley axes in crystal frame
        { 1,  1,  1},
        {-1,  1,  1},
        {-1, -1,  1},
        { 1, -1,  1}
    };

    // --- Acoustic Deformation Potentials & Intervalley Scattering ---
    double ac_deform_e   = 11.0 * units::eV;     // Electron acoustic deformation potential (eV)
    double ac_deform_h   =  3.4 * units::eV;     // Hole acoustic deformation potential (eV)
    double iv_min_energy = 10.3 * units::meV;    // Threshold kinetic energy for IV scattering
    IVRateModel iv_model = IVRateModel::Quadratic;
    double iv_quad_rate  = 6.72e-2;              // Hz
    double iv_quad_field = 217.0;                // V/m
    double iv_quad_power = 3.24;
    double iv_lin_rate0  = 1.5e6;                // Hz
    double iv_lin_rate1  = 1.5;                  // Hz / (V/cm)^power
    double iv_lin_power  = 4.0;

    // --- Crystal Lattice to Detector Coordinate Rotation ---
    Mat3 orient          = Mat3::identity();
};

// ============================================================================
// SuperconductorFilmConfig: Thin-Film Quasiparticle Cascade Parameters
// Reference: S.B. Kaplan et al., Phys. Rev. B 14 (1976)
// ============================================================================
struct SuperconductorFilmConfig {
    double film_absorption       = 0.20;                 // Geometric sensor coverage / absorption prob
    double film_thickness        = 600.0 * units::nm;    // Film thickness d
    double v_sound               = 3260.0;               // Speed of sound in film (m/s)
    double gap_energy            = 173.715 * units::ueV; // Superconducting gap Delta_0 (eV)
    double phonon_lifetime       = 242.0 * units::ps;    // Pair-breaking lifetime tau_0 at 2*Delta
    double phonon_lifetime_slope = 0.29;                 // Lifetime slope vs (E/Delta - 2)
    double low_qp_limit          = 3.0;                  // Min QP energy (in units of Delta) to radiate phonons
    double high_qp_limit         = 0.0;                  // Max QP energy split limit (0 = no split)
    double subgap_absorption     = 0.10;                 // Direct absorption probability below 2*Delta
    double absorber_gap          = 0.0;                  // Secondary TES trap gap (e.g. W = 15 ueV)
    double absorber_eff          = 1.0;                  // Quasiparticle collection efficiency
    double absorber_eff_slope    = 0.0;                  // Efficiency slope vs energy
    double temperature           = 0.05e-3;              // Film operating temperature (K)
    bool keep_all_phonons        = true;
};

} // namespace cmpmc
