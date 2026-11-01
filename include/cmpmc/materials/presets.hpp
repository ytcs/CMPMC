#pragma once

#include "cmpmc/materials/crystal_config.hpp"

namespace cmpmc::materials {

// ============================================================================
// Standard Material Datasheets
// Add new crystals or superconducting films here by returning a filled config.
// ============================================================================

// Germanium (Ge) - Matches G4CMP CrystalMaps/Ge/config.txt
// References: S. Tamura, PRB 31 (1985); Cabrera et al., arXiv:1004.1233 (2010)
inline CrystalConfig Germanium() {
    return CrystalConfig{}; // Default values in CrystalConfig are Ge
}

// Silicon (Si) - Matches G4CMP CrystalMaps/Si/config.txt
// References: S. Tamura, PRB 31 (1985); Jacoboni & Reggiani, Rev. Mod. Phys. 55 (1983)
inline CrystalConfig Silicon() {
    CrystalConfig c;
    c.name          = "Si";
    c.atomic_Z      = 14.0;
    c.atomic_A      = 28.085;
    c.density       = 2329.0; // kg/m^3
    c.C11           = 165.6 * units::GPa;
    c.C12           =  63.9 * units::GPa;
    c.C44           =  79.5 * units::GPa;
    c.beta          = -42.9 * units::GPa;
    c.gamma         = -94.5 * units::GPa;
    c.lambda        =  52.4 * units::GPa;
    c.mu            =  68.0 * units::GPa;
    c.B_iso         = 2.43e-42;
    c.A_anh         = 7.41e-56;
    c.fTT_frac      = 0.74;
    c.fLDOS         = 0.093;
    c.fSTDOS        = 0.531;
    c.fFTDOS        = 0.376;
    c.debye_freq    = 15.0 * units::THz;
    c.v_sound       = 9000.0;
    c.v_trans       = 5400.0;
    c.bandgap       = 1.17 * units::eV;
    c.pair_energy   = 3.81 * units::eV;
    c.fano_factor   = 0.15;
    c.h_mass_ratio  = 0.50;
    c.e_mass_tensor = {0.91, 0.19, 0.19};
    c.valley_dirs   = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    c.ac_deform_e   = 9.0 * units::eV;
    c.ac_deform_h   = 5.0 * units::eV;
    c.iv_min_energy = 12.0 * units::meV;
    c.iv_model      = IVRateModel::Linear;
    c.iv_lin_rate0  = 1.5e6;
    c.iv_lin_rate1  = 1.5;
    c.iv_lin_power  = 4.0;
    return c;
}

// Gallium Arsenide (GaAs) - Matches G4CMP CrystalMaps/GaAs/config.txt
inline CrystalConfig GalliumArsenide() {
    CrystalConfig c = Germanium();
    c.name          = "GaAs";
    c.atomic_Z      = 31.5;
    c.atomic_A      = 72.32;
    c.density       = 5317.0;
    c.C11           = 122.1 * units::GPa;
    c.C12           =  56.6 * units::GPa;
    c.C44           =  60.0 * units::GPa;
    c.B_iso         = 7.30e-42;
    c.A_anh         = 2.66e-55;
    c.fTT_frac      = 0.58;
    c.bandgap       = 1.519 * units::eV;
    c.pair_energy   = 4.20 * units::eV;
    c.fano_factor   = 0.12;
    c.e_mass_tensor = {0.067, 0.067, 0.067}; // Isotropic Gamma-valley
    c.valley_dirs   = {{1, 0, 0}};
    return c;
}

// Standard 600 nm Aluminum QET / KID Superconducting Film (S.B. Kaplan et al., PRB 14, 1976)
inline constexpr SuperconductorFilmConfig AluminumFilm() noexcept {
    return SuperconductorFilmConfig{};
}

// Thin Tungsten (W) Transition-Edge Sensor (TES) Film
inline constexpr SuperconductorFilmConfig TungstenTES() noexcept {
    SuperconductorFilmConfig f{};
    f.film_absorption       = 1.0;
    f.film_thickness        = 40.0 * units::nm;
    f.v_sound               = 2890.0;
    f.gap_energy            = 15.0 * units::ueV;
    f.phonon_lifetime       = 100.0 * units::ps;
    f.phonon_lifetime_slope = 0.29;
    f.subgap_absorption     = 0.50;
    return f;
}

} // namespace cmpmc::materials
