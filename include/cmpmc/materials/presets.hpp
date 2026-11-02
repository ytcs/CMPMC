#pragma once

// Built-in material presets: twelve crystals and two superconducting sensor films.
//
// Unless a comment states otherwise, phonon and charge constants are the values
// distributed with G4CMP (CrystalMaps/<material>/config.txt), so that results are
// directly comparable. Deliberate deviations are documented at the constant and
// cite their source.

#include "cmpmc/materials/crystal_config.hpp"

namespace cmpmc::materials {

// 1. Germanium (Ge) - Cubic
inline CrystalConfig Germanium() {
    CrystalConfig c{};
    // Primary (Debye-energy) phonon frequency. G4CMP CrystalMaps uses 2 THz, which is
    // below the Ge acoustic band edge; Leman, arXiv:1109.1193, Table I gives
    // nu_D = 8.64 THz (arXiv:2408.04732, Table 5: maximum optical frequency 9.02 THz).
    c.debye_freq = 8.64 * units::THz;
    return c;
}

// 2. Silicon (Si) - Cubic
inline CrystalConfig Silicon() {
    CrystalConfig c;
    c.name              = "Si";
    c.group             = BravaisGroup::Cubic;
    c.atomic_Z          = 14.0;
    c.atomic_A          = 28.085;
    c.density           = 2329.0;
    c.C11               = 165.6 * units::GPa;
    c.C12               =  63.9 * units::GPa;
    c.C44               =  79.5 * units::GPa;
    c.beta              = -42.9 * units::GPa;
    c.gamma             = -94.5 * units::GPa;
    c.lambda            =  52.4 * units::GPa;
    c.mu                =  68.0 * units::GPa;
    c.B_iso             = 2.43e-42;
    c.A_anh             = 7.41e-56;
    c.fTT_frac          = 0.74;
    c.fLDOS             = 0.093;
    c.fSTDOS            = 0.531;
    c.fFTDOS            = 0.376;
    c.debye_freq        = 15.0 * units::THz;
    c.v_sound           = 9000.0;
    c.v_trans           = 5400.0;
    c.bandgap           = 1.17 * units::eV;
    // eps_eh and F_eh from arXiv:2004.11499: 3.749 +- 0.001 eV and 0.119 +- 0.002
    // at 5.9 keV, 123 K (G4CMP CrystalMaps uses 3.81 eV and 0.15).
    c.pair_energy       = 3.749 * units::eV;
    c.fano_factor       = 0.119;
    c.epsilon_rel       = 11.68;
    c.alpha_nonparabolic= 0.5;
    c.h_mass_ratio      = 0.50;
    c.e_mass_tensor     = {0.91, 0.19, 0.19};
    c.valley_dirs       = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    c.ac_deform_e       = 9.0 * units::eV;
    c.ac_deform_h       = 5.0 * units::eV;
    c.iv_deform_eV_m    = {0.5e10, 0.8e10, 11.0e10, 0.3e10, 2.0e10, 2.0e10};
    c.iv_energy_eV      = {12.0e-3, 18.4e-3, 61.8e-3, 18.9e-3, 47.2e-3, 58.8e-3};
    c.iv_min_energy     = 12.0 * units::meV;
    c.iv_model          = IVRateModel::Linear;
    c.iv_lin_rate0      = 1.5e6;
    c.iv_lin_rate1      = 1.5;
    c.iv_lin_power      = 4.0;
    c.iv_quad_rate      = 3.5e-20;
    c.iv_quad_field     = 3395.0;
    c.iv_quad_power     = 7.47;
    return c;
}

// 3. Sapphire (Al2O3) - Rhombohedral
inline CrystalConfig Sapphire() {
    CrystalConfig c;
    c.name       = "Al2O3";
    c.group      = BravaisGroup::Rhombohedral;
    c.atomic_Z   = 10.0;
    c.atomic_A   = 20.39;
    c.density    = 3980.0;
    c.C11        = 495.0 * units::GPa;
    c.C12        = 161.0 * units::GPa;
    c.C13        = 111.1 * units::GPa;
    c.C14        = -22.7 * units::GPa;
    c.C33        = 499.9 * units::GPa;
    c.C44        = 147.7 * units::GPa;
    c.C66        = 166.7 * units::GPa;
    c.beta       = -27.02 * units::GPa;
    c.gamma      = -152.8 * units::GPa;
    c.lambda     =  139.8 * units::GPa;
    c.mu         = 166.24 * units::GPa;
    c.B_iso      = 0.025e-42;
    c.A_anh      = 12.7e-55;
    c.fTT_frac   = 0.67;
    c.fLDOS      = 0.137;
    c.fSTDOS     = 0.511;
    c.fFTDOS     = 0.351;
    c.debye_freq = 8.41 * units::THz;
    c.v_sound    = 11200.0;
    c.v_trans    = 6040.0;
    c.bandgap    = 8.8 * units::eV;
    c.pair_energy= 25.0 * units::eV;
    return c;
}

// 4. Calcium Tungstate (CaWO4) - Tetragonal
inline CrystalConfig CalciumTungstate() {
    CrystalConfig c;
    c.name       = "CaWO4";
    c.group      = BravaisGroup::Tetragonal;
    c.atomic_Z   = 32.0;
    c.atomic_A   = 72.0;
    c.density    = 6060.0;
    c.C11        = 151.5 * units::GPa;
    c.C12        =  65.6 * units::GPa;
    c.C13        =  45.0 * units::GPa;
    c.C16        =  18.8 * units::GPa;
    c.C33        = 134.0 * units::GPa;
    c.C44        =  35.4 * units::GPa;
    c.C66        =  40.0 * units::GPa;
    c.beta       =   0.0;
    c.gamma      =   0.0;
    c.lambda     =  57.94 * units::GPa;
    c.mu         =  40.78 * units::GPa;
    c.B_iso      = 15.0e-42;
    c.A_anh      = 7.9e-56;
    c.fTT_frac   = 0.74;
    c.fLDOS      = 0.115;
    c.fSTDOS     = 0.475;
    c.fFTDOS     = 0.408;
    c.debye_freq = 2.72 * units::THz;
    c.v_sound    = 4598.0;
    c.v_trans    = 2569.0;
    return c;
}

// 5. Alpha-Quartz (SiO2_alpha) - Rhombohedral
inline CrystalConfig AlphaQuartz() {
    CrystalConfig c;
    c.name       = "SiO2_alpha";
    c.group      = BravaisGroup::Rhombohedral;
    c.atomic_Z   = 10.0;
    c.atomic_A   = 20.03;
    c.density    = 2648.0;
    c.C11        = 86.103 * units::GPa;
    c.C12        = 10.430 * units::GPa;
    c.C13        = 12.894 * units::GPa;
    c.C14        = -17.206 * units::GPa;
    c.C33        = 108.552 * units::GPa;
    c.C44        = 59.215 * units::GPa;
    c.C66        = 44.000 * units::GPa;
    c.beta       =  51.28 * units::GPa;
    c.gamma      = -138.93 * units::GPa;
    c.lambda     =  5.466 * units::GPa;
    c.mu         = 48.646 * units::GPa;
    c.B_iso      = 4.05e-42;
    c.A_anh      = 36.69e-55;
    c.fTT_frac   = 0.71;
    c.fLDOS      = 0.10001;
    c.fSTDOS     = 0.53333;
    c.fFTDOS     = 0.36666;
    c.debye_freq = 4.2 * units::THz;
    c.v_sound    = 5780.0;
    c.v_trans    = 3340.0;
    return c;
}

// 6. Amorphous Silica (SiO2_amorph)
inline CrystalConfig AmorphousSilica() {
    CrystalConfig c = AlphaQuartz();
    c.name       = "SiO2_amorph";
    c.group      = BravaisGroup::Amorphous;
    c.density    = 2200.0;
    c.C11        = 78.0 * units::GPa;
    c.C12        = 16.0 * units::GPa;
    c.beta       = 54.18 * units::GPa;
    c.gamma      = -13.14 * units::GPa;
    c.lambda     = 16.10 * units::GPa;
    c.mu         = 31.20 * units::GPa;
    c.B_iso      = 4.65e-42;
    c.A_anh      = 7.99e-55;
    c.fTT_frac   = 0.54;
    c.v_sound    = 5970.0;
    c.v_trans    = 3760.0;
    return c;
}

// 7. Gallium Arsenide (GaAs) - Cubic
inline CrystalConfig GalliumArsenide() {
    CrystalConfig c = Germanium();
    c.name          = "GaAs";
    c.group         = BravaisGroup::Cubic;
    c.atomic_Z      = 31.5;
    c.atomic_A      = 72.32;
    c.density       = 5317.0;
    c.C11           = 121.57 * units::GPa;
    c.C12           =  57.40 * units::GPa;
    c.C44           =  59.50 * units::GPa;
    c.beta          = -54.71 * units::GPa;
    c.gamma         = -67.51 * units::GPa;
    c.lambda        =  47.20 * units::GPa;
    c.mu            =  44.20 * units::GPa;
    // Phonon constants as in arXiv:2408.04732, Tables 3-5 (A = 7.77e-55 s^4,
    // F_TT = 0.77, DOS 0.584:0.337:0.078), which agree with G4CMP CrystalMaps.
    c.B_iso         = 7.22e-42;
    c.A_anh         = 7.778e-55;
    c.fTT_frac      = 0.778;
    c.fLDOS         = 0.078;
    c.fSTDOS        = 0.584;
    c.fFTDOS        = 0.337;
    c.debye_freq    = 8.58 * units::THz;   // maximum optical frequency, arXiv:2408.04732, Table 5 (G4CMP CrystalMaps: 6.56 THz, maximum acoustic)
    c.v_sound       = 5130.0;
    c.v_trans       = 3020.0;
    c.bandgap       = 1.519 * units::eV;
    c.pair_energy   = 4.20 * units::eV;
    c.fano_factor   = 0.12;
    c.e_mass_tensor = {0.067, 0.067, 0.067};
    c.valley_dirs   = {{1, 0, 0}};
    return c;
}

// 8. Calcium Fluoride (CaF2) - Cubic
inline CrystalConfig CalciumFluoride() {
    CrystalConfig c;
    c.name       = "CaF2";
    c.group      = BravaisGroup::Cubic;
    c.atomic_Z   = 12.67;
    c.atomic_A   = 26.03;
    c.density    = 3180.0;
    c.C11        = 174.0 * units::GPa;
    c.C12        =  56.0 * units::GPa;
    c.C44        = 33.93 * units::GPa;
    c.beta       = -95.85 * units::GPa;
    c.gamma      = -43.69 * units::GPa;
    c.lambda     =  65.20 * units::GPa;
    c.mu         =  45.15 * units::GPa;
    c.B_iso      = 4.83e-42;
    c.A_anh      = 6.14e-55;
    c.fTT_frac   = 0.75;
    c.fLDOS      = 0.0819;
    c.fSTDOS     = 0.5990;
    c.fFTDOS     = 0.3180;
    c.debye_freq = 9.8 * units::THz;
    c.v_sound    = 6500.0;
    c.v_trans    = 3980.0;
    return c;
}

// 9. Lithium Fluoride (LiF) - Cubic
inline CrystalConfig LithiumFluoride() {
    CrystalConfig c;
    c.name       = "LiF";
    c.group      = BravaisGroup::Cubic;
    c.atomic_Z   = 6.0;
    c.atomic_A   = 12.97;
    c.density    = 2635.0;
    c.C11        = 111.2 * units::GPa;
    c.C12        =  42.0 * units::GPa;
    c.C44        =  62.8 * units::GPa;
    c.beta       = -83.94 * units::GPa;
    c.gamma      = -87.54 * units::GPa;
    c.lambda     =  30.72 * units::GPa;
    c.mu         =  51.50 * units::GPa;
    c.B_iso      = 1.17e-42;
    c.A_anh      = 5.50e-55;
    c.fTT_frac   = 0.68;
    c.fLDOS      = 0.0983;
    c.fSTDOS     = 0.5156;
    c.fFTDOS     = 0.3860;
    c.debye_freq = 10.73 * units::THz;
    c.v_sound    = 7290.0;
    c.v_trans    = 4450.0;
    return c;
}

// 10. Aluminum Bulk Crystal / 3D Superconductor (Al)
inline CrystalConfig AluminumCrystal() {
    CrystalConfig c;
    c.name         = "Al";
    c.group        = BravaisGroup::Cubic;
    c.atomic_Z     = 13.0;
    c.atomic_A     = 26.98;
    c.density      = 2700.0;
    c.C11          = 106.85 * units::GPa;
    c.C12          =  60.74 * units::GPa;
    c.C44          =  28.22 * units::GPa;
    c.beta         = -114.08 * units::GPa;
    c.gamma        =  -75.48 * units::GPa;
    c.lambda       =   58.63 * units::GPa;
    c.mu           =   26.16 * units::GPa;
    c.B_iso        = 0.0;
    c.A_anh        = 4.21e-54;
    c.fTT_frac     = 0.79;
    c.fLDOS        = 0.055;
    c.fSTDOS       = 0.6204;
    c.fFTDOS       = 0.3339;
    c.debye_freq   = 8.89 * units::THz;
    c.v_sound      = 6420.0;
    c.v_trans      = 3040.0;
    c.sc_delta0_eV = 173.715 * units::ueV;
    c.sc_tau0_qp_s = 438.0 * units::ns;
    c.sc_tau0_ph_s = 0.242 * units::ns;
    c.sc_dn_m2_s   = 0.0060; // 60 cm^2/s typical Al film diffusivity
    c.sc_teff_K    = 0.05;
    return c;
}

// 11. Niobium Bulk Crystal / 3D Superconductor (Nb)
inline CrystalConfig NiobiumCrystal() {
    CrystalConfig c;
    c.name         = "Nb";
    c.group        = BravaisGroup::Cubic;
    c.atomic_Z     = 41.0;
    c.atomic_A     = 92.91;
    c.density      = 8570.0;
    c.C11          = 246.5 * units::GPa;
    c.C12          = 134.8 * units::GPa;
    c.C44          =  28.7 * units::GPa;
    c.beta         = -379.0 * units::GPa;
    c.gamma        =   78.4 * units::GPa;
    c.lambda       =  145.0 * units::GPa;
    c.mu           =   39.7 * units::GPa;
    c.B_iso        = 0.0;
    c.A_anh        = 5.77e-54;
    c.fTT_frac     = 0.90;
    c.fLDOS        = 0.055;
    c.fSTDOS       = 0.6204;
    c.fFTDOS       = 0.3300;
    c.debye_freq   = 15.0 * units::THz;
    c.v_sound      = 4920.0;
    c.v_trans      = 2100.0;
    c.sc_delta0_eV = 1.50 * units::meV;
    c.sc_tau0_qp_s = 0.149 * units::ns;
    c.sc_tau0_ph_s = 0.00417 * units::ns;
    c.sc_dn_m2_s   = 0.0002;
    c.sc_teff_K    = 0.05;
    return c;
}

// 12. Copper (Cu) - Amorphous
inline CrystalConfig Copper() {
    CrystalConfig c = AmorphousSilica();
    c.name       = "Cu";
    c.group      = BravaisGroup::Amorphous;
    c.atomic_Z   = 29.0;
    c.atomic_A   = 63.55;
    c.density    = 8960.0;
    c.C11        = 168.3 * units::GPa;
    c.C12        = 121.4 * units::GPa;
    c.v_sound    = 4760.0;
    c.v_trans    = 2325.0;
    return c;
}

// Thin-film superconducting sensors.
// Aluminium QET film with the default SuperconductorFilmConfig parameters.
inline constexpr SuperconductorFilmConfig AluminumFilm() noexcept {
    return SuperconductorFilmConfig{};
}

// Tungsten transition-edge-sensor film.
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
