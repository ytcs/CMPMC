#pragma once

#include <numbers>

namespace cmpmc::units {

// ============================================================================
// Explicit Physical Units & Fundamental Constants
//
// Internal Base Units:
//   Length   : meter (m)
//   Time     : second (s)
//   Energy   : electron-volt (eV)
//   Mass     : kilogram (kg) for density/elasticity; eV/c^2 for particle kinematics
//   Voltage  : Volt (V)
//
// Usage in physics/material files:
//   double C11       = 126.0 * units::GPa;
//   double gap       = 173.7 * units::ueV;
//   double lifetime  = 242.0 * units::ps;
// ============================================================================

inline constexpr double pi          = std::numbers::pi;
inline constexpr double two_pi      = 2.0 * std::numbers::pi;
inline constexpr double half_pi     = 0.5 * std::numbers::pi;

// --- Length (meters) ---
inline constexpr double m           = 1.0;
inline constexpr double cm          = 1e-2;
inline constexpr double mm          = 1e-3;
inline constexpr double um          = 1e-6;
inline constexpr double nm          = 1e-9;
inline constexpr double Ang         = 1e-10;

// --- Time (seconds) ---
inline constexpr double s           = 1.0;
inline constexpr double ms          = 1e-3;
inline constexpr double us          = 1e-6;
inline constexpr double ns          = 1e-9;
inline constexpr double ps          = 1e-12;

// --- Frequency (Hz) ---
inline constexpr double Hz          = 1.0;
inline constexpr double MHz         = 1e6;
inline constexpr double GHz         = 1e9;
inline constexpr double THz         = 1e12;

// --- Energy (eV) & Potential (Volt) ---
inline constexpr double eV          = 1.0;
inline constexpr double keV         = 1e3;
inline constexpr double MeV         = 1e6;
inline constexpr double meV         = 1e-3;
inline constexpr double ueV         = 1e-6;
inline constexpr double Volt        = 1.0;

// --- Pressure / Elastic Stiffness (Pa = kg / (m * s^2)) ---
inline constexpr double Pa          = 1.0;
inline constexpr double GPa         = 1e9;

// --- Fundamental Physical Constants ---
inline constexpr double c_light     = 299792458.0;          // Speed of light (m/s)
inline constexpr double c_squared   = c_light * c_light;    // c^2 (m^2/s^2)
inline constexpr double h_Planck    = 4.135667696e-15;      // Planck constant h (eV*s)
inline constexpr double hbar_Planck = 6.582119569e-16;      // Reduced Planck hbar (eV*s)
inline constexpr double hbarc       = hbar_Planck * c_light;// hbar * c (eV*m)
inline constexpr double k_Boltzmann = 8.617333262e-5;       // Boltzmann constant k_B (eV/K)
inline constexpr double e_SI        = 1.602176634e-19;      // Elementary charge (C or J/eV)
inline constexpr double m_e_eV      = 510998.950;           // Electron rest mass m_e*c^2 (eV)
inline constexpr double m_e_kg      = 9.1093837015e-31;     // Electron mass (kg)

} // namespace cmpmc::units
