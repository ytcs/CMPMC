// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

// Material description: plain data structs for a crystal (CrystalConfig) and a
// superconducting sensor film (SuperconductorFilmConfig). All quantities are in
// SI units unless noted (energies in eV; see cmpmc/core/units.hpp).

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/units.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace cmpmc {

// Crystal symmetry class; selects which Voigt stiffness components are independent.
enum class BravaisGroup : uint8_t {
    Amorphous    = 0,
    Cubic        = 1,
    Tetragonal   = 2,
    Orthorhombic = 3,
    Hexagonal    = 4,
    Rhombohedral = 5,
    Monoclinic   = 6,
    Triclinic    = 7
};

// Intervalley scattering rate model for conduction electrons.
enum class IVRateModel : uint8_t {
    Quadratic = 0, // EDELWEISS parametrisation: rate * (E0^2 + E^2)^(power/2)
    Linear    = 1, // Power-law parametrization: rate0 + rate1 * E^power
    Matrix    = 2  // First-principles deformation potential matrix-element model
};

// CrystalConfig: physical parameters of a crystal. Covers all eight symmetry
// classes of BravaisGroup, Miller (hkl) orientation, grain-boundary scattering,
// conduction-band valleys, charge trapping, and bulk superconductor parameters.
struct CrystalConfig {
    std::string name     = "Ge";
    BravaisGroup group   = BravaisGroup::Cubic;
    double atomic_Z      = 32.0;                 // Atomic number Z (for NIEL)
    double atomic_A      = 72.63;                // Atomic mass A in g/mol (for NIEL)
    double density       = 5323.0;               // Mass density rho (kg / m^3)

    // --- 2nd-Order Elastic Stiffness Tensor (Pa) ---
    // For cubic crystals, simply set C11, C12, C44.
    // For non-cubic crystals (Rhombohedral Al2O3/SiO2, Tetragonal CaWO4, etc.),
    // set extra components (C13, C14, C16, C33, C66) or call `set_stiffness(i, j, val)`.
    double C11           = 126.0 * units::GPa;
    double C12           =  44.0 * units::GPa;
    double C13           =   0.0;
    double C14           =   0.0;
    double C16           =   0.0;
    double C33           =   0.0;
    double C44           =  67.0 * units::GPa;
    double C66           =   0.0;
    std::array<std::array<double, 6>, 6> C_custom{}; // Optional full 6x6 override
    bool has_custom_C    = false;

    // --- 3rd-Order Anharmonic & Isotope Parameters (S. Tamura, PRB 31, 1985) ---
    double beta          = -73.2 * units::GPa;
    double gamma         = -70.8 * units::GPa;
    double lambda        =  37.6 * units::GPa;
    double mu            =  56.1 * units::GPa;

    double B_iso         = 3.67e-41;             // s^3 (Isotope scattering: Gamma = B_iso * nu^4)
    double A_anh         = 6.43e-55;             // s^4 (L-phonon decay: Gamma = A_anh * nu^5)
    double fTT_frac      = 0.74;                 // Branching fraction for L -> T + T (vs L -> L + T)
    double polycrystal_mfp_m = 0.0;              // Grain-boundary elastic MFP (m); <= 0 disables

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
    double epsilon_rel   = 16.2;                 // Relative permittivity (dielectric constant)
    double neutral_impurity_cm3 = 1.0e11;        // Neutral impurity density (1/cm^3)
    double alpha_nonparabolic   = 0.3;           // Conduction band non-parabolicity (1/eV)
    double h_mass_ratio  = 0.350;                // Effective hole mass m_h / m_e
    Vec3 e_mass_tensor   = {1.588, 0.081, 0.081};// Valley electron mass tensor (m_L, m_T, m_T) / m_e
    std::vector<Vec3> valley_dirs = {            // Conduction band valley axes in crystal frame
        { 1,  1,  1},
        {-1,  1,  1},
        {-1, -1,  1},
        { 1, -1,  1}
    };

    // --- Acoustic & Intervalley Deformation Potentials ---
    double ac_deform_e   = 11.0 * units::eV;     // Electron acoustic deformation potential (eV)
    double ac_deform_h   =  3.4 * units::eV;     // Hole acoustic deformation potential (eV)
    std::vector<double> iv_deform_eV_m = {3.0e10, 0.2e10}; // Optical IV deformation potentials (eV/m)
    std::vector<double> iv_energy_eV   = {27.3 * units::meV, 10.3 * units::meV}; // IV phonon thresholds
    double iv_min_energy = 10.3 * units::meV;    // Minimum kinetic energy for IV scattering
    IVRateModel iv_model = IVRateModel::Quadratic;
    double iv_quad_rate  = 6.72e-2;              // Hz
    double iv_quad_field = 217.0;                // V/m
    double iv_quad_power = 3.24;
    double iv_lin_rate0  = 1.5e6;                // Hz
    double iv_lin_rate1  = 1.5;                  // Hz / (V/cm)^power
    double iv_lin_power  = 4.0;

    // --- Charge Trapping & Impact Ionization Parameters ---
    double e_trapping_mfp_m     = 0.0;           // Electron trapping MFP (m); <= 0 disables
    double h_trapping_mfp_m     = 0.0;           // Hole trapping MFP (m); <= 0 disables
    double e_trap_ion_d0_mfp_m  = 0.0;           // e- on D0 trap ionization MFP (m)
    double e_trap_ion_a_mfp_m   = 0.0;           // e- on A- trap ionization MFP (m)
    double h_trap_ion_d_mfp_m   = 0.0;           // h+ on D+ trap ionization MFP (m)
    double h_trap_ion_a0_mfp_m  = 0.0;           // h+ on A0 trap ionization MFP (m)
    double donor_trap_energy_eV = 10.0 * units::meV; // Shallow donor binding energy (eV)
    double accep_trap_energy_eV = 10.0 * units::meV; // Shallow acceptor binding energy (eV)

    // --- 3D Bulk Superconductor Parameters (for Al, Nb, Ta crystals/films) ---
    double sc_delta0_eV         = 0.0;           // 0K BCS bandgap Delta_0 (eV)
    double sc_tau0_qp_s         = 0.0;           // Characteristic QP phonon-emission lifetime (s)
    double sc_tau0_ph_s         = 0.0;           // Characteristic 2*Delta phonon pair-breaking time (s)
    double sc_dn_m2_s           = 0.0;           // Normal-state electron diffusivity D_n (m^2/s)
    double sc_teff_K            = 0.0;           // Effective quasiparticle temperature T_eff (K)
    double sc_trap_tau_s        = 0.0;           // Local QP trapping time (s)

    // --- Crystal Lattice to Detector Coordinate Rotation ---
    Mat3 orient          = Mat3::identity();

    // Set 1-based Voigt stiffness component C_ij (1..6, 1..6)
    constexpr void set_stiffness(int i_1based, int j_1based, double val_Pa) noexcept {
        if (i_1based < 1 || i_1based > 6 || j_1based < 1 || j_1based > 6) return;
        has_custom_C = true;
        C_custom[i_1based - 1][j_1based - 1] = val_Pa;
        C_custom[j_1based - 1][i_1based - 1] = val_Pa;
    }

    // Align Miller normal vector (h, k, l) with the detector +Z axis
    inline void set_miller_orientation(int h, int k, int l, double rot_rad = 0.0) noexcept {
        const Vec3 norm = Vec3{static_cast<double>(h), static_cast<double>(k), static_cast<double>(l)}.unit();
        const double theta = std::acos(std::clamp(norm.z, -1.0, 1.0));
        const double phi   = std::atan2(norm.y, norm.x);

        const double ct = std::cos(theta), st = std::sin(theta);
        const double cp = std::cos(phi),   sp = std::sin(phi);
        const double cr = std::cos(rot_rad), sr = std::sin(rot_rad);

        // R = R_z(phi) * R_y(theta) * R_z(rot) maps crystal frame to detector frame
        orient = Mat3{{
            {cp * ct * cr - sp * sr, -cp * ct * sr - sp * cr, cp * st},
            {sp * ct * cr + cp * sr, -sp * ct * sr + cp * cr, sp * st},
            {-st * cr,                st * sr,                ct}
        }};
    }

    // Build the full symmetric 6x6 Voigt matrix from the independent components
    // allowed by the symmetry class (standard reduced-stiffness relations).
    inline std::array<std::array<double, 6>, 6> build_voigt_matrix() const noexcept {
        std::array<std::array<double, 6>, 6> V = C_custom;
        if (!has_custom_C) {
            V[0][0] = C11;
            V[0][1] = C12;
            V[0][2] = C13;
            V[0][3] = C14;
            V[0][5] = C16;
            V[2][2] = C33;
            V[3][3] = C44;
            V[5][5] = C66;
        }

        switch (group) {
            case BravaisGroup::Amorphous:
                V[3][3] = 0.5 * (V[0][0] - V[0][1]);
                [[fallthrough]];
            case BravaisGroup::Cubic: {
                const double c11 = V[0][0], c12 = V[0][1], c44 = V[3][3];
                for (int i = 0; i < 6; ++i) {
                    for (int j = i; j < 6; ++j) {
                        if (i < 3 && j < 3) V[i][j] = (i == j) ? c11 : c12;
                        else if (i == j && i >= 3) V[i][i] = c44;
                        else V[i][j] = 0.0;
                    }
                }
                break;
            }
            case BravaisGroup::Hexagonal:
                V[0][5] = 0.0;
                V[5][5] = 0.5 * (V[0][0] - V[0][1]);
                [[fallthrough]];
            case BravaisGroup::Tetragonal: {
                V[1][1] = V[0][0];
                V[1][2] = V[0][2];
                V[1][5] = -V[0][5];
                V[4][4] = V[3][3];
                break;
            }
            case BravaisGroup::Rhombohedral: {
                const double c11 = V[0][0], c12 = V[0][1], c13 = V[0][2], c14 = V[0][3];
                const double c44 = V[3][3], c66 = 0.5 * (c11 - c12);
                V[0][4] = 0.0; V[0][5] = 0.0;
                V[1][1] = c11; V[1][2] = c13; V[1][3] = -c14; V[1][4] = 0.0; V[1][5] = 0.0;
                V[2][3] = 0.0; V[2][4] = 0.0; V[2][5] = 0.0;
                V[3][4] = 0.0; V[3][5] = 0.0;
                V[4][4] = c44; V[4][5] = c14;
                V[5][5] = c66;
                break;
            }
            case BravaisGroup::Orthorhombic:
            case BravaisGroup::Monoclinic:
            case BravaisGroup::Triclinic:
                break;
        }

        for (int i = 0; i < 6; ++i)
            for (int j = i + 1; j < 6; ++j)
                V[j][i] = V[i][j];

        return V;
    }
};

// Film phonon-loss model of the Kaplan cascade (KaplanFilm):
//   MeanFreePath:          escape probability P_esc = exp(-path/lambda) with a linear lifetime law
//                          lambda = v_s tau0 / (1 + slope (Omega/D - 2)) (arXiv:2302.05998,
//                          Eqs. (17)-(18)). Default.
//   KaplanRateCompetition: a film phonon of energy Omega >= 2D breaks a pair with probability
//                          1/(1 + tau_pb(Omega)/tau_l), otherwise it is lost to the substrate, with the exact
//                          T = 0 Kaplan rate 1/tau_pb(Omega) = I(Omega)/(pi D tau0) and the loss time tau_l
//                          (arXiv:1401.1937, Eqs. (2)-(5); tau_l/tau_0^phi is the film-loss parameter).
enum class FilmLossModel : uint8_t { MeanFreePath = 0, KaplanRateCompetition = 1 };

// SuperconductorFilmConfig: parameters of the thin-film quasiparticle cascade
// (S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976)).
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
    FilmLossModel loss_model     = FilmLossModel::MeanFreePath;
    double phonon_loss_time      = 0.0;                  // tau_l [s] for KaplanRateCompetition (<= 0: no loss)
};

} // namespace cmpmc
