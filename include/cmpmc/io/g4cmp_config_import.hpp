#pragma once

#include "cmpmc/materials/crystal_config.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cmpmc::io {

// Import of crystal parameters from the G4CMP `config.txt` lattice format
// (the files of a G4CMP `CrystalMaps/<material>/` directory) into a CrystalConfig.
//
// Format: one `key values... [unit]` entry per line; `#` starts a comment.
// Unknown keys are ignored. Supported keys:
//   amorphous | cubic | tetragonal | orthorhombic | hexagonal | rhombohedral |
//   monoclinic | triclinic                 crystal symmetry group
//   stiffness i j value unit              elastic constant C_ij (Voigt indices)
//   dyn beta gamma lambda mu unit         third-order (anharmonic) constants
//   scat B                                isotope scattering constant (s^3)
//   decay A                               anharmonic decay constant (s^4)
//   decayTT f                             L -> TT branching fraction
//   LDOS | STDOS | FTDOS f                mode density-of-states fractions
//   Debye value unit                      Debye frequency
//   vsound | vtrans value unit            longitudinal / transverse sound speed
//   bandgap | pairEnergy value unit       band gap, electron-hole pair energy
//   fanoFactor F                          Fano factor
//   hmass m                               hole mass ratio m_h / m_e
//   emass m1 m2 m3                        electron mass tensor (valley frame, / m_e)
//   valleyDir x y z                       conduction-valley axis (repeatable)
//   acDeform_e | acDeform_h value unit    acoustic deformation potentials
//   sc_tau0_qp | sc_tau0_ph value unit    superconductor characteristic times
// Recognised units: GPa, Pa, THz, GHz, Hz, eV, meV, ueV, m/s, km/s, ns, ps, s,
// V/m, V/cm, eV/cm, /cm3. An unrecognised unit is treated as 1 (SI).
//
// The format carries no density, so it is passed in. If no valleyDir is given, a
// single valley along [100] is assumed. Throws std::runtime_error if the file
// cannot be opened.
inline CrystalConfig import_g4cmp_config(const std::string& filepath,
                                            double density_kg_m3 = 5323.0,
                                            const std::string& name = "Custom") {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        throw std::runtime_error("cmpmc::io::import_g4cmp_config: cannot open " + filepath);
    }

    CrystalConfig cfg{};
    cfg.name    = name;
    cfg.density = density_kg_m3;
    cfg.valley_dirs.clear();

    auto parse_unit = [](const std::string& u) noexcept -> double {
        if (u == "GPa")       return units::GPa;
        if (u == "Pa")        return units::Pa;
        if (u == "THz")       return units::THz;
        if (u == "GHz")       return units::GHz;
        if (u == "Hz")        return units::Hz;
        if (u == "eV")        return units::eV;
        if (u == "meV")       return units::meV;
        if (u == "ueV")       return units::ueV;
        if (u == "m/s")       return 1.0;
        if (u == "km/s")      return 1e3;
        if (u == "ns")        return units::ns;
        if (u == "ps")        return units::ps;
        if (u == "s")         return units::s;
        if (u == "V/m")       return 1.0;
        if (u == "V/cm")      return 100.0;
        if (u == "eV/cm")     return 100.0 * units::eV;
        if (u == "/cm3")      return 1.0;
        return 1.0;
    };

    std::string line;
    while (std::getline(in, line)) {
        const auto hash_pos = line.find('#');
        if (hash_pos != std::string::npos) line = line.substr(0, hash_pos);
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "amorphous")      cfg.group = BravaisGroup::Amorphous;
        else if (key == "cubic")     cfg.group = BravaisGroup::Cubic;
        else if (key == "tetragonal")cfg.group = BravaisGroup::Tetragonal;
        else if (key == "orthorhombic") cfg.group = BravaisGroup::Orthorhombic;
        else if (key == "hexagonal") cfg.group = BravaisGroup::Hexagonal;
        else if (key == "rhombohedral") cfg.group = BravaisGroup::Rhombohedral;
        else if (key == "monoclinic")cfg.group = BravaisGroup::Monoclinic;
        else if (key == "triclinic") cfg.group = BravaisGroup::Triclinic;
        else if (key == "stiffness") {
            int i = 0, j = 0;
            double val = 0.0;
            std::string u;
            iss >> i >> j >> val >> u;
            cfg.set_stiffness(i, j, val * parse_unit(u));
        } else if (key == "dyn") {
            double b = 0, g = 0, l = 0, m = 0;
            std::string u;
            iss >> b >> g >> l >> m >> u;
            const double scale = parse_unit(u);
            cfg.beta   = b * scale;
            cfg.gamma  = g * scale;
            cfg.lambda = l * scale;
            cfg.mu     = m * scale;
        } else if (key == "scat") {
            iss >> cfg.B_iso;
        } else if (key == "decay") {
            iss >> cfg.A_anh;
        } else if (key == "decayTT") {
            iss >> cfg.fTT_frac;
        } else if (key == "LDOS") {
            iss >> cfg.fLDOS;
        } else if (key == "STDOS") {
            iss >> cfg.fSTDOS;
        } else if (key == "FTDOS") {
            iss >> cfg.fFTDOS;
        } else if (key == "Debye") {
            double val = 0.0;
            std::string u;
            iss >> val >> u;
            cfg.debye_freq = val * parse_unit(u);
        } else if (key == "vsound") {
            double val = 0.0;
            std::string u;
            iss >> val >> u;
            cfg.v_sound = val * parse_unit(u);
        } else if (key == "vtrans") {
            double val = 0.0;
            std::string u;
            iss >> val >> u;
            cfg.v_trans = val * parse_unit(u);
        } else if (key == "bandgap") {
            double val = 0.0;
            std::string u;
            iss >> val >> u;
            cfg.bandgap = val * parse_unit(u);
        } else if (key == "pairEnergy") {
            double val = 0.0;
            std::string u;
            iss >> val >> u;
            cfg.pair_energy = val * parse_unit(u);
        } else if (key == "fanoFactor") {
            iss >> cfg.fano_factor;
        } else if (key == "hmass") {
            iss >> cfg.h_mass_ratio;
        } else if (key == "emass") {
            iss >> cfg.e_mass_tensor.x >> cfg.e_mass_tensor.y >> cfg.e_mass_tensor.z;
        } else if (key == "valleyDir") {
            Vec3 v{};
            iss >> v.x >> v.y >> v.z;
            cfg.valley_dirs.push_back(v);
        } else if (key == "acDeform_e") {
            double val = 0.0; std::string u;
            iss >> val >> u;
            cfg.ac_deform_e = val * parse_unit(u);
        } else if (key == "acDeform_h") {
            double val = 0.0; std::string u;
            iss >> val >> u;
            cfg.ac_deform_h = val * parse_unit(u);
        } else if (key == "sc_tau0_qp") {
            double val = 0.0; std::string u;
            iss >> val >> u;
            cfg.sc_tau0_qp_s = val * parse_unit(u);
        } else if (key == "sc_tau0_ph") {
            double val = 0.0; std::string u;
            iss >> val >> u;
            cfg.sc_tau0_ph_s = val * parse_unit(u);
        }
    }

    if (cfg.valley_dirs.empty()) {
        cfg.valley_dirs.push_back({1.0, 0.0, 0.0});
    }
    return cfg;
}

} // namespace cmpmc::io
