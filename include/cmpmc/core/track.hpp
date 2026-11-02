// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include "cmpmc/core/math.hpp"
#include <cstdint>

namespace cmpmc {

// ============================================================================
// Condensed-Matter Excitation Types & Cache-Aligned Track States
// ============================================================================
enum class ExcitationType : uint8_t {
    Phonon         = 0, // Acoustic / optical lattice vibrations (L, FT, ST)
    Electron       = 1, // Conduction-band drift electron (with Herring-Vogt valley)
    Hole           = 2, // Valence-band drift hole
    BogoliubovQP   = 3, // Superconducting Bogoliubov quasiparticle
    Exciton        = 4, // Bound electron-hole pair (singlet/triplet)
    Magnon         = 5, // Spin-wave excitation (magnetic calorimeters / YIG)
    Plasmon        = 6, // Collective charge oscillation
    Custom         = 7  // User-defined excitation species
};

inline constexpr const char* excitation_name(ExcitationType t) noexcept {
    switch (t) {
        case ExcitationType::Phonon:       return "Phonon";
        case ExcitationType::Electron:     return "Electron";
        case ExcitationType::Hole:         return "Hole";
        case ExcitationType::BogoliubovQP: return "BogoliubovQP";
        case ExcitationType::Exciton:      return "Exciton";
        case ExcitationType::Magnon:       return "Magnon";
        case ExcitationType::Plasmon:      return "Plasmon";
        default:                           return "Custom";
    }
}

enum class PhononMode : uint8_t {
    L     = 0, // Longitudinal
    FT    = 1, // Fast Transverse
    ST    = 2, // Slow Transverse
    Count = 3
};

inline constexpr const char* mode_name(PhononMode m) noexcept {
    switch (m) {
        case PhononMode::L:  return "L";
        case PhononMode::FT: return "FT";
        case PhononMode::ST: return "ST";
        default:             return "?";
    }
}

enum class SurfaceId : uint8_t {
    TopFace    = 0,
    BottomFace = 1,
    SideWall   = 2,
    Count      = 3
};

// 64-byte Cache-Aligned Phonon Track (natively stores both wavevector k_dir and group velocity vg_dir)
struct alignas(64) PhononTrack {
    Vec3 pos;                        // 24B: Position in detector frame (m)
    Vec3 k_dir;                      // 24B: Wavevector unit direction
    Vec3 vg_dir;                     // 24B: Group velocity unit direction
    double energy   = 0.0;           //  8B: Phonon energy (eV)
    double vg_mag   = 0.0;           //  8B: Group velocity magnitude (m/s)
    double time     = 0.0;           //  8B: Global track time (s)
    double weight   = 1.0;           //  8B: Statistical weight
    uint32_t bounces= 0;             //  4B: Boundary reflection counter
    PhononMode mode = PhononMode::L; //  1B: Polarization mode (L, FT, ST)
};

// 64-byte Cache-Aligned Charge Carrier Track (Electron with valley index or Hole)
struct alignas(64) ChargeTrack {
    Vec3 pos;                                       // 24B: Position in detector frame (m)
    Vec3 p_eV;                                      // 24B: Crystal momentum hbar*c*k (eV)
    double ekin         = 0.0;                      //  8B: Kinetic energy (eV)
    double time         = 0.0;                      //  8B: Global time (s)
    double weight       = 1.0;                      //  8B: Statistical weight
    uint16_t valley     = 0;                        //  2B: Conduction-band valley index
    ExcitationType type = ExcitationType::Electron; //  1B: Electron or Hole
    bool alive          = true;                     //  1B: Active status
};

// 64-byte Cache-Aligned Generic Excitation Track (for BogoliubovQP, Exciton, Magnon, etc.)
struct alignas(64) GenericExcitationTrack {
    Vec3 pos;
    Vec3 mom_dir;
    double energy        = 0.0;
    double velocity      = 0.0;
    double time          = 0.0;
    double weight        = 1.0;
    uint32_t bounces     = 0;
    ExcitationType type  = ExcitationType::BogoliubovQP;
    uint8_t spin_or_band = 0;
};

struct SensorHit {
    Vec3 pos;
    double time;
    double edep;             // Weighted energy deposited (eV)
    double incident_energy;
    SurfaceId surface_id;
    ExcitationType excitation = ExcitationType::Phonon;
    PhononMode mode           = PhononMode::L;
    int32_t channel_id        = -1; // Patterned electrode channel ID (-1 if unpatterned)
};

} // namespace cmpmc
