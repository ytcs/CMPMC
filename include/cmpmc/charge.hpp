#pragma once

#include "cmpmc/core.hpp"
#include "cmpmc/lattice.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace cmpmc {

// ============================================================================
// ElectricField: Bias Voltage & 3D Field Evaluator
// Replaces: G4CMPFieldUtils & G4CMPMeshElectricField
// ============================================================================
struct UniformElectricField {
    Vec3 field_Vm = {0.0, 0.0, 0.0}; // Electric field E in V/m

    static constexpr UniformElectricField from_bias(double bias_volts, double thickness_m) noexcept {
        return UniformElectricField{{0.0, 0.0, -bias_volts / thickness_m}};
    }

    constexpr Vec3 at(const Vec3& /*pos*/) const noexcept {
        return field_Vm;
    }

    constexpr double potential(const Vec3& pos) const noexcept {
        return -field_Vm.dot(pos);
    }
};

struct LukePhononEmit {
    bool emitted = false;
    Vec3 q_dir;          // Emitted phonon wavevector direction in detector frame
    double phonon_eV;    // Emitted phonon energy (eV)
};

// ============================================================================
// ChargePhysics: Herring-Vogt Drift, Neganov-Trofimov-Luke Phonon Emission,
// and Intervalley Scattering
// Replaces: G4CMPVDriftProcess, G4CMPTimeStepper, G4CMPLukeScattering,
//           G4CMPLukeEmissionRate, G4CMPInterValleyScattering
// ============================================================================
class ChargePhysics {
public:
    explicit constexpr ChargePhysics(const CrystalLattice& lattice) noexcept
        : lattice_(lattice) {}

    inline ChargeTrack make_charge(ExcitationType type, const Vec3& pos, const Vec3& dir,
                                   double ekin_eV, uint16_t valley_idx, double time_s = 0.0) const noexcept {
        const Vec3 u = dir.unit();
        Vec3 p_eV{0.0, 0.0, 0.0};
        if (type == ExcitationType::Electron) {
            p_eV = lattice_.map_ekin_to_p_electron(valley_idx, u, ekin_eV);
        } else {
            const double m_h = lattice_.hole_mass_eV();
            const double p_mag = std::sqrt(std::max(0.0, ekin_eV * ekin_eV + 2.0 * ekin_eV * m_h));
            p_eV = u * p_mag;
        }
        return ChargeTrack{
            .pos    = pos,
            .p_eV   = p_eV,
            .ekin   = ekin_eV,
            .time   = time_s,
            .weight = 1.0,
            .valley = valley_idx,
            .type   = type,
            .alive  = true
        };
    }

    inline Vec3 velocity(const ChargeTrack& trk) const noexcept {
        if (trk.type == ExcitationType::Electron) {
            return lattice_.map_p_to_vel_electron(trk.valley, trk.p_eV);
        }
        return lattice_.map_p_to_vel_hole(trk.p_eV);
    }

    // Neganov-Trofimov-Luke acoustic phonon emission rate (Hz)
    // Reference: Cabrera et al., arXiv:1004.1233 (2010)
    //   Gamma_Luke(k) = (v_sound / (3 * l_0)) * (k_sound / k) * (k / k_sound - 1)^3
    inline double luke_emission_rate(const ChargeTrack& trk) const noexcept {
        const double v_sound = lattice_.config().v_sound;
        double k_mag   = 0.0;
        double k_sound = 0.0;
        double l0      = 0.0;

        if (trk.type == ExcitationType::Electron) {
            const Vec3 k_ell = lattice_.map_p_to_k_electron(trk.valley, trk.p_eV);
            const Vec3 k_sph = lattice_.ellipsoidal_to_spherical(trk.valley, k_ell);
            k_mag   = k_sph.norm();
            k_sound = lattice_.k_sound_electron();
            l0      = lattice_.l0_electron();
        } else {
            k_mag   = trk.p_eV.norm() / units::hbarc;
            k_sound = lattice_.k_sound_hole();
            l0      = lattice_.l0_hole();
        }

        if (k_mag <= k_sound) return 0.0;
        const double mach  = k_mag / k_sound;
        const double dm    = mach - 1.0;
        const double tstep = (3.0 * l0) / v_sound;
        const double tau   = tstep * mach / (dm * dm * dm);
        return (tau > 0.0) ? (1.0 / tau) : 0.0;
    }

    // Sample Neganov-Trofimov-Luke Cherenkov phonon emission and update recoil momentum in-place
    inline LukePhononEmit emit_luke_phonon(ChargeTrack& trk, Rng& rng) const noexcept {
        Vec3 k_trk{};
        double k_sound = 0.0;

        if (trk.type == ExcitationType::Electron) {
            const Vec3 k_ell = lattice_.map_p_to_k_electron(trk.valley, trk.p_eV);
            k_trk   = lattice_.ellipsoidal_to_spherical(trk.valley, k_ell);
            k_sound = lattice_.k_sound_electron();
        } else {
            k_trk   = trk.p_eV / units::hbarc;
            k_sound = lattice_.k_sound_hole();
        }

        const double k_mag = k_trk.norm();
        if (k_mag <= k_sound * 1.0001) return {};

        const Vec3 k_dir = k_trk / k_mag;
        const double v   = k_sound / k_mag;

        // Analytic inverse-CDF for Cherenkov angle: cos(theta) in (v, 1)
        const double u = rng.uniform();
        const double cos_th = std::clamp(v + (1.0 - v) * std::cbrt(1.0 - u), v, 1.0);
        const double sin_th = std::sqrt(std::max(0.0, 1.0 - cos_th * cos_th));
        const double phi    = units::two_pi * rng.uniform();

        const double q = 2.0 * (k_mag * cos_th - k_sound);
        if (q <= 0.0) return {};

        const Vec3 perp = k_dir.orthogonal_unit();
        Vec3 q_vec = (k_dir * q)
                         .rotate_around(perp, cos_th, sin_th)
                         .rotate_around(k_dir, std::cos(phi), std::sin(phi));

        const Vec3 k_recoil = k_trk - q_vec;

        if (trk.type == ExcitationType::Hole) {
            const double q_mag = q_vec.norm();
            const double e_ph  = q_mag * lattice_.config().v_sound * units::hbar_Planck;
            if (e_ph >= trk.ekin || e_ph <= 0.0) return {};

            trk.p_eV = k_recoil * units::hbarc;
            trk.ekin = lattice_.map_p_to_ekin_hole(trk.p_eV);
            return LukePhononEmit{true, q_vec.unit(), e_ph};
        } else {
            const double q_sph_mag = q_vec.norm();
            const double e_ph = q_sph_mag * lattice_.config().v_sound * units::hbar_Planck;
            if (e_ph >= trk.ekin || e_ph <= 0.0) return {};

            const Vec3 q_cryst   = lattice_.spherical_to_ellipsoidal(trk.valley, q_vec).unit();
            const Vec3 k_rec_ell = lattice_.spherical_to_ellipsoidal(trk.valley, k_recoil);
            const Vec3 p_rec_dir = lattice_.map_k_to_p_electron(trk.valley, k_rec_ell);

            const double e_recoil = trk.ekin - e_ph;
            trk.p_eV = lattice_.map_ekin_to_p_electron(trk.valley, p_rec_dir, e_recoil);
            trk.ekin = e_recoil;
            return LukePhononEmit{true, q_cryst, e_ph};
        }
    }

    // Intervalley scattering of a conduction-band electron
    inline void scatter_intervalley(ChargeTrack& trk, Rng& rng) const noexcept {
        const size_t nv = lattice_.num_valleys();
        if (trk.type != ExcitationType::Electron || nv < 2) return;

        const Vec3 k_cryst  = lattice_.map_p_to_k_electron(trk.valley, trk.p_eV);
        const Vec3 k_valley = lattice_.rotate_to_valley(trk.valley, k_cryst);

        const size_t offset = 1 + static_cast<size_t>(rng.uniform() * (nv - 1));
        const uint16_t new_valley = static_cast<uint16_t>((trk.valley + offset) % nv);

        const Vec3 k_new_cryst = lattice_.rotate_from_valley(new_valley, k_valley);
        const Vec3 p_new_dir   = lattice_.map_k_to_p_electron(new_valley, k_new_cryst);

        trk.valley = new_valley;
        trk.p_eV   = lattice_.map_ekin_to_p_electron(new_valley, p_new_dir, trk.ekin);
    }

private:
    const CrystalLattice& lattice_;
};

} // namespace cmpmc
