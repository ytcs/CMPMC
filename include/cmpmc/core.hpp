#pragma once

#include <cmath>
#include <cstdint>
#include <numbers>

namespace cmpmc {

// ============================================================================
// Physical Constants & Internal Units
// Base units:
//   Length : meter (m)
//   Time   : second (s)
//   Energy : electron-volt (eV)
//   Mass   : kg (for density kg/m^3 and stiffness Pa = kg/(m*s^2))
// ============================================================================
namespace units {
    inline constexpr double pi          = std::numbers::pi;
    inline constexpr double two_pi      = 2.0 * std::numbers::pi;
    inline constexpr double half_pi     = 0.5 * std::numbers::pi;

    // Length (meters)
    inline constexpr double m           = 1.0;
    inline constexpr double cm          = 1e-2;
    inline constexpr double mm          = 1e-3;
    inline constexpr double um          = 1e-6;
    inline constexpr double nm          = 1e-9;
    inline constexpr double Ang         = 1e-10;

    // Time (seconds)
    inline constexpr double s           = 1.0;
    inline constexpr double ms          = 1e-3;
    inline constexpr double us          = 1e-6;
    inline constexpr double ns          = 1e-9;
    inline constexpr double ps          = 1e-12;

    // Frequency (Hz)
    inline constexpr double Hz          = 1.0;
    inline constexpr double GHz         = 1e9;
    inline constexpr double THz         = 1e12;

    // Energy (eV) & Voltage (V)
    inline constexpr double eV          = 1.0;
    inline constexpr double keV         = 1e3;
    inline constexpr double meV         = 1e-3;
    inline constexpr double ueV         = 1e-6;
    inline constexpr double Volt        = 1.0;

    // Pressure / Stiffness (Pa)
    inline constexpr double Pa          = 1.0;
    inline constexpr double GPa         = 1e9;

    // Fundamental constants
    inline constexpr double c_light     = 299792458.0;     // m / s
    inline constexpr double c_squared   = c_light * c_light;
    inline constexpr double h_Planck    = 4.135667696e-15; // eV * s
    inline constexpr double hbar_Planck = 6.582119569e-16; // eV * s
    inline constexpr double hbarc       = hbar_Planck * c_light; // eV * m
    inline constexpr double k_Boltzmann = 8.617333262e-5;  // eV / K
    inline constexpr double e_SI        = 1.602176634e-19; // J / eV (or C)
    inline constexpr double m_e_eV      = 510998.950;      // Electron rest mass E = m_e c^2 (eV)
    inline constexpr double m_e_kg      = 9.1093837015e-31;// Electron mass (kg)
} // namespace units

// ============================================================================
// Generalized Condensed-Matter Excitation Types
// Designed for Phonons, Charge Carriers, Superconducting QPs, and future
// collective excitations (Excitons, Magnons, Plasmons, Triplets).
// ============================================================================
enum class ExcitationType : uint8_t {
    Phonon         = 0, // Acoustic / optical lattice vibrations (L, FT, ST)
    Electron       = 1, // Conduction-band drift electron (with Herring-Vogt valley)
    Hole           = 2, // Valence-band drift hole
    BogoliubovQP   = 3, // Superconducting Bogoliubov quasiparticle
    Exciton        = 4, // Bound electron-hole pair (singlet/triplet, e.g. scintillators/LHe)
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

// Phonon Polarization Modes (matches G4PhononPolarization)
enum class PhononMode : uint8_t {
    L  = 0, // Longitudinal
    FT = 1, // Fast Transverse
    ST = 2, // Slow Transverse
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

// ============================================================================
// 3D Vector & 3x3 Matrix (Zero-Overhead Value Types)
// ============================================================================
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vec3() noexcept = default;
    constexpr Vec3(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_) {}

    constexpr double operator[](int i) const noexcept {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    constexpr double& operator[](int i) noexcept {
        return i == 0 ? x : (i == 1 ? y : z);
    }

    constexpr Vec3 operator+(const Vec3& r) const noexcept { return {x + r.x, y + r.y, z + r.z}; }
    constexpr Vec3 operator-(const Vec3& r) const noexcept { return {x - r.x, y - r.y, z - r.z}; }
    constexpr Vec3 operator-() const noexcept { return {-x, -y, -z}; }
    constexpr Vec3 operator*(double s) const noexcept { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const noexcept {
        const double inv = 1.0 / s;
        return {x * inv, y * inv, z * inv};
    }

    constexpr Vec3& operator+=(const Vec3& r) noexcept { x += r.x; y += r.y; z += r.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& r) noexcept { x -= r.x; y -= r.y; z -= r.z; return *this; }
    constexpr Vec3& operator*=(double s) noexcept { x *= s; y *= s; z *= s; return *this; }

    constexpr double dot(const Vec3& r) const noexcept { return x * r.x + y * r.y + z * r.z; }
    constexpr Vec3 cross(const Vec3& r) const noexcept {
        return {
            y * r.z - z * r.y,
            z * r.x - x * r.z,
            x * r.y - y * r.x
        };
    }

    constexpr double norm2() const noexcept { return x * x + y * y + z * z; }
    inline double norm() const noexcept { return std::sqrt(norm2()); }

    inline Vec3 unit() const noexcept {
        const double n2 = norm2();
        if (n2 <= 0.0) return {0.0, 0.0, 1.0};
        const double inv = 1.0 / std::sqrt(n2);
        return {x * inv, y * inv, z * inv};
    }

    // Rotate vector around unit axis `k` by angle `phi` (Rodrigues' rotation formula)
    inline Vec3 rotate_around(const Vec3& k, double cos_phi, double sin_phi) const noexcept {
        return (*this) * cos_phi + k.cross(*this) * sin_phi + k * (k.dot(*this) * (1.0 - cos_phi));
    }

    // Construct any unit vector orthogonal to *this
    inline Vec3 orthogonal_unit() const noexcept {
        Vec3 a = std::abs(x) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        return cross(a).unit();
    }
};

constexpr Vec3 operator*(double s, const Vec3& v) noexcept { return v * s; }

struct Mat3 {
    double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    static constexpr Mat3 identity() noexcept {
        return Mat3{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    }

    static constexpr Mat3 diagonal(double xx, double yy, double zz) noexcept {
        return Mat3{{{xx, 0, 0}, {0, yy, 0}, {0, 0, zz}}};
    }

    static constexpr Mat3 from_rows(const Vec3& r0, const Vec3& r1, const Vec3& r2) noexcept {
        return Mat3{{{r0.x, r0.y, r0.z}, {r1.x, r1.y, r1.z}, {r2.x, r2.y, r2.z}}};
    }

    constexpr Vec3 operator*(const Vec3& v) const noexcept {
        return {
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        };
    }

    // Transpose-multiply (inverse rotation for orthogonal matrices)
    constexpr Vec3 mul_transpose(const Vec3& v) const noexcept {
        return {
            m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z,
            m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z,
            m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z
        };
    }
};

// ============================================================================
// Fast Stateless / Per-Thread PRNG: Xoshiro256**
// Replaces CLHEP HepRandom::getTheEngine()->flat()
// ============================================================================
struct Rng {
    uint64_t s[4];

    explicit constexpr Rng(uint64_t seed = 0x853c49e6748fea9bULL) noexcept {
        // SplitMix64 initialization
        uint64_t z = seed;
        for (int i = 0; i < 4; ++i) {
            z += 0x9e3779b97f4a7c15ULL;
            uint64_t x = z;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            s[i] = x ^ (x >> 31);
        }
    }

    static constexpr uint64_t rotl(uint64_t x, int k) noexcept {
        return (x << k) | (x >> (64 - k));
    }

    constexpr uint64_t next_u64() noexcept {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    // Uniform double in (0, 1) strictly positive so -std::log(uniform()) never overflows
    constexpr double uniform() noexcept {
        return (static_cast<double>((next_u64() >> 11) | 1ULL)) * 0x1.0p-53;
    }

    // Uniform double in [a, b)
    constexpr double uniform(double a, double b) noexcept {
        return a + (b - a) * uniform();
    }

    // Standard normal distribution N(0, 1) via Box-Muller
    inline double normal() noexcept {
        const double u1 = uniform();
        const double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(units::two_pi * u2);
    }

    // Isotropic 3D unit vector (replaces G4RandomDirection)
    inline Vec3 isotropic_dir() noexcept {
        const double cos_theta = 2.0 * uniform() - 1.0;
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = units::two_pi * uniform();
        return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
    }

    // Lambertian (cosine-weighted) hemisphere reflection around inward surface normal `n`
    // Matches G4CMP::LambertReflection
    inline Vec3 lambertian_dir(const Vec3& inward_normal) noexcept {
        const double phi = units::two_pi * uniform();
        const double sin_theta = std::sqrt(uniform());
        const double cos_theta = std::sqrt(std::max(0.0, 1.0 - sin_theta * sin_theta));

        const Vec3 u = inward_normal.orthogonal_unit();
        const Vec3 v = inward_normal.cross(u);
        return (u * (sin_theta * std::cos(phi)) +
                v * (sin_theta * std::sin(phi)) +
                inward_normal * cos_theta).unit();
    }
};

} // namespace cmpmc
