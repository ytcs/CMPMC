// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include "cmpmc/core/math.hpp"
#include "cmpmc/core/units.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace cmpmc {

// ============================================================================
// Stateless / Per-Thread PRNG (Xoshiro256**) & Angular Samplers
// ============================================================================
struct Rng {
    uint64_t s[4];

    explicit constexpr Rng(uint64_t seed = 0x853c49e6748fea9bULL) noexcept {
        uint64_t z = seed;
        for (int i = 0; i < 4; ++i) {
            z += 0x9e3779b97f4a7c15ULL;
            uint64_t x = z;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            s[i] = x ^ (x >> 31);
        }
    }

    // Independent generator for stream `stream` of a run seeded with `seed` (for
    // example one stream per event). Depends only on (seed, stream), so results
    // do not depend on how events are distributed over threads.
    static constexpr Rng for_stream(uint64_t seed, uint64_t stream) noexcept {
        return Rng(mix64(seed) ^ mix64(stream + 0x632be59bd9b4e019ULL));
    }

    static constexpr uint64_t mix64(uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
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

    // Uniform double in (0, 1) strictly positive so -std::log(uniform()) is always finite
    constexpr double uniform() noexcept {
        return (static_cast<double>((next_u64() >> 11) | 1ULL)) * 0x1.0p-53;
    }

    constexpr double uniform(double a, double b) noexcept {
        return a + (b - a) * uniform();
    }

    // Standard Gaussian N(0, 1)
    inline double normal() noexcept {
        const double u1 = uniform();
        const double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(units::two_pi * u2);
    }

    // Uniformly distributed 3D unit vector on the sphere
    inline Vec3 isotropic_dir() noexcept {
        const double cos_theta = 2.0 * uniform() - 1.0;
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = units::two_pi * uniform();
        return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
    }

    // Cosine-weighted (Lambertian) hemisphere direction around `inward_normal`
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
