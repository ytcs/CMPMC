// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace cmpmc {

// ============================================================================
// 3D Vector & 3x3 Matrix Value Types
// ============================================================================
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vec3() noexcept = default;
    constexpr Vec3(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_) {}

    constexpr double operator[](int i) const noexcept { return i == 0 ? x : (i == 1 ? y : z); }
    constexpr double& operator[](int i) noexcept { return i == 0 ? x : (i == 1 ? y : z); }

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
        return {y * r.z - z * r.y, z * r.x - x * r.z, x * r.y - y * r.x};
    }

    constexpr double norm2() const noexcept { return x * x + y * y + z * z; }
    inline double norm() const noexcept { return std::sqrt(norm2()); }

    inline Vec3 unit() const noexcept {
        const double n2 = norm2();
        if (n2 <= 0.0) return {0.0, 0.0, 1.0};
        const double inv = 1.0 / std::sqrt(n2);
        return {x * inv, y * inv, z * inv};
    }

    // Rodrigues' rotation around unit axis `k` by angle `phi`
    inline Vec3 rotate_around(const Vec3& k, double cos_phi, double sin_phi) const noexcept {
        return (*this) * cos_phi + k.cross(*this) * sin_phi + k * (k.dot(*this) * (1.0 - cos_phi));
    }

    // Any unit vector orthogonal to *this
    inline Vec3 orthogonal_unit() const noexcept {
        const Vec3 a = (std::abs(x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
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

// Fixed-capacity stack buffer to avoid std::vector heap allocations in physics cascades
template <typename T, size_t Capacity = 64>
struct SmallVec {
    std::array<T, Capacity> data{};
    size_t count = 0;

    constexpr void clear() noexcept { count = 0; }
    constexpr bool empty() const noexcept { return count == 0; }
    constexpr size_t size() const noexcept { return count; }
    constexpr void push_back(const T& v) noexcept {
        if (count < Capacity) data[count++] = v;
    }
    constexpr void pop_back() noexcept {
        if (count > 0) --count;
    }
    constexpr T& back() noexcept { return data[count - 1]; }
    constexpr const T& back() const noexcept { return data[count - 1]; }
    constexpr T* begin() noexcept { return data.data(); }
    constexpr T* end() noexcept { return data.data() + count; }
    constexpr const T* begin() const noexcept { return data.data(); }
    constexpr const T* end() const noexcept { return data.data() + count; }
    constexpr T& operator[](size_t i) noexcept { return data[i]; }
    constexpr const T& operator[](size_t i) const noexcept { return data[i]; }
};

// Buffer of the Kaplan film cascade (QP / phonon generations and re-emitted phonons).
// Capacity 256: the cascade of one phonon of energy Omega re-emits at most ~Omega/(2 Delta)
// phonons, so the buffer holds the full cascade of Debye-scale phonons (several hundred Delta).
// A smaller buffer would drop entries at push_back and lose energy; test_pair_breaking checks
// energy closure up to 350 Delta.
using KaplanBuffer = SmallVec<double, 256>;

} // namespace cmpmc
