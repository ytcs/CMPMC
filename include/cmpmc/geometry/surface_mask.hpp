#pragma once

#include "cmpmc/core.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cmpmc::geometry {

// ============================================================================
// ElectrodePatternMask: Fast 2D Surface Sensor Masking
// Replaces: G4CMPVElectrodePattern & complex Geant4 daughter volume arrays
//
// Physics & Detector Design Rationale:
//   Instead of placing thousands of tiny Geant4 daughter volumes (G4Box/G4Tubs)
//   for QET fins, KID resonators, or qubit pads on a crystal surface (which
//   causes O(N) G4Navigator voxel slowdown), CMPMC evaluates an analytic 2D
//   surface mask at the moment a phonon, charge, or quasiparticle hits a face.
//
//   If `find_channel(pos, rng)` returns `channel_id >= 0`, the track hits active
//   sensor metal for that readout channel; if it returns `-1`, the track hits
//   bare crystal between pads and undergoes bare-surface reflection.
// ============================================================================
struct ElectrodePad {
    enum class Shape : uint8_t {
        Disk          = 0, // Centered at (cx, cy) with radius r_max
        AnnularSector = 1, // r in [r_min, r_max], phi in [phi_min, phi_max]
        Rectangle     = 2  // [cx - hx, cx + hx] x [cy - hy, cy + hy]
    };

    Shape shape          = Shape::Disk;
    int32_t channel_id   = 0;
    double cx            = 0.0; // Center x (m) for Disk / Rectangle
    double cy            = 0.0; // Center y (m) for Disk / Rectangle
    double r_min         = 0.0; // Inner radius (m) or half-width hx (m) for Rectangle
    double r_max         = 0.0; // Outer radius (m) or half-height hy (m) for Rectangle
    double phi_min       = 0.0; // Start azimuth [0, 2*pi] for AnnularSector
    double phi_max       = units::two_pi; // End azimuth [0, 2*pi] for AnnularSector
    double fill_fraction = 1.0; // Sub-pad active metal coverage fraction in (0, 1]

    inline bool contains(double x, double y) const noexcept {
        switch (shape) {
            case Shape::Disk: {
                const double dx = x - cx;
                const double dy = y - cy;
                return (dx * dx + dy * dy) <= (r_max * r_max);
            }
            case Shape::AnnularSector: {
                const double dx = x - cx;
                const double dy = y - cy;
                const double r2 = dx * dx + dy * dy;
                if (r2 < r_min * r_min || r2 > r_max * r_max) return false;
                double phi = std::atan2(dy, dx);
                if (phi < 0.0) phi += units::two_pi;
                if (phi_min <= phi_max) {
                    return (phi >= phi_min && phi <= phi_max);
                } else {
                    // Wraps across phi = 0
                    return (phi >= phi_min || phi <= phi_max);
                }
            }
            case Shape::Rectangle: {
                return (std::abs(x - cx) <= r_min) && (std::abs(y - cy) <= r_max);
            }
        }
        return false;
    }
};

class ElectrodePatternMask {
public:
    ElectrodePatternMask() = default;
    explicit ElectrodePatternMask(std::vector<ElectrodePad> pads)
        : pads_(std::move(pads)) {}

    bool empty() const noexcept { return pads_.empty(); }
    size_t num_pads() const noexcept { return pads_.size(); }
    const std::vector<ElectrodePad>& pads() const noexcept { return pads_; }

    void add_pad(const ElectrodePad& pad) { pads_.push_back(pad); }

    // Evaluate which electrode channel (if any) is hit at `pos`.
    // Returns channel_id >= 0 if absorbed/active pad, or -1 if bare crystal gap.
    inline int32_t find_channel(const Vec3& pos, Rng& rng) const noexcept {
        for (const ElectrodePad& pad : pads_) {
            if (pad.contains(pos.x, pos.y)) {
                if (pad.fill_fraction >= 1.0 || rng.uniform() < pad.fill_fraction) {
                    return pad.channel_id;
                }
                return -1;
            }
        }
        return -1;
    }

    // Deterministic query (ignores stochastic sub-pad fill_fraction)
    inline int32_t find_geometric_channel(double x, double y) const noexcept {
        for (const ElectrodePad& pad : pads_) {
            if (pad.contains(x, y)) return pad.channel_id;
        }
        return -1;
    }

    // Preset 1: SuperCDMS-style 4-channel Phonon QET layout:
    //   Channel 0: Inner core disk (r <= 0.55 * R)
    //   Channels 1, 2, 3: Three 120-degree outer annular sectors (0.57 * R <= r <= 0.90 * R)
    //   Channel 4: Outer guard ring (0.92 * R <= r <= 0.98 * R)
    static ElectrodePatternMask supercdms_qet_mask(double crystal_radius_m,
                                                   double fill_fraction = 1.0) {
        std::vector<ElectrodePad> pads;
        // Inner core disk (Channel 0)
        pads.push_back(ElectrodePad{
            .shape         = ElectrodePad::Shape::Disk,
            .channel_id    = 0,
            .cx            = 0.0,
            .cy            = 0.0,
            .r_min         = 0.0,
            .r_max         = 0.55 * crystal_radius_m,
            .phi_min       = 0.0,
            .phi_max       = units::two_pi,
            .fill_fraction = fill_fraction
        });

        // 3 Outer 120-deg sectors (Channels 1, 2, 3) with 2-deg azimuthal gaps
        const double gap_rad = 2.0 * (units::pi / 180.0);
        const double sec_span = units::two_pi / 3.0;
        for (int s = 0; s < 3; ++s) {
            pads.push_back(ElectrodePad{
                .shape         = ElectrodePad::Shape::AnnularSector,
                .channel_id    = s + 1,
                .cx            = 0.0,
                .cy            = 0.0,
                .r_min         = 0.57 * crystal_radius_m,
                .r_max         = 0.90 * crystal_radius_m,
                .phi_min       = s * sec_span + 0.5 * gap_rad,
                .phi_max       = (s + 1) * sec_span - 0.5 * gap_rad,
                .fill_fraction = fill_fraction
            });
        }

        // Outer guard ring (Channel 4)
        pads.push_back(ElectrodePad{
            .shape         = ElectrodePad::Shape::AnnularSector,
            .channel_id    = 4,
            .cx            = 0.0,
            .cy            = 0.0,
            .r_min         = 0.92 * crystal_radius_m,
            .r_max         = 0.98 * crystal_radius_m,
            .phi_min       = 0.0,
            .phi_max       = units::two_pi,
            .fill_fraction = fill_fraction
        });

        return ElectrodePatternMask(std::move(pads));
    }

    // Preset 2: Rectangular NX x NY grid array of KID / Qubit / pixel pads
    static ElectrodePatternMask rectangular_grid_mask(int nx, int ny,
                                                      double pitch_x, double pitch_y,
                                                      double pad_half_w, double pad_half_h,
                                                      double fill_fraction = 1.0) {
        std::vector<ElectrodePad> pads;
        pads.reserve(static_cast<size_t>(nx * ny));
        int32_t ch = 0;
        for (int ix = 0; ix < nx; ++ix) {
            const double cx = (ix - 0.5 * (nx - 1)) * pitch_x;
            for (int iy = 0; iy < ny; ++iy) {
                const double cy = (iy - 0.5 * (ny - 1)) * pitch_y;
                pads.push_back(ElectrodePad{
                    .shape         = ElectrodePad::Shape::Rectangle,
                    .channel_id    = ch++,
                    .cx            = cx,
                    .cy            = cy,
                    .r_min         = pad_half_w,
                    .r_max         = pad_half_h,
                    .phi_min       = 0.0,
                    .phi_max       = 0.0,
                    .fill_fraction = fill_fraction
                });
            }
        }
        return ElectrodePatternMask(std::move(pads));
    }

private:
    std::vector<ElectrodePad> pads_;
};

} // namespace cmpmc::geometry
