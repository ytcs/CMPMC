// Audit test: Christoffel solver and group-velocity lookup tables.
//  (a) phase velocities along [100], [110], [111] vs the closed cubic forms;
//  (b) group velocity vs a central finite difference of omega(k) = |k| v_p(k_hat);
//  (c) error of the bilinear (cos theta, phi) LUT (map_k_to_vg) and of the FP32
//      octahedral LUT (map_k_to_vg_oct) against the exact solver over random k_hat,
//      reported as quantiles of the direction error (mrad) and |v_g| error.
#include "audit_common.hpp"
#include "cmpmc/lattice.hpp"
#include <algorithm>
#include <vector>
#include <chrono>

using namespace cmpmc;

static void stiff(const CrystalConfig& cfg, double C[3][3][3][3]) {
    const auto V = cfg.build_voigt_matrix();
    auto vi = [](int i, int j) { if (i == j) return i; if (i + j == 3) return 3; if (i + j == 2) return 4; return 5; };
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m)
        C[i][j][l][m] = V[vi(i, j)][vi(l, m)];
}

static void exact(const double C[3][3][3][3], double rho, const Vec3& k, Vec3 vg[3], double vp[3]) {
    physics::solve_christoffel_exact(C, rho, k.unit(), vg, vp);
}

static void run_cubic(const CrystalConfig& cfg) {
    double C[3][3][3][3]; stiff(cfg, C);
    const double r = cfg.density, c11 = cfg.C11, c12 = cfg.C12, c44 = cfg.C44;
    Vec3 vg[3]; double vp[3];
    const std::string n = cfg.name;
    exact(C, r, {1, 0, 0}, vg, vp);
    audit::check(n + " [100] v_L", vp[0], std::sqrt(c11 / r), 1e-10);
    audit::check(n + " [100] v_T (FT)", vp[1], std::sqrt(c44 / r), 1e-10);
    audit::check(n + " [100] v_T (ST)", vp[2], std::sqrt(c44 / r), 1e-10);
    exact(C, r, {1, 1, 0}, vg, vp);
    audit::check(n + " [110] v_L", vp[0], std::sqrt((c11 + c12 + 2 * c44) / (2 * r)), 1e-10);
    const double a = std::sqrt(c44 / r), b = std::sqrt((c11 - c12) / (2 * r));
    audit::check(n + " [110] v_FT", vp[1], std::max(a, b), 1e-10);
    audit::check(n + " [110] v_ST", vp[2], std::min(a, b), 1e-10);
    exact(C, r, {1, 1, 1}, vg, vp);
    audit::check(n + " [111] v_L", vp[0], std::sqrt((c11 + 2 * c12 + 4 * c44) / (3 * r)), 1e-10);
    audit::check(n + " [111] v_T", vp[1], std::sqrt((c11 - c12 + c44) / (3 * r)), 1e-10);
}

static void run_lut(const CrystalConfig& cfg) {
    double C[3][3][3][3]; stiff(cfg, C);
    const CrystalLattice lat(cfg);
    Rng rng(7);
    const std::string n = cfg.name;
    // (b) group velocity = grad_k omega (finite differences, h = 1e-6)
    double worst_fd = 0;
    for (int t = 0; t < 2000; ++t) {
        const Vec3 k = rng.isotropic_dir();
        Vec3 vg[3]; double vp[3]; exact(C, cfg.density, k, vg, vp);
        for (int m = 0; m < 3; ++m) {
            Vec3 g{};
            for (int dmn = 0; dmn < 3; ++dmn) {
                Vec3 kp = k, km = k; const double h = 1e-6;
                if (dmn == 0) { kp.x += h; km.x -= h; } else if (dmn == 1) { kp.y += h; km.y -= h; } else { kp.z += h; km.z -= h; }
                Vec3 a[3], b[3]; double pa[3], pb[3]; exact(C, cfg.density, kp, a, pa); exact(C, cfg.density, km, b, pb);
                const double w = (kp.norm() * pa[m] - km.norm() * pb[m]) / (2 * h);
                if (dmn == 0) g.x = w; else if (dmn == 1) g.y = w; else g.z = w;
            }
            // skip points within 1e-3 of a degeneracy (gradient undefined there)
            if (std::abs(vp[1] - vp[2]) < 1e-3 * vp[1] && m > 0) continue;
            worst_fd = std::max(worst_fd, (g - vg[m]).norm() / vg[m].norm());
        }
    }
    audit::check(n + " v_g vs finite-difference grad omega (max rel)", worst_fd, 0.0, 0.0, 1e-5);

    // (c) LUT errors
    for (int which = 0; which < 2; ++which) {
        for (int m = 0; m < 3; ++m) {
            std::vector<double> ang, mag;
            int nbad = 0; const int NT = 200000;
            for (int t = 0; t < NT; ++t) {
                const Vec3 k = rng.isotropic_dir();
                Vec3 vg[3]; double vp[3]; exact(C, cfg.density, k, vg, vp);
                const auto kin = which ? lat.map_k_to_vg_oct(static_cast<PhononMode>(m), k)
                                       : lat.map_k_to_vg(static_cast<PhononMode>(m), k);
                const double c = std::clamp(kin.vg_dir.dot(vg[m].unit()), -1.0, 1.0);
                ang.push_back(std::acos(c) * 1e3);
                mag.push_back(std::abs(kin.vg_mag - vg[m].norm()) / vg[m].norm());
                if (ang.back() > 17.45) ++nbad; // > 1 degree
            }
            std::sort(ang.begin(), ang.end()); std::sort(mag.begin(), mag.end());
            const std::string tag = n + (which ? " octLUT " : " map_k_to_vg ") + mode_name(static_cast<PhononMode>(m));
            audit::info(tag + " dir err median (mrad)", ang[NT / 2]);
            audit::info(tag + " dir err q99 (mrad)", ang[NT * 99 / 100]);
            audit::info(tag + " dir err max (mrad)", ang.back());
            audit::info(tag + " |vg| err q99 (rel)", mag[NT * 99 / 100]);
            audit::info(tag + " fraction dir err > 1 deg", double(nbad) / NT);
            if (m == 0) audit::check(tag + " q99 dir err < 5 mrad", ang[NT * 99 / 100], 0.0, 0.0, 5.0);
        }
    }
}

static void timing() {
    const CrystalLattice lat(materials::Germanium());
    Rng rng(3); std::vector<Vec3> ks(1 << 20); for (auto& k : ks) k = rng.isotropic_dir();
    double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (const auto& k : ks) sink += lat.map_k_to_vg(PhononMode::FT, k).vg_mag;
    auto t1 = std::chrono::steady_clock::now();
    for (const auto& k : ks) sink += lat.map_k_to_vg_oct(PhononMode::FT, k).vg_mag;
    auto t2 = std::chrono::steady_clock::now();
    audit::info("ns per map_k_to_vg call", std::chrono::duration<double, std::nano>(t1 - t0).count() / ks.size());
    audit::info("ns per map_k_to_vg_oct call", std::chrono::duration<double, std::nano>(t2 - t1).count() / ks.size());
    if (sink == 42.0) std::printf("x");
}

#ifndef CMPMC_AUDIT_REF
static void exact_mode() {
    CrystalLattice lat(materials::Germanium());
    lat.set_kinematics_mode(CrystalLattice::KinematicsMode::Exact);
    double C[3][3][3][3]; stiff(lat.config(), C);
    Rng rng(11); double worst = 0;
    for (int t = 0; t < 20000; ++t) { const Vec3 k = rng.isotropic_dir(); Vec3 vg[3]; double vp[3];
        exact(C, lat.config().density, k, vg, vp);
        worst = std::max(worst, (lat.map_k_to_vg(PhononMode::L, k).vg_dir - vg[0].unit()).norm()); }
    audit::check("Ge Exact kinematics mode max dir err (rad)", worst, 0.0, 0.0, 1e-9);
}
#endif

int main() {
#ifndef CMPMC_AUDIT_REF
    exact_mode();
#endif
    timing();
    std::printf("# test_christoffel\n");
    run_cubic(materials::Germanium());
    run_cubic(materials::Silicon());
    run_cubic(materials::GalliumArsenide());
    run_lut(materials::Germanium());
    run_lut(materials::Silicon());
    run_lut(materials::Sapphire());
    return audit::result();
}
