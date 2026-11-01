// R16 test: closed-form (trigonometric Cardano + cross-product) Christoffel solve against the
// Jacobi solve (R15 T1-validated) and the [100]/[110]/[111] closed forms, for cubic (Si, Ge, GaAs),
// trigonal (sapphire, alpha-quartz) and tetragonal (CaWO4) stiffness tensors.
//  (a) phase speeds of all modes, 2e5 random directions: max rel. diff <= 1e-10;
//  (b) group velocities: max |dv_g|/|v_g| over directions whose eigenvalue gap is > 1e-3
//      (well separated): <= 1e-8; near-degenerate directions fall back to Jacobi (identical);
//  (c) Si [100]: v_L = sqrt(C11/rho), v_T = sqrt(C44/rho) (closed form, Leman 1109.1193 Table I
//      inputs via the preset), [111] v_L = sqrt((C11 + 2C12 + 4C44)/(3 rho));
//  (d) cost per call (reported) and that CrystalLattice::map_k_to_vg_exact uses it.
#include "audit_common.hpp"
#include "cmpmc/lattice.hpp"
#include <chrono>
#include <string>

using namespace cmpmc;

static void build_C(const CrystalConfig& cfg, double C[3][3][3][3]) {
    auto vi = [](int i, int j) { if (i == j) return i; if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
                                 if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4; return 5; };
    const auto V = cfg.build_voigt_matrix();
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m)
        C[i][j][l][m] = V[vi(i, j)][vi(l, m)];
}

int main() {
    std::printf("# test_christoffel_cardano\n");
    for (const auto& cfg : {materials::Silicon(), materials::Germanium(), materials::GalliumArsenide(),
                            materials::Sapphire(), materials::AlphaQuartz(), materials::CalciumTungstate()}) {
        double C[3][3][3][3]; build_C(cfg, C);
        Rng rng(77);
        double dvp = 0, dvg = 0; int nfar = 0;
        for (int n = 0; n < 200000; ++n) {
            const Vec3 k = rng.isotropic_dir();
            Vec3 g1[3], g2[3]; double p1[3], p2[3];
            physics::solve_christoffel_exact(C, cfg.density, k, g1, p1);
            physics::solve_christoffel_cardano(C, cfg.density, k, g2, p2);
            double gapmin = 1e300;
            for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) gapmin = std::min(gapmin, std::abs(p1[a] * p1[a] - p1[b] * p1[b]) / (p1[0] * p1[0]));
            for (int m = 0; m < 3; ++m) {
                dvp = std::max(dvp, std::abs(p2[m] / p1[m] - 1.0));
                if (gapmin > 1e-3) { dvg = std::max(dvg, (g2[m] - g1[m]).norm() / g1[m].norm()); }
            }
            if (gapmin > 1e-3) ++nfar;
        }
        audit::check(cfg.name + " max rel |v_p cardano - jacobi|", dvp, 0.0, 0.0, 1e-10);
        audit::check(cfg.name + " max rel |v_g cardano - jacobi| (gap > 1e-3)", dvg, 0.0, 0.0, 1e-8);
        audit::info(cfg.name + " fraction of directions with gap > 1e-3", nfar / 200000.0);
        // timing
        Rng r2(5); Vec3 g[3]; double p[3]; double sink = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (int n = 0; n < 200000; ++n) { physics::solve_christoffel_exact(C, cfg.density, r2.isotropic_dir(), g, p); sink += p[0]; }
        auto t1 = std::chrono::steady_clock::now();
        for (int n = 0; n < 200000; ++n) { physics::solve_christoffel_cardano(C, cfg.density, r2.isotropic_dir(), g, p); sink += p[0]; }
        auto t2 = std::chrono::steady_clock::now();
        audit::info(cfg.name + " ns/call Jacobi", std::chrono::duration<double, std::nano>(t1 - t0).count() / 200000);
        audit::info(cfg.name + " ns/call Cardano", std::chrono::duration<double, std::nano>(t2 - t1).count() / 200000 + 0 * sink);
    }
    // (c) closed forms, Si
    const auto si = materials::Silicon();
    double C[3][3][3][3]; build_C(si, C);
    const auto V = si.build_voigt_matrix();
    Vec3 g[3]; double p[3];
    physics::solve_christoffel_cardano(C, si.density, {1, 0, 0}, g, p);   // degenerate T: Jacobi fallback
    audit::check("Si [100] v_L vs sqrt(C11/rho)", p[0], std::sqrt(V[0][0] / si.density), 1e-12);
    audit::check("Si [100] v_T vs sqrt(C44/rho)", p[1], std::sqrt(V[3][3] / si.density), 1e-12);
    const double s3 = 1.0 / std::sqrt(3.0);
    physics::solve_christoffel_cardano(C, si.density, {s3, s3, s3}, g, p);
    audit::check("Si [111] v_L vs sqrt((C11+2C12+4C44)/(3 rho))", p[0], std::sqrt((V[0][0] + 2 * V[0][1] + 4 * V[3][3]) / (3 * si.density)), 1e-12);
    const double s2 = 1.0 / std::sqrt(2.0);
    physics::solve_christoffel_cardano(C, si.density, {s2, s2, 0}, g, p);   // [110]: non-degenerate
    audit::check("Si [110] v_L vs sqrt((C11+C12+2C44)/(2 rho))", p[0], std::sqrt((V[0][0] + V[0][1] + 2 * V[3][3]) / (2 * si.density)), 1e-12);
    audit::check("Si [110] v_FT vs sqrt(C44/rho)", p[1], std::sqrt(V[3][3] / si.density), 1e-12);
    audit::check("Si [110] v_ST vs sqrt((C11-C12)/(2 rho))", p[2], std::sqrt((V[0][0] - V[0][1]) / (2 * si.density)), 1e-12);
    // (d) lattice Exact mode goes through the Cardano path
    CrystalLattice lat(si); lat.set_kinematics_mode(CrystalLattice::KinematicsMode::Exact);
    Rng r3(8); const Vec3 k = r3.isotropic_dir();
    const auto kin = lat.map_k_to_vg_exact(PhononMode::ST, k);
    physics::solve_christoffel_cardano(C, si.density, k, g, p);
    audit::check("lattice Exact ST |v_g| == cardano", kin.vg_mag, g[2].norm(), 1e-12);
    return audit::result();
}
