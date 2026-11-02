// Charge-carrier kinematics, valley orientation and intervalley rates.
//  (a) Herring-Vogt valley axes follow the crystal orientation: with orient = R the
//      electron velocity for detector-frame momentum R p must equal R v(p) of the
//      unrotated crystal;
//  (b) HV kinematics: E_kin = sum_i p_i^2/(2 m_i) and v = M^-1 p in the valley frame;
//  (c) intervalley rates: empirical Quadratic and Linear forms,
//      Matrix optical rate vs Jacoboni-Reggiani with Z_f = 2N-1 (as implemented) and
//      the per-phonon final-valley count (reported, not asserted).
#include "test_support.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"

using namespace cmpmc;

int main() {
    std::printf("# test_charge_misc\n");
    CrystalConfig c0 = materials::Germanium();
    CrystalConfig c1 = c0; c1.set_miller_orientation(1, 1, 0, 0.3);
    const CrystalLattice L0(c0), L1(c1);
    Rng rng(4); double worst = 0;
    for (int t = 0; t < 1000; ++t) {
        const Vec3 p = rng.isotropic_dir() * 3.0;      // eV (hbar c k)
        const size_t iv = t % L0.num_valleys();
        const Vec3 v0 = L0.map_p_to_vel_electron(iv, p);
        const Vec3 v1 = L1.map_p_to_vel_electron(iv, c1.orient * p);
        worst = std::max(worst, (v1 - c1.orient * v0).norm() / v0.norm());
    }
    cmpmc_test::check("Ge e- velocity covariant under crystal rotation (max rel)", worst, 0.0, 0.0, 1e-12);

    // (b) HV kinematics in valley 0 of Ge ([111] axis): p along the valley axis
    {
        const Vec3 ax = Vec3{1, 1, 1}.unit();
        const double p = 2.0; const double mL = c0.e_mass_tensor.x * units::m_e_eV;
        cmpmc_test::check("Ge e- E_kin(p || [111]) = p^2/(2 m_l)", L0.map_p_to_ekin_electron(0, ax * p), p * p / (2 * mL), 1e-6);
        cmpmc_test::check("Ge e- |v|(p || [111]) = p c/m_l", L0.map_p_to_vel_electron(0, ax * p).norm(), p * units::c_light / mL, 1e-6);
    }
    // (c) IV rates
    {
        const CrystalConfig& g = c0; const Vec3 E{0, 0, 150.0};
        const double q = g.iv_quad_rate * std::pow(g.iv_quad_field * g.iv_quad_field + 150.0 * 150.0, 0.5 * g.iv_quad_power);
        cmpmc_test::check("Ge IV Quadratic rate at 150 V/m (closed form)", physics::intervalley_scattering_rate(0.02, E, g), q, 1e-12);
        const CrystalConfig s = materials::Silicon(); const Vec3 Es{0, 0, 2.0e4};
        const double l = s.iv_lin_rate0 + s.iv_lin_rate1 * std::pow(200.0, s.iv_lin_power);
        cmpmc_test::check("Si IV Linear rate at 200 V/cm (closed form)", physics::intervalley_scattering_rate(0.02, Es, s), l, 1e-12);
        // Matrix optical part: J&R W = Z_f m_d^{3/2} D^2 sqrt(E - hw)(..)/(sqrt2 pi rho hbar^2 hw) (T = 0 emission)
        CrystalConfig sm = s; sm.iv_model = IVRateModel::Matrix; sm.neutral_impurity_cm3 = 0.0;
        const double Ek = 0.08;  // eV, above all six Si thresholds
        const double impl = physics::intervalley_matrix_rate(Ek, sm);
        // Si: g-type phonons (first three) Z_f = 1, f-type (last three) Z_f = 4 (Jacoboni & Reggiani 1983)
        CrystalConfig one = sm; one.valley_dirs = {{1,0,0},{0,1,0}};   // 2N-1 = 3 -> rescale below
        double jr = 0;
        for (size_t i = 0; i < sm.iv_energy_eV.size(); ++i) {
            CrystalConfig ci = sm; ci.iv_energy_eV = {sm.iv_energy_eV[i]}; ci.iv_deform_eV_m = {sm.iv_deform_eV_m[i]};
            const double per = physics::intervalley_matrix_rate(Ek, ci) / 5.0;   // per final valley (2N-1 = 5 for N=3)
            jr += per * (i < 3 ? 1.0 : 4.0);
        }
        cmpmc_test::info("Si IV Matrix rate impl (1/s) at 80 meV", impl);
        cmpmc_test::info("Si IV Matrix rate with J&R Z_f (g:1, f:4)", jr);
        cmpmc_test::info("ratio impl / J&R", impl / jr);
        (void)one;
    }
    return cmpmc_test::result();
}
