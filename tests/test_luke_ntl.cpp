// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Charge drift and Neganov-Trofimov-Luke (NTL) emission.
//  (a) l0 = pi hbar^4 rho / (2 m^3 Xi^2) vs Leman, arXiv:1109.1193, Table II;
//  (b) implemented emission rate vs the theta-integral of Leman Eq. (22)
//      P(theta) dtheta = (v/l0)(k/k_s)^2 (cos theta - k_s/k)^2 sin theta dtheta;
//  (c) sampled emission angle and phonon energy vs the same density;
//  (d) energy closure per carrier in a uniform field (NTL heat = W_bias = e V_b per collected pair; per carrier: sum E_Luke +
//      E_kin(final) - E_kin(initial) = e E . (r_final - r_initial));
//  (e) step-size convergence of the NTL phonon spectrum (mean, q90 of E_Luke).
#include "test_support.hpp"
#include "cmpmc/transport/engine.hpp"
#include <algorithm>
#include <vector>

using namespace cmpmc;

static void l0_table() {
    struct R { CrystalConfig c; double le, lh; };
    for (const auto& r : {R{materials::Silicon(), 16.9e-6, 7.5e-6}, R{materials::Germanium(), 257e-6, 108e-6}}) {
        const CrystalLattice lat(r.c);
        cmpmc_test::check(r.c.name + " l0 electron vs Leman Table II", lat.l0_electron(), r.le, 0.05);
        cmpmc_test::check(r.c.name + " l0 hole vs Leman Table II", lat.l0_hole(), r.lh, 0.05);
    }
}

static void rate_and_angles(const CrystalConfig& cfg) {
    const CrystalLattice lat(cfg);
    const ChargePhysics cp(lat);
    const double ks = lat.k_sound_hole(), l0 = lat.l0_hole(), v = cfg.v_sound;
    Rng rng(99);
    for (double M : {1.5, 2.0, 5.0, 20.0}) {
        ChargeTrack h = cp.make_charge(ExcitationType::Hole, {0, 0, 0}, {0, 0, 1}, 0.0, 0);
        h.p_eV = Vec3{0, 0, M * ks * units::hbarc};
        h.ekin = lat.map_p_to_ekin_hole(h.p_eV);
        // theta integral of Leman (22)
        const int N = 200000; double s = 0; const double th_max = std::acos(1.0 / M);
        for (int i = 0; i < N; ++i) { const double th = (i + 0.5) * th_max / N;
            const double c = std::cos(th) - 1.0 / M; s += (v / l0) * M * M * c * c * std::sin(th) * th_max / N; }
        const std::string n = cfg.name + " hole M=" + std::to_string(M).substr(0, 4);
        cmpmc_test::check(n + " rate impl vs int Leman (22)", cp.luke_emission_rate(h), s, 1e-4);
        cmpmc_test::info(n + " ratio to Leman (21) [(M-1)^3/M^2 form]", cp.luke_emission_rate(h) / ((v / (3 * l0)) * std::pow(M - 1, 3) / (M * M)));
        // sampled cos(theta) and phonon energy
        double sc = 0, se = 0; const int NS = 400000;
        for (int i = 0; i < NS; ++i) {
            ChargeTrack t = h;
            const auto e = cp.emit_luke_phonon(t, rng);
            if (!e.emitted) continue;
            sc += e.q_dir.dot(Vec3{0, 0, 1});
            se += e.phonon_eV;
        }
        // reference moments of p(c) ~ (c - 1/M)^2 on [1/M, 1]: <c> = 1/M + (3/4)(1 - 1/M)
        const double a = 1.0 / M, mc = a + 0.75 * (1 - a);
        // E_ph = hbar v q, q = 2 k (c - 1/M): <E> = hbar v 2 k (3/4)(1 - 1/M)
        const double me = units::hbar_Planck * v * 2.0 * M * ks * 0.75 * (1 - a);
        cmpmc_test::check(n + " <cos theta_q> vs Leman (22)", sc / NS, mc, 2e-3);
        cmpmc_test::check(n + " <E_phonon> vs Leman (22)", se / NS, me, 2e-3);
    }
}

struct DriftResult { double closure_worst = 0, e_mean = 0, e_q90 = 0, n_ph_per_carrier = 0; };

static DriftResult drift(const CrystalConfig& cfg, ExcitationType type, double E_Vm, double thick, int ncar,
                         double dk_frac) {
    const CrystalLattice lat(cfg);
    DetectorModel det;
    det.solid.radius = 5.0 * units::cm;
    det.solid.half_z = 0.5 * thick;
    det.e_field = UniformElectricField{{0.0, 0.0, -E_Vm}};
    det.surfaces[0].prob_charge_abs = 1.0; det.surfaces[1].prob_charge_abs = 1.0; det.surfaces[2].prob_charge_abs = 1.0;
    det.charge_max_dk_frac = dk_frac;
    const SimulationEngine eng(lat, det);
    const ChargePhysics cp(lat);
    Rng rng(2024);
    DriftResult R; std::vector<double> eph;
    const double z0 = (type == ExcitationType::Electron) ? -0.49 * thick : 0.49 * thick;
    for (int c = 0; c < ncar; ++c) {
        std::vector<ChargeTrack> cs{cp.make_charge(type, {0, 0, z0}, rng.isotropic_dir(), 1e-6, c % std::max<size_t>(1, lat.num_valleys()))};
        std::vector<PhononTrack> ph; std::vector<SensorHit> hits; SimStats st;
        eng.transport_charges(cs, ph, hits, st, rng);
        if (hits.empty()) continue;
        const double ek_f = hits[0].incident_energy - 0.5 * cfg.bandgap;
        const double sgn = (type == ExcitationType::Electron) ? -1.0 : 1.0;
        const double work = sgn * (det.e_field.field_Vm.dot(hits[0].pos - Vec3{0, 0, z0}));
        R.closure_worst = std::max(R.closure_worst, std::abs(st.total_luke_phonon_eV + ek_f - 1e-6 - work) / work);
        for (const auto& p : ph) eph.push_back(p.energy);
        R.n_ph_per_carrier += double(ph.size()) / ncar;
    }
    std::sort(eph.begin(), eph.end());
    double s = 0; for (double e : eph) s += e;
    R.e_mean = eph.empty() ? 0 : s / eph.size();
    R.e_q90 = eph.empty() ? 0 : eph[eph.size() * 9 / 10];
    return R;
}

int main() {
    std::printf("# test_luke_ntl\n");
    l0_table();
    rate_and_angles(materials::Germanium());
    rate_and_angles(materials::Silicon());
    for (const auto& cfg : {materials::Germanium(), materials::Silicon()}) {
        for (auto type : {ExcitationType::Electron, ExcitationType::Hole}) {
            for (double E : {100.0, 1000.0, 1.0e4}) {   // V/m (1, 10, 100 V/cm)
                const std::string n = cfg.name + (type == ExcitationType::Electron ? " e- " : " h+ ") +
                                      "E=" + std::to_string(int(E)) + "V/m";
                const double thick = 1e-3;
                const DriftResult coarse = drift(cfg, type, E, thick, 20, 0.0);
                cmpmc_test::info(n + " energy closure error, fixed-time step (charge_max_dk_frac = 0)", coarse.closure_worst);
                cmpmc_test::info(n + " <E_Luke> meV, fixed-time step", coarse.e_mean * 1e3);
                cmpmc_test::info(n + " q90 E_Luke meV, fixed-time step", coarse.e_q90 * 1e3);
                const DriftResult fine = drift(cfg, type, E, thick, 20, 0.05);
                const DriftResult finer = drift(cfg, type, E, thick, 20, 0.02);
                cmpmc_test::check(n + " energy closure, dk-limited step", fine.closure_worst, 0.0, 0.0, 1e-3);
                cmpmc_test::info(n + " <E_Luke> meV, dk <= 0.05 k_s", fine.e_mean * 1e3);
                cmpmc_test::info(n + " q90 E_Luke meV, dk <= 0.05 k_s", fine.e_q90 * 1e3);
                cmpmc_test::check(n + " <E_Luke> converged (0.05 vs 0.02 k_s, 5%)", fine.e_mean, finer.e_mean, 0.05);
            }
        }
    }
    return cmpmc_test::result();
}
