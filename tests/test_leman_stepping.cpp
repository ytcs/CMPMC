// R16 test: Leman second-order charge time stepping (DetectorModel::ChargeStepping::LemanSecondOrder;
// 1109.1193 Sec. V.E Eqs. (26)-(30) [FT]) against a converged first-order reference.
//  Reference: first-order stepping with the wavevector change per step <= 0.01 max(k_s, |k|)
//  (R15 P5 at 5x finer cap; P5's own default is 0.05). Second order: sampling step <= 0.25 max(k_s, |k|).
//  For e-/h+ in Ge and Si at 1, 10, 100 V/cm over 1 mm:
//  (a) per-carrier energy closure (campaign WP7-E-v3 (E11)): sum E_Luke + dE_kin = e E . dr, <= 1e-3;
//  (b) NTL phonon mean energy <E_Luke> and phonons per carrier vs the reference, 3%;
//  (c) cost: charge steps per carrier relative to P5 (reported; check: fewer steps than P5);
//  (d) transport_charges_fast_chord now runs the Leman path (identical output for one seed).
#include "audit_common.hpp"
#include "cmpmc/engine.hpp"
#include <string>
#include <vector>

using namespace cmpmc;

struct R { double closure = 0, emean = 0, nph = 0, steps = 0; };

static R drift(const CrystalConfig& cfg, ExcitationType type, double E_Vm, int ncar, bool second, double frac) {
    const double thick = 1e-3;
    const CrystalLattice lat(cfg);
    DetectorModel det;
    det.solid.radius = 5.0 * units::cm; det.solid.half_z = 0.5 * thick;
    det.e_field = UniformElectricField{{0.0, 0.0, -E_Vm}};
    if (second) { det.charge_stepping = DetectorModel::ChargeStepping::LemanSecondOrder; det.charge_sampling_dk_frac = frac; }
    else det.charge_max_dk_frac = frac;
    const SimulationEngine eng(lat, det);
    const ChargePhysics cp(lat);
    Rng rng(4242);
    R r; double s = 0, n = 0;
    const double z0 = (type == ExcitationType::Electron) ? -0.49 * thick : 0.49 * thick;
    for (int c = 0; c < ncar; ++c) {
        std::vector<ChargeTrack> cs{cp.make_charge(type, {0, 0, z0}, rng.isotropic_dir(), 1e-6, c % std::max<size_t>(1, lat.num_valleys()))};
        std::vector<PhononTrack> ph; std::vector<SensorHit> hits; SimStats st;
        eng.transport_charges(cs, ph, hits, st, rng);
        if (hits.empty()) continue;
        const double ekf = hits[0].incident_energy - 0.5 * cfg.bandgap;
        const double sg = (type == ExcitationType::Electron) ? -1.0 : 1.0;
        const double work = sg * det.e_field.field_Vm.dot(hits[0].pos - Vec3{0, 0, z0});
        r.closure = std::max(r.closure, std::abs(st.total_luke_phonon_eV + ekf - 1e-6 - work) / work);
        for (const auto& p : ph) { s += p.energy; n += 1; }
        r.steps += double(st.charge_steps) / ncar;
    }
    r.emean = s / n; r.nph = n / ncar;
    return r;
}

int main() {
    std::printf("# test_leman_stepping\n");
    const int NC = 60;
    for (const auto& cfg : {materials::Germanium(), materials::Silicon()})
        for (auto type : {ExcitationType::Electron, ExcitationType::Hole})
            for (double E : {100.0, 1000.0, 1.0e4}) {
                const std::string n = cfg.name + (type == ExcitationType::Electron ? " e- " : " h+ ") + std::to_string(int(E)) + " V/m";
                const R ref = drift(cfg, type, E, NC, false, 0.01);
                const R p5 = drift(cfg, type, E, NC, false, 0.05);
                const R lem = drift(cfg, type, E, NC, true, 0.25);
                audit::check(n + " Leman (E11) closure", lem.closure, 0.0, 0.0, 1e-3);
                audit::check(n + " Leman <E_Luke> vs fine first-order ref", lem.emean, ref.emean, 0.03);
                audit::check(n + " Leman phonons/carrier vs ref", lem.nph, ref.nph, 0.03);
                audit::info(n + " P5 <E_Luke> rel. deviation from ref", p5.emean / ref.emean - 1.0);
                audit::info(n + " steps/carrier: ref, P5, Leman", ref.steps);
                audit::info("   ", p5.steps);
                audit::info("   ", lem.steps);
                audit::check(n + " Leman uses fewer steps than P5", lem.steps < p5.steps ? 1.0 : 0.0, 1.0, 0.0);
            }
    // (d) transport_charges_fast_chord is the Leman path: identical output for the same seed
    {
        const CrystalLattice si(materials::Silicon());
        DetectorModel det; det.solid.radius = 5.0 * units::cm; det.solid.half_z = 0.5e-3;
        det.e_field = UniformElectricField{{0.0, 0.0, -1.0e4}};
        DetectorModel dl = det; dl.charge_stepping = DetectorModel::ChargeStepping::LemanSecondOrder;
        const SimulationEngine e1(si, det), e2(si, dl);
        const ChargePhysics cp(si);
        double s1 = 0, s2 = 0; size_t n1 = 0, n2 = 0;
        for (int pass = 0; pass < 2; ++pass) {
            Rng rng(99);
            std::vector<ChargeTrack> cs;
            for (int c = 0; c < 10; ++c) cs.push_back(cp.make_charge(ExcitationType::Hole, {0, 0, 0.49e-3}, {0, 0, -1}, 1e-6, 0));
            std::vector<PhononTrack> ph; std::vector<SensorHit> hits; SimStats st;
            if (pass == 0) e1.transport_charges_fast_chord(cs, ph, hits, st, rng); else e2.transport_charges(cs, ph, hits, st, rng);
            double s = 0; for (const auto& p : ph) s += p.energy;
            (pass == 0 ? s1 : s2) = s; (pass == 0 ? n1 : n2) = ph.size();
        }
        audit::check("fast_chord == Leman transport_charges: sum E_Luke", s1, s2, 0.0, 0.0);
        audit::check("fast_chord == Leman transport_charges: N phonons", double(n1), double(n2), 0.0, 0.0);
    }
    return audit::result();
}
