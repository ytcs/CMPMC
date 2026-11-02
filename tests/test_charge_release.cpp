// Energy released at carrier collection (DetectorModel::charge_release).
//  Energy bookkeeping at collection: the pair energy is converted to heat at the contact,
//  E_rel = E_kin + (E_c - mu_F) for e-, E_kin + (mu_F - E_v) for h+.
//  (a) closed form: e- + h+ at rest release E_g exactly, for any mu_F; mid-gap default;
//  (b) per-mode energy closure in the engine: sum(hit edep) + sum(emitted phonons) + lost
//      = sum(E_rel), to 1e-12 (RecordOnly, BulkDebye, ElectrodeFilm);
//  (c) BulkDebye emits phonons of h nu_D (Debye frequency from the material preset);
//  (d) ElectrodeFilm: the re-injected population is soft (born below a few Delta_Al; Leman,
//      arXiv:1109.1193, p. 17: "low frequency, ballistic"): mean energy and the fraction below 2 Delta_Al and
//      4 Delta_Al are reported, the mean is checked to be < 4 Delta_Al, i.e. >= 100x softer than the
//      bulk release (h nu_D(Ge) = 35.7 meV); the absorbed fraction is the film's eta_pb;
//  (e) numerical seed-quantum size (80 vs 160 Delta; the hot-electron cascade feeds phonons up to
//      h nu_D(Al) ~ 200 Delta): absorbed fraction and mean re-emitted energy agree to < 2%.
#include "test_support.hpp"
#include "cmpmc/transport/engine.hpp"
#include <string>

using namespace cmpmc;

struct Out { double erel = 0, hit = 0, ph = 0, lost = 0, n_ph = 0, below2 = 0, below4 = 0, emean = 0; };

static Out run(const CrystalLattice& lat, DetectorModel det, ExcitationType type, int N, uint64_t seed) {
    det.e_field = UniformElectricField::from_bias(0.0, 2.0 * det.solid.half_z);
    const SimulationEngine eng(lat, det);
    const ChargePhysics cp(lat);
    Rng rng(seed);
    std::vector<ChargeTrack> cs;
    for (int i = 0; i < N; ++i)
        cs.push_back(cp.make_charge(type, {0, 0, det.solid.half_z - 1e-6}, {0, 0, 1}, 1e-6, 0, 0.0));
    std::vector<PhononTrack> ph; std::vector<SensorHit> hits; SimStats S;
    eng.transport_charges(cs, ph, hits, S, rng);
    Out o;
    for (const auto& h : hits) { o.erel += h.incident_energy; o.hit += h.edep; }
    const double D = det.surfaces[0].kaplan_sensor->gap_energy;
    for (const auto& p : ph) { o.ph += p.energy * p.weight; o.n_ph += 1; if (p.energy < 2 * D) o.below2++; if (p.energy < 4 * D) o.below4++; }
    o.lost = S.escaped_or_killed_eV;
    o.emean = o.n_ph > 0 ? o.ph / o.n_ph : 0.0;
    o.below2 /= std::max(1.0, o.n_ph); o.below4 /= std::max(1.0, o.n_ph);
    return o;
}

int main() {
    std::printf("# test_charge_release\n");
    const CrystalLattice ge(materials::Germanium());
    const double Eg = ge.config().bandgap, hnuD = ge.debye_energy();
    DetectorModel det;
    det.solid.half_z = 1e-3; det.solid.radius = 1e-2;
    const double D = det.surfaces[0].kaplan_sensor->gap_energy;
    const int N = 1000;

    // (a) + (b) RecordOnly, with mu_F mid-gap and at E_v + 0.2 eV
    for (double muF : {-1.0, 0.2}) {
        DetectorModel d = det; d.fermi_level_above_Ev_eV = muF;
        const Out e = run(ge, d, ExcitationType::Electron, N, 1), h = run(ge, d, ExcitationType::Hole, N, 2);
        const double mu = muF < 0 ? 0.5 * Eg : muF;
        cmpmc_test::check("RecordOnly mu_F=" + std::to_string(mu).substr(0, 5) + ": (E_rel,e + E_rel,h)/N vs E_g", (e.erel + h.erel) / N, Eg, 1e-4);
        cmpmc_test::check("  e- E_rel/N vs E_c - mu_F", e.erel / N, Eg - mu, 1e-4);
        cmpmc_test::check("  closure hit = E_rel (no phonons)", e.hit + h.hit, e.erel + h.erel, 1e-12);
    }
    // (b) + (c) BulkDebye
    {
        DetectorModel d = det; d.charge_release = DetectorModel::ChargeRelease::BulkDebye;
        const Out o = run(ge, d, ExcitationType::Electron, N, 3);
        cmpmc_test::check("BulkDebye closure hit + phonons + lost = E_rel", o.hit + o.ph + o.lost, o.erel, 1e-12);
        cmpmc_test::check("BulkDebye phonon energy vs E_rel/round(E_rel/h nu_D)", o.emean,
                     (o.erel / N) / std::max(1.0, std::round((o.erel / N) / hnuD)), 1e-9);
        cmpmc_test::info("BulkDebye phonons per carrier", o.n_ph / N);
    }
    // (b) + (d) + (e) ElectrodeFilm, mean-free-path and rate-competition film loss
    for (int lm = 0; lm < 2; ++lm) {
        DetectorModel d = det; d.charge_release = DetectorModel::ChargeRelease::ElectrodeFilm;
        if (lm == 1) for (auto& sp : d.surfaces) if (sp.kaplan_sensor) {
            sp.kaplan_sensor->loss_model = FilmLossModel::KaplanRateCompetition;
            sp.kaplan_sensor->phonon_loss_time = KaplanFilm::loss_time_from_geometry(
                sp.kaplan_sensor->film_thickness, sp.kaplan_sensor->v_sound, 0.94);   // Ge/Al AMM normal-incidence L transmission (test_interface_amm)
        }
        const std::string tag = lm ? "ElectrodeFilm [Kaplan rate competition]" : "ElectrodeFilm [mean-free-path loss law]";
        Out o[2];
        for (int k = 0; k < 2; ++k) {
            DetectorModel dk = d; dk.electrode_seed_max_gap = k ? 160.0 : 80.0;
            for (auto t : {ExcitationType::Electron, ExcitationType::Hole}) {
                const Out r = run(ge, dk, t, N, 10 + k + (t == ExcitationType::Hole ? 5 : 0));
                o[k].erel += r.erel; o[k].hit += r.hit; o[k].ph += r.ph; o[k].lost += r.lost; o[k].n_ph += r.n_ph;
                o[k].below2 += r.below2 * r.n_ph; o[k].below4 += r.below4 * r.n_ph;
            }
            o[k].emean = o[k].ph / o[k].n_ph; o[k].below2 /= o[k].n_ph; o[k].below4 /= o[k].n_ph;
        }
        cmpmc_test::check(tag + " closure hit + phonons + lost = E_rel", o[0].hit + o[0].ph + o[0].lost, o[0].erel, 1e-12);
        cmpmc_test::info(tag + " absorbed fraction (film eta)", o[0].hit / o[0].erel);
        cmpmc_test::info(tag + " re-injected fraction", o[0].ph / o[0].erel);
        cmpmc_test::info(tag + " mean re-injected phonon energy / Delta_Al", o[0].emean / D);
        cmpmc_test::info(tag + " fraction of phonons < 2 Delta_Al", o[0].below2);
        cmpmc_test::info(tag + " fraction of phonons < 4 Delta_Al", o[0].below4);
        cmpmc_test::check(tag + " soft: mean re-injected energy < 4 Delta_Al", o[0].emean < 4 * D ? 1.0 : 0.0, 1.0, 0.0);
        cmpmc_test::info(tag + " h nu_D(Ge) / mean re-injected energy", hnuD / o[0].emean);
        cmpmc_test::check(tag + " seed 160D vs 80D: absorbed fraction", o[1].hit / o[1].erel, o[0].hit / o[0].erel, 0.02);
        cmpmc_test::check(tag + " seed 160D vs 80D: mean re-injected energy", o[1].emean, o[0].emean, 0.02);
    }
    return cmpmc_test::result();
}
