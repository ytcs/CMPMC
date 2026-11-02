// Audit test: phonon surface reflection branching vs G4CMP
// (G4CMPSurfaceProperty::{Specular,Diffuse,Anharmonic}ReflProb, normalised in
// G4CMPPhononBoundaryProcess::DoReflection), using the examples/phonon coefficients
// (which G4CMP itself labels "not well-motivated ... only to demonstrate").
//  (a) the three probabilities at 100, 300, 400, 500, 600, 1000 GHz;
//  (b) engine-level: fraction of surface events that are surface down-conversions,
//      and energy conservation of that channel;
//  (c) R16: default surfaces are diffuse with p_anh = 0, and the demo set is guarded.
#include "test_support.hpp"
#include "cmpmc/transport/engine.hpp"

using namespace cmpmc;

static double poly(const std::vector<double>& c, double f) { double a = 0; for (size_t i = c.size(); i > 0; --i) a = a * f + c[i - 1]; return a; }

struct G4Split { double s, d, a; };
static G4Split g4cmp(const SurfaceProperty& sp, double specProb, double f) {
    const double anh = (f > sp.anh_cutoff) ? 0.0 : poly(sp.anh_coeffs, f);
    const double diff = (f > sp.anh_cutoff) ? 1.0 - specProb : poly(sp.diff_coeffs, std::min(f, sp.refl_cutoff));
    const double spec = (f > sp.refl_cutoff) ? 1.0 - diff - anh : poly(sp.spec_coeffs, f);
    const double n = spec + diff + anh; return {spec / n, diff / n, anh / n};
}

int main() {
    std::printf("# test_surface\n");
    const SurfaceProperty sp = SurfaceProperty::demo_polynomial_bare_surface();
    for (double f : {100.0, 300.0, 400.0, 500.0, 600.0, 1000.0}) {
        const double E = f * units::GHz * units::h_Planck;
        const G4Split g = g4cmp(sp, 0.0, f);
        const std::string t = "f=" + std::to_string(int(f)) + " GHz ";
#ifdef CMPMC_AUDIT_REF
        const double ps = sp.specular_prob(E), pa = 0.0;   // legacy: spec = P_spec(min(f, cutoff)), no anharmonic branch
#else
        const auto sl = sp.reflection_split(E); const double ps = sl.p_spec, pa = sl.p_anh;
#endif
        cmpmc_test::check(t + "p_spec vs G4CMP", ps, g.s, 1e-9, 1e-9);
        cmpmc_test::check(t + "p_anh  vs G4CMP", pa, g.a, 1e-9, 1e-9);
    }
#ifndef CMPMC_AUDIT_REF
    // engine-level: bare box, no bulk scattering, all faces bare walls
    CrystalConfig c = materials::Silicon(); c.B_iso = 0.0; c.A_anh = 0.0;
    const CrystalLattice lat(c);
    DetectorModel det; det.max_bounces = 1; det.min_phonon_eV = 0.0;
    det.allow_unvalidated = true;   // the demo coefficients are guarded (R16)
    for (auto& s : det.surfaces) s = SurfaceProperty::demo_polynomial_bare_surface();
    const SimulationEngine eng(lat, det);
    Rng rng(9);
    const double f = 450.0, E = f * units::GHz * units::h_Planck;
    std::vector<PhononTrack> st; std::vector<SensorHit> hits; SimStats S;
    const int N = 200000;
    for (int i = 0; i < N; ++i) st.push_back(eng.make_track({0, 0, 0}, rng.isotropic_dir(), PhononMode::FT, E));
    eng.transport(st, hits, S, rng);
    cmpmc_test::check("engine: surface down-conversion fraction at 450 GHz", double(S.surface_downconversions) / N,
                 g4cmp(det.surfaces[0], 0.0, f).a, 0.02);
    cmpmc_test::check("engine: energy accounted (killed/escaped == injected)", S.escaped_or_killed_eV, N * E, 1e-9);
    // (c) R16 campaign defaults: fully diffuse, no surface down-conversion (1109.1193 Sec. IV.G [FT]),
    //     and the demo coefficients are refused without opt-in
    const DetectorModel dflt;
    for (double fq : {100.0, 300.0, 450.0, 600.0, 2000.0}) {
        for (size_t k = 0; k < dflt.surfaces.size(); ++k) {
            const auto sl = dflt.surfaces[k].reflection_split(fq * units::GHz * units::h_Planck);
            cmpmc_test::check("default surface " + std::to_string(k) + " f=" + std::to_string(int(fq)) + " p_diff", sl.p_diff, 1.0, 0.0, 1e-15);
            cmpmc_test::check("default surface " + std::to_string(k) + " f=" + std::to_string(int(fq)) + " p_anh", sl.p_anh, 0.0, 0.0, 0.0);
        }
    }
    bool threw = false;
    try { DetectorModel d2; d2.surfaces[0] = SurfaceProperty::demo_polynomial_sensor_surface(); const SimulationEngine e2(lat, d2); }
    catch (const UnvalidatedModuleError&) { threw = true; }
    cmpmc_test::check("demo surface refused without allow_unvalidated", threw ? 1.0 : 0.0, 1.0, 0.0);
#endif
    return cmpmc_test::result();
}
