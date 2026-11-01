// Audit test: phonon surface reflection branching vs G4CMP
// (G4CMPSurfaceProperty::{Specular,Diffuse,Anharmonic}ReflProb, normalised in
// G4CMPPhononBoundaryProcess::DoReflection), using the examples/phonon coefficients
// (which G4CMP itself labels "not well-motivated ... only to demonstrate").
//  (a) the three probabilities at 100, 300, 400, 500, 600, 1000 GHz;
//  (b) engine-level: fraction of surface events that are surface down-conversions,
//      and energy conservation of that channel.
#include "audit_common.hpp"
#include "cmpmc/engine.hpp"

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
    const SurfaceProperty sp = SurfaceProperty::g4cmp_bare_wall_surface();
    for (double f : {100.0, 300.0, 400.0, 500.0, 600.0, 1000.0}) {
        const double E = f * units::GHz * units::h_Planck;
        const G4Split g = g4cmp(sp, 0.0, f);
        const std::string t = "f=" + std::to_string(int(f)) + " GHz ";
#ifdef CMPMC_AUDIT_REF
        const double ps = sp.specular_prob(E), pa = 0.0;   // legacy: spec = P_spec(min(f, cutoff)), no anharmonic branch
#else
        const auto sl = sp.reflection_split(E); const double ps = sl.p_spec, pa = sl.p_anh;
#endif
        audit::check(t + "p_spec vs G4CMP", ps, g.s, 1e-9, 1e-9);
        audit::check(t + "p_anh  vs G4CMP", pa, g.a, 1e-9, 1e-9);
    }
#ifndef CMPMC_AUDIT_REF
    // engine-level: bare box, no bulk scattering, all faces bare walls
    CrystalConfig c = materials::Silicon(); c.B_iso = 0.0; c.A_anh = 0.0;
    const CrystalLattice lat(c);
    DetectorModel det; det.max_bounces = 1; det.min_phonon_eV = 0.0;
    for (auto& s : det.surfaces) s = SurfaceProperty::g4cmp_bare_wall_surface();
    const SimulationEngine eng(lat, det);
    Rng rng(9);
    const double f = 450.0, E = f * units::GHz * units::h_Planck;
    std::vector<PhononTrack> st; std::vector<SensorHit> hits; SimStats S;
    const int N = 200000;
    for (int i = 0; i < N; ++i) st.push_back(eng.make_track({0, 0, 0}, rng.isotropic_dir(), PhononMode::FT, E));
    eng.transport(st, hits, S, rng);
    audit::check("engine: surface down-conversion fraction at 450 GHz", double(S.surface_downconversions) / N,
                 g4cmp(det.surfaces[0], 0.0, f).a, 0.02);
    audit::check("engine: energy accounted (killed/escaped == injected)", S.escaped_or_killed_eV, N * E, 1e-9);
#endif
    return audit::result();
}
