// Build check: every public header (including the SoA / wavefront engines that use
// KaplanSurrogateLUT) compiles together, and a short phonon run conserves energy.
#include "audit_common.hpp"
#include "cmpmc/engine.hpp"
#include "cmpmc/wavefront_soa.hpp"
#include "cmpmc/gpu/wavefront_pipeline.hpp"

using namespace cmpmc;

int main() {
    std::printf("# test_headers_compile\n");
    const CrystalLattice lat(materials::Germanium());
    DetectorModel det;
    const SimulationEngine eng(lat, det);
    Rng rng(1);
    std::vector<PhononTrack> st; std::vector<SensorHit> hits; SimStats S;
    double Ein = 0;
    for (int i = 0; i < 2000; ++i) { st.push_back(eng.make_track({0, 0, 0}, rng.isotropic_dir(), PhononMode::L, 0.004)); Ein += 0.004; }
    eng.transport(st, hits, S, rng);
    audit::check("Ge phonon run: deposited + escaped/killed = injected", S.total_deposited_eV + S.escaped_or_killed_eV, Ein, 1e-9);
    return audit::result();
}
