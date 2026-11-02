// R16 test: runtime guards on modules without a validated reference path (core/validation.hpp).
// Each guarded entry point must throw UnvalidatedModuleError by default and run when the caller
// opts in (DetectorModel::allow_unvalidated / the NEE argument). Validated paths must not throw.
#include "test_support.hpp"
#include "cmpmc/transport/engine.hpp"
#include "cmpmc/experimental/next_event_estimation.hpp"
#include "cmpmc/experimental/wavefront_soa.hpp"
#include "cmpmc/experimental/wavefront_pipeline.hpp"
#include <functional>
#include <string>

using namespace cmpmc;

static void expect(const std::string& name, bool should_throw, const std::function<void()>& f) {
    bool threw = false;
    try { f(); } catch (const UnvalidatedModuleError&) { threw = true; }
    cmpmc_test::check(name + (should_throw ? " throws by default" : " runs"), threw ? 1.0 : 0.0, should_throw ? 1.0 : 0.0, 0.0);
}

int main() {
    std::printf("# test_guards\n");
    const CrystalLattice ge(materials::Germanium());
    for (bool allow : {false, true}) {
        const std::string sfx = allow ? " [opt-in]" : "";
        DetectorModel det; det.allow_unvalidated = allow;
        Rng rng(3);
        expect("jump_to_ballistic" + sfx, !allow, [&] {
            DetectorModel d = det; d.jump_to_ballistic.enabled = true;
            SimulationEngine e(ge, d); std::vector<PhononTrack> st{e.make_track({0, 0, 0}, {0, 0, 1}, PhononMode::L, 1e-3)};
            std::vector<SensorHit> h; SimStats S; e.transport(st, h, S, rng); });
        expect("transport_charges_fast_chord (now Leman, validated)" + sfx, false, [&] {
            SimulationEngine e(ge, det); std::vector<ChargeTrack> c; std::vector<PhononTrack> p;
            std::vector<SensorHit> h; SimStats S; e.transport_charges_fast_chord(c, p, h, S, rng); });
        expect("transport_superconductor_coupled" + sfx, !allow, [&] {
            SimulationEngine e(ge, det); std::vector<PhononTrack> p; std::vector<GenericExcitationTrack> q;
            std::vector<SensorHit> h; SimStats S; e.transport_superconductor_coupled(p, q, h, S, rng); });
        expect("IV Matrix in transport_charges" + sfx, !allow, [&] {
            CrystalConfig c = materials::Silicon(); c.iv_model = IVRateModel::Matrix;
            const CrystalLattice si(c); SimulationEngine e(si, det); std::vector<ChargeTrack> cs; std::vector<PhononTrack> p;
            std::vector<SensorHit> h; SimStats S; e.transport_charges(cs, p, h, S, rng); });
        expect("experimental::WavefrontSoAEngine" + sfx, !allow, [&] { experimental::WavefrontSoAEngine e(ge, det); });
        expect("experimental::WavefrontPipelineEngine" + sfx, !allow, [&] { experimental::WavefrontPipelineEngine e(ge, det); });
        expect("evaluate_sparse_sensor_nee" + sfx, !allow, [&] {
            double t = 0; experimental::evaluate_sparse_sensor_nee({0, 0, -0.01}, {0, 0, 1}, {0, 0, 0.01}, {0, 0, -1},
                                                              1e-4, 1e-3, 1.0, t, 3500.0, allow); });
    }
    // validated paths do not throw
    DetectorModel det; Rng rng(4);
    expect("transport (default surfaces, no shortcuts)", false, [&] {
        SimulationEngine e(ge, det); std::vector<PhononTrack> st{e.make_track({0, 0, 0}, {0, 0, 1}, PhononMode::L, 1e-3)};
        std::vector<SensorHit> h; SimStats S; e.transport(st, h, S, rng); });
    expect("transport_charges (Linear IV)", false, [&] {
        SimulationEngine e(ge, det); std::vector<ChargeTrack> c; std::vector<PhononTrack> p;
        std::vector<SensorHit> h; SimStats S; e.transport_charges(c, p, h, S, rng); });
    return cmpmc_test::result();
}
