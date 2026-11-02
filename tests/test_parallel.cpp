// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// run_events(): reproducibility across thread counts, agreement with a serial
// loop, energy closure, hit bookkeeping, and propagation of worker exceptions.
#include "test_support.hpp"
#include "cmpmc/transport/parallel.hpp"

#include <cstring>
#include <stdexcept>

using namespace cmpmc;

namespace {

// Even events: one 7.5 meV phonon. Odd events: one e-/h+ pair in the field.
struct MixedGenerator {
    const SimulationEngine* engine;
    const ChargePhysics* charges;
    size_t n_valleys;
    void operator()(EventInput& in) const {
        if (in.event_id % 2 == 0) {
            in.phonons.push_back(engine->make_track({0, 0, 0}, in.rng.isotropic_dir(), PhononMode::FT, 7.5e-3));
            in.injected_eV += 7.5e-3;
        } else {
            const auto v = static_cast<uint16_t>(in.rng.uniform() * n_valleys) % n_valleys;
            in.charges.push_back(charges->make_charge(ExcitationType::Electron, {0, 0, 0}, in.rng.isotropic_dir(), 1e-6, v));
            in.charges.push_back(charges->make_charge(ExcitationType::Hole, {0, 0, 0}, in.rng.isotropic_dir(), 1e-6, 0));
        }
    }
};

bool same_stats(const SimStats& a, const SimStats& b) {
    return std::memcmp(&a, &b, sizeof(SimStats)) == 0;
}

bool same_hits(const RunResult& a, const RunResult& b) {
    if (a.hits.size() != b.hits.size() || a.hit_offsets != b.hit_offsets) return false;
    for (size_t i = 0; i < a.hits.size(); ++i) {
        const SensorHit &x = a.hits[i], &y = b.hits[i];
        if (x.edep != y.edep || x.time != y.time || x.pos.x != y.pos.x || x.pos.y != y.pos.y ||
            x.pos.z != y.pos.z || x.surface_id != y.surface_id)
            return false;
    }
    return true;
}

} // namespace

int main() {
    std::printf("# test_parallel\n");
    static_assert(sizeof(SimStats) % sizeof(double) == 0, "SimStats is compared bytewise");

    const CrystalLattice lattice(materials::Germanium());
    DetectorModel det;
    det.e_field = UniformElectricField::from_bias(2.0, 2.54 * units::cm);
    det.luke_sample_fraction = 0.02;   // weighted Luke phonons keep the test fast
    const SimulationEngine engine(lattice, det);
    const ChargePhysics charges(lattice);
    const MixedGenerator gen{&engine, &charges, lattice.config().valley_dirs.size()};

    const uint64_t N = 203;  // not a multiple of the chunk size
    RunOptions opt;
    opt.seed = 77;

    opt.threads = 1;
    const RunResult r1 = run_events(engine, N, gen, opt);
    cmpmc_test::check("threads used (requested 1)", r1.threads_used, 1, 0.0);

    for (unsigned t : {2u, 3u, 4u}) {
        opt.threads = t;
        const RunResult rt = run_events(engine, N, gen, opt);
        const std::string n = std::to_string(t) + " threads: ";
        cmpmc_test::check(n + "threads used", rt.threads_used, t, 0.0);
        cmpmc_test::check(n + "statistics bit-identical to 1 thread", same_stats(r1.stats, rt.stats) ? 1 : 0, 1, 0.0);
        cmpmc_test::check(n + "hits identical to 1 thread (order, values)", same_hits(r1, rt) ? 1 : 0, 1, 0.0);
    }

    // Same physics as a plain serial loop over the same per-event streams.
    SimStats serial;
    std::vector<SensorHit> serial_hits;
    {
        std::vector<PhononTrack> ph;
        std::vector<ChargeTrack> ch;
        for (uint64_t ev = 0; ev < N; ++ev) {
            Rng rng = Rng::for_stream(opt.seed, ev);
            ph.clear(); ch.clear();
            EventInput in{ev, rng, ph, ch, 0.0};
            gen(in);
            serial.total_injected_eV += in.injected_eV;
            if (!ch.empty()) engine.transport_charges(ch, ph, serial_hits, serial, rng);
            if (!ph.empty()) engine.transport(ph, serial_hits, serial, rng);
        }
    }
    cmpmc_test::check("vs serial loop: total steps", static_cast<double>(r1.stats.total_steps),
                      static_cast<double>(serial.total_steps), 0.0);
    cmpmc_test::check("vs serial loop: charge steps", static_cast<double>(r1.stats.charge_steps),
                      static_cast<double>(serial.charge_steps), 0.0);
    cmpmc_test::check("vs serial loop: hit count", static_cast<double>(r1.hits.size()),
                      static_cast<double>(serial_hits.size()), 0.0);
    cmpmc_test::check("vs serial loop: deposited energy", r1.stats.total_deposited_eV, serial.total_deposited_eV, 1e-12);

    // Bookkeeping.
    cmpmc_test::check("hit_offsets has n + 1 entries", static_cast<double>(r1.hit_offsets.size()),
                      static_cast<double>(N + 1), 0.0);
    cmpmc_test::check("last hit offset = number of hits", static_cast<double>(r1.hit_offsets.back()),
                      static_cast<double>(r1.hits.size()), 0.0);

    // Phonon-only energy closure through run_events.
    {
        RunOptions o; o.threads = 4; o.seed = 5;
        auto phonon_gen = [&](EventInput& in) {
            in.phonons.push_back(engine.make_track({0, 0, 0}, in.rng.isotropic_dir(), PhononMode::L, 4e-3));
            in.injected_eV += 4e-3;
        };
        const RunResult r = run_events(engine, 500, phonon_gen, o);
        cmpmc_test::check("phonon run: deposited + escaped = injected",
                          r.stats.total_deposited_eV + r.stats.escaped_or_killed_eV, r.stats.total_injected_eV, 1e-9);
        o.keep_hits = false;
        const RunResult rn = run_events(engine, 500, phonon_gen, o);
        cmpmc_test::check("keep_hits = false: no hits stored", static_cast<double>(rn.hits.size()), 0.0, 0.0);
        cmpmc_test::check("keep_hits = false: same statistics", same_stats(r.stats, rn.stats) ? 1 : 0, 1, 0.0);
    }

    // An exception in a worker thread reaches the caller.
    {
        CrystalConfig cfg = materials::Germanium();
        cfg.iv_model = IVRateModel::Matrix;   // guarded: throws without allow_unvalidated
        const CrystalLattice guarded_lattice(cfg);
        const SimulationEngine guarded(guarded_lattice, det);
        const ChargePhysics gc(guarded_lattice);
        const MixedGenerator g{&guarded, &gc, cfg.valley_dirs.size()};
        RunOptions o; o.threads = 4;
        bool caught = false;
        try { (void)run_events(guarded, 64, g, o); } catch (const UnvalidatedModuleError&) { caught = true; }
        cmpmc_test::check("worker exception rethrown on caller", caught ? 1 : 0, 1, 0.0);
    }
    return cmpmc_test::result();
}
