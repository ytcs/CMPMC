// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

// Multi-threaded event loop on top of SimulationEngine.
//
// run_events() simulates events 0 .. n-1. For each event it calls a user
// generator that fills the initial phonon and charge stacks, then transports the
// charges (their Luke phonons join the phonon stack) and the phonons.
//
// Results are reproducible and independent of the thread count:
//   - event i draws from its own stream Rng::for_stream(seed, i);
//   - events are processed in fixed-size chunks of consecutive events, and chunk
//     results are merged in chunk order, so statistics and the hit list are
//     bit-identical for 1, 2, ... N threads.
//
// The engine, lattice and detector are shared read-only. The generator is called
// concurrently from several threads, so it must not modify shared state without
// synchronisation; it should take all randomness from EventInput::rng.

#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include "cmpmc/transport/engine.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace cmpmc {

// What the generator fills for one event.
struct EventInput {
    uint64_t event_id;
    Rng& rng;                            // the event's random stream
    std::vector<PhononTrack>& phonons;   // initial phonons (empty on entry)
    std::vector<ChargeTrack>& charges;   // initial charge carriers (empty on entry)
    double injected_eV = 0.0;            // add the energy of every track you create
};

struct RunOptions {
    unsigned threads          = 0;    // 0: std::thread::hardware_concurrency()
    uint64_t seed             = 1;    // run seed; event i uses Rng::for_stream(seed, i)
    uint32_t events_per_chunk = 16;   // unit of work; part of the reproducibility contract
    bool keep_hits            = true; // false: accumulate statistics only (saves memory)
};

struct RunResult {
    SimStats stats;                      // summed over all events
    std::vector<SensorHit> hits;         // all hits, ordered by event (if keep_hits)
    std::vector<uint64_t> hit_offsets;   // hits of event i: [hit_offsets[i], hit_offsets[i+1])
    unsigned threads_used = 0;
    double wall_seconds   = 0.0;
};

namespace detail {

struct ChunkResult {
    SimStats stats;
    std::vector<SensorHit> hits;
    std::vector<uint64_t> hits_per_event;
};

inline unsigned resolve_threads(unsigned requested, uint64_t n_chunks) noexcept {
    unsigned n = requested ? requested : std::max(1u, std::thread::hardware_concurrency());
    return static_cast<unsigned>(std::max<uint64_t>(1, std::min<uint64_t>(n, n_chunks)));
}

} // namespace detail

template <class Generator>
RunResult run_events(const SimulationEngine& engine, uint64_t n_events, Generator&& generate,
                     const RunOptions& options = {}) {
    const auto t0 = std::chrono::steady_clock::now();
    const uint64_t chunk_size = std::max<uint32_t>(1, options.events_per_chunk);
    const uint64_t n_chunks = (n_events + chunk_size - 1) / chunk_size;
    const unsigned n_threads = detail::resolve_threads(options.threads, n_chunks);

    std::vector<detail::ChunkResult> chunks(n_chunks);
    std::atomic<uint64_t> next_chunk{0};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;

    auto worker = [&]() {
        std::vector<PhononTrack> phonons;
        std::vector<ChargeTrack> charges;
        try {
            for (;;) {
                const uint64_t c = next_chunk.fetch_add(1, std::memory_order_relaxed);
                if (c >= n_chunks || failed.load(std::memory_order_relaxed)) return;
                detail::ChunkResult& out = chunks[c];
                const uint64_t first = c * chunk_size;
                const uint64_t last = std::min(n_events, first + chunk_size);
                out.hits_per_event.reserve(last - first);
                std::vector<SensorHit> event_hits;
                for (uint64_t ev = first; ev < last; ++ev) {
                    Rng rng = Rng::for_stream(options.seed, ev);
                    phonons.clear();
                    charges.clear();
                    event_hits.clear();
                    EventInput input{ev, rng, phonons, charges, 0.0};
                    generate(input);
                    out.stats.total_injected_eV += input.injected_eV;
                    if (!charges.empty()) engine.transport_charges(charges, phonons, event_hits, out.stats, rng);
                    if (!phonons.empty()) engine.transport(phonons, event_hits, out.stats, rng);
                    out.hits_per_event.push_back(options.keep_hits ? event_hits.size() : 0);
                    if (options.keep_hits) out.hits.insert(out.hits.end(), event_hits.begin(), event_hits.end());
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!error) error = std::current_exception();
            failed.store(true, std::memory_order_relaxed);
        }
    };

    if (n_threads == 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_threads);
        for (unsigned t = 0; t < n_threads; ++t) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
    }
    if (error) std::rethrow_exception(error);

    RunResult result;
    result.threads_used = n_threads;
    size_t total_hits = 0;
    for (const auto& c : chunks) total_hits += c.hits.size();
    result.hits.reserve(total_hits);
    result.hit_offsets.reserve(n_events + 1);
    result.hit_offsets.push_back(0);
    for (auto& c : chunks) {
        result.stats += c.stats;
        for (uint64_t n : c.hits_per_event) result.hit_offsets.push_back(result.hit_offsets.back() + n);
        result.hits.insert(result.hits.end(), c.hits.begin(), c.hits.end());
        std::vector<SensorHit>().swap(c.hits);
    }
    result.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return result;
}

} // namespace cmpmc
