#include "cmpmc/engine.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

using namespace cmpmc;

// ============================================================================
// 1. Phonon Caustics Benchmark (Replicates G4CMP examples/caustics)
// ============================================================================
static void run_caustics_demo(const CrystalLattice& ge_lattice) {
    std::cout << "========================================================================\n";
    std::cout << " [1] Anisotropic Phonon Focusing Caustics (Ge [100] Crystal)\n";
    std::cout << "========================================================================\n";

    constexpr int NX = 61;
    constexpr int NY = 31;
    std::array<std::array<int, NX>, NY> grid{};

    constexpr uint64_t N_RAYS = 2'000'000;
    Rng rng(42);

    const auto t0 = std::chrono::steady_clock::now();
    for (uint64_t i = 0; i < N_RAYS; ++i) {
        // Sample wavevector in upper hemisphere
        Vec3 k = rng.isotropic_dir();
        if (k.z < 0.0) k.z = -k.z;

        const PhononMode mode = (i % 2 == 0) ? PhononMode::ST : PhononMode::FT;
        const auto kin = ge_lattice.map_k_to_vg(mode, k);
        if (kin.vg_dir.z <= 0.2) continue;

        // Project onto z = 1 plane, x/z and y/z in [-1.0, 1.0]
        const double px = kin.vg_dir.x / kin.vg_dir.z;
        const double py = kin.vg_dir.y / kin.vg_dir.z;
        if (std::abs(px) < 1.0 && std::abs(py) < 1.0) {
            const int ix = static_cast<int>((px + 1.0) * 0.5 * (NX - 1));
            const int iy = static_cast<int>((py + 1.0) * 0.5 * (NY - 1));
            grid[iy][ix]++;
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double dt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    int max_count = 1;
    for (const auto& row : grid)
        for (int c : row)
            max_count = std::max(max_count, c);

    const char* palette = " .:-=+*#%@";
    std::cout << "  Projected Group Velocity Intensity on (001) Face (FT + ST modes):\n";
    for (int y = 0; y < NY; ++y) {
        std::cout << "   ";
        for (int x = 0; x < NX; ++x) {
            const double norm = std::sqrt(static_cast<double>(grid[y][x]) / max_count);
            const int idx = std::clamp(static_cast<int>(norm * 9.99), 0, 9);
            std::cout << palette[idx];
        }
        std::cout << "\n";
    }
    std::cout << std::fixed << std::setprecision(2)
              << "  -> Mapped " << (N_RAYS / 1e6) << "M anisotropic wavevectors in "
              << dt_ms << " ms (" << (N_RAYS / (dt_ms * 1e3)) << " M queries/s)\n\n";
}

// ============================================================================
// 2. Full Detector Phonon + Kaplan QP Transport (Replicates G4CMP examples/phonon)
// ============================================================================
static void run_detector_benchmark(const CrystalLattice& ge_lattice, uint32_t num_primaries, unsigned num_threads) {
    std::cout << "========================================================================\n";
    std::cout << " [2] Full Ge Cylinder + Al QET Sensor Cascade (G4CMP examples/phonon)\n";
    std::cout << "========================================================================\n";

    DetectorModel detector{};
    SimulationEngine engine(ge_lattice, detector);

    // Primary Debye-energy (2 THz = 8.27 meV) Longitudinal phonons at center (0, 0, 0)
    // Each 2 THz L-phonon undergoes rapid anharmonic downconversion into a cascade of ST/FT phonons.
    const double primary_energy = ge_lattice.debye_energy();

    std::cout << "  Crystal          : Germanium [100] (R = 3.81 cm, H = 2.54 cm)\n"
              << "  Sensors          : Top & Bottom Al QETs (600 nm, 2*Delta = 347.4 ueV, 20% cov)\n"
              << "  Primary Phonons  : " << num_primaries << " L-phonons @ "
              << (primary_energy / units::meV) << " meV (2.0 THz Debye energy)\n"
              << "  Worker Threads   : " << num_threads << "\n\n";

    std::vector<SimStats> thread_stats(num_threads);
    std::vector<std::vector<SensorHit>> thread_hits(num_threads);
    std::vector<std::thread> workers;

    const auto t0 = std::chrono::steady_clock::now();

    const uint32_t base_count = num_primaries / num_threads;
    const uint32_t remainder  = num_primaries % num_threads;

    for (unsigned tid = 0; tid < num_threads; ++tid) {
        const uint32_t count = base_count + (tid < remainder ? 1 : 0);
        workers.emplace_back([&, tid, count]() {
            Rng rng(10000ULL + tid * 7919ULL);
            std::vector<PhononTrack> stack;
            stack.reserve(256);
            thread_hits[tid].reserve(count * 8);

            for (uint32_t i = 0; i < count; ++i) {
                const Vec3 k_dir = rng.isotropic_dir();
                stack.push_back(engine.make_track({0.0, 0.0, 0.0}, k_dir, PhononMode::L, primary_energy, 0.0));
                thread_stats[tid].total_injected_eV += primary_energy;
                engine.transport(stack, thread_hits[tid], thread_stats[tid], rng);
            }
        });
    }

    for (auto& w : workers) w.join();

    const auto t1 = std::chrono::steady_clock::now();
    const double dt_s = std::chrono::duration<double>(t1 - t0).count();

    SimStats total{};
    std::vector<SensorHit> all_hits;
    for (unsigned tid = 0; tid < num_threads; ++tid) {
        total += thread_stats[tid];
        all_hits.insert(all_hits.end(), thread_hits[tid].begin(), thread_hits[tid].end());
    }

    const uint64_t total_decays = total.decays_L_to_TT + total.decays_L_to_LT;
    const double energy_balance_err = std::abs(
        total.total_injected_eV - (total.total_deposited_eV + total.escaped_or_killed_eV)
    ) / total.total_injected_eV;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  --- Performance Summary ---\n"
              << "  Wall-clock time          : " << (dt_s * 1e3) << " ms\n"
              << "  Total Ray Steps          : " << total.total_steps
              << " (" << (total.total_steps / (dt_s * 1e6)) << " M steps/s)\n"
              << "  Primary + Decay Tracks   : " << (num_primaries + total_decays)
              << " (" << ((num_primaries + total_decays) / (dt_s * 1e3)) << " k tracks/s)\n\n";

    std::cout << "  --- Condensed Matter Physics Summary ---\n"
              << "  Anharmonic Decays (L->TT): " << total.decays_L_to_TT
              << " (" << (100.0 * total.decays_L_to_TT / std::max<uint64_t>(1, total_decays)) << "%)\n"
              << "  Anharmonic Decays (L->LT): " << total.decays_L_to_LT
              << " (" << (100.0 * total.decays_L_to_LT / std::max<uint64_t>(1, total_decays)) << "%)\n"
              << "  Isotope Scatterings      : " << total.isotope_scatters << "\n"
              << "  Boundary Reflections     : " << total.surface_refls << "\n"
              << "  Al QET Sensor Hits       : " << total.sensor_hits << "\n"
              << "  Total Injected Energy    : " << (total.total_injected_eV / units::eV) << " eV\n"
              << "  Total QET Deposited Edep : " << (total.total_deposited_eV / units::eV) << " eV ("
              << (100.0 * total.total_deposited_eV / total.total_injected_eV) << "%)\n"
              << "  Sub-gap / Escaped Energy : " << (total.escaped_or_killed_eV / units::eV) << " eV ("
              << (100.0 * total.escaped_or_killed_eV / total.total_injected_eV) << "%)\n"
              << "  Energy Conservation Diff : " << std::scientific << energy_balance_err << std::fixed << "\n\n";

    // Build arrival time histogram (0 to 80 us in 20 bins of 4 us)
    constexpr int N_BINS = 20;
    constexpr double MAX_TIME_US = 80.0;
    std::array<double, N_BINS> time_hist{};
    uint64_t mode_counts[3] = {0, 0, 0};

    for (const auto& h : all_hits) {
        mode_counts[static_cast<int>(h.mode)]++;
        const double t_us = h.time / units::us;
        const int bin = static_cast<int>((t_us / MAX_TIME_US) * N_BINS);
        if (bin >= 0 && bin < N_BINS) {
            time_hist[bin] += h.edep;
        }
    }

    double max_bin_e = 1e-30;
    for (double e : time_hist) max_bin_e = std::max(max_bin_e, e);

    std::cout << "  --- Sensor Hit Polarization Breakdown ---\n"
              << "  L  : " << (100.0 * mode_counts[0] / std::max<size_t>(1, all_hits.size())) << "%  (Target LDOS  ~ 9.8%)\n"
              << "  FT : " << (100.0 * mode_counts[1] / std::max<size_t>(1, all_hits.size())) << "%  (Target FTDOS ~ 36.7%)\n"
              << "  ST : " << (100.0 * mode_counts[2] / std::max<size_t>(1, all_hits.size())) << "%  (Target STDOS ~ 53.5%)\n\n";

    std::cout << "  --- QET Sensor Arrival Time Pulse Profile (0 - 80 us) ---\n";
    for (int b = 0; b < N_BINS; ++b) {
        const double t_start = b * (MAX_TIME_US / N_BINS);
        const double t_end   = (b + 1) * (MAX_TIME_US / N_BINS);
        const int bar_len = static_cast<int>(40.0 * (time_hist[b] / max_bin_e));
        std::cout << "  " << std::setw(4) << static_cast<int>(t_start) << "-"
                  << std::setw(2) << static_cast<int>(t_end) << " us | "
                  << std::string(bar_len, '#') << "\n";
    }
    std::cout << "========================================================================\n";
}

int main(int argc, char** argv) {
    uint32_t num_primaries = 5000;
    if (argc > 1) {
        num_primaries = static_cast<uint32_t>(std::stoul(argv[1]));
    }
    const unsigned num_threads = std::max(1u, std::thread::hardware_concurrency());

    const auto t_init0 = std::chrono::steady_clock::now();
    CrystalLattice ge_lattice(CrystalLattice::germanium_config());
    const auto t_init1 = std::chrono::steady_clock::now();
    const double init_ms = std::chrono::duration<double, std::milli>(t_init1 - t_init0).count();

    std::cout << "Initialized Germanium Christoffel 2D LUT (181x361x3) + Decay Inverse-CDF in "
              << std::fixed << std::setprecision(2) << init_ms << " ms.\n\n";

    run_caustics_demo(ge_lattice);
    run_detector_benchmark(ge_lattice, num_primaries, num_threads);

    return 0;
}
