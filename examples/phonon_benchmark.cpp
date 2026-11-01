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
        Vec3 k = rng.isotropic_dir();
        if (k.z < 0.0) k.z = -k.z;

        const PhononMode mode = (i % 2 == 0) ? PhononMode::ST : PhononMode::FT;
        const auto kin = ge_lattice.map_k_to_vg(mode, k);
        if (kin.vg_dir.z <= 0.2) continue;

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
}

// ============================================================================
// 3. End-to-End Nuclear Recoil Partition + Charge Drift + NTL Phonon Cascade
//    (Replicates G4CMPEnergyPartition + G4CMPVDriftProcess + G4CMPLukeScattering)
// ============================================================================
static void run_recoil_charge_luke_demo(const CrystalLattice& ge_lattice) {
    std::cout << "========================================================================\n";
    std::cout << " [3] End-to-End Nuclear Recoil + Herring-Vogt Charge Drift + Luke Gain\n";
    std::cout << "========================================================================\n";

    constexpr double BIAS_VOLTS = 4.0 * units::Volt;
    constexpr double RECOIL_EV  = 100.0 * units::eV;

    DetectorModel detector{};
    detector.e_field = UniformElectricField::from_bias(BIAS_VOLTS, 2.0 * detector.solid.half_z);
    detector.luke_sample_fraction = 0.05; // Track 5% of Luke phonons with weight = 20x

    SimulationEngine engine(ge_lattice, detector);
    EnergyPartition partitioner(ge_lattice);

    Rng rng(1337);
    std::vector<ChargeTrack> charges;
    std::vector<PhononTrack> phonons;
    std::vector<SensorHit> hits;
    SimStats stats{};

    const auto t0 = std::chrono::steady_clock::now();

    // 1. Partition a 100 eV Nuclear Recoil at the center of the Ge crystal
    PartitionOptions part_opts{};
    part_opts.max_sampled_phonons = 1000; // Downsample prompt phonon cloud to 1000 tracks (with exact weight)
    const auto part = partitioner.partition_deposit(
        {0.0, 0.0, 0.0}, 0.0, RECOIL_EV, /*is_nuclear_recoil=*/true, charges, phonons, rng, part_opts);

    // 2. Drift e-/h+ pairs in the 4V bias field (Herring-Vogt valleys + IV scatter + Luke emission)
    engine.transport_charges(charges, phonons, hits, stats, rng);

    // 3. Transport all prompt + Luke-emitted phonons through the crystal to the Al QET sensors
    engine.transport(phonons, hits, stats, rng);

    const auto t1 = std::chrono::steady_clock::now();
    const double dt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Expected total calorimeter energy = E_recoil + N_pairs * e * V_bias
    const double expected_luke_work_eV = static_cast<double>(part.n_pairs_true) * BIAS_VOLTS;
    const double expected_total_eV     = RECOIL_EV + expected_luke_work_eV;
    const double measured_total_eV     = stats.total_deposited_eV + stats.escaped_or_killed_eV;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Recoil Event             : " << RECOIL_EV << " eV Nuclear Recoil in Ge @ "
              << BIAS_VOLTS << " V Bias\n"
              << "  Lindhard Ionization Yield: " << (100.0 * part.lindhard_yield) << "% ("
              << part.ionizing_eV << " eV ionizing, " << part.niel_eV << " eV NIEL)\n"
              << "  Generated e-/h+ Pairs    : " << part.n_pairs_true << " pairs ("
              << (part.n_pairs_true * 2) << " carriers across 4 L-valleys)\n"
              << "  Prompt Recoil Phonons    : " << part.n_phonons_true << " phonons ("
              << part.prompt_phonon_eV << " eV)\n"
              << "  Charge Drift Steps       : " << stats.charge_steps << " ("
              << stats.intervalley_scatters << " intervalley scatters)\n"
              << "  Luke Phonons Emitted     : " << stats.luke_phonons_emitted
              << " (Total Luke energy = " << stats.total_luke_phonon_eV << " eV)\n"
              << "  Expected NTL Work (N*eV) : " << expected_luke_work_eV << " eV\n"
              << "  Electrode Charge Recomb  : " << stats.charge_recomb_edep_eV << " eV\n"
              << "  Total Calorimeter Energy : " << measured_total_eV << " eV (Expected "
              << expected_total_eV << " eV)\n"
              << "  End-to-End Runtime       : " << dt_ms << " ms\n";
    std::cout << "========================================================================\n\n";
}

// ============================================================================
// 4. Multi-Symmetry Crystal Presets & Polycrystalline Scattering (Burst 1)
// ============================================================================
static void run_crystal_symmetries_demo() {
    std::cout << "========================================================================\n";
    std::cout << " [4] All-Bravais Crystal Symmetries & Polycrystal Scattering (Burst 1)\n";
    std::cout << "========================================================================\n";

    CrystalConfig al2o3_cfg = materials::Sapphire();
    al2o3_cfg.set_miller_orientation(1, 0, 0);
    CrystalLattice al2o3(al2o3_cfg);

    CrystalConfig cawo4_cfg = materials::CalciumTungstate();
    CrystalLattice cawo4(cawo4_cfg);

    CrystalConfig nb_cfg = materials::NiobiumCrystal();
    nb_cfg.polycrystal_mfp_m = 150.0 * units::um; // Polycrystalline Nb grain size
    CrystalLattice nb_poly(nb_cfg);

    const auto kin_al2o3 = al2o3.map_k_to_vg(PhononMode::L, Vec3{0.0, 0.0, 1.0});
    const auto kin_cawo4 = cawo4.map_k_to_vg(PhononMode::ST, Vec3{1.0, 1.0, 0.0}.unit());
    const auto kin_nb    = nb_poly.map_k_to_vg(PhononMode::L, Vec3{1.0, 0.0, 0.0});

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Sapphire (Al2O3, Rhombohedral, [100] orient) : v_g(L, +z)  = "
              << kin_al2o3.vg_mag << " m/s\n"
              << "  CaWO4 (Tetragonal C16=-17 GPa)               : v_g(ST,110) = "
              << kin_cawo4.vg_mag << " m/s\n"
              << "  Polycrystalline Nb (150 um grain MFP)        : v_g(L, +x)  = "
              << kin_nb.vg_mag << " m/s, Gamma_poly = "
              << std::scientific << std::setprecision(2) << nb_poly.polycrystal_rate(kin_nb.vg_mag)
              << " s^-1\n";
    std::cout << "========================================================================\n\n";
}

// ============================================================================
// 5. Matrix-Element IV Scattering, Charge Trapping & Trap Ionization (Burst 2)
// ============================================================================
static void run_charge_trapping_demo() {
    std::cout << "========================================================================\n";
    std::cout << " [5] Matrix-Element IV Rate, Charge Trapping & Trap Ionization (Burst 2)\n";
    std::cout << "========================================================================\n";

    CrystalConfig si_cfg = materials::Silicon();
    si_cfg.iv_model            = IVRateModel::Matrix;
    si_cfg.e_trapping_mfp_m    = 5.0 * units::cm;   // 5 cm electron trapping MFP
    si_cfg.e_trap_ion_d0_mfp_m = 1.5 * units::cm;   // 1.5 cm donor trap impact ionization MFP
    CrystalLattice si_lat(si_cfg);

    DetectorModel det{};
    det.e_field = UniformElectricField::from_bias(6.0 * units::Volt, 2.0 * det.solid.half_z);
    det.luke_sample_fraction = 0.0; // Focus on charge transport only
    SimulationEngine engine(si_lat, det);

    Rng rng(2026);
    std::vector<ChargeTrack> charges;
    std::vector<PhononTrack> phonons;
    std::vector<SensorHit> hits;
    SimStats stats{};

    EnergyPartition part(si_lat);
    part.partition_deposit({0.0, 0.0, 0.0}, 0.0, 200.0 * units::eV, /*is_nuclear_recoil=*/false, charges, phonons, rng);
    engine.transport_charges(charges, phonons, hits, stats, rng);

    const double iv_mat_50meV = physics::intervalley_matrix_rate(50.0 * units::meV, si_cfg);
    const double r_cloud_nm   = physics::ChargeCloudGenerator::compute_cloud_radius(50, si_cfg) / units::nm;

    std::cout << std::scientific << std::setprecision(3)
              << "  Si Matrix-Element IV Rate @ 50 meV : " << iv_mat_50meV << " s^-1\n"
              << std::fixed << std::setprecision(2)
              << "  3D Charge Cloud Radius (50 pairs)  : " << r_cloud_nm << " nm\n"
              << "  Intervalley Scatters (Matrix mode) : " << stats.intervalley_scatters << "\n"
              << "  Trap Impact Ionizations (e->D+)    : " << stats.trap_ionizations << "\n"
              << "  Bulk Trapped Carriers              : " << stats.charges_trapped << "\n"
              << "  Collected Electrode Carriers       : " << stats.charges_collected << "\n";
    std::cout << "========================================================================\n\n";
}

// ============================================================================
// 6. All 5 NIEL Partition Models & Exact Fano Binomial Verification (Burst 3)
// ============================================================================
static void run_niel_fano_demo() {
    std::cout << "========================================================================\n";
    std::cout << " [6] All 5 NIEL Ionization Models & Exact Fano Binomial (Burst 3)\n";
    std::cout << "========================================================================\n";

    const CrystalConfig si_cfg = materials::Silicon();
    const CrystalConfig ge_cfg = materials::Germanium();
    constexpr double E_REC = 250.0 * units::eV;

    const double y_lind   = physics::compute_ionization_yield(physics::NIELModel::Lindhard,   E_REC, si_cfg);
    const double y_lewin  = physics::compute_ionization_yield(physics::NIELModel::LewinSmith, E_REC, si_cfg);
    const double y_sarkis = physics::compute_ionization_yield(physics::NIELModel::Sarkis,     E_REC, si_cfg);
    const double y_impact = physics::compute_ionization_yield(physics::NIELModel::ImpactTUNL, E_REC, si_cfg);
    const double y_emp_ge = physics::compute_ionization_yield(physics::NIELModel::Empirical,  E_REC, ge_cfg);

    Rng rng(999);
    double sum = 0.0, sum2 = 0.0;
    constexpr int N_TRIALS = 20000;
    constexpr double MEAN_PAIRS = 4.5;
    constexpr double FANO_SI    = 0.115;
    for (int i = 0; i < N_TRIALS; ++i) {
        const double n = static_cast<double>(physics::sample_fano_binomial(MEAN_PAIRS, FANO_SI, rng));
        sum  += n;
        sum2 += n * n;
    }
    const double mean_meas = sum / N_TRIALS;
    const double var_meas  = (sum2 / N_TRIALS) - mean_meas * mean_meas;
    const double fano_meas = var_meas / mean_meas;

    std::cout << std::fixed << std::setprecision(3)
              << "  250 eV Si Nuclear Recoil Yield (Lindhard)    : " << (100.0 * y_lind)   << "%\n"
              << "  250 eV Si Nuclear Recoil Yield (Lewin-Smith) : " << (100.0 * y_lewin)  << "%\n"
              << "  250 eV Si Nuclear Recoil Yield (Sarkis 2023) : " << (100.0 * y_sarkis) << "%\n"
              << "  250 eV Si Nuclear Recoil Yield (IMPACT@TUNL) : " << (100.0 * y_impact) << "%\n"
              << "  250 eV Ge Nuclear Recoil Yield (CDMSlite Emp): " << (100.0 * y_emp_ge) << "%\n"
              << "  Exact Fano Binomial (mean=4.5, F=0.115)      : measured mean=" << mean_meas
              << ", measured F=" << fano_meas << "\n";
    std::cout << "========================================================================\n\n";
}

// ============================================================================
// 7. 3D Superconducting Bogoliubov Quasiparticle Coupled Cascade (Burst 4)
// ============================================================================
static void run_superconductor_3d_demo() {
    std::cout << "========================================================================\n";
    std::cout << " [7] 3D Superconducting Bogoliubov QP + Phonon Coupled Cascade (Burst 4)\n";
    std::cout << "========================================================================\n";

    CrystalConfig al_cfg = materials::AluminumCrystal();
    al_cfg.sc_teff_K     = 0.15;               // 150 mK effective bath temperature
    al_cfg.sc_trap_tau_s = 25.0 * units::us;   // 25 us vortex/local trapping lifetime
    CrystalLattice al_sc(al_cfg);

    DetectorModel sc_det{};
    sc_det.solid.radius = 2.0 * units::mm;
    sc_det.solid.half_z = 0.5 * units::um;     // 1 um thick 3D superconducting Al film
    sc_det.surfaces[0].prob_abs = 0.05;        // 5% tunnel junction / collector pad absorption
    sc_det.surfaces[1].prob_abs = 0.05;
    sc_det.surfaces[2].prob_abs = 0.0;
    SimulationEngine sc_engine(al_sc, sc_det);

    Rng rng(777);
    std::vector<PhononTrack> phonons;
    std::vector<GenericExcitationTrack> qps;
    std::vector<SensorHit> hits;
    SimStats stats{};

    // Inject 200 high-energy pair-breaking phonons (E = 10 * Delta_0 = 1.8 meV) inside 3D Al film
    const double e_init = 10.0 * al_cfg.sc_delta0_eV;
    for (int i = 0; i < 200; ++i) {
        phonons.push_back(sc_engine.make_track({0.0, 0.0, 0.0}, rng.isotropic_dir(), PhononMode::L, e_init, 0.0));
        stats.total_injected_eV += e_init;
    }

    const auto t0 = std::chrono::steady_clock::now();
    sc_engine.transport_superconductor_coupled(phonons, qps, hits, stats, rng, /*planar_2d_diffusion=*/false);
    const auto t1 = std::chrono::steady_clock::now();
    const double dt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    const auto& sc = al_sc.superconductor();
    const double e_cons_diff = std::abs(stats.total_injected_eV - (stats.total_deposited_eV + stats.escaped_or_killed_eV))
                               / stats.total_injected_eV;

    std::cout << std::fixed << std::setprecision(3)
              << "  3D Al Superconductor Gap Delta(150 mK) : " << (sc.delta_T_eV / units::ueV) << " ueV (T_c = "
              << sc.t_crit_K << " K)\n"
              << "  Cooper-Pair Breakings (Phonon -> 2 QP) : " << stats.sc_pair_breakings << "\n"
              << "  QP Phonon Radiations  (QP -> QP' + Ph) : " << stats.qp_phonons_radiated << "\n"
              << "  QP Thermal Recombinations (QP -> 2D Ph): " << stats.qp_recombinations << "\n"
              << "  QP Local Trapping Events               : " << stats.qp_local_traps << "\n"
              << "  3D QP Brownian Diffusion Steps         : " << stats.qp_diffusion_steps << "\n"
              << "  Energy Conservation Relative Diff      : " << std::scientific << e_cons_diff
              << std::fixed << " (Runtime: " << dt_ms << " ms)\n";
    std::cout << "========================================================================\n\n";
}

// ============================================================================
// 8. 3D/2D Tetrahedral FEM Mesh Field + Patterned 2D Electrode Mask (Burst 5)
// ============================================================================
static void run_mesh_field_and_mask_demo(const CrystalLattice& ge_lattice) {
    std::cout << "========================================================================\n";
    std::cout << " [8] 3D/2D FEM Mesh Electric Field + Patterned 2D QET Mask (Burst 5)\n";
    std::cout << "========================================================================\n";

    // (a) Build a 3D Kuhn-tetrahedralized FEM mesh (12x12x12 hex cells = 10,368 tetrahedra)
    //     for an iZIP-style non-uniform potential V(x,y,z) = -2000 * z + 8000 * (x^2 + y^2)
    const double R = 3.81 * units::cm;
    const double H = 1.27 * units::cm;
    const auto mesh_3d = geometry::TriLinearMesh3D::build_cartesian_grid(
        R, R, H, 12, 12, 12,
        [](const Vec3& p) {
            return -2000.0 * p.z + 8000.0 * (p.x * p.x + p.y * p.y);
        }
    );

    // Build a 2D axisymmetric (r, z) triangular FEM mesh (24x24 cells = 1,152 triangles)
    const auto mesh_2d = geometry::BiLinearMesh2D::build_grid(
        0.0, R, -H, H, 24, 24,
        [](double r, double z) {
            return -2000.0 * z + 8000.0 * (r * r);
        },
        /*axisymmetric_rz=*/true
    );

    const Vec3 test_pt{1.0 * units::cm, -0.5 * units::cm, 0.3 * units::cm};
    const Vec3 E_3d = mesh_3d.at(test_pt);
    const Vec3 E_2d = mesh_2d.at(test_pt);

    // (b) Configure detector with 3D Tetrahedral FEM field and 5-channel SuperCDMS QET mask on TopFace
    DetectorModel det{};
    det.solid.radius = R;
    det.solid.half_z = H;
    det.mesh_field_3d = mesh_3d;
    det.luke_sample_fraction = 0.05;
    det.surfaces[static_cast<size_t>(SurfaceId::TopFace)].electrode_mask =
        geometry::ElectrodePatternMask::supercdms_qet_mask(R, /*fill_fraction=*/0.85);

    SimulationEngine engine(ge_lattice, det);
    EnergyPartition partitioner(ge_lattice);

    Rng rng(2026);
    std::vector<PhononTrack> phonons;
    std::vector<ChargeTrack> charges;
    std::vector<SensorHit> hits;
    SimStats stats{};

    PartitionOptions part_opts{};
    part_opts.max_sampled_phonons = 500;
    partitioner.partition_deposit({0.5 * units::cm, 0.2 * units::cm, 0.0}, 0.0,
                                  50.0 * units::eV, /*is_nuclear_recoil=*/false,
                                  charges, phonons, rng, part_opts);

    engine.transport_charges(charges, phonons, hits, stats, rng);
    engine.transport(phonons, hits, stats, rng);

    std::array<uint32_t, 5> ch_hits{0, 0, 0, 0, 0};
    for (const SensorHit& h : hits) {
        if (h.surface_id == SurfaceId::TopFace && h.channel_id >= 0 && h.channel_id < 5) {
            ch_hits[h.channel_id]++;
        }
    }

    // (c) Compare Standard `map_k_to_vg` vs Zero-Trig FP32 Octahedral `map_k_to_vg_oct`
    double max_cos_err = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Vec3 kd = rng.isotropic_dir();
        const auto k_std = ge_lattice.map_k_to_vg(PhononMode::ST, kd);
        const auto k_oct = ge_lattice.map_k_to_vg_oct(PhononMode::ST, kd);
        max_cos_err = std::max(max_cos_err, 1.0 - k_std.vg_dir.dot(k_oct.vg_dir));
    }

    KaplanSurrogateLUT kaplan_lut(KaplanFilm::aluminum_sensor());
    SmallVec<double, 64> refl;
    const double fast_edep = kaplan_lut.absorb_phonon_fast(2.0 * units::meV, refl, rng);

    std::cout << std::fixed << std::setprecision(2)
              << "  3D Tetrahedral FEM Mesh  : " << mesh_3d.num_tetrahedra() << " tets, E(test) = ("
              << E_3d.x << ", " << E_3d.y << ", " << E_3d.z << ") V/m\n"
              << "  2D Axisym (r,z) FEM Mesh : " << mesh_2d.num_triangles() << " tris, E(test) = ("
              << E_2d.x << ", " << E_2d.y << ", " << E_2d.z << ") V/m\n"
              << "  TopFace QET Channel Hits : Ch0(Inner)=" << ch_hits[0]
              << ", Ch1=" << ch_hits[1] << ", Ch2=" << ch_hits[2]
              << ", Ch3=" << ch_hits[3] << ", Ch4(Guard)=" << ch_hits[4] << "\n"
              << "  Octahedral FP32 LUT Err  : max(1 - cos dtheta) = " << std::scientific << max_cos_err
              << std::fixed << ", KaplanSurrogateLUT sample = " << (fast_edep / units::meV) << " meV\n";
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

    std::cout << "Initialized Germanium Lattice (Christoffel LUT + Decay Inverse-CDF + 4 HV L-Valleys)\n"
              << "  Init time : " << std::fixed << std::setprecision(2) << init_ms << " ms\n"
              << "  m_cond_e  : " << (ge_lattice.electron_mass_eV() / units::m_e_eV) << " m_e, "
              << "l0_e = " << (ge_lattice.l0_electron() / units::nm) << " nm, "
              << "l0_h = " << (ge_lattice.l0_hole() / units::nm) << " nm\n\n";

    run_caustics_demo(ge_lattice);
    run_detector_benchmark(ge_lattice, num_primaries, num_threads);
    run_recoil_charge_luke_demo(ge_lattice);
    run_crystal_symmetries_demo();
    run_charge_trapping_demo();
    run_niel_fano_demo();
    run_superconductor_3d_demo();
    run_mesh_field_and_mask_demo(ge_lattice);

    return 0;
}
