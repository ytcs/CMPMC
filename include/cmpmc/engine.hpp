#pragma once

#include "cmpmc/charge.hpp"
#include "cmpmc/core.hpp"
#include "cmpmc/geometry.hpp"
#include "cmpmc/kaplan.hpp"
#include "cmpmc/lattice.hpp"
#include "cmpmc/partition.hpp"
#include "cmpmc/physics/custom_excitation_example.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace cmpmc {

struct SimStats {
    // Phonon statistics
    uint64_t total_steps          = 0;
    uint64_t isotope_scatters     = 0;
    uint64_t polycrystal_scatters = 0;
    uint64_t decays_L_to_TT       = 0;
    uint64_t decays_L_to_LT       = 0;
    uint64_t surface_refls        = 0;
    uint64_t sensor_hits          = 0;
    uint64_t killed_tracks        = 0;

    // Charge carrier & Luke statistics
    uint64_t charge_steps         = 0;
    uint64_t luke_phonons_emitted = 0;
    uint64_t intervalley_scatters = 0;
    uint64_t charges_trapped      = 0;
    uint64_t trap_ionizations     = 0;
    uint64_t charge_surface_refls = 0;
    uint64_t charges_collected    = 0;
    double total_luke_phonon_eV   = 0.0;
    double charge_recomb_edep_eV  = 0.0;

    // 3D Superconducting Bogoliubov Quasiparticle statistics
    uint64_t sc_pair_breakings    = 0;
    uint64_t qp_diffusion_steps   = 0;
    uint64_t qp_phonons_radiated  = 0;
    uint64_t qp_recombinations    = 0;
    uint64_t qp_local_traps       = 0;

    // Energy accounting
    double total_injected_eV      = 0.0;
    double total_deposited_eV     = 0.0;
    double escaped_or_killed_eV   = 0.0;

    SimStats& operator+=(const SimStats& r) noexcept {
        total_steps          += r.total_steps;
        isotope_scatters     += r.isotope_scatters;
        polycrystal_scatters += r.polycrystal_scatters;
        decays_L_to_TT       += r.decays_L_to_TT;
        decays_L_to_LT       += r.decays_L_to_LT;
        surface_refls        += r.surface_refls;
        sensor_hits          += r.sensor_hits;
        killed_tracks        += r.killed_tracks;
        charge_steps         += r.charge_steps;
        luke_phonons_emitted += r.luke_phonons_emitted;
        intervalley_scatters += r.intervalley_scatters;
        charges_trapped      += r.charges_trapped;
        trap_ionizations     += r.trap_ionizations;
        charge_surface_refls += r.charge_surface_refls;
        charges_collected    += r.charges_collected;
        total_luke_phonon_eV += r.total_luke_phonon_eV;
        charge_recomb_edep_eV+= r.charge_recomb_edep_eV;
        sc_pair_breakings    += r.sc_pair_breakings;
        qp_diffusion_steps   += r.qp_diffusion_steps;
        qp_phonons_radiated  += r.qp_phonons_radiated;
        qp_recombinations    += r.qp_recombinations;
        qp_local_traps       += r.qp_local_traps;
        total_injected_eV    += r.total_injected_eV;
        total_deposited_eV   += r.total_deposited_eV;
        escaped_or_killed_eV += r.escaped_or_killed_eV;
        return *this;
    }
};

// ============================================================================
// DetectorModel & Unified Multi-Excitation SimulationEngine
// ============================================================================
struct DetectorModel {
    CylinderSolid solid{};
    UniformElectricField e_field = UniformElectricField::from_bias(0.0, 2.54 * units::cm);
    std::optional<geometry::TriLinearMesh3D> mesh_field_3d = std::nullopt;
    std::optional<geometry::BiLinearMesh2D>  mesh_field_2d = std::nullopt;
    std::array<SurfaceProperty, static_cast<size_t>(SurfaceId::Count)> surfaces{
        SurfaceProperty::g4cmp_al_sensor_surface(), // TopFace
        SurfaceProperty::g4cmp_al_sensor_surface(), // BottomFace
        SurfaceProperty::g4cmp_bare_wall_surface()  // SideWall
    };
    uint32_t max_bounces       = 5000;
    double min_phonon_eV       = 10.0 * units::ueV; // Cutoff below 2*Delta_W
    double luke_sample_fraction= 1.0;               // Fraction of Luke phonons tracked (with weight = 1/frac)

    inline Vec3 electric_field_at(const Vec3& pos) const noexcept {
        if (mesh_field_3d.has_value()) return mesh_field_3d->at(pos);
        if (mesh_field_2d.has_value()) return mesh_field_2d->at(pos);
        return e_field.at(pos);
    }
};

class SimulationEngine {
public:
    SimulationEngine(const CrystalLattice& lattice, const DetectorModel& detector)
        : lattice_(lattice), detector_(detector), charge_phys_(lattice) {}

    inline PhononTrack make_track(const Vec3& pos, const Vec3& k_dir,
                                  PhononMode mode, double energy_eV,
                                  double time_s = 0.0, double weight = 1.0) const noexcept {
        const Vec3 k_unit = k_dir.unit();
        const auto kin = lattice_.map_k_to_vg(mode, k_unit);
        return PhononTrack{
            .pos     = pos,
            .k_dir   = k_unit,
            .vg_dir  = kin.vg_dir,
            .energy  = energy_eV,
            .vg_mag  = kin.vg_mag,
            .time    = time_s,
            .weight  = weight,
            .bounces = 0,
            .mode    = mode
        };
    }

    // ========================================================================
    // 1. Charge Carrier Transport (Electrons & Holes + Luke Phonon Emission)
    // ========================================================================
    void transport_charges(std::vector<ChargeTrack>& charge_stack,
                           std::vector<PhononTrack>& phonon_stack,
                           std::vector<SensorHit>& hits,
                           SimStats& stats,
                           Rng& rng) const {
        const double e_gap = lattice_.config().bandgap;
        const double luke_prob = std::clamp(detector_.luke_sample_fraction, 0.0, 1.0);
        const double luke_wt_mult = (luke_prob > 0.0) ? (1.0 / luke_prob) : 0.0;

        while (!charge_stack.empty()) {
            ChargeTrack trk = charge_stack.back();
            charge_stack.pop_back();

            const double q_sign = (trk.type == ExcitationType::Electron) ? -1.0 : +1.0;

            for (uint32_t step = 0; step < 500'000; ++step) {
                stats.charge_steps++;

                const Vec3 E_Vm = detector_.electric_field_at(trk.pos);
                const Vec3 v_cur = charge_phys_.velocity(trk);
                const double v_cur_mag = v_cur.norm();

                const double r_luke = charge_phys_.luke_emission_rate(trk);
                const double r_iv   = (trk.type == ExcitationType::Electron)
                    ? lattice_.intervalley_rate(trk.ekin, E_Vm)
                    : 0.0;
                const double r_trap  = physics::charge_trapping_rate(trk.type, v_cur_mag, lattice_.config());
                const double r_ion_e = physics::trap_ionization_e_rate(trk.type, v_cur_mag, lattice_.config());
                const double r_ion_h = physics::trap_ionization_h_rate(trk.type, v_cur_mag, lattice_.config());
                const double r_tot   = r_luke + r_iv + r_trap + r_ion_e + r_ion_h;

                const double dt_max_accel = 100.0 * units::ps;
                const double dt_phys = (r_tot > 1.0)
                    ? (-std::log(rng.uniform()) / r_tot)
                    : std::numeric_limits<double>::infinity();
                const double dt = std::min(dt_phys, dt_max_accel);

                const Vec3 dp_eV = E_Vm * (q_sign * units::c_light * dt);
                const Vec3 p_mid = trk.p_eV + dp_eV * 0.5;
                const Vec3 v_mid = (trk.type == ExcitationType::Electron)
                    ? lattice_.map_p_to_vel_electron(trk.valley, p_mid)
                    : lattice_.map_p_to_vel_hole(p_mid);
                const double v_mid_mag = v_mid.norm();
                const double d_step    = v_mid_mag * dt;

                if (v_mid_mag > 1e-6) {
                    const Vec3 v_dir = v_mid / v_mid_mag;
                    const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, v_dir);

                    if (bhit.distance <= d_step) {
                        const double dt_bound = bhit.distance / v_mid_mag;
                        trk.pos  += v_dir * bhit.distance;
                        trk.time += dt_bound;
                        trk.p_eV += E_Vm * (q_sign * units::c_light * dt_bound);
                        trk.ekin  = (trk.type == ExcitationType::Electron)
                            ? lattice_.map_p_to_ekin_electron(trk.valley, trk.p_eV)
                            : lattice_.map_p_to_ekin_hole(trk.p_eV);

                        const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];
                        int32_t ch_id = -1;
                        bool hit_active_pad = true;
                        if (surf.electrode_mask.has_value()) {
                            ch_id = surf.electrode_mask->find_channel(trk.pos, rng);
                            hit_active_pad = (ch_id >= 0);
                        }
                        if (!hit_active_pad ||
                            !charge_phys_.should_absorb_boundary(trk, bhit.inward_normal,
                                                                 surf.prob_charge_abs, surf.min_k_elec, surf.min_k_hole, rng)) {
                            stats.charge_surface_refls++;
                            trk.pos += bhit.inward_normal * 1e-11;
                            charge_phys_.reflect_at_boundary(trk, bhit.inward_normal, surf.prob_charge_spec, rng);
                            continue;
                        }

                        const double recomb_eV = (trk.ekin + 0.5 * e_gap) * trk.weight;
                        stats.charges_collected++;
                        stats.charge_recomb_edep_eV += recomb_eV;
                        stats.total_deposited_eV    += recomb_eV;
                        hits.push_back(SensorHit{
                            .pos             = trk.pos,
                            .time            = trk.time,
                            .edep            = recomb_eV,
                            .incident_energy = trk.ekin + 0.5 * e_gap,
                            .surface_id      = bhit.surface_id,
                            .excitation      = trk.type,
                            .mode            = PhononMode::L,
                            .channel_id      = ch_id
                        });
                        break;
                    }
                }

                trk.pos  += v_mid * dt;
                trk.time += dt;
                trk.p_eV += dp_eV;
                trk.ekin  = (trk.type == ExcitationType::Electron)
                    ? lattice_.map_p_to_ekin_electron(trk.valley, trk.p_eV)
                    : lattice_.map_p_to_ekin_hole(trk.p_eV);

                if (dt_phys <= dt_max_accel) {
                    const double u_proc = rng.uniform() * r_tot;
                    if (u_proc < r_trap) {
                        stats.charges_trapped++;
                        stats.escaped_or_killed_eV += (trk.ekin + 0.5 * e_gap) * trk.weight;
                        break;
                    } else if (u_proc < r_trap + r_ion_e + r_ion_h) {
                        stats.trap_ionizations++;
                        const bool lib_elec = (u_proc < r_trap + r_ion_e);
                        const double e_knock = 1.0e-6 * units::eV; // G4CMPDriftTrapIonization knock-on energy
                        const size_t n_val = std::max<size_t>(1, lattice_.num_valleys());
                        const uint16_t kv = lib_elec ? static_cast<uint16_t>(rng.next_u64() % n_val) : 0;
                        ChargeTrack knockon = charge_phys_.make_charge(
                            lib_elec ? ExcitationType::Electron : ExcitationType::Hole,
                            trk.pos, rng.isotropic_dir(), e_knock, kv, trk.time);
                        knockon.weight = trk.weight;
                        charge_stack.push_back(knockon);
                    } else if (u_proc < r_trap + r_ion_e + r_ion_h + r_iv) {
                        charge_phys_.scatter_intervalley(trk, rng);
                        stats.intervalley_scatters++;
                    } else if (r_luke > 0.0) {
                        const LukePhononEmit emit = charge_phys_.emit_luke_phonon(trk, rng);
                        if (emit.emitted) {
                            stats.luke_phonons_emitted++;
                            stats.total_luke_phonon_eV += emit.phonon_eV * trk.weight;

                            if (luke_prob > 0.0 && (luke_prob >= 1.0 || rng.uniform() < luke_prob)) {
                                const PhononMode m = lattice_.sample_dos_mode(rng);
                                phonon_stack.push_back(make_track(
                                    trk.pos, emit.q_dir, m, emit.phonon_eV,
                                    trk.time, trk.weight * luke_wt_mult
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    // ========================================================================
    // 2. Phonon Transport (Ballistic Anisotropic Ray-Tracing + Kaplan QP Cascade)
    // ========================================================================
    void transport(std::vector<PhononTrack>& stack,
                   std::vector<SensorHit>& hits,
                   SimStats& stats,
                   Rng& rng) const {
        SmallVec<double, 64> kaplan_reflected;

        while (!stack.empty()) {
            PhononTrack trk = stack.back();
            stack.pop_back();

            for (;;) {
                if (trk.energy < detector_.min_phonon_eV || trk.bounces >= detector_.max_bounces) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += trk.energy * trk.weight;
                    break;
                }

                stats.total_steps++;

                const double rate_iso  = lattice_.isotope_rate(trk.energy);
                const double rate_poly = lattice_.polycrystal_rate(trk.vg_mag);
                const double rate_down = lattice_.downconversion_rate(trk.mode, trk.energy);
                const double rate_tot  = rate_iso + rate_poly + rate_down;

                const double d_phys = (rate_tot > 1e-12)
                    ? (-std::log(rng.uniform()) * (trk.vg_mag / rate_tot))
                    : std::numeric_limits<double>::infinity();

                const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, trk.vg_dir);

                if (d_phys < bhit.distance) {
                    trk.pos  += trk.vg_dir * d_phys;
                    trk.time += d_phys / trk.vg_mag;

                    const double u_proc = rng.uniform() * rate_tot;
                    if (rate_down > 0.0 && u_proc < rate_down) {
                        perform_anharmonic_decay(trk, stack, stats, rng);
                    } else {
                        if (u_proc < rate_down + rate_iso) {
                            stats.isotope_scatters++;
                        } else {
                            stats.polycrystal_scatters++;
                        }
                        trk.mode  = lattice_.sample_dos_mode(rng);
                        trk.k_dir = rng.isotropic_dir();
                        const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                        trk.vg_dir = kin.vg_dir;
                        trk.vg_mag = kin.vg_mag;
                    }
                    continue;
                }

                trk.pos  += trk.vg_dir * bhit.distance;
                trk.time += bhit.distance / trk.vg_mag;
                trk.bounces++;

                const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];
                int32_t ch_id = -1;
                bool hit_active_pad = true;
                if (surf.electrode_mask.has_value()) {
                    ch_id = surf.electrode_mask->find_channel(trk.pos, rng);
                    hit_active_pad = (ch_id >= 0);
                }

                if (hit_active_pad && surf.kaplan_sensor.has_value() && rng.uniform() < surf.prob_abs) {
                    const KaplanFilm& film = *surf.kaplan_sensor;
                    if (rng.uniform() < film.film_absorption) {
                        const double edep = film.absorb_phonon(trk.energy, kaplan_reflected, rng);
                        if (edep > 0.0) {
                            const double weighted_edep = edep * trk.weight;
                            stats.sensor_hits++;
                            stats.total_deposited_eV += weighted_edep;
                            hits.push_back(SensorHit{
                                .pos             = trk.pos,
                                .time            = trk.time,
                                .edep            = weighted_edep,
                                .incident_energy = trk.energy,
                                .surface_id      = bhit.surface_id,
                                .excitation      = ExcitationType::Phonon,
                                .mode            = trk.mode,
                                .channel_id      = ch_id
                            });
                        }
                        double sum_refl = 0.0;
                        for (double re : kaplan_reflected) sum_refl += re;
                        stats.escaped_or_killed_eV += std::max(0.0, trk.energy - edep - sum_refl) * trk.weight;

                        if (kaplan_reflected.empty()) {
                            break;
                        }
                        const Vec3 reemit_pos = trk.pos + bhit.inward_normal * 1e-11;
                        for (size_t i = 1; i < kaplan_reflected.size(); ++i) {
                            const PhononMode m = lattice_.sample_dos_mode(rng);
                            const Vec3 kd = sample_inward_wavevector(m, bhit.inward_normal, rng);
                            stack.push_back(make_track(reemit_pos, kd, m, kaplan_reflected[i], trk.time, trk.weight));
                        }
                        trk.pos    = reemit_pos;
                        trk.energy = kaplan_reflected[0];
                        trk.mode   = lattice_.sample_dos_mode(rng);
                        trk.k_dir  = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                        const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                        trk.vg_dir = kin.vg_dir;
                        trk.vg_mag = kin.vg_mag;
                        continue;
                    }
                }

                if (rng.uniform() > surf.prob_refl) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += trk.energy * trk.weight;
                    break;
                }

                stats.surface_refls++;
                trk.pos += bhit.inward_normal * 1e-11;

                const double p_spec = surf.specular_prob(trk.energy);
                if (rng.uniform() < p_spec) {
                    Vec3 k_refl = trk.k_dir - bhit.inward_normal * (2.0 * trk.k_dir.dot(bhit.inward_normal));
                    auto kin = lattice_.map_k_to_vg(trk.mode, k_refl);
                    if (kin.vg_dir.dot(bhit.inward_normal) <= 1e-6) {
                        k_refl = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                        kin = lattice_.map_k_to_vg(trk.mode, k_refl);
                    }
                    trk.k_dir  = k_refl;
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                } else {
                    trk.k_dir  = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                    const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                }
            }
        }
    }

    // ========================================================================
    // 3. Extensible Custom Excitation Transport (BogoliubovQP, Exciton, Magnon)
    //    Uses C++20 `physics::ExcitationProcessConcept` for zero-overhead inlining!
    // ========================================================================
    template <physics::ExcitationProcessConcept Model>
    void transport_custom_excitation(const Model& model,
                                     std::vector<GenericExcitationTrack>& stack,
                                     std::vector<SensorHit>& hits,
                                     SimStats& stats,
                                     Rng& rng) const {
        while (!stack.empty()) {
            GenericExcitationTrack trk = stack.back();
            stack.pop_back();

            for (;;) {
                if (trk.energy < model.min_energy_eV() || trk.bounces >= detector_.max_bounces) {
                    stats.killed_tracks++;
                    stats.escaped_or_killed_eV += trk.energy * trk.weight;
                    break;
                }

                stats.total_steps++;
                const double rate = model.bulk_rate(trk);
                const double d_phys = (rate > 1e-12 && trk.velocity > 0.0)
                    ? (-std::log(rng.uniform()) * (trk.velocity / rate))
                    : std::numeric_limits<double>::infinity();

                const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, trk.mom_dir);

                if (d_phys < bhit.distance) {
                    trk.pos  += trk.mom_dir * d_phys;
                    trk.time += d_phys / std::max(trk.velocity, 1e-6);
                    model.perform_bulk_interaction(trk, stack, rng);
                    continue;
                }

                trk.pos  += trk.mom_dir * bhit.distance;
                trk.time += bhit.distance / std::max(trk.velocity, 1e-6);
                trk.bounces++;

                const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];
                if (rng.uniform() < surf.prob_abs) {
                    const double weighted_edep = trk.energy * trk.weight;
                    stats.sensor_hits++;
                    stats.total_deposited_eV += weighted_edep;
                    hits.push_back(SensorHit{
                        .pos             = trk.pos,
                        .time            = trk.time,
                        .edep            = weighted_edep,
                        .incident_energy = trk.energy,
                        .surface_id      = bhit.surface_id,
                        .excitation      = trk.type,
                        .mode            = PhononMode::L
                    });
                    break;
                }

                stats.surface_refls++;
                trk.pos     += bhit.inward_normal * 1e-11;
                trk.mom_dir  = rng.lambertian_dir(bhit.inward_normal);
            }
        }
    }

    // ========================================================================
    // 4. Full 3D Coupled Superconducting Bogoliubov Quasiparticle + Phonon
    //    Transport (G4CMPBogoliubovQP, G4CMPSCPairBreakingProcess,
    //    G4CMPQPRadiatesPhononProcess, G4CMPQPRecombinationProcess,
    //    G4CMPQPLocalTrappingProcess, G4CMPQPDiffusion)
    // ========================================================================
    void transport_superconductor_coupled(std::vector<PhononTrack>& phonon_stack,
                                          std::vector<GenericExcitationTrack>& qp_stack,
                                          std::vector<SensorHit>& hits,
                                          SimStats& stats,
                                          Rng& rng,
                                          bool planar_2d_diffusion = false) const {
        const auto& sc = lattice_.superconductor();
        if (!sc.is_enabled()) {
            transport(phonon_stack, hits, stats, rng);
            return;
        }

        while (!phonon_stack.empty() || !qp_stack.empty()) {
            // (A) Transport 3D Bulk Phonons (with Cooper-pair breaking into 2 Bogoliubov QPs)
            while (!phonon_stack.empty()) {
                PhononTrack trk = phonon_stack.back();
                phonon_stack.pop_back();

                for (;;) {
                    if (trk.energy < detector_.min_phonon_eV || trk.bounces >= detector_.max_bounces) {
                        stats.killed_tracks++;
                        stats.escaped_or_killed_eV += trk.energy * trk.weight;
                        break;
                    }
                    stats.total_steps++;

                    const double rate_iso  = lattice_.isotope_rate(trk.energy);
                    const double rate_poly = lattice_.polycrystal_rate(trk.vg_mag);
                    const double rate_down = lattice_.downconversion_rate(trk.mode, trk.energy);
                    const double rate_pb   = sc.phonon_pair_breaking_rate(trk.energy);
                    const double rate_tot  = rate_iso + rate_poly + rate_down + rate_pb;

                    const double d_phys = (rate_tot > 1e-12)
                        ? (-std::log(rng.uniform()) * (trk.vg_mag / rate_tot))
                        : std::numeric_limits<double>::infinity();

                    const BoundaryHit bhit = detector_.solid.intersect_exit(trk.pos, trk.vg_dir);

                    if (d_phys < bhit.distance) {
                        trk.pos  += trk.vg_dir * d_phys;
                        trk.time += d_phys / trk.vg_mag;

                        const double u_proc = rng.uniform() * rate_tot;
                        if (rate_pb > 0.0 && u_proc < rate_pb) {
                            // Cooper-pair breaking: Phonon (E >= 2*Delta) -> QP1(E1) + QP2(E - E1)
                            stats.sc_pair_breakings++;
                            const double e1 = sc.sample_pair_breaking_qp_energy(trk.energy, rng);
                            const double e2 = trk.energy - e1;
                            qp_stack.push_back(GenericExcitationTrack{
                                .pos          = trk.pos,
                                .mom_dir      = rng.isotropic_dir(),
                                .energy       = e1,
                                .velocity     = sc.qp_velocity(e1),
                                .time         = trk.time,
                                .weight       = trk.weight,
                                .bounces      = 0,
                                .type         = ExcitationType::BogoliubovQP,
                                .spin_or_band = 0
                            });
                            qp_stack.push_back(GenericExcitationTrack{
                                .pos          = trk.pos,
                                .mom_dir      = rng.isotropic_dir(),
                                .energy       = e2,
                                .velocity     = sc.qp_velocity(e2),
                                .time         = trk.time,
                                .weight       = trk.weight,
                                .bounces      = 0,
                                .type         = ExcitationType::BogoliubovQP,
                                .spin_or_band = 0
                            });
                            break;
                        } else if (rate_down > 0.0 && u_proc < rate_pb + rate_down) {
                            perform_anharmonic_decay(trk, phonon_stack, stats, rng);
                        } else {
                            if (u_proc < rate_pb + rate_down + rate_iso) stats.isotope_scatters++;
                            else stats.polycrystal_scatters++;
                            trk.mode  = lattice_.sample_dos_mode(rng);
                            trk.k_dir = rng.isotropic_dir();
                            const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                            trk.vg_dir = kin.vg_dir;
                            trk.vg_mag = kin.vg_mag;
                        }
                        continue;
                    }

                    // Boundary reflection / absorption
                    trk.pos  += trk.vg_dir * bhit.distance;
                    trk.time += bhit.distance / trk.vg_mag;
                    trk.bounces++;
                    const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];
                    if (rng.uniform() < surf.prob_abs) {
                        const double weighted_edep = trk.energy * trk.weight;
                        stats.sensor_hits++;
                        stats.total_deposited_eV += weighted_edep;
                        hits.push_back(SensorHit{
                            .pos             = trk.pos,
                            .time            = trk.time,
                            .edep            = weighted_edep,
                            .incident_energy = trk.energy,
                            .surface_id      = bhit.surface_id,
                            .excitation      = ExcitationType::Phonon,
                            .mode            = trk.mode
                        });
                        break;
                    }
                    stats.surface_refls++;
                    trk.pos    += bhit.inward_normal * 1e-11;
                    trk.k_dir   = sample_inward_wavevector(trk.mode, bhit.inward_normal, rng);
                    const auto kin = lattice_.map_k_to_vg(trk.mode, trk.k_dir);
                    trk.vg_dir = kin.vg_dir;
                    trk.vg_mag = kin.vg_mag;
                }
            }

            // (B) Transport 3D Diffusing Bogoliubov Quasiparticles
            while (!qp_stack.empty()) {
                GenericExcitationTrack qp = qp_stack.back();
                qp_stack.pop_back();

                for (uint32_t step = 0; step < 100'000; ++step) {
                    if (qp.energy < sc.delta_T_eV || qp.bounces >= detector_.max_bounces) {
                        stats.killed_tracks++;
                        stats.escaped_or_killed_eV += qp.energy * qp.weight;
                        break;
                    }
                    stats.qp_diffusion_steps++;

                    const double r_rad  = sc.qp_radiates_phonon_rate(qp.energy);
                    const double r_rec  = sc.qp_recombination_rate(qp.energy);
                    const double r_trap = sc.qp_local_trapping_rate();
                    const double r_tot  = r_rad + r_rec + r_trap;

                    const double dt_max  = 1.0 * units::ns;
                    const double dt_phys = (r_tot > 1e-6)
                        ? (-std::log(rng.uniform()) / r_tot)
                        : std::numeric_limits<double>::infinity();
                    const double dt = std::min(dt_phys, dt_max);

                    // Sample diffusive displacement step (G4CMPQPDiffusion)
                    const Vec3 dr = sc.sample_diffusive_step(qp.energy, dt, rng, planar_2d_diffusion);
                    const double dr_mag = dr.norm();

                    if (dr_mag > 1e-15) {
                        const Vec3 dir = dr / dr_mag;
                        const BoundaryHit bhit = detector_.solid.intersect_exit(qp.pos, dir);
                        if (bhit.distance <= dr_mag) {
                            const double frac = bhit.distance / dr_mag;
                            qp.pos  += dir * bhit.distance;
                            qp.time += dt * frac;
                            qp.bounces++;

                            const SurfaceProperty& surf = detector_.surfaces[static_cast<size_t>(bhit.surface_id)];
                            if (rng.uniform() < surf.prob_abs) {
                                const double weighted_edep = qp.energy * qp.weight;
                                stats.sensor_hits++;
                                stats.total_deposited_eV += weighted_edep;
                                hits.push_back(SensorHit{
                                    .pos             = qp.pos,
                                    .time            = qp.time,
                                    .edep            = weighted_edep,
                                    .incident_energy = qp.energy,
                                    .surface_id      = bhit.surface_id,
                                    .excitation      = ExcitationType::BogoliubovQP,
                                    .mode            = PhononMode::L
                                });
                                break;
                            }
                            stats.surface_refls++;
                            qp.pos += bhit.inward_normal * 1e-11;
                            qp.mom_dir = rng.lambertian_dir(bhit.inward_normal);
                            continue;
                        }
                    }

                    qp.pos  += dr;
                    qp.time += dt;

                    if (dt_phys <= dt_max) {
                        const double u_proc = rng.uniform() * r_tot;
                        if (r_trap > 0.0 && u_proc < r_trap) {
                            // Local vortex/impurity trapping (G4CMPQPLocalTrappingProcess)
                            stats.qp_local_traps++;
                            stats.escaped_or_killed_eV += qp.energy * qp.weight;
                            break;
                        } else if (r_rec > 0.0 && u_proc < r_trap + r_rec) {
                            // Recombination with thermal bath QP -> 2*Delta phonon
                            stats.qp_recombinations++;
                            const double e_ph = qp.energy + sc.delta_T_eV;
                            // Account for the extra thermal bath QP energy sc.delta_T_eV injected by bath
                            stats.total_injected_eV += sc.delta_T_eV * qp.weight;
                            const PhononMode m = lattice_.sample_dos_mode(rng);
                            phonon_stack.push_back(make_track(
                                qp.pos, rng.isotropic_dir(), m, e_ph, qp.time, qp.weight
                            ));
                            break;
                        } else if (r_rad > 0.0) {
                            // QP radiates a phonon: QP(E) -> QP(E - Omega) + Phonon(Omega)
                            const double omega = sc.sample_radiated_phonon_energy(qp.energy, rng);
                            if (omega > 0.0 && qp.energy - omega >= sc.delta_T_eV) {
                                stats.qp_phonons_radiated++;
                                qp.energy  -= omega;
                                qp.velocity = sc.qp_velocity(qp.energy);
                                const PhononMode m = lattice_.sample_dos_mode(rng);
                                phonon_stack.push_back(make_track(
                                    qp.pos, rng.isotropic_dir(), m, omega, qp.time, qp.weight
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

private:
    const CrystalLattice& lattice_;
    const DetectorModel& detector_;
    ChargePhysics charge_phys_;

    inline Vec3 sample_inward_wavevector(PhononMode mode, const Vec3& inward_normal, Rng& rng) const noexcept {
        for (int trial = 0; trial < 16; ++trial) {
            const Vec3 kd = rng.lambertian_dir(inward_normal);
            const auto kin = lattice_.map_k_to_vg(mode, kd);
            if (kin.vg_dir.dot(inward_normal) > 1e-6) {
                return kd;
            }
        }
        return inward_normal;
    }

    inline void perform_anharmonic_decay(PhononTrack& trk,
                                         std::vector<PhononTrack>& stack,
                                         SimStats& stats,
                                         Rng& rng) const {
        const Vec3 dir1 = trk.k_dir;
        const Vec3 perp = dir1.orthogonal_unit();
        const double phi = units::two_pi * rng.uniform();
        const double c_phi = std::cos(phi);
        const double s_phi = std::sin(phi);

        if (rng.uniform() < lattice_.config().fTT_frac) {
            stats.decays_L_to_TT++;
            const double x   = lattice_.sample_tt_energy_fraction(rng);
            const double th1 = lattice_.theta_TT(x);
            const double th2 = lattice_.theta_TT(1.0 - x);

            const Vec3 k1 = dir1.rotate_around(perp, std::cos(th1), std::sin(th1)).rotate_around(dir1, c_phi, s_phi).unit();
            const Vec3 k2 = dir1.rotate_around(perp, std::cos(-th2), std::sin(-th2)).rotate_around(dir1, c_phi, s_phi).unit();

            const PhononMode m1 = lattice_.sample_transverse_mode(rng);
            const PhononMode m2 = lattice_.sample_transverse_mode(rng);

            stack.push_back(make_track(trk.pos, k2, m2, (1.0 - x) * trk.energy, trk.time, trk.weight));

            trk.energy = x * trk.energy;
            trk.mode   = m1;
            trk.k_dir  = k1;
            const auto kin1 = lattice_.map_k_to_vg(m1, k1);
            trk.vg_dir = kin1.vg_dir;
            trk.vg_mag = kin1.vg_mag;
        } else {
            stats.decays_L_to_LT++;
            const double x   = lattice_.sample_lt_energy_fraction(rng);
            const double thL = lattice_.theta_L_to_LT(x);
            const double thT = lattice_.theta_T_from_LT(x);

            const Vec3 kL = dir1.rotate_around(perp, std::cos(thL), std::sin(thL)).rotate_around(dir1, c_phi, s_phi).unit();
            const Vec3 kT = dir1.rotate_around(perp, std::cos(-thT), std::sin(-thT)).rotate_around(dir1, c_phi, s_phi).unit();

            const PhononMode mT = lattice_.sample_transverse_mode(rng);

            stack.push_back(make_track(trk.pos, kT, mT, (1.0 - x) * trk.energy, trk.time, trk.weight));

            trk.energy = x * trk.energy;
            trk.mode   = PhononMode::L;
            trk.k_dir  = kL;
            const auto kinL = lattice_.map_k_to_vg(PhononMode::L, kL);
            trk.vg_dir = kinL.vg_dir;
            trk.vg_mag = kinL.vg_mag;
        }
    }
};

} // namespace cmpmc
