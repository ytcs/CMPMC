// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

// Energy partition of a deposit into charge carriers and prompt phonons, and the
// step accumulator / hit merger for primary-particle input.

#include "cmpmc/charge/charge_transport.hpp"
#include "cmpmc/core/core.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"
#include "cmpmc/charge/niel_models.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cmpmc {

// EnergyPartition: splits an energy deposit into electron-hole pairs and prompt
// phonons, using one of the nuclear-recoil ionisation-yield models of
// cmpmc/charge/niel_models.hpp and the Fano-binomial pair-number sampler.
struct PartitionOptions {
    physics::NIELModel niel_model = physics::NIELModel::Lindhard;
    physics::EmpiricalNIELParams empirical_niel{};
    bool use_fano_fluctuations   = true;
    bool use_charge_cloud        = true;  // Distribute e-/h+ pairs in a 3D cloud
    uint32_t max_sampled_pairs   = 1000;  // Downsample e-/h+ pairs above this count (with weight)
    uint32_t max_sampled_phonons = 5000;  // Downsample prompt phonons above this count (with weight)
};

class EnergyPartition {
public:
    using Options = PartitionOptions;

    struct Result {
        double total_edep_eV     = 0.0;
        double ionizing_eV       = 0.0;
        double niel_eV           = 0.0;
        double lindhard_yield    = 1.0;
        uint64_t n_pairs_true    = 0;
        uint64_t n_pairs_sampled = 0;
        uint64_t n_phonons_true  = 0;
        uint64_t n_phonons_sampled = 0;
        double charge_weight     = 1.0;
        double phonon_weight     = 1.0;
        double prompt_phonon_eV  = 0.0;
    };

    explicit constexpr EnergyPartition(const CrystalLattice& lattice) noexcept
        : lattice_(lattice), charge_phys_(lattice) {}

    // Lindhard ionisation yield Y(E) (Lindhard 1963, in the Lewin & Smith 1996 form;
    // see cmpmc/charge/niel_models.hpp)
    double lindhard_ionization_yield(double energy_eV, double z1 = 0.0, double a1 = 0.0) const noexcept {
        if (energy_eV <= 0.0) return 0.0;
        const double z2 = lattice_.config().atomic_Z;
        const double a2 = lattice_.config().atomic_A;
        if (z1 <= 0.0) z1 = z2; // Self-recoil default (e.g. neutron scatter on Ge/Si nucleus)
        if (a1 <= 0.0) a1 = a2;

        const double z1_23 = std::cbrt(z1 * z1);
        const double z2_23 = std::cbrt(z2 * z2);
        const double zpow  = z1_23 + z2_23;
        const double asum  = a1 + a2;

        const double el = 30.724 * units::eV * z1 * z2 * std::sqrt(zpow) * asum / a2;
        const double fl = (0.0793 * z1_23 * std::sqrt(z2 * asum * asum * asum / (a1 * a1 * a1 * a2)))
                          / std::pow(zpow, 0.75);
        const double eps = energy_eV / el;

        const double numer = fl * (eps + 3.4008 * std::pow(eps, 0.16667) + 0.40244 * std::pow(eps, 0.75));
        return numer / (1.0 + numer);
    }

    // Partition an energy deposit (Electronic Recoil or Nuclear Recoil) into
    // conduction electrons, valence holes, and prompt Debye phonons.
    template <typename PhononContainer>
    Result partition_deposit(const Vec3& pos,
                             double time_s,
                             double total_edep_eV,
                             bool is_nuclear_recoil,
                             std::vector<ChargeTrack>& out_charges,
                             PhononContainer& out_phonons,
                             Rng& rng,
                             const Options& opts = Options{}) const {
        Result res{};
        res.total_edep_eV  = total_edep_eV;
        res.lindhard_yield = is_nuclear_recoil
            ? physics::compute_ionization_yield(opts.niel_model, total_edep_eV, lattice_.config(), opts.empirical_niel)
            : 1.0;
        res.ionizing_eV    = total_edep_eV * res.lindhard_yield;
        res.niel_eV        = total_edep_eV - res.ionizing_eV;

        const double e_gap  = lattice_.config().bandgap;
        const double e_pair = lattice_.config().pair_energy;
        const double fano   = lattice_.config().fano_factor;

        // 1. Sample number of e-/h+ pairs with the interpolated Fano binomial sampler
        if (res.ionizing_eV >= e_gap && e_pair > 0.0) {
            const double mean_pairs = res.ionizing_eV / e_pair;
            if (!opts.use_fano_fluctuations || mean_pairs < 1.0) {
                res.n_pairs_true = static_cast<uint64_t>(std::floor(mean_pairs + rng.uniform()));
            } else {
                const uint64_t max_pairs = static_cast<uint64_t>(std::floor(res.ionizing_eV / e_gap));
                res.n_pairs_true = std::min(max_pairs, physics::sample_fano_binomial(mean_pairs, fano, rng));
            }
        }

        // Each e-/h+ pair locks up `e_gap` in bandgap potential energy;
        // also give each carrier a small initial thermal/hot kinetic energy (e.g. 5 meV)
        const double init_carrier_ekin = 5.0 * units::meV;
        const double total_pair_energy = static_cast<double>(res.n_pairs_true) * (e_gap + 2.0 * init_carrier_ekin);
        res.prompt_phonon_eV = std::max(0.0, total_edep_eV - total_pair_energy);

        // 2. Spawn sampled e-/h+ pairs (with optional 3D charge-cloud spatial distribution)
        if (res.n_pairs_true > 0) {
            res.n_pairs_sampled = std::min<uint64_t>(res.n_pairs_true, opts.max_sampled_pairs);
            res.charge_weight   = static_cast<double>(res.n_pairs_true) / res.n_pairs_sampled;
            const size_t n_valleys = std::max<size_t>(1, lattice_.num_valleys());
            const double r_cloud = opts.use_charge_cloud
                ? physics::ChargeCloudGenerator::compute_cloud_radius(res.n_pairs_true, lattice_.config())
                : 0.0;

            for (uint64_t i = 0; i < res.n_pairs_sampled; ++i) {
                const Vec3 pair_pos = pos + physics::ChargeCloudGenerator::sample_offset(r_cloud, rng);
                const uint16_t valley = static_cast<uint16_t>(rng.next_u64() % n_valleys);
                ChargeTrack elec = charge_phys_.make_charge(
                    ExcitationType::Electron, pair_pos, rng.isotropic_dir(), init_carrier_ekin, valley, time_s);
                elec.weight = res.charge_weight;

                ChargeTrack hole = charge_phys_.make_charge(
                    ExcitationType::Hole, pair_pos, rng.isotropic_dir(), init_carrier_ekin, 0, time_s);
                hole.weight = res.charge_weight;

                out_charges.push_back(elec);
                out_charges.push_back(hole);
            }
        }

        // 3. Spawn prompt Debye-energy phonons from remaining recoil + relaxation energy
        const double e_debye = lattice_.debye_energy();
        if (res.prompt_phonon_eV > 0.0 && e_debye > 0.0) {
            res.n_phonons_true = std::max<uint64_t>(1, static_cast<uint64_t>(std::round(res.prompt_phonon_eV / e_debye)));
            res.n_phonons_sampled = std::min<uint64_t>(res.n_phonons_true, opts.max_sampled_phonons);
            res.phonon_weight = static_cast<double>(res.n_phonons_true) / res.n_phonons_sampled;
            const double single_phonon_eV = res.prompt_phonon_eV / res.n_phonons_true;

            for (uint64_t i = 0; i < res.n_phonons_sampled; ++i) {
                const PhononMode mode = lattice_.sample_dos_mode(rng);
                const Vec3 k_dir      = rng.isotropic_dir();
                const auto kin        = lattice_.map_k_to_vg(mode, k_dir);
                out_phonons.push_back({
                    .pos     = pos,
                    .k_dir   = k_dir,
                    .vg_dir  = kin.vg_dir,
                    .energy  = single_phonon_eV,
                    .vg_mag  = kin.vg_mag,
                    .time    = time_s,
                    .weight  = res.phonon_weight,
                    .bounces = 0,
                    .mode    = mode
                });
            }
        }

        return res;
    }

private:
    const CrystalLattice& lattice_;
    ChargePhysics charge_phys_;
};

// PrimaryStepInput, StepAccumulator and HitMerger: interface for feeding primary-
// particle steps (e.g. from an external particle-transport code) into CMPMC.
//
// Consolidates adjacent continuous dE/dx or nuclear recoil steps along a
// primary particle track within `combining_step_length` (default 1 mm),
// partitions the merged deposit into e-/h+ pairs and prompt phonons, and
// distributes secondaries along the trajectory chord with Gaussian jitter
// (dl +- dl/6) for charged tracks or at the endpoint for neutral tracks.
struct PrimaryStepInput {
    int32_t track_id       = 1;
    Vec3 start_pos         = {0.0, 0.0, 0.0};
    Vec3 end_pos           = {0.0, 0.0, 0.0};
    double time_s          = 0.0;
    double edep_eV         = 0.0; // Electronic / total energy deposit (eV)
    double eniel_eV        = 0.0; // Non-ionizing nuclear recoil deposit (eV)
    bool is_charged        = true;
    bool hit_boundary      = false;
    bool track_stopped     = false;
};

struct StepAccumulator {
    int32_t track_id   = -1;
    uint32_t n_steps   = 0;
    Vec3 start_pos     = {0.0, 0.0, 0.0};
    Vec3 end_pos       = {0.0, 0.0, 0.0};
    double time_s      = 0.0;
    double edep_eV     = 0.0;
    double eniel_eV    = 0.0;
    bool is_charged    = true;

    constexpr void clear() noexcept {
        n_steps    = 0;
        start_pos  = {0.0, 0.0, 0.0};
        end_pos    = {0.0, 0.0, 0.0};
        time_s     = 0.0;
        edep_eV    = 0.0;
        eniel_eV   = 0.0;
        is_charged = true;
    }

    inline void add(const PrimaryStepInput& step) noexcept {
        if (n_steps == 0) {
            track_id   = step.track_id;
            start_pos  = step.start_pos;
            is_charged = step.is_charged;
        }
        end_pos   = step.end_pos;
        time_s    = std::max(time_s, step.time_s);
        edep_eV  += step.edep_eV;
        eniel_eV += step.eniel_eV;
        n_steps++;
    }
};

class HitMerger {
public:
    explicit HitMerger(const CrystalLattice& lattice,
                       double combining_step_length_m = 1.0 * units::mm,
                       const PartitionOptions& opts = PartitionOptions{})
        : partitioner_(lattice),
          combining_step_length_(combining_step_length_m),
          opts_(opts) {}

    // Process a single primary track step; flushes accumulated steps into
    // `out_charges` and `out_phonons` whenever distance >= combining_step_length,
    // boundary is reached, or track stops.
    bool process_step(const PrimaryStepInput& step,
                      std::vector<ChargeTrack>& out_charges,
                      std::vector<PhononTrack>& out_phonons,
                      Rng& rng) {
        const double dist = (step.end_pos - accum_.end_pos).norm();
        const bool do_add = (accum_.n_steps == 0 ||
                             dist < combining_step_length_ ||
                             step.hit_boundary ||
                             step.track_stopped);

        const bool has_energy = (step.edep_eV > 0.0 || step.eniel_eV > 0.0);
        if (do_add && has_energy) {
            accum_.add(step);
        }

        const bool ready = (accum_.n_steps > 0) &&
                           ((step.end_pos - accum_.end_pos).norm() >= combining_step_length_ ||
                            !do_add ||
                            step.hit_boundary ||
                            step.track_stopped);

        if (ready) {
            flush(out_charges, out_phonons, rng);
        }

        if (!do_add && has_energy) {
            accum_.add(step);
        }

        return ready;
    }

    // Flush any remaining accumulated energy deposits at the end of an event
    EnergyPartition::Result flush(std::vector<ChargeTrack>& out_charges,
                                  std::vector<PhononTrack>& out_phonons,
                                  Rng& rng) {
        EnergyPartition::Result res{};
        if (accum_.n_steps == 0 || (accum_.edep_eV <= 0.0 && accum_.eniel_eV <= 0.0)) {
            accum_.clear();
            return res;
        }

        const size_t c_start = out_charges.size();
        const size_t p_start = out_phonons.size();

        const double total_e = accum_.edep_eV + accum_.eniel_eV;
        const bool is_nr     = (accum_.eniel_eV > 0.5 * total_e);
        const Vec3 mid_pos   = (accum_.start_pos + accum_.end_pos) * 0.5;

        res = partitioner_.partition_deposit(
            mid_pos, accum_.time_s, total_e, is_nr, out_charges, out_phonons, rng, opts_);

        // For charged tracks with non-zero chord length, spread generated secondaries
        // along [start_pos, end_pos] with Gaussian step jitter dl +- dl/6
        const Vec3 traj   = accum_.end_pos - accum_.start_pos;
        const double len  = traj.norm();
        const size_t nsec = (out_charges.size() - c_start) + (out_phonons.size() - p_start);

        if (accum_.is_charged && len > 1e-7 && nsec > 1) {
            const Vec3 tdir   = traj / len;
            const double dl   = len / static_cast<double>(nsec);
            const double sigl = dl / 6.0;
            Vec3 cur_pos      = accum_.start_pos;

            for (size_t i = c_start; i < out_charges.size(); ++i) {
                cur_pos += tdir * std::max(0.0, dl + sigl * rng.normal());
                out_charges[i].pos = cur_pos;
            }
            for (size_t i = p_start; i < out_phonons.size(); ++i) {
                cur_pos += tdir * std::max(0.0, dl + sigl * rng.normal());
                out_phonons[i].pos = cur_pos;
            }
        }

        accum_.clear();
        return res;
    }

private:
    EnergyPartition partitioner_;
    double combining_step_length_;
    PartitionOptions opts_;
    StepAccumulator accum_{};
};

} // namespace cmpmc
