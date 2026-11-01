#pragma once

#include "cmpmc/charge.hpp"
#include "cmpmc/core.hpp"
#include "cmpmc/lattice.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cmpmc {

// ============================================================================
// LindhardNIEL & EnergyPartition
// Replaces: G4CMPEnergyPartition, G4CMPLindhardNIEL, G4CMPFanoBinomial
// ============================================================================
struct PartitionOptions {
    bool use_fano_fluctuations   = true;
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

    // Lewin & Smith (1996) / Lindhard & Robinson Ionization Yield Y(E)
    // Matches G4CMPLindhardNIEL::PartitionNIEL
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
        res.lindhard_yield = is_nuclear_recoil ? lindhard_ionization_yield(total_edep_eV) : 1.0;
        res.ionizing_eV    = total_edep_eV * res.lindhard_yield;
        res.niel_eV        = total_edep_eV - res.ionizing_eV;

        const double e_gap  = lattice_.config().bandgap;
        const double e_pair = lattice_.config().pair_energy;
        const double fano   = lattice_.config().fano_factor;

        // 1. Sample number of e-/h+ pairs with Fano fluctuations
        if (res.ionizing_eV >= e_gap && e_pair > 0.0) {
            const double mean_pairs = res.ionizing_eV / e_pair;
            if (!opts.use_fano_fluctuations || mean_pairs < 1.0) {
                res.n_pairs_true = static_cast<uint64_t>(std::floor(mean_pairs + rng.uniform()));
            } else {
                const double sigma = std::sqrt(fano * mean_pairs);
                const double sampled = std::round(mean_pairs + sigma * rng.normal());
                res.n_pairs_true = static_cast<uint64_t>(std::clamp(sampled, 1.0, std::floor(res.ionizing_eV / e_gap)));
            }
        }

        // Each e-/h+ pair locks up `e_gap` in bandgap potential energy;
        // also give each carrier a small initial thermal/hot kinetic energy (e.g. 5 meV)
        const double init_carrier_ekin = 5.0 * units::meV;
        const double total_pair_energy = static_cast<double>(res.n_pairs_true) * (e_gap + 2.0 * init_carrier_ekin);
        res.prompt_phonon_eV = std::max(0.0, total_edep_eV - total_pair_energy);

        // 2. Spawn sampled e-/h+ pairs
        if (res.n_pairs_true > 0) {
            res.n_pairs_sampled = std::min<uint64_t>(res.n_pairs_true, opts.max_sampled_pairs);
            res.charge_weight   = static_cast<double>(res.n_pairs_true) / res.n_pairs_sampled;
            const size_t n_valleys = std::max<size_t>(1, lattice_.num_valleys());

            for (uint64_t i = 0; i < res.n_pairs_sampled; ++i) {
                const uint16_t valley = static_cast<uint16_t>(rng.next_u64() % n_valleys);
                ChargeTrack elec = charge_phys_.make_charge(
                    ExcitationType::Electron, pos, rng.isotropic_dir(), init_carrier_ekin, valley, time_s);
                elec.weight = res.charge_weight;

                ChargeTrack hole = charge_phys_.make_charge(
                    ExcitationType::Hole, pos, rng.isotropic_dir(), init_carrier_ekin, 0, time_s);
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

} // namespace cmpmc
