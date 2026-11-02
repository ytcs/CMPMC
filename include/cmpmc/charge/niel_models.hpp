#pragma once

// Nuclear-recoil ionisation-yield models and Fano pair-number sampling.

#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/units.hpp"
#include "cmpmc/materials/crystal_config.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace cmpmc::physics {

// Nuclear-recoil ionisation-yield (NIEL partition) models and the Fano binomial
// pair-number sampler:
//   1. Lindhard & Robinson (1963)
//   2. Lewin & Smith (1996)
//   3. Sarkis et al., PRA 107, 062811 (2023)   (silicon, E >= 47.8 eV)
//   4. IMPACT@TUNL, PRL 131, 091801 (2023)     (silicon, E <= 10 keV)
//   5. CDMSlite empirical k(E), PRD 105 (2022)
// plus the interpolated Fano binomial sampler (CDMS HVeV Run 1 data release, 2019).
// The model set and parameterisations follow G4CMP.

enum class NIELModel : uint8_t {
    Lindhard   = 0,
    LewinSmith = 1,
    Sarkis     = 2,
    ImpactTUNL = 3,
    Empirical  = 4
};

struct EmpiricalNIELParams {
    bool energy_dependent_k = true;
    double k_fixed          = 0.157;
    double k_low            = 0.080;
    double k_high           = 0.157;
    double e_low_eV         = 100.0 * units::eV;
    double e_high_eV        = 2000.0 * units::eV;
};

// --- 1. Lindhard & Robinson (1963) Ionization Yield ---
inline double lindhard_ionization_yield(double energy_eV,
                                        const CrystalConfig& cfg,
                                        double z1 = 0.0,
                                        double a1 = 0.0) noexcept {
    if (energy_eV <= 0.0) return 0.0;
    const double z2 = cfg.atomic_Z;
    const double a2 = cfg.atomic_A;
    if (z1 <= 0.0) z1 = z2;
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

// --- 2. Lewin & Smith (1996) Ionization Yield ---
inline double lewin_smith_ionization_yield(double energy_eV, const CrystalConfig& cfg) noexcept {
    if (energy_eV <= 0.0) return 0.0;
    const double Z = cfg.atomic_Z;
    const double A = cfg.atomic_A;
    const double z23 = std::cbrt(Z * Z);
    const double epsilon = (11.5 / units::keV) * energy_eV / (Z * z23 * z23); // * Z^(-7/3)
    const double k = 0.133 * z23 / std::sqrt(A);
    const double h = 0.7 * std::pow(epsilon, 0.6) + 3.0 * std::pow(epsilon, 0.15) + epsilon;
    return (k * h) / (1.0 + k * h);
}

// --- 3. IMPACT@TUNL (2023) Silicon Power-Law Yield ---
// Reference: PRL 131, 091801 (2023) — Y(E_r) = Y10 * (E_r / 10 keV)^B for Si <= 10 keV
inline double impact_tunl_ionization_yield(double energy_eV, const CrystalConfig& cfg) noexcept {
    if (energy_eV <= 0.0) return 0.0;
    const bool is_silicon = (std::abs(cfg.atomic_Z - 14.0) < 0.5) || (std::abs(cfg.atomic_A - 28.09) < 0.5);
    if (is_silicon && energy_eV <= 10.0 * units::keV) {
        constexpr double Y10 = 0.302;
        constexpr double B   = 0.261;
        return Y10 * std::pow(energy_eV / (10.0 * units::keV), B);
    }
    return lewin_smith_ionization_yield(energy_eV, cfg);
}

// --- 4. Sarkis et al. (2023) Silicon Ionization Yield ---
// Reference: Phys. Rev. A 107, 062811 (2023) — Silicon nuclear recoil yield down to 47.83 eV
inline double sarkis_ionization_yield(double energy_eV, const CrystalConfig& cfg) noexcept {
    if (energy_eV <= 0.0) return 0.0;
    const bool is_silicon = (std::abs(cfg.atomic_Z - 14.0) < 0.5) && (std::abs(cfg.atomic_A - 28.09) < 1.0);
    if (!is_silicon || energy_eV > 3.0 * units::MeV) {
        return lewin_smith_ionization_yield(energy_eV, cfg);
    }

    struct Knot { double e_eV; double yield; };
    static constexpr std::array<Knot, 64> knots = {{
        {47.83270, 0.000000}, {52.76735, 0.009114}, {58.75143, 0.016981}, {65.97224, 0.023622},
        {74.67257, 0.029413}, {85.16158, 0.034885}, {97.80973, 0.040365}, {113.05138, 0.045962},
        {131.41991, 0.051702}, {153.56956, 0.057597}, {180.26945, 0.063652}, {212.44431, 0.069870},
        {251.21921, 0.076252}, {297.95152, 0.082797}, {354.26532, 0.089510}, {422.11910, 0.096396},
        {503.86504, 0.103464}, {602.35199, 0.110722}, {720.98828, 0.118187}, {863.90755, 0.125872},
        {1036.03980, 0.133800}, {1243.36450, 0.141991}, {1493.04550, 0.150478}, {1793.74340, 0.159282},
        {2155.85770, 0.168439}, {2591.89340, 0.177980}, {3116.95740, 0.187961}, {3749.18370, 0.198399},
        {4510.42640, 0.209326}, {5427.00620, 0.220779}, {6530.63220, 0.232847}, {7859.32390, 0.245567},
        {9445.41000, 0.258878}, {11368.94900, 0.273042}, {13685.03100, 0.288012}, {16473.85900, 0.303838},
        {19832.02600, 0.320553}, {23876.07900, 0.338204}, {28746.29200, 0.356826}, {34611.68400, 0.376451},
        {41676.09800, 0.397077}, {50184.87400, 0.418697}, {60434.10300, 0.441291}, {72780.28500, 0.464810},
        {87653.00800, 0.489194}, {105570.14000, 0.514350}, {127155.61000, 0.540142}, {153161.31000, 0.566445},
        {184493.27000, 0.593129}, {222243.38000, 0.620053}, {267727.15000, 0.647094}, {322529.80000, 0.674139},
        {388561.59000, 0.701059}, {468124.30000, 0.727705}, {563991.37000, 0.753920}, {679504.88000, 0.779502},
        {818691.76000, 0.804197}, {986404.40000, 0.827705}, {1188489.80000, 0.849716}, {1431993.30000, 0.869965},
        {1725404.10000, 0.888262}, {2078951.70000, 0.904513}, {2504962.10000, 0.918716}, {3018287.70000, 0.930944}
    }};

    if (energy_eV <= knots.front().e_eV) return 0.0;
    if (energy_eV >= knots.back().e_eV)  return knots.back().yield;

    auto it = std::upper_bound(knots.begin(), knots.end(), energy_eV,
                               [](double e, const Knot& k) { return e < k.e_eV; });
    const size_t idx = static_cast<size_t>(std::distance(knots.begin(), it));
    const Knot& k0 = knots[idx - 1];
    const Knot& k1 = knots[idx];
    const double t = (energy_eV - k0.e_eV) / (k1.e_eV - k0.e_eV);
    return k0.yield + t * (k1.yield - k0.yield);
}

// --- 5. CDMSlite Empirical k(E) Ionization Yield ---
// Reference: Phys. Rev. D 105, 122002 (2022)
inline double empirical_ionization_yield(double energy_eV,
                                         const CrystalConfig& cfg,
                                         const EmpiricalNIELParams& emp = {}) noexcept {
    if (energy_eV <= 0.0) return 0.0;
    if (emp.energy_dependent_k && (energy_eV > emp.e_high_eV || emp.e_low_eV == emp.e_high_eV)) {
        return lewin_smith_ionization_yield(energy_eV, cfg);
    }

    double k = emp.k_fixed;
    if (emp.energy_dependent_k) {
        const double dk = emp.k_high - emp.k_low;
        const double dE = emp.e_high_eV - emp.e_low_eV;
        k = std::max(0.0, emp.k_low + (energy_eV - emp.e_low_eV) * (dk / dE));
    }

    const double Z = cfg.atomic_Z;
    const double z23 = std::cbrt(Z * Z);
    const double epsilon = (11.5 / units::keV) * energy_eV / (Z * z23 * z23);
    const double h = 0.7 * std::pow(epsilon, 0.6) + 3.0 * std::pow(epsilon, 0.15) + epsilon;
    return (k * h) / (1.0 + k * h);
}

// --- Unified NIEL Dispatcher ---
inline double compute_ionization_yield(NIELModel model,
                                       double energy_eV,
                                       const CrystalConfig& cfg,
                                       const EmpiricalNIELParams& emp = {}) noexcept {
    switch (model) {
        case NIELModel::Lindhard:   return lindhard_ionization_yield(energy_eV, cfg);
        case NIELModel::LewinSmith: return lewin_smith_ionization_yield(energy_eV, cfg);
        case NIELModel::Sarkis:     return sarkis_ionization_yield(energy_eV, cfg);
        case NIELModel::ImpactTUNL: return impact_tunl_ionization_yield(energy_eV, cfg);
        case NIELModel::Empirical:  return empirical_ionization_yield(energy_eV, cfg, emp);
    }
    return lindhard_ionization_yield(energy_eV, cfg);
}

// --- 6. Interpolated Fano Binomial Sampler ---
// Reference: CDMS HVeV Run 1 data release (2019)
namespace detail {

inline double choose_binom(int64_t n, int64_t x) noexcept {
    if (x > n / 2) x = n - x;
    if (x < 0) return 0.0;
    if (x == 0) return 1.0;
    if (x == 1) return static_cast<double>(n);
    double c = 1.0;
    for (int64_t k = 1; k <= x; ++k) {
        c *= static_cast<double>(n - k + 1) / static_cast<double>(k);
    }
    return c;
}

inline double pdf_binomial(int64_t x, int64_t n, double p) noexcept {
    if (x > n || x < 0 || p <= 0.0 || p > 1.0) return 0.0;
    if (p == 1.0) return (x == n) ? 1.0 : 0.0;
    const double lnbin = std::log(p) * static_cast<double>(x) +
                         std::log(1.0 - p) * static_cast<double>(n - x);
    if (lnbin < -708.3) return 0.0;
    if (lnbin > 709.7)  return 1e300;
    return choose_binom(n, x) * std::exp(lnbin);
}

} // namespace detail

inline uint64_t sample_fano_binomial(double mean, double fano, Rng& rng) noexcept {
    if (mean <= 0.0) return 0;
    if (mean <= 1.0) return 1;

    // Fano == 1 -> Poisson distribution
    if (std::abs(fano - 1.0) < 1e-9) {
        if (mean > 30.0) {
            return static_cast<uint64_t>(std::max(0.0, std::round(mean + std::sqrt(mean) * rng.normal())));
        }
        const double L = std::exp(-mean);
        uint64_t k = 0;
        double p = 1.0;
        do {
            ++k;
            p *= rng.uniform();
        } while (p > L);
        return k - 1;
    }

    const double prob = std::clamp(1.0 - fano, 1e-6, 1.0 - 1e-6);
    const double ntry = mean / prob;

    // Gaussian limit for large N (same threshold as G4CMP)
    if (mean > 9.0 * (fano / prob) && mean > 9.0 * (prob / fano)) {
        const double sampled = std::round(mean + std::sqrt(mean * fano) * rng.normal());
        return static_cast<uint64_t>(std::max(1.0, sampled));
    }

    // Interpolated Binomial distribution for small N
    const int64_t nlo = static_cast<int64_t>(std::floor(ntry));
    const int64_t nhi = static_cast<int64_t>(std::ceil(ntry));
    if (nlo == nhi) {
        uint64_t hits = 0;
        for (int64_t i = 0; i < nlo; ++i) {
            if (rng.uniform() < prob) ++hits;
        }
        return hits;
    }

    const double p_lo = std::min(mean / static_cast<double>(nlo), 1.0);
    const double p_hi = std::min(mean / static_cast<double>(nhi), 1.0);
    const double dp   = ( std::abs(p_hi - p_lo) > 1e-15 ) ? ((prob - p_lo) / (p_hi - p_lo)) : 0.0;

    const double pdf_max = std::max(
        detail::pdf_binomial(static_cast<int64_t>(std::floor(mean)), nlo, p_lo),
        detail::pdf_binomial(static_cast<int64_t>(std::ceil(mean)),  nlo, p_lo)
    );
    if (pdf_max <= 0.0) {
        return static_cast<uint64_t>(std::round(mean));
    }

    for (int iter = 0; iter < 1000; ++iter) {
        const int64_t pick = static_cast<int64_t>(rng.uniform() * static_cast<double>(nhi + 1));
        const double pick_p = (detail::pdf_binomial(pick, nlo, p_lo) * (1.0 - dp) +
                               detail::pdf_binomial(pick, nhi, p_hi) * dp) / pdf_max;
        if (rng.uniform() <= pick_p) {
            return static_cast<uint64_t>(std::max<int64_t>(0, pick));
        }
    }
    return static_cast<uint64_t>(std::round(mean));
}

} // namespace cmpmc::physics
