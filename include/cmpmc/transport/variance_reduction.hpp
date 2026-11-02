#pragma once

// Variance-reduction utilities used by the transport engine.

#include "cmpmc/core/rng.hpp"
#include <cstdint>

namespace cmpmc {

// Russian roulette for high-bounce tracks. Every 32 bounces beyond
// `threshold_bounces` a track survives with probability `survival_prob` and its
// weight is divided by it (unbiased). Returns false when the track is terminated.
// threshold_bounces == 0 disables the game.
inline bool apply_russian_roulette(double& weight,
                                   uint32_t bounces,
                                   uint32_t threshold_bounces,
                                   double survival_prob,
                                   Rng& rng) noexcept {
    if (threshold_bounces == 0 || bounces <= threshold_bounces) return true;
    if ((bounces - threshold_bounces) % 32 != 0) return true;
    if (rng.uniform() < survival_prob) {
        weight /= survival_prob;
        return true;
    }
    return false;
}

} // namespace cmpmc
