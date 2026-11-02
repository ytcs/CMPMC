#pragma once

#include "cmpmc/core/rng.hpp"
#include "cmpmc/core/track.hpp"
#include <concepts>
#include <vector>

namespace cmpmc {

// Requirements on a user-defined excitation species for
// SimulationEngine::transport_custom_excitation<Model>(). The engine owns free
// flight, boundary intersection and sensor absorption; the model supplies the
// bulk physics. Because the model is a template parameter, its methods are
// inlined into the transport loop (no virtual dispatch).
//
// See examples/custom_excitation.cpp and docs/extending.md.
template <typename Model>
concept ExcitationModel = requires(const Model& m,
                                   GenericExcitationTrack& trk,
                                   std::vector<GenericExcitationTrack>& secondaries,
                                   Rng& rng) {
    // 1. Bulk interaction rate Gamma(E) in [1/s]
    { m.bulk_rate(trk) } -> std::convertible_to<double>;

    // 2. Execute bulk scattering / decay / pair-breaking in-place
    { m.perform_bulk_interaction(trk, secondaries, rng) } -> std::same_as<void>;

    // 3. Minimum survival energy cutoff (eV)
    { m.min_energy_eV() } -> std::convertible_to<double>;
};

} // namespace cmpmc
