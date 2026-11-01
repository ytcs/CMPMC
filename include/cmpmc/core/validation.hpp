#pragma once

#include <stdexcept>
#include <string>

namespace cmpmc {

// ============================================================================
// Runtime guard for modules that are NOT validated against a reference
// (VALIDATION.md, "unsupported" rows; R15 CMPMC-AUD Sec. 4 DO-NOT-USE verdicts).
// A guarded entry point throws std::logic_error unless the caller opts in with
// DetectorModel::allow_unvalidated = true (or the explicit argument of a free
// function). Opting in is for benchmarking and development only; results of a
// guarded module must not feed a campaign number.
// ============================================================================
struct UnvalidatedModuleError : std::logic_error {
    using std::logic_error::logic_error;
};

inline void require_validated(bool allow_unvalidated, const char* module, const char* reason) {
    if (allow_unvalidated) return;
    throw UnvalidatedModuleError(std::string("CMPMC: '") + module + "' is not validated (" + reason +
                                 "; see VALIDATION.md). Set allow_unvalidated = true to run it anyway.");
}

} // namespace cmpmc
