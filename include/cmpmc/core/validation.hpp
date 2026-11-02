#pragma once

#include <stdexcept>
#include <string>

namespace cmpmc {

// Runtime guard for experimental modules, i.e. modules without a validated
// reference path (listed as "experimental" in docs/validation.md).
//
// A guarded entry point throws UnvalidatedModuleError unless the caller opts in
// with DetectorModel::allow_unvalidated = true (or the explicit argument of a
// free function). Opting in is intended for benchmarking and development; results
// of a guarded module should not be used for physics conclusions.

// Thrown by require_validated() when an experimental module is used without opt-in.
struct UnvalidatedModuleError : std::logic_error {
    using std::logic_error::logic_error;
};

// Throws UnvalidatedModuleError naming `module` and `reason` unless allow_unvalidated is true.
inline void require_validated(bool allow_unvalidated, const char* module, const char* reason) {
    if (allow_unvalidated) return;
    throw UnvalidatedModuleError(std::string("CMPMC: '") + module + "' is not validated (" + reason +
                                 "; see docs/validation.md). Set allow_unvalidated = true to run it anyway.");
}

} // namespace cmpmc
