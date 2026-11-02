// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once
// Minimal check/report helpers shared by the test programs.
//
// Each check prints one line "name value reference rel_err PASS|FAIL"; the program
// returns nonzero from result() if any check failed. A test that needs external
// data which is not available calls skip(), which returns the CTest skip code.
#include <cmath>
#include <cstdio>
#include <string>

namespace cmpmc_test {

inline constexpr int kSkipReturnCode = 77;

inline int g_fail = 0;

inline void check(const std::string& name, double val, double ref, double tol_rel, double tol_abs = 0.0) {
    const double err = std::abs(val - ref);
    const double rel = (ref != 0.0) ? err / std::abs(ref) : err;
    const bool ok = (err <= tol_abs) || (rel <= tol_rel);
    std::printf("%-58s %14.6g %14.6g %10.3g %s\n", name.c_str(), val, ref, rel, ok ? "PASS" : "FAIL");
    if (!ok) ++g_fail;
}

inline void info(const std::string& name, double val) {
    std::printf("%-58s %14.6g\n", name.c_str(), val);
}

inline int result() {
    std::printf("# %d check(s) failed\n", g_fail);
    return g_fail ? 1 : 0;
}

inline int skip(const std::string& reason) {
    std::printf("# skipped: %s\n", reason.c_str());
    return kSkipReturnCode;
}

} // namespace cmpmc_test
