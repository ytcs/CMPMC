#pragma once
// Shared helpers for the CMPMC physics-audit tests (R15 CMPMC-AUD).
// Each test prints "name value reference rel_err status" lines and returns
// nonzero if any hard check fails. Define CMPMC_AUDIT_REPORT_ONLY to never fail
// (used when the same source is compiled against the unpatched 921db93 tree).
#include <cmath>
#include <cstdio>
#include <string>

namespace audit {
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
#ifdef CMPMC_AUDIT_REPORT_ONLY
    std::printf("# report-only build: %d check(s) failed\n", g_fail);
    return 0;
#else
    std::printf("# %d check(s) failed\n", g_fail);
    return g_fail ? 1 : 0;
#endif
}
} // namespace audit
