// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

// Every public header compiles in one translation unit: the umbrella header
// (supported API) together with the experimental headers, and the version macros
// agree with the constexpr version object.
#include "test_support.hpp"
#include "cmpmc/cmpmc.hpp"
#include "cmpmc/experimental/jump_to_ballistic.hpp"
#include "cmpmc/experimental/next_event_estimation.hpp"
#include "cmpmc/experimental/wavefront_soa.hpp"
#include "cmpmc/experimental/wavefront_pipeline.hpp"

#include <string>

int main() {
    std::printf("# test_umbrella_header\n");
    cmpmc_test::check("version major", cmpmc::version.major, CMPMC_VERSION_MAJOR, 0.0);
    cmpmc_test::check("version minor", cmpmc::version.minor, CMPMC_VERSION_MINOR, 0.0);
    cmpmc_test::check("version patch", cmpmc::version.patch, CMPMC_VERSION_PATCH, 0.0);
    const std::string expected = std::to_string(CMPMC_VERSION_MAJOR) + "." + std::to_string(CMPMC_VERSION_MINOR) +
                                 "." + std::to_string(CMPMC_VERSION_PATCH);
    cmpmc_test::check("version string matches macros", expected == cmpmc::version_string ? 1.0 : 0.0, 1.0, 0.0);
    return cmpmc_test::result();
}
