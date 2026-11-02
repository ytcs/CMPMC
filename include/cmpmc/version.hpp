#pragma once

// Library version. CMakeLists.txt reads the project version from this file.
#define CMPMC_VERSION_MAJOR 0
#define CMPMC_VERSION_MINOR 1
#define CMPMC_VERSION_PATCH 0
#define CMPMC_VERSION_STRING "0.1.0"

namespace cmpmc {

struct Version {
    int major;
    int minor;
    int patch;
};

inline constexpr Version version{CMPMC_VERSION_MAJOR, CMPMC_VERSION_MINOR, CMPMC_VERSION_PATCH};
inline constexpr const char* version_string = CMPMC_VERSION_STRING;

} // namespace cmpmc
