cmake_minimum_required(VERSION 3.21)

# Single source of truth for the version: include/cmpmc/version.hpp
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/include/cmpmc/version.hpp" _cmpmc_version_hpp)
string(REGEX MATCH "CMPMC_VERSION_MAJOR ([0-9]+)" _ "${_cmpmc_version_hpp}")
set(_cmpmc_major ${CMAKE_MATCH_1})
string(REGEX MATCH "CMPMC_VERSION_MINOR ([0-9]+)" _ "${_cmpmc_version_hpp}")
set(_cmpmc_minor ${CMAKE_MATCH_1})
string(REGEX MATCH "CMPMC_VERSION_PATCH ([0-9]+)" _ "${_cmpmc_version_hpp}")
set(_cmpmc_patch ${CMAKE_MATCH_1})

project(cmpmc
    VERSION ${_cmpmc_major}.${_cmpmc_minor}.${_cmpmc_patch}
    DESCRIPTION "Monte Carlo transport of phonons, charge carriers and quasiparticles in cryogenic crystals"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(CMPMC_BUILD_TESTS        "Build the test suite"                              ${PROJECT_IS_TOP_LEVEL})
option(CMPMC_BUILD_EXAMPLES     "Build the example programs"                        ${PROJECT_IS_TOP_LEVEL})
option(CMPMC_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (tests/examples)" OFF)
option(CMPMC_NATIVE_ARCH        "Compile tests/examples with -march=native"         OFF)
option(CMPMC_INSTALL            "Generate install and package-config rules"         ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

# ---- Library (header-only) ---------------------------------------------------
add_library(cmpmc INTERFACE)
add_library(cmpmc::cmpmc ALIAS cmpmc)
target_compile_features(cmpmc INTERFACE cxx_std_20)
target_include_directories(cmpmc INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(cmpmc INTERFACE Threads::Threads)

# ---- Build settings for in-tree executables (not propagated to consumers) ----
add_library(cmpmc_build_options INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cmpmc_build_options INTERFACE -Wall -Wextra -Wpedantic)
    if(CMPMC_WARNINGS_AS_ERRORS)
        target_compile_options(cmpmc_build_options INTERFACE -Werror)
    endif()
    if(CMPMC_NATIVE_ARCH)
        target_compile_options(cmpmc_build_options INTERFACE -march=native)
    endif()
elseif(MSVC)
    target_compile_options(cmpmc_build_options INTERFACE /W4 /permissive-)
    if(CMPMC_WARNINGS_AS_ERRORS)
        target_compile_options(cmpmc_build_options INTERFACE /WX)
    endif()
endif()

if(CMPMC_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(CMPMC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---- Install / package config ------------------------------------------------
if(CMPMC_INSTALL)
    install(TARGETS cmpmc EXPORT cmpmcTargets)
    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(EXPORT cmpmcTargets
        NAMESPACE cmpmc::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cmpmc)
    configure_package_config_file(cmake/cmpmcConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/cmpmcConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cmpmc)
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cmpmcConfigVersion.cmake
        COMPATIBILITY SameMinorVersion
        ARCH_INDEPENDENT)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmpmcConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmpmcConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cmpmc)
endif()
