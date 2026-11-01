# Physics-audit tests (R15 CMPMC-AUD). Each test prints value / reference / rel. error
# lines and exits nonzero if a check fails.
set(CMPMC_AUDIT_TESTS
    test_anharmonic_decay
    test_christoffel
    test_presets
    test_rates
    test_luke_ntl
    test_charge_misc
    test_fano_partition
    test_kaplan
    test_surface
    test_headers_compile
    test_pair_breaking
    test_interface_amm
    test_guards
    test_charge_release
    test_leman_stepping
)
foreach(t ${CMPMC_AUDIT_TESTS})
    add_executable(${t} ${t}.cpp)
    target_link_libraries(${t} PRIVATE cmpmc)
    target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${t} PRIVATE -O2 -Wall -Wextra)
    endif()
    add_test(NAME ${t} COMMAND ${t})
    set_tests_properties(${t} PROPERTIES TIMEOUT 1800)
endforeach()
