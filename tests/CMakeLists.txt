# Each test program prints "name value reference rel_err PASS|FAIL" lines and
# exits nonzero if a check fails. See docs/validation.md for what each test
# checks and against which reference.

set(CMPMC_G4CMP_CRYSTALMAPS_DIR "" CACHE PATH
    "Optional: path to a G4CMP CrystalMaps directory, used to cross-check the material presets")

# name|comma-separated labels
set(CMPMC_TESTS
    "test_umbrella_header|unit"
    "test_christoffel|phonon,validation"
    "test_christoffel_cardano|phonon,validation"
    "test_anharmonic_decay|phonon,validation"
    "test_rates|phonon,validation"
    "test_surface|phonon,validation"
    "test_interface_amm|phonon,validation"
    "test_presets|materials,validation"
    "test_config_import|materials,unit"
    "test_g4cmp_import|materials,external-data"
    "test_kaplan|superconductor,validation"
    "test_pair_breaking|superconductor,validation"
    "test_luke_ntl|charge,validation"
    "test_leman_stepping|charge,validation"
    "test_charge_misc|charge,validation"
    "test_charge_release|charge,validation"
    "test_fano_partition|charge,validation"
    "test_engine_energy_closure|transport,unit"
    "test_guards|transport,unit"
)

foreach(entry IN LISTS CMPMC_TESTS)
    string(REPLACE "|" ";" parts "${entry}")
    list(GET parts 0 name)
    list(GET parts 1 labels)
    string(REPLACE "," ";" labels "${labels}")
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE cmpmc::cmpmc cmpmc_build_options)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        TIMEOUT 1800
        SKIP_RETURN_CODE 77)
endforeach()

set_tests_properties(test_g4cmp_import PROPERTIES
    ENVIRONMENT "CMPMC_G4CMP_CRYSTALMAPS_DIR=${CMPMC_G4CMP_CRYSTALMAPS_DIR}")

# Every public header must compile on its own (one translation unit per header).
file(GLOB_RECURSE CMPMC_PUBLIC_HEADERS RELATIVE ${PROJECT_SOURCE_DIR}/include
     CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/cmpmc/*.hpp)
set(CMPMC_HEADER_CHECK_SOURCES "")
foreach(header IN LISTS CMPMC_PUBLIC_HEADERS)
    string(MAKE_C_IDENTIFIER "${header}" stem)
    set(src ${CMAKE_CURRENT_BINARY_DIR}/header_check/${stem}.cpp)
    file(CONFIGURE OUTPUT ${src} CONTENT "#include \"${header}\"\n")
    list(APPEND CMPMC_HEADER_CHECK_SOURCES ${src})
endforeach()
add_library(cmpmc_header_check OBJECT ${CMPMC_HEADER_CHECK_SOURCES})
target_link_libraries(cmpmc_header_check PRIVATE cmpmc::cmpmc cmpmc_build_options)
