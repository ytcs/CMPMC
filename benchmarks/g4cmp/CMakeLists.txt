# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 TC Yu and the CMPMC contributors
#
# G4CMP side of the comparison benchmark. Requires an installed Geant4 and a
# built G4CMP library:
#
#   cmake -S benchmarks/g4cmp -B build-g4cmp-bench \
#         -DGeant4_DIR=<geant4>/lib/cmake/Geant4 -DG4CMP_DIR=<g4cmp-install>
#   cmake --build build-g4cmp-bench

cmake_minimum_required(VERSION 3.16)
project(g4cmp_bench CXX)

find_package(Geant4 REQUIRED)
set(G4CMP_DIR "" CACHE PATH "G4CMP install prefix (contains include/G4CMP and lib/)")

find_path(G4CMP_INCLUDE_DIR G4CMPPhysicsList.hh
          HINTS ${G4CMP_DIR}/include/G4CMP ${G4CMP_DIR}/include REQUIRED)
find_library(G4CMP_LIBRARY NAMES G4cmp HINTS ${G4CMP_DIR}/lib ${G4CMP_DIR}/lib64 REQUIRED)
find_library(QHULL_LIBRARY NAMES qhullcpp HINTS ${G4CMP_DIR}/lib ${G4CMP_DIR}/lib64)
find_library(QHULL_R_LIBRARY NAMES qhull_p qhullstatic_p HINTS ${G4CMP_DIR}/lib ${G4CMP_DIR}/lib64)

add_executable(g4cmp_bench g4cmp_bench.cc)
target_include_directories(g4cmp_bench PRIVATE ${G4CMP_INCLUDE_DIR})
target_link_libraries(g4cmp_bench PRIVATE ${G4CMP_LIBRARY} ${Geant4_LIBRARIES})
if(QHULL_LIBRARY)
    target_link_libraries(g4cmp_bench PRIVATE ${QHULL_LIBRARY} ${QHULL_R_LIBRARY})
endif()
