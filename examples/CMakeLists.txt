# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 TC Yu and the CMPMC contributors

foreach(example quickstart custom_excitation phonon_benchmark)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE cmpmc::cmpmc cmpmc_build_options)
endforeach()
