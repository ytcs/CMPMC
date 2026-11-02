// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors

#pragma once

// Umbrella header for the supported CMPMC API.
//
// Experimental modules (cmpmc/experimental/*) are deliberately not included;
// include them explicitly. See docs/experimental.md.

#include "cmpmc/version.hpp"

#include "cmpmc/core/core.hpp"

#include "cmpmc/materials/crystal_config.hpp"
#include "cmpmc/materials/presets.hpp"
#include "cmpmc/io/g4cmp_config_import.hpp"
#include "cmpmc/lattice/crystal_lattice.hpp"

#include "cmpmc/phonon/phonon_kernels.hpp"
#include "cmpmc/phonon/interface_transmission.hpp"

#include "cmpmc/charge/charge_kernels.hpp"
#include "cmpmc/charge/charge_transport.hpp"
#include "cmpmc/charge/niel_models.hpp"
#include "cmpmc/charge/energy_partition.hpp"

#include "cmpmc/superconductor/kaplan_film.hpp"
#include "cmpmc/superconductor/bulk_superconductor.hpp"

#include "cmpmc/geometry/solids.hpp"
#include "cmpmc/geometry/surface_property.hpp"
#include "cmpmc/geometry/mesh_field.hpp"
#include "cmpmc/geometry/electrode_pattern.hpp"

#include "cmpmc/transport/excitation_model.hpp"
#include "cmpmc/transport/variance_reduction.hpp"
#include "cmpmc/transport/engine.hpp"
#include "cmpmc/transport/parallel.hpp"
