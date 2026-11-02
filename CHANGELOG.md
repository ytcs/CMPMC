# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions follow
[Semantic Versioning](https://semver.org/); before 1.0, minor versions may contain
breaking changes.

## [Unreleased]

### Changed
- Headers are reorganised into module directories:

  | Previous path | New path |
  |---|---|
  | `cmpmc/core.hpp` | `cmpmc/core/core.hpp` |
  | `cmpmc/lattice.hpp` | `cmpmc/lattice/crystal_lattice.hpp` |
  | `cmpmc/materials/config_reader.hpp` | `cmpmc/io/g4cmp_config_import.hpp` |
  | `cmpmc/physics/phonon_physics.hpp` | `cmpmc/phonon/phonon_kernels.hpp` |
  | `cmpmc/physics/interface_transmission.hpp` | `cmpmc/phonon/interface_transmission.hpp` |
  | `cmpmc/charge.hpp` | `cmpmc/charge/charge_transport.hpp` |
  | `cmpmc/physics/charge_physics.hpp` | `cmpmc/charge/charge_kernels.hpp` |
  | `cmpmc/partition.hpp` | `cmpmc/charge/energy_partition.hpp` |
  | `cmpmc/physics/niel_models.hpp` | `cmpmc/charge/niel_models.hpp` |
  | `cmpmc/kaplan.hpp` | `cmpmc/superconductor/kaplan_film.hpp` |
  | `cmpmc/physics/superconductor_3d.hpp` | `cmpmc/superconductor/bulk_superconductor.hpp` |
  | `cmpmc/geometry.hpp` | `cmpmc/geometry/solids.hpp` + `cmpmc/geometry/surface_property.hpp` |
  | `cmpmc/geometry/surface_mask.hpp` | `cmpmc/geometry/electrode_pattern.hpp` |
  | `cmpmc/engine.hpp` | `cmpmc/transport/engine.hpp` |
  | `cmpmc/physics/custom_excitation_example.hpp` | `cmpmc/transport/excitation_model.hpp` (concept only) |
  | `cmpmc/physics/fast_shortcuts.hpp` | `cmpmc/transport/variance_reduction.hpp`, `cmpmc/experimental/jump_to_ballistic.hpp`, `cmpmc/experimental/next_event_estimation.hpp` |
  | `cmpmc/wavefront_soa.hpp` | `cmpmc/experimental/wavefront_soa.hpp` |
  | `cmpmc/gpu/wavefront_pipeline.hpp` | `cmpmc/experimental/wavefront_pipeline.hpp` |

- Renamed identifiers:

  | Previous | New |
  |---|---|
  | `materials::load_g4cmp_config_file` | `io::import_g4cmp_config` |
  | `FilmLossModel::G4CMPMeanFreePath` | `FilmLossModel::MeanFreePath` (same value, still the default) |
  | `SurfaceProperty::g4cmp_demo_al_sensor_surface` | `SurfaceProperty::demo_polynomial_sensor_surface` |
  | `SurfaceProperty::g4cmp_demo_bare_wall_surface` | `SurfaceProperty::demo_polynomial_bare_surface` |
  | `physics::Superconductor3D` | `physics::BulkSuperconductor` |
  | `physics::ExcitationProcessConcept` | `ExcitationModel` |
  | `physics::apply_russian_roulette` | `apply_russian_roulette` |
  | `physics::JumpToBallisticConfig`, `physics::jump_to_ballistic_cascade` | `experimental::…` |
  | `physics::evaluate_sparse_sensor_nee` | `experimental::evaluate_sparse_sensor_nee` |
  | `WavefrontSoAEngine`, `PhononBatchSoA16` | `experimental::…` |
  | `gpu::WavefrontPipelineEngine` and related types | `experimental::…` |
  | `ElectrodePatternMask::supercdms_qet_mask` | `ElectrodePatternMask::five_channel_qet_mask` |

- The example models `BogoliubovQPModel` and `MagnonModel` moved out of the library
  into `examples/custom_excitation.cpp`.
- CMake: the project is now named `cmpmc` and exports the `cmpmc::cmpmc` target. It
  provides install rules and a package config (`find_package(cmpmc)`), plus the options
  `CMPMC_BUILD_TESTS`, `CMPMC_BUILD_EXAMPLES`, `CMPMC_INSTALL`,
  `CMPMC_WARNINGS_AS_ERRORS` and `CMPMC_NATIVE_ARCH`. `-march=native` is no longer the
  default.
- The `UnvalidatedModuleError` message now points to `docs/validation.md`.

### Added
- `cmpmc/cmpmc.hpp` umbrella header and `cmpmc/version.hpp`.
- `examples/quickstart.cpp` and `examples/custom_excitation.cpp`.
- Documentation set under `docs/`, plus `CONTRIBUTING.md` and this changelog.
- CTest labels, a skip return code, and `test_umbrella_header`.
- `test_config_import` (importer on a synthetic file) and a build-time check that every
  public header compiles on its own.
- `test_g4cmp_import`, which checks the presets against a G4CMP `CrystalMaps` tree when
  `CMPMC_G4CMP_CRYSTALMAPS_DIR` is set and reports SKIPPED otherwise.

### Fixed
- `materials/crystal_config.hpp` did not compile on its own because it was missing
  `#include <algorithm>`.

## [0.1.0]

First version.
