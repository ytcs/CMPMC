# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo

**CMPMC** is a fast, standalone C++20 condensed-matter Monte Carlo engine that ports and optimizes the physics of [G4CMP](https://github.com/kelseymh/G4CMP) without any Geant4 or CLHEP dependencies. It is architected from the ground up as a **multi-excitation transport engine** supporting **Phonons**, **Conduction-Band Drift Electrons** (with anisotropic Herring-Vogt valleys), **Valence-Band Drift Holes**, **Superconducting Bogoliubov Quasiparticles**, and extensible collective excitations (**Excitons, Magnons, Plasmons**).

See:
- **[`PLAN.md`](PLAN.md)**: Architectural breakdown, G4CMP-to-CMPMC module mapping, and the **10,000x Speedup Roadmap** (Octahedral LUTs, SIMD SoA wavefronts, GPU hardware texture/RT cores, and algorithmic physics shortcuts).
- **[`CONTRIBUTING_PHYSICS.md`](CONTRIBUTING_PHYSICS.md)**: Physicist's cookbook for adding new crystals/films, modifying scattering/decay formulas, or plugging in new excitation species in ~30 lines of C++20.

## Complete Feature Set (100% G4CMP Physics Port + Optimizations)

- **Zero External Dependencies (`include/cmpmc/core.hpp`):** Header-only C++20 library requiring only the C++ standard library and `std::thread`.
- **All 8 Bravais Symmetries, 12 `CrystalMaps` Presets & `config.txt` Reader (`include/cmpmc/materials/`):**
  - Supports `Amorphous`, `Cubic`, `Tetragonal`, `Orthorhombic`, `Hexagonal`, `Rhombohedral`, `Monoclinic`, and `Triclinic` crystals (`CrystalConfig::build_voigt_matrix()`) + Miller $(h, k, l)$ orientation.
  - Built-in presets (`presets.hpp`) for `Ge`, `Si`, `Al2O3` (Sapphire), `CaWO4`, `SiO2_alpha` (Quartz), `SiO2_amorph`, `GaAs`, `CaF2`, `LiF`, `Al`, `Nb`, and `Cu`, plus a standalone `G4LatticeReader` parser (`load_g4cmp_config_file()`).
- **Anisotropic Phonon Kinematics & Scattering (`include/cmpmc/lattice.hpp`, `include/cmpmc/physics/phonon_physics.hpp`):**
  - Exact $3 \times 3$ symmetric Christoffel matrix solver ($D_{il} = \rho^{-1} C_{ijlm} \hat{k}_j \hat{k}_m$) and group velocity tensor contraction.
  - Pre-baked 2D $(\cos\theta, \phi)$ bilinear lookup table (`181 x 361` bins) and zero-trig FP32 `256 x 256` Octahedral lookup table (`map_k_to_vg_oct`).
  - Pre-baked 1024-bin Inverse-CDF tables for Tamura anharmonic decay (`FV4` $L \to L'+T$ and `FTV5` $L \to T+T$), isotopic $\omega^4$ scattering, and polycrystalline grain-boundary elastic scattering.
- **Complete Charge Carrier Transport (`include/cmpmc/charge.hpp`, `include/cmpmc/physics/charge_physics.hpp`):**
  - Herring-Vogt anisotropic electron effective-mass valleys, effective-mass hole drift, and Neganov-Trofimov-Luke Cherenkov phonon emission (analytic inverse-CDF cone sampler).
  - First-principles deformation-potential matrix-element optical phonon + Erginsoy neutral impurity intervalley scattering (`G4CMPInterValleyRate`), bulk charge trapping (`G4CMPDriftTrappingProcess`), trap impact ionization (`G4CMPDriftTrapIonization`), boundary reflection/recombination (`G4CMPDriftBoundaryProcess`), and 3D initial charge clouds (`G4CMPChargeCloud`).
- **All 5 NIEL Partition Models & Exact Fano Binomial Sampler (`include/cmpmc/partition.hpp`, `include/cmpmc/physics/niel_models.hpp`):**
  - Complete ports of `Lindhard`, `LewinSmith`, `Sarkis` (64-knot low-energy Si table), `ImpactTUNL`, and CDMSlite `Empirical` $k(E)$ models, plus the exact interpolated Fano binomial sampler (`G4CMPFanoBinomial`).
- **Lumped Kaplan Film & Full 3D Superconducting Bogoliubov Quasiparticles (`include/cmpmc/kaplan.hpp`, `include/cmpmc/physics/superconductor_3d.hpp`):**
  - Singularity-free `KaplanFilm` cascade + $O(1)$ branchless `KaplanSurrogateLUT`.
  - Full 3D Bogoliubov quasiparticle Brownian diffusion ($D(E) = D_0 \sqrt{1 - (\Delta/E)^2}$), Cooper-pair breaking, QP phonon radiation, thermal recombination, and local vortex trapping (`transport_superconductor_coupled`).
- **2D/3D Tetrahedral FEM Mesh Electric Fields & Patterned 2D Sensor Masks (`include/cmpmc/geometry/`):**
  - `TriLinearMesh3D` and `BiLinearMesh2D` barycentric facet-walking FEM field interpolators (`G4CMPTriLinearInterp`, `G4CMPBiLinearInterp`, `G4CMPMeshElectricField`).
  - `ElectrodePatternMask` (`G4CMPVElectrodePattern`) for multi-channel QET/KID/qubit surface pad layouts.
- **Extensible Collective Excitations (`include/cmpmc/physics/custom_excitation_example.hpp`):**
  - C++20 `ExcitationProcessConcept` (`transport_custom_excitation`) with ready-to-use `SingletExcitonProcess` and `FerromagneticMagnonProcess` models.

## Quick Start

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/phonon_benchmark 5000
```
