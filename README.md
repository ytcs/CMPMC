# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo

**CMPMC** is a fast, standalone C++20 condensed-matter Monte Carlo engine that ports and optimizes the physics of [G4CMP](https://github.com/kelseymh/G4CMP) without any Geant4 or CLHEP dependencies. It is architected from the ground up as a **multi-excitation transport engine** supporting **Phonons**, **Conduction-Band Drift Electrons** (with anisotropic Herring-Vogt valleys), **Valence-Band Drift Holes**, **Superconducting Bogoliubov Quasiparticles**, and extensible collective excitations (**Excitons, Magnons, Plasmons**).

See [`PLAN.md`](PLAN.md) for the full architectural breakdown, G4CMP-to-CMPMC module mapping, mathematical optimizations, and porting roadmap.

## Features in the Prototype

- **Zero External Dependencies (`include/cmpmc/core.hpp`):** Header-only C++20 library requiring only the C++ standard library and `std::thread`.
- **Anisotropic Crystal Kinematics & Herring-Vogt Valleys (`include/cmpmc/lattice.hpp`):**
  - Exact $3 \times 3$ symmetric Christoffel matrix solver ($D_{il} = \rho^{-1} C_{ijlm} \hat{k}_j \hat{k}_m$) and group velocity tensor contraction.
  - Pre-baked 2D $(\cos\theta, \phi)$ bilinear lookup table (`181 x 361` bins per polarization mode) + Miller orientation rotation.
  - Pre-baked 1024-bin Inverse-CDF tables for Tamura anharmonic decay (`FV4` $L \to L'+T$ and `FTV5` $L \to T+T$), replacing G4CMP's 30-iteration rejection sampling loops with a single branchless table lookup.
  - Full Herring-Vogt anisotropic electron effective-mass tensors, valley rotation matrices ($[111]$ L-valleys in Ge, $[100]$ X-valleys in Si), and acoustic deformation scattering lengths ($l_{0,e}, l_{0,h}$).
- **Charge Carrier Drift, Intervalley Scattering & Luke Phonon Emission (`include/cmpmc/charge.hpp`):**
  - Full port of `G4CMPVDriftProcess`, `G4CMPLukeScattering`, `G4CMPLukeEmissionRate`, and `G4CMPInterValleyScattering` (Quadratic & Linear rate models).
  - Analytic inverse-CDF Cherenkov angle sampling $\cos\theta_{\text{ph}} = (k_s/k) + (1 - k_s/k)(1 - u)^{1/3}$ that eliminates G4CMP's 1,000-trial rejection loop.
- **Primary Energy Partition (`include/cmpmc/partition.hpp`):**
  - Port of `G4CMPEnergyPartition`, `G4CMPLindhardNIEL` (Lewin & Smith nuclear recoil ionization yield), and `G4CMPFanoBinomial` (Fano-fluctuated $e^-/h^+$ pair generation + prompt Debye phonon cloud with statistical weighting).
- **Lumped Kaplan Superconducting Thin-Film Cascade (`include/cmpmc/kaplan.hpp`):**
  - Full port of `G4CMPKaplanQP` with analytic trigonometric variable substitutions that identically cancel the $1/\sqrt{x^2 - \Delta^2}$ endpoint singularities, using stack-allocated `SmallVec` buffers (zero heap allocations).
- **Analytic CSG Ray-Tracing & Multi-Excitation Engine (`include/cmpmc/geometry.hpp`, `include/cmpmc/engine.hpp`):**
  - 64-byte cache-aligned `PhononTrack`, `ChargeTrack`, and `GenericExcitationTrack` structs.
  - Analytic `CylinderSolid` and `BoxSolid` exit intersection with zero Geant4 boundary turnaround steps.

## Quick Start

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/phonon_benchmark 5000
```
