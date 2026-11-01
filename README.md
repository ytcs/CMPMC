# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo

**CMPMC** is a fast, standalone C++20 condensed-matter Monte Carlo engine that ports and optimizes the core physics of [G4CMP](https://github.com/kelseymh/G4CMP) without any Geant4 or CLHEP dependencies.

See [`PLAN.md`](PLAN.md) for the full architectural breakdown, G4CMP-to-CMPMC module mapping, mathematical optimizations, and porting roadmap.

## Features in the Prototype

- **Zero External Dependencies:** Header-only C++20 library (`include/cmpmc/`) requiring only the C++ standard library and `std::thread`.
- **Anisotropic Crystal Kinematics (`include/cmpmc/lattice.hpp`):**
  - Exact $3 \times 3$ symmetric Christoffel matrix solver ($D_{il} = \rho^{-1} C_{ijlm} \hat{k}_j \hat{k}_m$) and group velocity tensor contraction.
  - Pre-baked 2D $(\cos\theta, \phi)$ bilinear lookup table (`181 x 361` bins per polarization mode) + Miller orientation rotation (~65M queries/sec).
  - Pre-baked 1024-bin Inverse-CDF tables for Tamura anharmonic decay (`FV4` $L \to L'+T$ and `FTV5` $L \to T+T$), replacing G4CMP's 30-iteration rejection sampling loops with a single branchless table lookup.
- **Lumped Kaplan Superconducting Thin-Film Cascade (`include/cmpmc/kaplan.hpp`):**
  - Full port of `G4CMPKaplanQP` (Cooper-pair breaking, quasiparticle phonon emission, escape probability, direct/subgap absorption).
  - Analytic trigonometric change-of-variables $x = \Delta + (E - 2\Delta)\sin^2(\frac{\pi}{2} u)$ and $x = \Delta + (E - \Delta)u^2$ that identically cancel the $1/\sqrt{x^2 - \Delta^2}$ endpoint singularities.
  - Zero heap allocations via fixed-capacity stack buffers (`SmallVec`).
- **Analytic CSG Ray-Tracing & Event-Based Transport (`include/cmpmc/geometry.hpp`, `include/cmpmc/engine.hpp`):**
  - 64-byte cache-aligned `PhononTrack` natively storing both wavevector `k_dir` and group velocity `vg_dir`.
  - Analytic `CylinderSolid` and `BoxSolid` exit intersection with zero Geant4 boundary turnaround steps.

## Quick Start

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/phonon_benchmark 10000
```
