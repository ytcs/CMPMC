# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo
## Architecture & G4CMP Porting Plan

**CMPMC** is a standalone, zero-dependency C++20 Monte Carlo transport engine for low-temperature condensed matter detectors (superconducting qubits, TES/KID dark matter & neutrino calorimeters, and cryogenic semiconductor crystals). It replaces the Geant4 + G4CMP stack downstream of primary energy deposition with a flat, data-oriented ray-tracing and particle-transport kernel.

---

## 1. Why Replace Geant4 in G4CMP?

G4CMP's physics (~15,000 LOC in `G4CMP/library/src/`) is self-contained, but its coupling to Geant4 introduces severe performance and architectural bottlenecks:

1. **The $\vec{k} \nparallel \vec{v}_g$ Impedance Mismatch:**
   - Geant4's `G4Track` assumes momentum $\vec{p}$ and velocity $\vec{v}$ are parallel.
   - In anisotropic crystals, the wavevector $\vec{k}$ and group velocity $\vec{v}_g = \nabla_{\vec{k}}\omega(\vec{k})$ point in different directions (producing phonon caustics).
   - G4CMP works around this by overwriting `G4Track` momentum with $\hat{v}_g$ and storing $\vec{k}$ in a heap-allocated `G4CMPPhononTrackInfo` sidecar queried via `dynamic_cast` on every step.
   - **CMPMC Fix:** `PhononTrack` natively stores both `k_dir` and `vg_dir` in a single 64-byte cache-aligned struct.

2. **Boundary "Turnaround Step" Overhead:**
   - When a phonon reflects off a crystal boundary in Geant4, `G4Navigator` has already associated the boundary state with entering the next volume. G4CMP must execute zero-length turnaround steps (`stepLength <= 1e-13 m` in `G4CMPBoundaryUtils.cc`) and re-query the volume hierarchy.
   - **CMPMC Fix:** Analytic ray-primitive intersection returns `(d_geom, normal, surface_id)`. Reflections update `(k_dir, vg_dir)` in-place without leaving the volume.

3. **Unoptimized Rejection Sampling Loops in Inner Loops:**
   - **Anharmonic Decay (`G4CMPPhononDownconversionProcess.cc`):** `FV4` ($L \to L'+T$, envelope `BUFF=26`) and `FTV5` ($L \to T+T$, envelope `BUFF=8.5`) run 10–30 rejection iterations with `sin`/`cos`/`acos` per decay.
   - **Kaplan QP Cascade (`G4CMPVKaplanQP.cc`):** `QPEnergyRand` and `PhononEnergyRand` have $1/\sqrt{x^2 - \Delta^2}$ endpoint singularities regulated by `BUFF=1000`, requiring ~15–30 rejection iterations per Cooper pair broken.
   - **CMPMC Fix:**
     - Precomputed 1024-bin monotonic Inverse-CDF lookup tables for `FV4` and `FTV5` (O(1) branchless sampling).
     - Analytic trigonometric substitution $x = \Delta + (E - 2\Delta)\sin^2(\frac{\pi}{2} u)$ in `KaplanFilm::sample_qp_energy()`, which cancels the inverse-square-root singularities at both endpoints ($x = \Delta$ and $x = E - \Delta$) and achieves >95% first-draw acceptance.

4. **Heap Allocations & `G4Navigator` Virtual Dispatch:**
   - Every secondary phonon in Geant4 allocates `G4Track`, `G4DynamicParticle`, `G4CMPPhononTrackInfo`, and `G4TouchableHistory` on the heap.
   - **CMPMC Fix:** Preallocated per-thread `std::vector<PhononTrack>` stack and fixed-capacity `SmallVec` buffers inside `KaplanFilm` result in **zero heap allocations** during transport.

---

## 2. G4CMP to CMPMC Module Mapping

| G4CMP Source Files | CMPMC Header / Module | Port Status | Notes |
| :--- | :--- | :---: | :--- |
| `G4ThreeVector`, `Randomize.hh`, `G4SystemOfUnits.hh` | `include/cmpmc/core.hpp` | **Done (Prototype)** | `Vec3`, `Mat3`, `Xoshiro256ss` PRNG, physical units (`eV`, `meV`, `THz`, `ns`, `mm`). |
| `G4LatticeLogical.cc`, `G4LatticePhysical.cc`, `G4CMPPhononKinematics.cc` | `include/cmpmc/lattice.hpp` | **Done (Prototype)** | Exact 3x3 Christoffel eigensolver, $\vec{v}_g$ tensor contraction, bilinear $(\cos\theta, \phi)$ LUT (`181x361`), inverse-CDF decay LUTs, built-in `Ge` & `Si` configs. |
| `G4CMPPhononScatteringProcess.cc`, `G4CMPPhononDownconversionProcess.cc` | `include/cmpmc/lattice.hpp` & `include/cmpmc/engine.hpp` | **Done (Prototype)** | Combined total MFP ($\Gamma_{\text{tot}} = \Gamma_{\text{iso}} + \Gamma_{\text{down}}$), single `std::log(u)` per step, exact isotropic $\vec{k}$ mode conversion & anharmonic decay kinematics. |
| `G4CMPKaplanQP.cc`, `G4CMPVKaplanQP.cc` | `include/cmpmc/kaplan.hpp` | **Done (Prototype)** | Full Kaplan superconducting thin-film cascade, Cooper-pair breaking, QP phonon radiation, escape probability, subgap/direct absorption, singularity-free sampler. |
| `G4CMPPhononBoundaryProcess.cc`, `G4CMPSurfaceProperty.cc`, `G4Tubs`, `G4Box` | `include/cmpmc/geometry.hpp` | **Done (Prototype)** | Analytic `Cylinder` and `Box` ray intersection, frequency-dependent polynomial `SurfaceProperty` (specular, diffuse, anharmonic), sensor attachment. |
| `G4SteppingManager`, `G4Track`, `G4CMPPhononTrackInfo` | `include/cmpmc/engine.hpp` | **Done (Prototype)** | 64-byte `PhononTrack`, zero-allocation event loop, OpenMP multi-threading. |
| `CrystalMaps/*/config.txt` parser | `src/config_loader.cpp` | **Phase 2 (TODO)** | Parse `CrystalMaps/` files (`Al2O3`, `CaWO4`, `GaAs`, `LiF`, `CaF2`, `SiO2`) including triclinic/hexagonal/tetragonal stiffness tensors. |
| `G4CMPBiLinearInterp.cc`, `G4CMPTriLinearInterp.cc`, `G4CMPMeshElectricField.cc` | `include/cmpmc/field.hpp` | **Phase 3 (TODO)** | Replace `G4TessellatedSolid` tetrahedral walk with flat barycentric adjacency walk + bounding grid for 2D/3D FEM electric fields. |
| `G4CMPVDriftProcess.cc`, `G4CMPLukeScattering.cc`, `G4CMPTimeStepper.cc` | `include/cmpmc/charge.hpp` | **Phase 3 (TODO)** | Herring-Vogt valley transform for $e^-$, effective-mass hole drift, Neganov-Trofimov-Luke Cherenkov phonon emission, intervalley scattering. |
| `G4CMPEnergyPartition.cc`, `G4CMPLindhardNIEL.cc`, `G4CMPFanoBinomial.cc` | `include/cmpmc/partition.hpp` | **Phase 4 (TODO)** | Lindhard nuclear recoil ionization yield, Fano binomial charge sampling, primary downconversion phonon cloud generation. |
| `G4CMPBogoliubovQP.cc`, `G4CMPSCPairBreakingProcess.cc` | `include/cmpmc/superconductor.hpp` | **Phase 5 (TODO)** | Explicit 3D Bogoliubov quasiparticle random-walk diffusion, pair-breaking, recombination, and trapping in 3D superconducting films. |

---

## 3. Roadmap to Complete the Port

### Phase 1: Core Phonon + Kaplan QP Engine (Completed in Prototype)
- [x] Zero-dependency `Vec3`, `Mat3`, and `Xoshiro256ss` PRNG (`include/cmpmc/core.hpp`).
- [x] Exact 3x3 symmetric Christoffel matrix solver ($D_{il} = \rho^{-1} C_{ijlm} \hat{k}_j \hat{k}_m$), mode polarization classifier ($L, FT, ST$), and group velocity tensor contraction (`include/cmpmc/lattice.hpp`).
- [x] Precomputed 2D $(\cos\theta, \phi) \to \vec{v}_g(\hat{k})$ bilinear LUT (`181 x 361` grid) + Miller orientation rotation.
- [x] Precomputed 1024-bin Inverse-CDF tables for Tamura anharmonic decay (`FV4` and `FTV5`).
- [x] Lumped Kaplan superconducting thin-film cascade (`include/cmpmc/kaplan.hpp`) with analytic trigonometric singularity cancellation.
- [x] Analytic `Cylinder` and `Box` CSG primitives + frequency-dependent `SurfaceProperty` (`include/cmpmc/geometry.hpp`).
- [x] Event-based transport engine (`include/cmpmc/engine.hpp`) and benchmark executable (`examples/phonon_benchmark.cpp`).

### Phase 2: Generalized Crystal Symmetries & Patterned 2D Surface Masks
- [ ] Add full $6 \times 6$ Voigt stiffness matrix loader for non-cubic crystals (`Al2O3` sapphire, `CaWO4`, `SiO2` quartz) from `CrystalMaps/*/config.txt`.
- [ ] Add `SurfaceMask2D` on planar crystal faces (`z = ±H/2`) to support complex patterned QET/TES fins, KID resonators, and qubit ground planes via 2D spatial lookup (bounding boxes + 2D polygons with holes) instead of 3D boolean CSG subtraction.

### Phase 3: Charge Carrier Transport & Luke Phonon Emission
- [ ] Port Herring-Vogt valley matrices (`valleyDir`, anisotropic `emass` tensor, relativistic/post-Newtonian corrections from `G4LatticeLogical.cc`).
- [ ] Implement flat tetrahedral mesh walker (`TetMeshField`) replacing `G4CMPBiLinearInterp` / `G4CMPTriLinearInterp` (store `Mat3 inv_edge_matrix` per tetrahedron for 9-multiply barycentric walking).
- [ ] Port `G4CMPLukeScattering` and `G4CMPInterValleyScattering` into `ChargeEngine`, pushing emitted Luke phonons directly onto the `PhononTrack` stack.

### Phase 4: Primary Energy Partition (`Lindhard` + `Fano`)
- [ ] Port `G4CMPEnergyPartition`, `G4CMPLindhardNIEL`, and `G4CMPFanoBinomial` so energy deposits `(pos, t, Edep, Eniel, Z, A)` from any external source (or HDF5/ROOT/CSV hit list) automatically partition into $e^-/h^+$ pairs and primary Debye-scale phonons.

### Phase 5: SIMD / GPU Kernel Acceleration
- [ ] Add Structure-of-Arrays (`PhononBatchSoA`) wrapper around `SimulationEngine::step_track`.
- [ ] Annotate `core.hpp`, `lattice.hpp`, `kaplan.hpp`, and `geometry.hpp` methods with `CMPMC_HOST_DEVICE` (`__host__ __device__`) to compile directly with `nvcc` (CUDA) or `hipcc` (ROCm).
