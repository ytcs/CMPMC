# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo
## Multi-Excitation Architecture & G4CMP Porting Plan

**CMPMC** is a standalone, zero-dependency C++20 Monte Carlo transport engine for low-temperature condensed matter detectors (superconducting qubits, TES/KID dark matter & neutrino calorimeters, and cryogenic semiconductor crystals). It replaces the Geant4 + G4CMP stack downstream of primary energy deposition with a flat, data-oriented multi-excitation transport kernel supporting **Phonons, Conduction Electrons, Valence Holes, Bogoliubov Quasiparticles**, and future collective excitations (**Excitons, Magnons, Plasmons**).

---

## 1. Why Replace Geant4 in G4CMP?

G4CMP's physics (~15,000 LOC in `G4CMP/library/src/`) is self-contained, but its coupling to Geant4 introduces severe performance and architectural bottlenecks:

1. **The $\vec{k} \nparallel \vec{v}_g$ & Anisotropic Valley Impedance Mismatch:**
   - Geant4's `G4Track` assumes momentum $\vec{p}$ and velocity $\vec{v}$ are parallel.
   - In anisotropic crystals, phonon wavevector $\vec{k}$ and group velocity $\vec{v}_g = \nabla_{\vec{k}}\omega(\vec{k})$ point in different directions (producing phonon caustics), and conduction-band electrons have ellipsoidal Herring-Vogt effective mass tensors $M = \text{diag}(m_\parallel, m_\perp, m_\perp)$ in tilted valley frames.
   - G4CMP works around this by overwriting `G4Track` momentum and storing true wavevectors/valley indices in heap-allocated `G4CMPPhononTrackInfo` and `G4CMPDriftTrackInfo` sidecars queried via `dynamic_cast` on every step.
   - **CMPMC Fix:** Flat 64-byte cache-aligned structs (`PhononTrack`, `ChargeTrack`, `GenericExcitationTrack`) natively store wavevectors, group velocities, transport momenta, and valley indices.

2. **Boundary "Turnaround Step" Overhead:**
   - When a particle reflects off a crystal boundary in Geant4, `G4Navigator` has already associated the boundary state with entering the next volume. G4CMP must execute zero-length turnaround steps (`stepLength <= 1e-13 m` in `G4CMPBoundaryUtils.cc`) and re-query the volume hierarchy.
   - **CMPMC Fix:** Analytic ray-primitive intersection returns `(d_geom, inward_normal, surface_id)`. Reflections update `(k_dir, vg_dir)` in-place without leaving the volume.

3. **Unoptimized Rejection Sampling Loops in Inner Loops:**
   - **Anharmonic Decay (`G4CMPPhononDownconversionProcess.cc`):** `FV4` ($L \to L'+T$, envelope `BUFF=26`) and `FTV5` ($L \to T+T$, envelope `BUFF=8.5`) run 10–30 rejection iterations with `sin`/`cos`/`acos` per decay.
   - **Kaplan QP Cascade (`G4CMPVKaplanQP.cc`):** `QPEnergyRand` and `PhononEnergyRand` have $1/\sqrt{x^2 - \Delta^2}$ endpoint singularities regulated by `BUFF=1000`, requiring ~15–30 rejection iterations per Cooper pair broken.
   - **Neganov-Trofimov-Luke Phonon Emission (`G4CMPLukeScattering.cc`):** Wraps Cherenkov angle sampling in a 1,000-trial `while (!goodThrow)` loop.
   - **CMPMC Fix:**
     - Precomputed 1024-bin monotonic Inverse-CDF lookup tables for `FV4` and `FTV5` (O(1) branchless sampling).
     - Analytic trigonometric substitution $x = \Delta + (E - 2\Delta)\sin^2(\frac{\pi}{2} u)$ in `KaplanFilm::sample_qp_energy()` and $x = \Delta + (E - \Delta)u^2$ in `KaplanFilm::sample_phonon_energy()`, which cancel the inverse-square-root endpoint singularities identically.
     - Analytic inverse-CDF $\cos\theta_{\text{ph}} = (k_s/k) + (1 - k_s/k)(1 - u)^{1/3}$ in `ChargePhysics::emit_luke_phonon()`, which is strictly inside the forward Cherenkov cone on the first draw.

---

## 2. G4CMP to CMPMC Module Mapping

| G4CMP Source Files | CMPMC Header / Module | Port Status | Notes |
| :--- | :--- | :---: | :--- |
| `G4ThreeVector`, `Randomize.hh`, `G4SystemOfUnits.hh` | `include/cmpmc/core.hpp` | **Done (Prototype)** | `Vec3`, `Mat3`, `Xoshiro256ss` PRNG, physical units, `ExcitationType` enum (`Phonon`, `Electron`, `Hole`, `BogoliubovQP`, `Exciton`, `Magnon`, `Plasmon`). |
| `G4LatticeLogical.cc`, `G4LatticePhysical.cc`, `G4CMPPhononKinematics.cc` | `include/cmpmc/lattice.hpp` | **Done (Prototype)** | Exact 3x3 Christoffel eigensolver, $\vec{v}_g$ tensor contraction, bilinear $(\cos\theta, \phi)$ LUT (`181x361`), inverse-CDF decay LUTs, Herring-Vogt valley transforms, built-in `Ge` & `Si` configs. |
| `G4CMPPhononScatteringProcess.cc`, `G4CMPPhononDownconversionProcess.cc` | `include/cmpmc/lattice.hpp` & `include/cmpmc/engine.hpp` | **Done (Prototype)** | Combined total MFP ($\Gamma_{\text{tot}} = \Gamma_{\text{iso}} + \Gamma_{\text{down}}$), single `std::log(u)` per step, exact isotropic $\vec{k}$ mode conversion & anharmonic decay kinematics. |
| `G4CMPKaplanQP.cc`, `G4CMPVKaplanQP.cc` | `include/cmpmc/kaplan.hpp` | **Done (Prototype)** | Full Kaplan superconducting thin-film cascade, Cooper-pair breaking, QP phonon radiation, escape probability, singularity-free sampler. |
| `G4CMPVDriftProcess.cc`, `G4CMPLukeScattering.cc`, `G4CMPInterValleyScattering.cc`, `G4CMPTimeStepper.cc` | `include/cmpmc/charge.hpp` & `include/cmpmc/engine.hpp` | **Done (Prototype)** | Herring-Vogt valley transform for $e^-$, effective-mass hole drift, Neganov-Trofimov-Luke Cherenkov phonon emission, Quadratic & Linear intervalley scattering. |
| `G4CMPEnergyPartition.cc`, `G4CMPLindhardNIEL.cc`, `G4CMPFanoBinomial.cc` | `include/cmpmc/partition.hpp` | **Done (Prototype)** | Lindhard nuclear recoil ionization yield, Fano charge-pair sampling, prompt Debye phonon cloud generation, statistical weighting. |
| `G4CMPPhononBoundaryProcess.cc`, `G4CMPSurfaceProperty.cc`, `G4Tubs`, `G4Box` | `include/cmpmc/geometry.hpp` | **Done (Prototype)** | Analytic `CylinderSolid` and `BoxSolid` ray intersection, frequency-dependent polynomial `SurfaceProperty` (specular, diffuse, anharmonic), sensor attachment. |
| `CrystalMaps/*/config.txt` parser | `src/config_loader.cpp` | **Phase 2 (TODO)** | Parse `CrystalMaps/` files (`Al2O3`, `CaWO4`, `GaAs`, `LiF`, `CaF2`, `SiO2`) including triclinic/hexagonal/tetragonal stiffness tensors. |
| `G4CMPBiLinearInterp.cc`, `G4CMPTriLinearInterp.cc`, `G4CMPMeshElectricField.cc` | `include/cmpmc/field.hpp` | **Phase 3 (TODO)** | Replace `G4TessellatedSolid` tetrahedral walk with flat barycentric adjacency walk + bounding grid for 2D/3D FEM electric fields. |
| `G4CMPBogoliubovQP.cc`, `G4CMPSCPairBreakingProcess.cc` | `include/cmpmc/superconductor.hpp` | **Phase 4 (TODO)** | Explicit 3D Bogoliubov quasiparticle random-walk diffusion, pair-breaking, recombination, and trapping in 3D superconducting films. |
| Future Excitations (`Exciton`, `Magnon`, `Plasmon`) | `include/cmpmc/excitations.hpp` | **Phase 5 (TODO)** | Singlet/triplet exciton diffusion & radiative decay (scintillators / LHe), spin-wave magnon dispersion & Gilbert damping, GPU/SIMD batching. |

---

## 3. Roadmap to Complete the Port

### Phase 1: Core Multi-Excitation Prototype (Completed)
- [x] Zero-dependency `Vec3`, `Mat3`, `ExcitationType`, and `Xoshiro256ss` PRNG (`include/cmpmc/core.hpp`).
- [x] Exact 3x3 symmetric Christoffel matrix solver, mode polarization classifier ($L, FT, ST$), group velocity LUT, Inverse-CDF anharmonic decay tables, and Herring-Vogt valley matrices (`include/cmpmc/lattice.hpp`).
- [x] Lumped Kaplan superconducting thin-film cascade (`include/cmpmc/kaplan.hpp`) with analytic trigonometric singularity cancellation.
- [x] Conduction-band electron & valence-band hole drift, Neganov-Trofimov-Luke Cherenkov phonon emission, and intervalley scattering (`include/cmpmc/charge.hpp`).
- [x] Lindhard nuclear recoil ionization partition + Fano charge-pair & prompt phonon generator (`include/cmpmc/partition.hpp`).
- [x] Analytic `CylinderSolid` and `BoxSolid` CSG primitives + frequency-dependent `SurfaceProperty` (`include/cmpmc/geometry.hpp`).
- [x] Unified multi-excitation transport engine (`include/cmpmc/engine.hpp`) and benchmark suite (`examples/phonon_benchmark.cpp`).

### Phase 2: Generalized Crystal Symmetries & Patterned 2D Surface Masks
- [ ] Add full $6 \times 6$ Voigt stiffness matrix loader for non-cubic crystals (`Al2O3` sapphire, `CaWO4`, `SiO2` quartz, `GaAs`, `LiF`) from `CrystalMaps/*/config.txt`.
- [ ] Add `SurfaceMask2D` on planar crystal faces (`z = ±H/2`) to support complex patterned QET/TES fins, KID resonators, and qubit ground planes via 2D spatial lookup.

### Phase 3: 3D FEM Mesh Electric Fields & Charge Trapping
- [ ] Implement flat tetrahedral mesh walker (`TetMeshField`) replacing `G4CMPBiLinearInterp` / `G4CMPTriLinearInterp` (store `Mat3 inv_edge_matrix` per tetrahedron for 9-multiply barycentric walking).
- [ ] Port `G4CMPChargeTrapping` and `G4CMPTrapIonization` (shallow/deep impurity trapping and impact ionization).

### Phase 4: 3D Superconducting Bogoliubov Quasiparticles & New Excitations
- [ ] Port `G4CMPBogoliubovQP`, `G4CMPSCPairBreakingProcess`, `G4CMPSCRecombinationProcess`, and `G4CMendozaQP` into a dedicated `BogoliubovEngine` using `GenericExcitationTrack`.
- [ ] Add `Exciton` (singlet/triplet diffusion + VUV photon/phonon emission) and `Magnon` (anisotropic spin-wave dispersion $\omega(\vec{k})$) transport handlers.

### Phase 5: SIMD / GPU Kernel Acceleration
- [ ] Add Structure-of-Arrays (`ExcitationBatchSoA`) wrapper and `__host__ __device__` annotations for CUDA/ROCm/SYCL compilation.
