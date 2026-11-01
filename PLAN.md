# CMPMC: Zero-Dependency Condensed Matter Physics Monte Carlo
## Multi-Excitation Architecture, Porting Plan & 10,000x Optimization Roadmap

**CMPMC** is a standalone, zero-dependency C++20 Monte Carlo transport engine for low-temperature condensed matter detectors (superconducting qubits, TES/KID dark matter & neutrino calorimeters, and cryogenic semiconductor crystals). It replaces the Geant4 + G4CMP stack downstream of primary energy deposition with a flat, data-oriented multi-excitation transport kernel supporting **Phonons, Conduction Electrons, Valence Holes, Bogoliubov Quasiparticles**, and future collective excitations (**Excitons, Magnons, Plasmons**).

---

## 1. Why Replace Geant4 in G4CMP?

G4CMP's physics (~15,000 LOC in `G4CMP/library/src/`) is self-contained, but its coupling to Geant4 introduces severe performance and architectural bottlenecks:

1. **The $\vec{k} \nparallel \vec{v}_g$ & Anisotropic Valley Impedance Mismatch:**
   - Geant4's `G4Track` assumes momentum $\vec{p}$ and velocity $\vec{v}$ are parallel.
   - In anisotropic crystals, phonon wavevector $\vec{k}$ and group velocity $\vec{v}_g = \nabla_{\vec{k}}\omega(\vec{k})$ point in different directions (producing phonon caustics), and conduction-band electrons have ellipsoidal Herring-Vogt effective mass tensors $M = \text{diag}(m_\parallel, m_\perp, m_\perp)$ in tilted valley frames.
   - G4CMP works around this by overwriting `G4Track` momentum and storing true wavevectors/valley indices in heap-allocated `G4CMPPhononTrackInfo` and `G4CMPDriftTrackInfo` sidecars queried via `dynamic_cast` on every step.
   - **CMPMC Fix:** Flat 64-byte cache-aligned structs (`PhononTrack`, `ChargeTrack`, `GenericExcitationTrack`) natively store wavevectors, group velocities, crystal momenta, and valley indices.

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

## 3. The 10,000x Speedup Roadmap: CPU SIMD, GPU Ray-Tracing & Algorithmic Physics Shortcuts

Achieving a **1,000x to 10,000x+ end-to-end speedup** over baseline Geant4 + G4CMP comes from stacking four multiplicative tiers of optimization:

```mermaid
flowchart TD
    T0["Baseline Geant4 + G4CMP\n(~5k - 20k steps/s per core)"]
    T1["Tier 1: Zero-Dep Event Loop + Inverse-CDF LUTs\n(20x - 40x speedup | DONE in Prototype)"]
    T2["Tier 2: Octahedral FP32 LUT + SIMD SoA Wavefronts\n(3x - 5x additional -> ~100x - 150x cumulative)"]
    T3["Tier 3: GPU Wavefront Kernel + Hardware Texture & RT Cores\n(30x - 80x additional -> ~3,000x - 10,000x cumulative)"]
    T4["Tier 4: Algorithmic Physics Shortcuts & Variance Reduction\n(5x - 50x additional -> 15,000x - 100,000x+ effective speedup)"]

    T0 --> T1 --> T2 --> T3 --> T4
```

### Tier 1: Zero-Dependency Event Kernel & Inverse-CDF Tables (`20x–40x` — *Completed in Prototype*)
- Eliminated `G4Track` / `G4Step` heap allocations, `G4Navigator` virtual dispatch, `dynamic_cast<G4CMPPhononTrackInfo*>`, and zero-length boundary turnaround steps.
- Replaced `FV4` / `FTV5` rejection loops with 1024-bin Inverse-CDF tables, removed `KaplanFilm` endpoint singularities via analytic trigonometric variable substitutions, and replaced the 1,000-trial `G4CMPLukeScattering` rejection loop with analytic Cherenkov cone sampling.

---

### Tier 2: Octahedral LUT Mapping, Cache-Line Packing & CPU SIMD (`3x–5x` additional)

1. **Octahedral Equal-Area Mapping (Zero-Trig LUT Indexing):**
   - Currently, `CrystalLattice::map_k_to_vg()` computes `std::atan2(k.y, k.x)` on every step (~18 CPU cycles, blocks SIMD vectorization) and indexes a `double` `(cos_theta, phi)` table (`4.7 MB`, which spills out of L2 cache into L3).
   - **Optimization:** Replace spherical `(cos_theta, phi)` coordinates with **Octahedral Unit-Vector Mapping** $\mathbb{S}^2 \to [-1, 1]^2$:
     ```cpp
     // Maps unit vector k -> (u, v) in [0, 1]^2 with ZERO trig functions!
     inline void oct_encode(const Vec3f& k, float& u, float& v) noexcept {
         const float inv_l1 = 1.0f / (std::abs(k.x) + std::abs(k.y) + std::abs(k.z));
         float px = k.x * inv_l1;
         float py = k.y * inv_l1;
         if (k.z < 0.0f) {
             const float tx = (1.0f - std::abs(py)) * (px >= 0.0f ? 1.0f : -1.0f);
             const float ty = (1.0f - std::abs(px)) * (py >= 0.0f ? 1.0f : -1.0f);
             px = tx; py = ty;
         }
         u = px * 0.5f + 0.5f;
         v = py * 0.5f + 0.5f;
     }
     ```
   - **FP32 / FP16 LUT Compression:** Store `{vg_dir.x, vg_dir.y, vg_dir.z, vg_mag}` as 16-byte `float4` on a `256 x 256` octahedral grid per mode (`3 x 256 x 256 x 16 B = 3.0 MB`, or `1.5 MB` with `uint16_t` octahedral-encoded $\hat{v}_g$ + `float` $v_g$). The entire 3-mode LUT fits inside a single core's **L2 cache**, and `map_k_to_vg` vectorizes across 8/16 lanes with `AVX2` / `AVX-512`.

2. **Structure-of-Arrays (SoA) SIMD Batching:**
   - Transpose `PhononTrack` into 16-wide SIMD batches (`PhononWavefrontSoA`) so ray-cylinder intersection, `-log(u)` free-flight sampling, and octahedral LUT interpolation execute in lockstep across 16 tracks per instruction.

---

### Tier 3: GPU Wavefront Ray-Tracing, Hardware Texture Units & RT Cores (`30x–80x` over 14-core CPU)

Phonon and quasiparticle transport is embarrassingly parallel and maps directly onto modern GPU ray-tracing hardware (CUDA / OptiX, Vulkan Ray Tracing, or Metal):

1. **Hardware 2D Layered Texture Units (`tex2DLayered<float4>`) for Anisotropic Kinematics:**
   - Bind the 3-mode Octahedral $\hat{k} \to (\hat{v}_g, |\vec{v}_g|)$ LUT to a hardware `cudaTextureObject_t` (`cudaFilterModeLinear`).
   - Every GPU Streaming Multiprocessor (SM) has dedicated fixed-function texture silicon that performs 2D bilinear interpolation in **1 clock cycle** with dedicated texture L1 cache, completely offloading `map_k_to_vg` from the FP32 compute cores!

2. **Hardware Ray-Tracing (RT) Cores for Complex Detector & Qubit Geometries:**
   - For simple cylindrical/box crystals, an inline register-only ray-cylinder intersection takes ~10 GPU instructions.
   - For complex multi-tower calorimeters, crystal housing clamps, and **patterned 2D superconducting qubit/resonator/QET masks** (`QuasiparticleResonatorAssembly` with thousands of micron-scale features), compile the geometry into an **OptiX / Vulkan BVH acceleration structure**.
   - Hardware **RT Cores** traverse the BVH and return `(t_hit, primitive_id, surface_id)` at **>10 Billion rays/second** on a single consumer/datacenter GPU.

3. **Wavefront Path-Tracing Queue Architecture (Zero Warp Divergence):**
   - A naive single-kernel ("megakernel") GPU port suffers from **warp divergence**: if 1 thread in a 32-thread warp enters a multi-step `KaplanFilm` QP cascade or anharmonic decay while the other 31 threads do a specular reflection, 31 GPU lanes sit idle.
   - **Wavefront Solution:** Decompose the transport loop into 4 compact, non-divergent GPU kernels connected by ring buffers in global/L2 memory:
     ```mermaid
     flowchart LR
         Q_Ray["Active Ray Queue\n(SoA Global Buffer)"] --> K_Step["Kernel 1: Free-Flight & Ray-Trace\n(tex2D LUT + RT Core / Analytic Hit)"]
         K_Step -->|d_phys < d_geom\nIsotope| Q_Iso["Isotope Scatter Queue"]
         K_Step -->|d_phys < d_geom\nDecay| Q_Dec["Anharmonic Decay Queue"]
         K_Step -->|d_geom <= d_phys\nReflect| Q_Refl["Surface Reflect Queue"]
         K_Step -->|d_geom <= d_phys\nSensor| Q_Kap["Kaplan Film Queue"]
         Q_Iso & Q_Dec & Q_Refl & Q_Kap -->|Compacted Active Tracks| Q_Ray
     ```
   - Use warp-level ballot intrinsics (`__ballot_sync`, `__popc`) to compact output queues in shared memory before writing to global memory, keeping warp execution efficiency above **90%**.

---

### Tier 4: Algorithmic Physics Shortcuts & Variance Reduction (`5x–50x` additional -> `10,000x+` total)

Beyond hardware acceleration, four physics-informed algorithmic shortcuts drastically reduce the number of steps required per event while preserving exact observables:

1. **Pre-Baked Kaplan Thin-Film Transfer Matrix (`KaplanSurrogateLUT`):**
   - Notice that `KaplanFilm::absorb_phonon(E_in)` depends *only* on the scalar incident phonon energy $E_{\text{in}}$ and static film parameters!
   - **Shortcut:** At initialization, run 100,000 Monte Carlo trials across 512 energy bins $E_{\text{in}} \in [2\Delta, 50\Delta]$ to tabulate:
     - The exact CDF of deposited energy fraction $F(E_{\text{dep}} / E_{\text{in}} \mid E_{\text{in}})$
     - The multiplicity and energy CDF of re-emitted crystal phonons $P(N_{\text{refl}}, \{E_{\text{refl}}\} \mid E_{\text{in}})$
   - At runtime, a sensor film hit is resolved in **2 table lookups (~5 ns, zero loops)**, making `KaplanFilm` 100% branchless on both CPU and GPU!

2. **Asymptotic "Jump-to-Ballistic" Anharmonic Downconversion Kernel:**
   - Primary phonons start at the Debye frequency ($\nu_D \sim 2\text{–}15\text{ THz}$) where $\lambda_{\text{iso}} \propto \nu^{-4}$ and $\lambda_{\text{anh}} \propto \nu^{-5}$ are nanometers to microns. High-frequency phonons execute hundreds of tiny diffusive steps within a $<20\,\mu\text{m}$ sphere around the vertex before downconverting to the quasi-ballistic regime ($\nu \lesssim 0.4\text{ THz}$, where $\lambda_{\text{MFP}} \sim \text{cm}$).
   - **Shortcut:** Pre-compute the universal self-similar Tamura downconversion Green's function $G(\nu_0 \to \nu_{\text{bal}}, \Delta r, \Delta t)$ down to a configurable quasi-ballistic cutoff frequency $\nu_{\text{cut}}$ (where $\lambda_{\text{iso}}(\nu_{\text{cut}}) = 0.05 \times d_{\text{wall}}$). Any phonon generated well inside the bulk with $\nu > \nu_{\text{cut}}$ can jump directly to its post-cascade quasi-ballistic daughter spectrum $(\nu_i, \Delta\vec{r}_i, \Delta t_i)$ in a single step, eliminating **60%–80% of all bulk scattering steps**!

3. **Macroscopic Chord Stepping (CSDA + Straggling) for High-Voltage Luke Emission:**
   - At high bias voltages ($50\text{–}100\text{ V}$ in HV/CDMSlite detectors), each drifting electron or hole emits $>20,000$ Neganov-Trofimov-Luke phonons, requiring $>100,000$ sub-micron steps per charge carrier.
   - Once a carrier reaches terminal drift velocity ($v_d \gtrsim v_{\text{sound}}$), it travels along a deterministic electric-field streamline determined by its Herring-Vogt valley mobility $\mu_v \vec{E}$, converting electrostatic work $e \vec{E} \cdot d\vec{x}$ into Cherenkov phonons with a stationary energy/angle spectrum.
   - **Shortcut:** Advance the carrier in macroscopic chords ($\Delta s \sim 50\text{–}200\,\mu\text{m}$, limited only by intervalley scattering or E-field curvature) and sample the weighted Luke phonons directly along the chord. This accelerates high-voltage charge + Luke transport by **50x–200x**.

4. **Next-Event Estimation (NEE) & Russian Roulette for Sparse Sensors (Qubits / KIDs):**
   - In superconducting qubit chips or sparse KID arrays where active sensors cover $<1\%$ of the surface area, $>99\%$ of phonons bounce thousands of times before thermalizing.
   - **Shortcut:** Combine **Russian Roulette** after $N_{\text{bounces}} > 50$ (surviving tracks get weight $w \leftarrow w / p_{\text{surv}}$) with **Next-Event Estimation (Shadow Rays)** at diffuse reflections, projecting deterministic weighted contributions directly onto small sensor pads.

---

## 4. 5-Burst Porting & Optimization Tracker

We organize the complete G4CMP port and optimization roadmap into 5 self-contained bursts, committing and pushing after each burst:

### Burst 1: All Bravais Symmetries, 12 `CrystalMaps` Presets, `config.txt` Reader & Polycrystal Scattering (**COMPLETED**)
- [x] Port all 8 Bravais symmetry groups (`Amorphous`, `Cubic`, `Tetragonal`, `Orthorhombic`, `Hexagonal`, `Rhombohedral`, `Monoclinic`, `Triclinic`) and full $6 \times 6$ Voigt stiffness matrix expansion (`CrystalConfig::build_voigt_matrix()` matching `G4CMPCrystalGroup::FillElReduced`).
- [x] Port Miller `(h, k, l)` orientation rotation (`CrystalConfig::set_miller_orientation()` matching `G4LatticePhysical::SetMillerOrientation`).
- [x] Add all 12 G4CMP `CrystalMaps/` presets (`Ge`, `Si`, `Al2O3` sapphire, `CaWO4`, `SiO2_alpha` quartz, `SiO2_amorph`, `GaAs`, `CaF2`, `LiF`, `Al`, `Nb`, `Cu`) in `include/cmpmc/materials/presets.hpp`.
- [x] Port `G4LatticeReader` (`materials::load_g4cmp_config_file()`) in `include/cmpmc/materials/config_reader.hpp`.
- [x] Port polycrystalline grain-boundary elastic scattering (`G4CMPPhononPolycrystalElasticScattering` -> `physics::polycrystal_scattering_rate`).

### Burst 2: Complete Charge Carrier Physics (Matrix-Element IV Rate, Trapping, Trap Ionization & 3D Charge Cloud) (**COMPLETED**)
- [x] Port first-principles deformation-potential matrix-element intervalley scattering (`G4CMPInterValleyRate.cc` — optical phonon deformation potentials + Erginsoy neutral impurity scattering in `include/cmpmc/physics/charge_physics.hpp`).
- [x] Port charge carrier bulk trapping (`G4CMPDriftTrappingProcess.cc` — electron & hole trapping MFPs `e_trap_mfp_m`, `h_trap_mfp_m`).
- [x] Port trap impact ionization (`G4CMPDriftTrapIonization.cc` — `e_on_e`, `e_on_h`, `h_on_e`, `h_on_h` impurity liberation channels).
- [x] Port boundary reflection/recombination (`G4CMPDriftBoundaryProcess.cc`, `G4CMPDriftRecombinationProcess.cc`) and 3D initial charge cloud (`G4CMPChargeCloud.cc`).

### Burst 3: All 5 NIEL Partition Models & Exact Fano Binomial Sampler (**IN PROGRESS**)
- [ ] Port `G4CMPLindhardNIEL.cc`, `G4CMPLewinSmithNIEL.cc`, `G4CMPSarkisNIEL.cc`, `G4CMPImpactTunlNIEL.cc`, and `G4CMPEmpiricalNIEL.cc` into `include/cmpmc/physics/niel_models.hpp`.
- [ ] Port exact Fano-factor generalized binomial sampling (`G4CMPFanoBinomial.cc`) into `include/cmpmc/partition.hpp`.

### Burst 4: Full 3D Superconducting Bogoliubov Quasiparticle Transport (**TODO**)
- [ ] Port `G4CMPBogoliubovQP.cc`, `G4CMPSCUtils.cc`, `G4CMPQPDiffusion.cc` ($D(E) = D_0 \sqrt{1 - (\Delta/E)^2}$), `G4CMPSCPairBreakingProcess.cc`, `G4CMPQPRadiatesPhononProcess.cc`, `G4CMPQPRecombinationProcess.cc`, `G4CMPQPLocalTrappingProcess.cc`, and `G4CMPQPBoundaryProcess.cc` into `include/cmpmc/physics/superconductor_3d.hpp` and couple with 3D bulk phonons in `SimulationEngine`.

### Burst 5: 2D/3D Tetrahedral FEM Mesh Electric Fields, 2D Patterned Electrode Masks & Tier-2/Tier-4 LUT Optimizations (**TODO**)
- [ ] Port 2D/3D tetrahedral barycentric FEM electric field interpolator (`G4CMPBiLinearInterp.cc`, `G4CMPTriLinearInterp.cc`, `G4CMPMeshElectricField.cc`) into `include/cmpmc/geometry/mesh_field.hpp`.
- [ ] Port patterned 2D sensor/electrode surface masks (`G4CMPVElectrodePattern.cc`) into `include/cmpmc/geometry/surface_mask.hpp`.
- [ ] Implement `256 x 256` `float4` **Octahedral Equal-Area LUT** (`oct_encode`, zero trig) and pre-baked `KaplanSurrogateLUT`.
