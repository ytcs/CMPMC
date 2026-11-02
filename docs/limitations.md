# Known limitations and planned work

Each item below gives the current behaviour, the target model and its source, and
the test a future implementation has to pass. Items marked *experimental* are guarded
at runtime ([experimental.md](experimental.md)).

## Superconducting films

### Pair-breaking efficiency beyond the lumped Kaplan cascade
- **Current:** `FilmLossModel::KaplanRateCompetition` uses the exact T = 0 Kaplan τ_pb(Ω)
  and a film-loss time τ_l. It reproduces arXiv:1401.1937 Fig. 4 to ≤ 0.011 for
  τ_l/τ₀^φ = 0.5–4. At τ_l/τ₀^φ = 8 it comes out 0.041 low, because the event cascade still
  loses 6–11 % of the 2Δ–4Δ phonons where the kinetic-equation result has already
  saturated.
- **Missing:**
  1. The film's Eliashberg α²F(Ω) in place of the Debye Ω² law.
  2. T > 0: thermal factors, stimulated emission, and recombination with thermal and
     driven QP populations (the Chang–Scalapino kinetic equations, arXiv:1401.1937
     Eqs. (1)–(5)).
  3. Thickness enters only through τ_l = 4d/(η_t v_s). A depth-resolved model with
     mode-resolved AMM transmission is missing.
- **Approach:** a deterministic solver for the coupled f(E), n(Ω) kinetic equations
  that writes an offline table η_pb(Ω/Δ, τ_l/τ₀^φ, T/T_c) for `KaplanFilm` to read.
- **Test:** all five curves of arXiv:1401.1937 Fig. 4 to 0.01, including τ_l/τ₀^φ = 8.

### Bulk 3D quasiparticle transport (*experimental*)
- **Current:** QPs recombine only with thermal QPs (no QP–QP / Rothwarf–Taylor term).
  The gap is fixed, the rate tables are clamped above 60Δ, boundaries are diffuse only,
  and the QP velocity is a placeholder (10⁶ m/s).
- **Approach:** Rothwarf–Taylor dN/dt = −RN² + 2βN_ph (N_ph: number of 2Δ phonons) with R from Chang–Scalapino. QP–QP
  recombination becomes pair events on a local density grid. The rate tables are
  replaced by `KaplanFilm::kaplan_pair_breaking_rate`, which is exact for any Ω. v_F and
  D_n come from a sourced Al film.
- **Test:** the Rothwarf–Taylor steady state (arXiv:1401.1937 Eq. (4)) and QP-number
  decay in a uniform box.

## Phonons

### Anisotropic acoustic mismatch and thin-film interference
- **Current:** the isotropic AMM/DMM in `phonon/interface_transmission.hpp`.
- **Approach:** for a given k_∥, solve the sextic for the six k_z roots in each
  anisotropic medium and select them by the sign of v_g,z. Then solve the 6×6
  continuity of displacement and traction, and compute energy-flux coefficients from
  the anisotropic Poynting vector. Film interference matters for λ ≳ d_film (sub-THz
  phonons in 40–60 nm films) and needs a transfer matrix through the film.
- **Test:** flux conservation, the isotropic limit equal to the current AMM, and
  reciprocity.

### Created-phonon spectrum of the prompt cascade
- **Current:** prompt phonons are created as a single line at hν_D with the DOS mode
  mix.
- **Approach:** a branch-resolved spectrum for nuclear recoils and an electron–phonon
  emission kernel for electron recoils, with molecular-dynamics tables above about 1 eV.
- **Test:** energy closure and the first moments of each population.

### Dispersive phonons and transverse-mode decay
- **Current:** the isotropic continuum approximation with linear dispersion. Only
  longitudinal phonons decay.
- **Approach:** tabulate first-principles ω_ν(q) on the octahedral grid with
  v_g = ∇_q ω, and use Tamura three-phonon rates with the full cubic tensor, including
  the T → T + L channels.
- **Test:** the low-q limit equal to the Christoffel solve, and the rebuilt A and B
  constants.

### Faster exact kinematics
- **Current:** the closed-form Christoffel solve is accurate to 1e-11 but only 1.3–1.65×
  faster than Jacobi. Most of the cost is in the two rank-4 tensor contractions.
- **Approach:** contract with the 6×6 Voigt matrix and use the sparsity of the cubic,
  trigonal and tetragonal systems (cubic: about 20 flops).
- **Test:** `test_christoffel_cardano`, unchanged.

## Charge carriers

### Trap kinetics
- **Current:** constant trapping mean free paths.
- **Approach:** occupation-dependent capture and emission,
  dn_j/dt = σ_j(u) u n_f (N_j − n_j) − Γ_em n_j − σ_rec u n_f′ n_j, on a per-detector voxel
  grid shared by all carriers of an event. Trapped charge adds a space-charge field.
  Every capture and release emits its binding energy as phonons, through the
  `ChargeRelease` machinery.
- **Test:** the Hecht limit at constant σ and N ≫ n; detailed balance at equilibrium;
  energy closure.

### Intervalley Matrix rate (*experimental*)
- **Defect:** `intervalley_matrix_rate` uses a final-valley count Z_f = 2N − 1 for every
  phonon type, which in Si is about 4× the per-type value. It also adds the Erginsoy
  neutral-impurity term to the intervalley rate, although that scattering is
  intravalley.
- **Approach:** a per-phonon final-valley set (g/f type), and Erginsoy as a separate
  intravalley elastic process. Z_f values must come from a primary source (Jacoboni and
  Reggiani, Rev. Mod. Phys. 55, 645 (1983)).
- **Test:** the closed-form rate at several energies; Si drift velocity vs. field at
  77 K against measured data.

### Non-parabolic Herring–Vogt kinematics
- ħ²k\*²/2m_c = E(1 + αE). For Luke-dominated transport at 1–100 V/cm the carrier
  kinetic energy is ≲ 10 meV, so αE ≤ 0.005. The correction is below the stepping
  tolerance and has not been implemented.

## Variance reduction and alternative engines (*experimental*)

### SoA and wavefront engines
- **Defects:** the SoA engine has only the L → TT branch, isotropic daughters and diffuse
  walls, and keeps only the first re-emitted phonon. The wavefront pipeline has
  isotropic daughters. Neither engine conserves energy to the standard of the scalar
  engine.
- **Approach:** route each lane through the scalar engine's decay, surface and film
  kernels.
- **Test:** the same arrival-time histogram (KS test) and collected energy (1e-9) as
  `SimulationEngine::transport` on a fixed seed set.

### Jump-to-ballistic
- **Defects:** a fixed-radius random displacement instead of a space–time distribution,
  rounded daughter counts, a cap of 48 daughters, and isotropic daughter wavevectors.
- **Approach:** tabulate the quasi-diffusion Green's function G(f₀ → f_b; r, t) from full
  infinite-medium transport runs and sample (r, t, f, mode) jointly from it.
- **Test:** arrival-time and position distributions against `transport` (KS) and
  energy closure.

### Next-event estimation
- **Defect:** the estimator is Lambertian in real space and ignores phonon focusing.
  It is biased whenever v_g is not parallel to k̂, and it is not wired into the engine.
- **Approach:** estimate in k-space by sampling k̂ in the pre-image of the pad under the
  v_g map.
- **Test:** unbiasedness against analogue transport on a sparse-pad geometry.
