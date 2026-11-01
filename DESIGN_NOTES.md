# CMPMC design notes for the items not implemented in R16

Numbering follows the R15 CMPMC-AUD "SOTA upgrade list" (§6), then the unsupported modules (§12–§14).
Each note states the target, the formula and its source, the reference the implementation must reproduce (the future CTest), and why it was not done in R16
(the brief's rule: an item needing more than about 2 h of work gets a design note). Tags as in VALIDATION.md.

## 1. η_pb(Ω, d_film, T) beyond the lumped Kaplan cascade (R16 did the reduced scope)
- **Done in R16:** `FilmLossModel::KaplanRateCompetition` (exact T = 0 Kaplan τ_pb(Ω), loss time τ_l), case-II photon seed, `pair_breaking_efficiency()` → η_pb, F_qp.
  It reproduces 1401.1937 Fig. 4 [FT] to ≤ 0.011 for τ_l/τ₀^φ = 0.5–4. At τ_l/τ₀^φ = 8 it is 0.041 low. There the kinetic-equation result sits at its own thick-film value, while the
  event cascade still loses 6–11% of the 2Δ–4Δ phonons.
- **Remaining:** (i) the Eliashberg α²F(Ω) of the film instead of the Debye Ω² law in the QP-emission kernel (Al: tabulated α²F, source to fetch); (ii) T > 0: thermal factors
  (1 − f − f′) and stimulated emission, and recombination with the thermal and driven QP populations (Chang–Scalapino kinetic equations, 1401.1937 Eqs. (1)–(5) [FT], CS (A9), (A10));
  (iii) the "E = Δ for the excess" convention of 1401.1937 (+4%, its §V) vs counting QPs at their energies; (iv) d_film enters only through τ_l = 4d/(η_t v_s). A d-resolved model
  (phonon mean free path vs depth, AMM η_t per mode from `interface_transmission.hpp`) belongs here too.
- **Implementation:** a deterministic solver of the coupled f(E), n(Ω) kinetic equations on a 1 µeV grid with γ-broadened ρ(E) (1401.1937 §II [FT]). The solver writes an offline table
  η_pb(Ω/Δ, τ_l/τ₀^φ, T/T_c), and `KaplanFilm` reads it. **Test:** all five Fig. 4 curves, including τ_l/τ₀^φ = 8, to 0.01, and Fig. 5 (probe on/off). **Effort:** 3 d for the table route.

## 2. Anisotropic AMM and thin-film interference (R16 did the isotropic AMM/DMM)
- **Remaining:** the substrate is anisotropic. Use the Christoffel-matched boundary problem: for a given k_∥, solve the sextic for the six k_z roots in each medium, select by the sign of
  v_g,z (not of k_z), and solve the 6×6 continuity of u and of the tractions σ·n̂. Energy-flux coefficients use the anisotropic Poynting vector. Film interference applies for
  λ ≳ d_film (sub-THz NTL phonons in 40–60 nm films, campaign WP8-R §3.2): a transfer matrix through the film onto a substrate.
- **Test:** conservation Σ = 1 (as in test_interface_amm), the isotropic limit equal to `amm_psv`, and reciprocity. **Effort:** 1–2 weeks.

## 5. Trap kinetics (campaign WP7-E-v3 (E9)–(E12)) instead of constant λ_tr
- **Formula:** dn^e_j/dt = σ^e_j(u)u n_e^f(𝓃_j − n^e_j) − Γ^em_{e,j}n^e_j − σ^rec_j u n_h^f n^e_j (and e ↔ h) (E9). Trapped charge has a space-charge field E_es, which enters
  `electric_field_at`. The ledger rows are E_trap (μ_F-referenced, Thm E-3) and release phonons E^tr → ℋ (E1).
- **Design:** a per-detector voxel grid of n_j(r, t) shared by all carriers of an event, or across events for the history. Capture rate σ(u)u(𝓃 − n) replaces v/λ. Emission is a
  Poisson process per occupied trap that creates a carrier with the time tag. Every capture and release emits E^tr as phonons, through the BulkDebye or ElectrodeFilm machinery of R16.
- **Test:** the Hecht limit (E10) at constant σ and 𝓃 ≫ n, two routes as in C4a; detailed balance of capture/emission at equilibrium; (E1) closure with E_trap. **Effort:** about 1 week.

## 6. IV Matrix rate with per-type Z_f and intra-valley impurity scattering (module guarded)
- **Defect:** `intervalley_matrix_rate` uses Z_f = 2N − 1 for every phonon. In Si that is ×3.96 the per-type value (R15 T6). The Erginsoy neutral-impurity term is added to
  the *inter*valley rate, but neutral-impurity scattering is intra-valley momentum relaxation.
- **Source search (R16):** Jacoboni & Reggiani, RMP 55, 645 (1983) is closed access (Unpaywall: `is_oa = false`). The arXiv searches (lit.py) found no Si g/f-phonon Z_f table. Three
  candidates were fetched and read without result: cond-mat/0506045, 0806.4347 and 2605.18581 (text kept in literature/papers-CMPMC-DEV/). The usual values, g: Z = 1 and f: Z = 4
  (Si), and Z = 3 (Ge L–L), are **[UNVERIFIED]**. They must not enter code as constants until a source is read.
- **Design:** a config vector `iv_final_valleys[i]` per optical phonon i, plus an `iv_type[i] ∈ {g, f}` that selects the final-valley set when sampling. Erginsoy moves to a separate
  intra-valley elastic process: it randomises k̂ within the valley and keeps the valley.
- **Test:** the J&R rate formula W = Z_f m_d^{3/2}D_f²√(E ∓ ħω)(1 + 2αE)[n + ½ ∓ ½]/(√2πρħ²ω) in closed form at several E. Drift velocity vs field in Si at 77 K against the
  measured Canali/Jacoboni data, read in full text. **Effort:** 1 d once the source is in hand.

## 7. Created-phonon spectrum of the prompt cascade
- **Formula:** branch-projected ρ^ph_κ (WP2-v3 (3.6)) for NR at x_R ≲ 1; the e–ph emission kernel of WP5-v4 §2 for ER; MD tables for NR above an eV. This replaces the mono-energetic hν_D line
  in `EnergyPartition` (campaign WP8-R IF-R1).
- **Test:** energy closure (W4) and the first moments of each population against the WP2/WP5 kernels. **Effort:** M (about 1 week); it depends on campaign tables not yet in hand.

## 8. Dispersive phonons and T-mode decay beyond the isotropic continuum approximation
- **Design:** DFPT ω_ν(q) (QE ph.x, available via run-heavy.sh) tabulated on the octahedral grid, with v_g = ∇_qω. Tamura three-phonon rates with the full cubic tensor, including the T → T + L
  channels (2408.04732 §2.1 refs [30–32] [FT]). **Test:** the low-q limit equal to the Christoffel solve (test_christoffel_cardano), and A, B rebuilt (T3). **Effort:** L (weeks).

## 9. Bulk 3D Bogoliubov QP transport (module guarded)
- **Defects (R15):** recombination only with thermal QPs (no QP–QP or Rothwarf–Taylor term), fixed gap, LUTs clamped above 60Δ, diffuse boundaries only, and `qp_velocity` = 10⁶ m/s as a
  placeholder.
- **Design:** Rothwarf–Taylor dN/dt = −RN² + 2βN_2Δ (1401.1937 Eqs. (2)–(3) [FT]), with R from Chang–Scalapino (A9) (R ≡ 2R_CS, 1401.1937 §II). QP–QP recombination becomes pair events in a
  local density grid. The τ_pb LUT is un-clamped by using `KaplanFilm::kaplan_pair_breaking_rate` (exact for any Ω, R16). v_F and D_n come from a sourced Al film.
- **Test:** the RT steady state N² = N_T² + Γ(βτ_l + 1)/R (1401.1937 Eq. (4)) and QP number decay in a uniform box. **Effort:** M.

## 10. Non-parabolic Herring–Vogt kinematics
- **Formula:** ħ²k*²/2m_c = E(1 + αE), v = ħk*/(m_c(1 + 2αE)) (in the HV frame), with the Luke threshold and rate re-derived for the modified density of states. α is from
  Jacoboni–Reggiani **[UNVERIFIED]** (Si ≈ 0.5 eV⁻¹, Ge ≈ 0.3 eV⁻¹).
- **Why not done:** NTL carriers at 1–100 V/cm have E_kin ≲ 10 meV, so αE ≤ 0.005 and the correction to v and to the Luke rate is ≤ 1% (below the 2–3% stepping tolerance). The
  value of α also needs the source of §6. Effort S once α is sourced.

## 11. Exact kinematics as default (R16 did the closed-form solve)
- The Cardano solve is accurate to 1e-11 but only 1.3–1.65× faster than Jacobi. The cost is in the two rank-4 contractions, D_il and v_g (≈ 0.5 µs at −O3), not in the eigen-solve.
  **Design:** contract with the 6×6 Voigt matrix and use the cubic, trigonal and tetragonal sparsity (cubic: D from 3 constants, about 20 flops). Then compare with the 52 ns of the
  octahedral LUT. **Test:** test_christoffel_cardano unchanged. **Effort:** S.

## 12. SoA / GPU wavefront engines (guarded)
- **Defects (R15):** SoA has L→TT only (no L′T branch), isotropic daughters, diffuse-only walls and keeps only the first re-emitted phonon. GPU has isotropic daughters. The benchmark's
  wavefront energy balance is 9.70 / 14.29 eV deposited/injected.
- **Design:** route each lane's decay through `perform_anharmonic_decay` (Leman angles, P1/P2), the surface through `reflection_split` and `amm_interface`, and push every re-emitted phonon
  (KaplanBuffer). **Test:** the same arrival-time histogram and collected energy as `SimulationEngine::transport` on a fixed seed set (KS test, energy to 1e-9). **Effort:** M.

## 13. Jump-to-ballistic (guarded)
- **Defects:** a fixed-radius random shell displacement instead of a space–time distribution, rounded N_iso, a cap of 48 daughters, isotropic k̂, and the L fraction set by fLDOS
  rather than the steady mix.
- **Design:** tabulate the quasi-diffusion Green's function G(f₀ → f_b; r, t) from full `transport` runs in an infinite medium (per material, per f₀). Sample (r, t, f, mode) jointly
  from it, and accept a jump only where the clearance exceeds the q99 of r. **Test:** the arrival-time and position distributions against `transport` (KS) and energy closure.
  **Effort:** M.

## 14. Next-event estimation (guarded)
- **Defect:** a Lambertian real-space estimator. It ignores phonon focusing, so the estimate is biased whenever v_g ≠ k̂ direction, and it is not wired into the engine.
- **Design:** use the estimator in k-space. For a diffuse re-emission at x, the contribution to pad P is ∫ over k̂ (Lambertian in k̂) of 1[ray along v_g(k̂) hits P]
  × attenuation. Evaluate it by sampling k̂ in the pre-image of P under the v_g map (octahedral LUT inverse, or rejection with a focusing Jacobian). **Test:** unbiasedness against analogue
  transport on a sparse-pad geometry (mean within 2σ, variance reduction reported). **Effort:** M.
