# CMPMC validation table

Every module, the CTest target that checks it (`tests/`, registered in `tests/CMakeLists.txt`), the reference it is checked against, and the tolerance.
Build with `cmake -B build && cmake --build build && (cd build && ctest)`. Each test prints `name value reference rel_err PASS/FAIL` lines.

Status: **validated** (checked against the reference in the row), **validated, opt-in** (validated but not the default), **port** (identical to G4CMP, not independently re-derived),
**unsupported** (no reference path; the entry point throws `UnvalidatedModuleError` unless `allow_unvalidated` is set; results must not feed a campaign number).
Tags: [FT] source read in full text; [DB] database value; [D] derived here; [E] estimate; [UNVERIFIED] quoted, not checked.

## Phonons
| Module | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| Christoffel solve (Jacobi), v_g | test_christoffel (T1) | [100]/[110]/[111] closed forms; finite-difference ∇_kω | 1e-10; 1e-9 | validated |
| Christoffel closed-form (Cardano), `KinematicsMode::Exact` | test_christoffel_cardano | Jacobi solve on 2e5 directions × 6 crystals (cubic, trigonal, tetragonal); [100]/[110]/[111] closed forms | v_p 5e-13, v_g 2e-11 (gap > 1e-3; Jacobi fallback at degeneracies); 1e-12 | validated (1.3–1.65× faster than Jacobi only; not default) |
| Octahedral FP32 v_g LUT (default kinematics) | test_christoffel (T1) | exact solve | q99 direction error 0.2 mrad (L), 2.5 mrad (T) | validated |
| Bilinear 181×361 v_g LUT | test_christoffel (T1) | exact solve | q99 14–16 mrad near z (reported) | legacy, not default |
| Isotope B, rates Bν⁴, Aν⁵ | test_presets (T3), test_rates (T4) | B rebuilt from isotope abundances (2408.04732 Eqs. (3)–(9) [FT]) | 10% (observed 0.5–2.7%) | validated (ICA, f ≪ f_D) |
| L→L′T, L→TT energy split and angles | test_anharmonic_decay (T2) | Leman 1109.1193 Eqs. (9)–(14) [FT] | moments 1e-3; Σk = k₀ to 1e-13 | validated (ICA, T = 0) |
| Presets (C_ij, A, B, DOS, ν_D, ε_eh, F) | test_presets (T3) | G4CMP config.txt (1e-6); 2408.04732 T3–T5 [FT] (2–5%); WP5-v4 §6.4 (ε_eh 0.5%, F 5%); ν_D(Ge) 1% | as listed | validated |
| Surface reflection, campaign default (diffuse, p_anh = 0) | test_surface (c) | Leman 1109.1193 §IV.G [FT] ("diffusive scattering ... a good model") | exact | validated (no sourced specular / down-conversion coefficients exist) |
| Surface reflection, G4CMP demo polynomials (`g4cmp_demo_*`) | test_surface (a), (b) | G4CMPSurfaceProperty / PhononBoundaryProcess formulas | 1e-9; branching 2% | port of demonstration values: **unsupported** (`demo_values`, guarded) |
| AMM/DMM interface transmission (isotropic), `SurfaceProperty::amm_interface` | test_interface_amm | energy-flux conservation R+T = 1 (not imposed); 4Z₁Z₂/(Z₁+Z₂)² normal limit; identical media; SH and P-SV reciprocity; campaign WP8-R §3.2 table (24 entries); DMM detailed balance | 1e-10; 1e-12; 1e-12; 1e-10; 6e-4 abs; 1e-14 | validated, opt-in (isotropic half-space, ideal interface = upper bound; measured Si→Al 0.3–0.55, 1805.02495 [FT] via WP8-R) |
| Prompt phonons (Debye line, DOS mode mix) | test_fano_partition (T7) | WP5-v4 (W4) bookkeeping | 1e-15 closure | validated as bookkeeping; no created-phonon spectrum (DESIGN_NOTES §7) |
| Engine energy closure (phonon run) | test_headers_compile (T10) | deposited + escaped = injected | 1e-9 | validated |

## Superconducting films
| Module | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| Kaplan pair-breaking / QP-emission samplers | test_kaplan (T8) | Kaplan densities (moments) | 2e-3 | validated (T = 0) |
| η_pb no-loss limit | test_kaplan (T8), test_pair_breaking (d) | 2Δ/Ω for Ω < 4Δ (WP5-v4 (W13) limits); 0.59 asymptote (1401.1937 [FT]) | 2e-3; 0.01 abs | validated |
| Exact T = 0 Kaplan τ_pb(Ω) | test_pair_breaking (a) | 1/τ_pb(2Δ) = 1/τ₀ closed form; brute-force quadrature | 1e-9; 1e-5 (observed 3e-11) | validated |
| η_pb(Ω, τ_l/τ₀^φ), F_qp: `FilmLossModel::KaplanRateCompetition`, `pair_breaking_efficiency()` | test_pair_breaking (b), (c), (e) | 2Δ/Ω below 4Δ, all τ_l; 1401.1937 Fig. 4 [FT] (digitised, ±0.005), photon (case II) seed | exact below 4Δ; 0.015 abs for τ_l/τ₀^φ = 0.5–4 (observed ≤ 0.011); 0.045 at τ_l/τ₀^φ = 8 (observed 0.041, documented model difference) | validated, opt-in. **Reduced scope:** lumped film, Kaplan/Debye α²F, T = 0, no probe; Chang–Scalapino / Eliashberg / T > 0 → DESIGN_NOTES §1 |
| Cascade energy closure at Debye-scale energies | test_pair_breaking (f) | absorbed + re-emitted = incident | 1e-12 per event at 50–350Δ | validated (buffer 256; was 64 with silent loss > ~130Δ) |
| G4CMP film escape law | test_kaplan (T8) | exact Kaplan τ_pb | within 6% on 2.5–40Δ (reported) | port (legacy default) |
| KaplanSurrogateLUT | test_kaplan (T8) | full cascade in a re-entry cavity | 3% (observed < 1%) | validated |
| Bulk 3D Bogoliubov QP transport (`transport_superconductor_coupled`) | test_guards | — | — | **unsupported** (guarded; DESIGN_NOTES §9) |

## Charges
| Module | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| ℓ₀, Luke rate, angle and energy | test_luke_ntl (T5) | Leman 1109.1193 Table II, Eq. (22) [FT] | 5%; 1e-4; 2e-3 | validated (parabolic, T = 0) |
| First-order Δk-capped stepping (default) | test_luke_ntl (T5) | WP7-E-v3 (E11) closure; step convergence 0.05 vs 0.02 k | 1e-3; 5% | validated |
| Leman second-order stepping (`ChargeStepping::LemanSecondOrder`) | test_leman_stepping | first-order with Δk ≤ 0.01 k (Ge, Si; e⁻, h⁺; 1–100 V/cm); (E11) | ⟨E_Luke⟩ and phonons/carrier 3% (observed ≤ 2.4%); closure 1e-3 | validated, opt-in (1.3–4.5× fewer steps than P5) |
| `transport_charges_fast_chord` | test_leman_stepping (d), test_guards | identical to the Leman path | exact | validated (now the Leman path; former fixed-Mach chord removed) |
| HV valleys, orientation | test_charge_misc (T6) | covariance under rotation; E_kin, |v| closed forms along the valley axis | 1e-12; 1e-6 | validated (parabolic) |
| IV Quadratic / Linear | test_charge_misc (T6) | G4CMP formulas | 1e-12 | port (empirical fits, fitted field range) |
| IV Matrix (`IVRateModel::Matrix`) in transport | test_guards; test_charge_misc reports ×3.96 vs Jacoboni–Reggiani Z_f | — | — | **unsupported** (guarded; DESIGN_NOTES §6) |
| Collection release E_rel = E_kin + (E_c − μ_F) / (μ_F − E_v) | test_charge_release (a), (b) | WP7-E-v3 §1.3 ledger ("collection"); e⁻+h⁺ at rest = E_g | 1e-4; closure 1e-12 | validated (μ_F an input; default mid-gap [E]) |
| `ChargeRelease::BulkDebye` | test_charge_release (c) | phonons at hν_D | 1e-9 | validated, opt-in (bookkeeping) |
| `ChargeRelease::ElectrodeFilm` (soft surface-recombination phonons) | test_charge_release (d), (e) | closure; seed-quantum convergence 80 vs 160Δ | 1e-12; 2% (observed ≤ 1.2%) | validated as numerics, opt-in; the soft spectrum (mean ≈ 1.0Δ_Al, 93% < 2Δ_Al for 600 nm Al on Ge) is a model prediction, not yet compared with data |
| Trapping, trap ionisation (constant λ) | test_charge_misc (T6) | G4CMP | — | port; (E9) kinetics → DESIGN_NOTES §5 |
| Fano sampler, (W4) closure | test_fano_partition (T7) | ⟨N⟩, Var/⟨N⟩ = F (18 cases); WP5-v4 (W4) | 1%, 5% (observed < 1%); closure | validated |

## Shortcuts and alternative engines
| Module | Test | Status |
|---|---|---|
| Jump-to-ballistic (`jump_to_ballistic.enabled`) | test_guards | **unsupported** (guarded; DESIGN_NOTES §13) |
| Next-event estimation (`evaluate_sparse_sensor_nee`) | test_guards | **unsupported** (guarded; DESIGN_NOTES §14) |
| Russian roulette | — | unbiased (w/p); usable |
| `WavefrontSoAEngine`, `gpu::WavefrontPipelineEngine` | test_guards | **unsupported** (guarded; DESIGN_NOTES §12) |

## Constants that are still G4CMP example values
`materials::AluminumFilm()` (the default sensor film): `film_absorption` 0.20, `subgap_absorption` 0.10, `gap_energy` 173.715 µeV, `phonon_lifetime_slope` 0.29.
`phonon_lifetime` 242 ps is Kaplan's τ₀^φ(Al) as used by G4CMP (1401.1937 [FT] uses 0.26 ns). Use `amm_interface` and `KaplanRateCompetition`, with `subgap_absorption = 0`, for campaign runs.
