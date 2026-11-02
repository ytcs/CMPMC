# Validation

Each model is listed below with the test that checks it, the reference it is compared
against, and the tolerance. The tests are in `tests/` and run with

```bash
cmake -B build && cmake --build build -j && ctest --test-dir build
```

Each test prints one line per check (`name value reference rel_err PASS|FAIL`) and exits
nonzero if any check fails. CTest labels let you select subsets: `ctest -L phonon`,
`-L charge`, `-L superconductor`, `-L materials`, `-L transport`, `-L validation`,
`-L unit`.

**Status values**

- **validated**: checked against the reference in the row. These models are the defaults.
- **validated, opt-in**: validated, but you have to enable it explicitly.
- **reference parameterisation**: the same formula and parameters as G4CMP, with no
  independent derivation here.
- **experimental**: no reference comparison exists yet. The entry point throws
  `cmpmc::UnvalidatedModuleError` unless you set `allow_unvalidated`. See
  [experimental.md](experimental.md).

Every model assumes the isotropic continuum approximation (ICA) and T = 0 unless the row
says otherwise.

## Phonons

| Model | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| Christoffel solve (Jacobi) and group velocity | `test_christoffel` | Closed forms along [100], [110], [111]; finite-difference ∇ₖω | 1e-10; 1e-9 | validated |
| Closed-form Christoffel solve (trigonometric Cardano), `KinematicsMode::Exact` | `test_christoffel_cardano` | Jacobi solve on 2·10⁵ directions for 6 crystals (cubic, trigonal, tetragonal); closed forms | v_p 5e-13, v_g 2e-11 (falls back to Jacobi at degeneracies) | validated, opt-in |
| Octahedral FP32 group-velocity table (default kinematics) | `test_christoffel` | Exact solve | 99th-percentile direction error 0.2 mrad (L), 2.5 mrad (T) | validated |
| Bilinear 181×361 (cos θ, φ) table, `KinematicsMode::Bilinear` | `test_christoffel` | Exact solve | 99th-percentile error 14–16 mrad near the poles (reported) | alternative, not default |
| Isotope scattering constant B and rates Bν⁴, Aν⁵ | `test_presets`, `test_rates` | B rebuilt from isotope abundances (Tamura; arXiv:2408.04732 Eqs. (3)–(9)) | 10 % (observed 0.5–2.7 %) | validated (f ≪ f_D) |
| L → L′T and L → TT energy split and daughter angles | `test_anharmonic_decay` | Leman, arXiv:1109.1193 Eqs. (9)–(14) | moments 1e-3; Σk = k₀ to 1e-13 | validated |
| Material presets (C_ij, A, B, DOS, ν_D, ε_eh, F) | `test_presets` | arXiv:2408.04732 Tables 3–5 (2–5 %); Si ε_eh and F from arXiv:2004.11499 | as listed in the test | validated |
| G4CMP `config.txt` import (`io::import_g4cmp_config`) | `test_config_import` | Synthetic file: keywords, units, comments, unknown keys, defaults, missing file | exact / 1e-12 | validated |
| Presets vs. G4CMP `CrystalMaps` | `test_g4cmp_import` | G4CMP `config.txt` read through `io::import_g4cmp_config` | 1e-6 | runs only when `CMPMC_G4CMP_CRYSTALMAPS_DIR` is set |
| Surface reflection, default (fully diffuse, no surface down-conversion) | `test_surface` | Leman, arXiv:1109.1193 §IV.G | exact | validated |
| Surface reflection, polynomial model (`SurfaceProperty::demo_polynomial_*`) | `test_surface` | Specular/diffuse/anharmonic three-way split formulas | 1e-9; branching 2 % | experimental (example coefficients) |
| Isotropic AMM / DMM interface transmission, `SurfaceProperty::amm_interface` | `test_interface_amm` | Energy-flux conservation R + T = 1 (not imposed by the solver); 4Z₁Z₂/(Z₁+Z₂)² at normal incidence; identical media; SH and P-SV reciprocity; independently computed reference table (24 entries); DMM detailed balance | 1e-10; 1e-12; 1e-12; 1e-10; 6e-4 abs; 1e-14 | validated, opt-in (ideal interface, an upper bound on transmission) |
| Prompt phonons (Debye line, DOS mode mix) | `test_fano_partition` | Energy bookkeeping of the partition | 1e-15 | validated as bookkeeping (no created-phonon spectrum, see [limitations.md](limitations.md)) |
| Engine energy closure, phonon run | `test_engine_energy_closure` | deposited + escaped = injected | 1e-9 | validated |

## Superconducting films

| Model | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| Kaplan pair-breaking and QP phonon-emission samplers | `test_kaplan` | Kaplan et al., PRB 14, 4854 (1976): moments of the densities | 2e-3 | validated |
| Pair-breaking efficiency η_pb, no-loss limit | `test_kaplan`, `test_pair_breaking` | 2Δ/Ω for Ω < 4Δ; asymptote ≈ 0.59 (arXiv:1401.1937) | 2e-3; 0.01 abs | validated |
| Exact Kaplan pair-breaking time τ_pb(Ω) | `test_pair_breaking` | 1/τ_pb(2Δ) = 1/τ₀ (closed form); brute-force quadrature | 1e-9; 1e-5 | validated |
| η_pb(Ω, τ_l/τ₀^φ) and F_qp: `FilmLossModel::KaplanRateCompetition`, `KaplanFilm::pair_breaking_efficiency()` | `test_pair_breaking` | 2Δ/Ω below 4Δ for all τ_l; arXiv:1401.1937 Fig. 4 (digitised, ±0.005); photon seed | exact below 4Δ; 0.015 abs for τ_l/τ₀^φ = 0.5–4; 0.045 at τ_l/τ₀^φ = 8 | validated, opt-in (lumped film, T = 0, see [limitations.md](limitations.md)) |
| Cascade energy closure at 50–350 Δ | `test_pair_breaking` | absorbed + re-emitted = incident | 1e-12 per event | validated |
| Mean-free-path film escape law, `FilmLossModel::MeanFreePath` (default) | `test_kaplan` | Exact Kaplan τ_pb | within 6 % over 2.5–40 Δ (reported) | reference parameterisation |
| `KaplanSurrogateLUT` | `test_kaplan` | Full cascade in a re-entry cavity | 3 % (observed < 1 %) | validated |
| Bulk 3D quasiparticle transport (`transport_superconductor_coupled`) | `test_guards` | none | n/a | experimental |

## Charge carriers

| Model | Test | Reference | Tolerance | Status |
|---|---|---|---|---|
| Luke scattering length ℓ₀, rate, emission angle and energy | `test_luke_ntl` | Leman, arXiv:1109.1193 Table II, Eq. (22) | 5 %; 1e-4; 2e-3 | validated (parabolic bands) |
| First-order Δk-capped stepping (default) | `test_luke_ntl` | Energy closure (field work = Luke phonon energy + Δ kinetic energy); step convergence 0.05 vs 0.02 k | 1e-3; 5 % | validated |
| Second-order stepping, `ChargeStepping::LemanSecondOrder` | `test_leman_stepping` | First-order stepping with Δk ≤ 0.01 k (Ge, Si; e⁻, h⁺; 1–100 V/cm); energy closure | ⟨E_Luke⟩ and phonons per carrier 3 %; closure 1e-3 | validated, opt-in (1.3–4.5× fewer steps) |
| `transport_charges_fast_chord` | `test_leman_stepping`, `test_guards` | Identical to the second-order path | exact | validated |
| Herring–Vogt valleys, orientation | `test_charge_misc` | Covariance under rotation; E_kin and \|v\| closed forms | 1e-12; 1e-6 | validated |
| Intervalley rate, Quadratic and Linear | `test_charge_misc` | Formula check | 1e-12 | reference parameterisation (empirical fits) |
| Intervalley rate, `IVRateModel::Matrix` | `test_guards`, `test_charge_misc` (reports the factor vs. Jacoboni–Reggiani) | none | n/a | experimental |
| Release energy at collection, E_rel = E_kin + (E_c − μ_F) for e⁻, E_kin + (μ_F − E_v) for h⁺ | `test_charge_release` | An e⁻ + h⁺ pair at rest releases E_g; closure | 1e-4; 1e-12 | validated (μ_F is an input, default mid-gap) |
| `ChargeRelease::BulkDebye` | `test_charge_release` | Phonons at hν_D | 1e-9 | validated, opt-in |
| `ChargeRelease::ElectrodeFilm` (soft surface-recombination phonons) | `test_charge_release` | Closure; convergence in the seed-quantum size (80 vs. 160 Δ) | 1e-12; 2 % | validated numerically, opt-in. The soft spectrum is a model prediction that has not yet been compared with data |
| Trapping and trap ionisation (constant mean free path) | `test_charge_misc` | Formula check | n/a | reference parameterisation |
| Fano sampler, partition closure | `test_fano_partition` | ⟨N⟩ and Var/⟨N⟩ = F for 18 cases | 1 %, 5 % (observed < 1 %) | validated |

## Variance reduction and alternative engines

| Module | Test | Status |
|---|---|---|
| Russian roulette (`DetectorModel::russian_roulette_bounces`) | none (unbiased by construction, w → w/p) | supported |
| Jump-to-ballistic (`DetectorModel::jump_to_ballistic`) | `test_guards` | experimental |
| Next-event estimation (`experimental::evaluate_sparse_sensor_nee`) | `test_guards` | experimental |
| `experimental::WavefrontSoAEngine`, `experimental::WavefrontPipelineEngine` | `test_guards` | experimental |

## Default parameters taken from G4CMP

`materials::AluminumFilm()`, the default sensor film, uses the G4CMP example values for
`film_absorption` (0.20), `subgap_absorption` (0.10), `gap_energy` (173.715 µeV) and
`phonon_lifetime_slope` (0.29). `phonon_lifetime` (242 ps) is Kaplan's τ₀^φ for Al as
used by G4CMP. arXiv:1401.1937 uses 0.26 ns.

For quantitative sensor-response studies, use `SurfaceProperty::amm_interface` and
`FilmLossModel::KaplanRateCompetition` with `subgap_absorption = 0`.
