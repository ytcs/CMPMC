# Physics models

This page summarises each model and its sources. For the check behind each model and
its tolerance, see [validation.md](validation.md). For simplifications and planned
extensions, see [limitations.md](limitations.md). Internal units are SI, except
energies, which are in eV. `core/units.hpp` defines the conversion constants.

## Crystal elasticity and phonon kinematics

*`materials/crystal_config.hpp`, `phonon/phonon_kernels.hpp`, `lattice/crystal_lattice.hpp`*

- **Stiffness tensor.** The independent elastic constants of each crystal system
  (amorphous, cubic, tetragonal, orthorhombic, hexagonal, rhombohedral, monoclinic,
  triclinic) are expanded to the full 6 × 6 Voigt matrix. A Miller (hkl) orientation
  rotates the crystal so that the (hkl) normal points along the detector +z axis.
- **Christoffel equation.** For a wavevector direction k̂, the phase velocities are the
  square roots of the eigenvalues of D_il = ρ⁻¹ C_ijlm k̂_j k̂_m. The group velocity
  v_g = ∇_k ω is obtained by contracting the stiffness tensor with the polarisation
  vectors. Two solvers are available: a Jacobi eigen-solve, and a closed-form
  trigonometric (Cardano) solve that falls back to Jacobi near degeneracies.
- **Lookup tables.** The default `KinematicsMode::Octahedral` stores (v̂_g, |v_g|) for the
  three modes on a 256 × 256 octahedral map of the sphere in FP32. Table lookup needs
  no trigonometric functions. `KinematicsMode::Exact` solves the Christoffel equation
  on every call, and `KinematicsMode::Bilinear` uses a 181 × 361 (cos θ, φ) table.

## Phonon scattering and decay

*`phonon/phonon_kernels.hpp`, `transport/engine.hpp`*

- **Isotope scattering.** Rate Γ_iso = B ν⁴ (Tamura). The default treats scattering as
  elastic with the daughter mode sampled from the density of states.
- **Anharmonic decay.** Rate Γ_anh = A ν⁵ for longitudinal phonons, split into
  L → L′ + T and L → T + T with the branching fraction `fTT_frac`. The daughter energy
  fractions follow Leman, arXiv:1109.1193, Eqs. (9) and (12), sampled from inverse-CDF
  tables. The daughter angles conserve crystal momentum exactly (Eqs. (10), (11),
  (13), (14)).
- **Polycrystal grain scattering.** An optional elastic scattering channel with
  constant mean free path `polycrystal_mfp_m`.
- **Free flight.** The engine draws one exponential path length from the total rate.
  If the boundary comes first, it handles the surface; otherwise it selects the
  channel in proportion to its rate.

## Surfaces and interfaces

*`geometry/surface_property.hpp`, `phonon/interface_transmission.hpp`*

- **Default reflection.** Fully diffuse (Lambertian in k̂, re-mapped to v_g) with no
  surface down-conversion. Leman, arXiv:1109.1193, §IV.G, reports that diffuse
  scattering describes phonon–surface reflection well.
- **Polynomial model.** Frequency-dependent specular, diffuse and anharmonic
  probabilities given as polynomials in f (GHz), with cut-offs. The coefficients
  shipped as `demo_polynomial_*` are example values and are guarded as experimental.
- **Acoustic mismatch (AMM).** `IsotropicAMMInterface` solves the elastodynamic
  boundary problem between two isotropic half-spaces: normal incidence, SH, and
  oblique P-SV with mode conversion. It returns energy-flux transmission by mode and
  angle. When it is set on a surface, it replaces the film's constant entry
  probability.
- **Diffuse mismatch (DMM).** t₁₂ = S₂/(S₁ + S₂) with S = v_L⁻² + 2v_T⁻² (Debye model).

## Superconducting sensor films

*`superconductor/kaplan_film.hpp`, `superconductor/bulk_superconductor.hpp`*

- **Kaplan cascade.** A phonon of energy Ω ≥ 2Δ that enters the film breaks a Cooper
  pair. The QP energy follows the Kaplan density [E₁E₂ + Δ²]/√((E₁² − Δ²)(E₂² − Δ²)).
  QPs relax by phonon emission with density (E − E′)²(E′ − Δ²/E)/√(E′² − Δ²) (Kaplan et
  al., PRB 14, 4854 (1976)). Re-emitted phonons either escape back into the crystal
  or break further pairs. Energy is conserved per event.
- **Film loss.** `FilmLossModel::MeanFreePath` (default) uses
  P_esc = exp(−path/λ) with λ = v_s τ₀/(1 + slope·(Ω/Δ − 2)) (arXiv:2302.05998,
  Eqs. (17)–(18)). `FilmLossModel::KaplanRateCompetition` lets the exact T = 0 Kaplan
  pair-breaking rate compete with a film-loss time τ_l = 4d/(η_t v_s) (arXiv:1401.1937).
- **Pair-breaking efficiency.** `KaplanFilm::pair_breaking_efficiency()` returns η_pb and
  the QP Fano factor for a phonon or photon seed.
- **Surrogate table.** `KaplanSurrogateLUT` tabulates complete cascade outcomes,
  including the individual re-emitted phonons, versus incident energy. Sampling an
  outcome costs O(1).
- **Bulk quasiparticles (experimental).** Energy-dependent diffusion
  D(E) = D_n √(1 − (Δ/E)²), plus pair-breaking, phonon radiation, recombination with
  thermal QPs, and local trapping, using the BCS gap Δ(T).

## Charge carriers

*`charge/charge_transport.hpp`, `charge/charge_kernels.hpp`, `transport/engine.hpp`*

- **Band structure.** Electrons occupy Herring–Vogt ellipsoidal valleys with mass
  tensor diag(m_∥, m_⊥, m_⊥) along each valley axis. The axes rotate with the crystal
  orientation. Holes use an isotropic effective mass. Bands are parabolic.
- **Luke (Neganov–Trofimov–Luke) emission.** Above the sound speed, carriers emit
  acoustic phonons inside the Cherenkov cone. The scattering length ℓ₀, the rate
  (v/3ℓ₀)(M − 1)³/M and the angle and energy distributions follow Leman,
  arXiv:1109.1193, Table II and Eq. (22). The emission angle is drawn by an analytic
  inverse CDF.
- **Time stepping.** By default (first-order) each step changes the carrier wavevector
  by at most `charge_max_dk_frac`·max(k_s, |k|), so the rate is effectively constant
  over the step. `ChargeStepping::LemanSecondOrder` lets the rate vary linearly over a
  longer sampling step and draws the event time by inverse CDF (arXiv:1109.1193
  §V.E, Eqs. (26)–(30)).
- **Intervalley scattering.** `IVRateModel::Quadratic` uses
  rate·(E₀² + E²)^(power/2) and `Linear` uses rate₀ + rate₁·E^power. Both are empirical
  fits over a limited field range. `Matrix` is a deformation-potential optical-phonon
  model plus Erginsoy neutral-impurity scattering, and is experimental.
- **Trapping and impact ionisation.** Constant mean free paths for e⁻/h⁺ trapping and
  for the four trap-ionisation channels.
- **Boundaries.** A carrier is absorbed (collected) with probability `prob_charge_abs`,
  subject to an optional minimum perpendicular wavevector. Otherwise it is reflected,
  specularly or diffusely, conserving kinetic energy.
- **Collection energy.** A collected carrier releases
  E_rel = E_kin + (E_c − μ_F) for electrons, or E_kin + (μ_F − E_v) for holes.
  `ChargeRelease` controls where that energy goes: it can be recorded only, emitted as
  bulk Debye phonons, or deposited into the electrode film. In the film case, the
  phonons that escape the film are re-injected as a soft, ballistic population.
- **Fields.** `UniformElectricField`, `geometry::TriLinearMesh3D` (tetrahedral) and
  `geometry::BiLinearMesh2D` (triangular, planar or axisymmetric) return
  barycentric-interpolated fields with cell walking.

## Energy partition

*`charge/energy_partition.hpp`, `charge/niel_models.hpp`*

- **Ionisation yield for nuclear recoils.** Five models are available:
  - Lindhard et al. (1963)
  - Lewin & Smith (1996)
  - Sarkis et al., Phys. Rev. A 107, 062811 (2023): Si, tabulated
  - IMPACT@TUNL, Phys. Rev. Lett. 131, 091801 (2023): Si power law
  - an empirical k(E) with user parameters
- **Pair statistics.** The number of e⁻/h⁺ pairs is drawn from an interpolated binomial
  whose mean is E_ion/ε_eh and whose variance-to-mean ratio is exactly F. The
  remainder becomes prompt phonons at hν_D with the density-of-states mode mix.
  Optional statistical down-sampling assigns weights.
- **Charge cloud.** Initial pairs can be spread over a sphere whose radius scales with N^(1/3) of the pair count, then folded back
  into the solid.
- **Step accumulation.** `StepAccumulator` and `HitMerger` merge the energy deposits of
  a primary track into one partition call per track segment.

## References

- J. Lindhard et al., Mat. Fys. Medd. Dan. Vid. Selsk. 33, 10 (1963).
- J. D. Lewin and P. F. Smith, Astropart. Phys. 6, 87 (1996).
- S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
- C. Jacoboni and L. Reggiani, Rev. Mod. Phys. 55, 645 (1983).
- S. W. Leman, *Invited review: Physics and Monte Carlo techniques as relevant to
  cryogenic, phonon, and ionization readout of CDMS radiation-detectors*,
  arXiv:1109.1193.
- arXiv:1401.1937: pair-breaking efficiency and quasiparticle kinetics in
  superconducting films (Chang–Scalapino / Rothwarf–Taylor treatment).
- arXiv:2302.05998: phonon escape from superconducting films (mean-free-path loss law).
- arXiv:2408.04732: phonon constants (isotope, anharmonic, density of states) for
  Si, Ge and GaAs.
- arXiv:2004.11499: electron–hole pair energy and Fano factor in Si.
- arXiv:1805.02495: measured Si → Al phonon transmission.
- M. Sarkis et al., Phys. Rev. A 107, 062811 (2023).
- IMPACT@TUNL collaboration, Phys. Rev. Lett. 131, 091801 (2023).
- G4CMP, https://github.com/kelseymh/G4CMP: reference implementation for several
  parameterisations and the source of the `config.txt` lattice format.
