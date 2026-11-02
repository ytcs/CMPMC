# Architecture

CMPMC is header-only and lives under `include/cmpmc/`. Code is in `namespace cmpmc`,
with pure-physics kernels in `cmpmc::physics`, material data in `cmpmc::materials`,
file import in `cmpmc::io`, field and electrode geometry in `cmpmc::geometry`, and
unvalidated modules in `cmpmc::experimental`.

## Module map

```
include/cmpmc/
├── cmpmc.hpp                 umbrella header (supported API only)
├── version.hpp               CMPMC_VERSION_* macros, cmpmc::version
├── core/                     units & constants, Vec3/Mat3, Rng (xoshiro256**),
│                             track records, UnvalidatedModuleError
├── materials/                CrystalConfig, SuperconductorFilmConfig, built-in presets
├── io/                       import of G4CMP config.txt lattice files
├── lattice/                  CrystalLattice: precomputed kinematics and samplers per crystal
├── phonon/                   Christoffel solver, scattering/decay rates and PDFs,
│                             AMM/DMM interface transmission
├── charge/                   Herring–Vogt drift, Luke emission, intervalley scattering,
│                             trapping, NIEL yield models, energy partition, hit merging
├── superconductor/           Kaplan thin-film cascade and surrogate table,
│                             bulk BCS quasiparticle kernels
├── geometry/                 analytic solids, surface properties, FEM field meshes,
│                             patterned electrode masks
├── transport/                SimulationEngine and DetectorModel, the ExcitationModel
│                             concept, Russian roulette
└── experimental/             SoA and wavefront engines, jump-to-ballistic,
                              next-event estimation (not in cmpmc.hpp)
```

The dependencies run in one direction:
`core → materials → phonon / charge / superconductor → lattice → geometry → transport`.
`experimental/` may depend on anything. Nothing outside it depends on it, except the
`DetectorModel::jump_to_ballistic` option (see [experimental.md](experimental.md)).

## Data flow of an event

```
energy deposit ──► EnergyPartition ──► e⁻/h⁺ pairs (ChargeTrack)
                                   └─► prompt phonons (PhononTrack)

ChargeTrack ──► SimulationEngine::transport_charges ──► Luke phonons, collection hits,
                                                        release energy
PhononTrack ──► SimulationEngine::transport ──► scattering / decay / surface reflection
                                            └─► sensor film (KaplanFilm) ──► SensorHit
```

1. **Material setup.** A `CrystalConfig` (a preset, user-filled, or imported) is turned
   into a `CrystalLattice`. The constructor solves the Christoffel equation over the
   sphere and builds the octahedral group-velocity table, the decay inverse-CDF tables
   and the charge-transport constants. Build it once and share it read-only.
2. **Detector setup.** `DetectorModel` holds the solid, the electric field (uniform or
   FEM mesh), one `SurfaceProperty` per surface (reflection model, optional Kaplan
   sensor film, optional AMM interface, optional electrode mask), cut-offs and the
   transport options.
3. **Transport.** `SimulationEngine` is a const object holding references to the
   lattice and detector. Each `transport*` call consumes a stack of tracks, appends
   `SensorHit`s, and accumulates `SimStats`. Secondaries are pushed onto the same stack,
   so a call returns only when every excitation has been absorbed, has escaped, or has
   fallen below the cut-off.

## Track records

`PhononTrack`, `ChargeTrack` and `GenericExcitationTrack` are flat, 64-byte-aligned
structs. A phonon stores its wavevector direction and its group-velocity direction
separately, because the two are not parallel in an anisotropic crystal. A charge
carrier stores its valley index and its crystal momentum. Each track carries a
statistical weight, so the engine can sample Luke phonons and apply Russian roulette
without bias.

## Energy accounting

`SimStats` holds `total_injected_eV` (the caller sets it), `total_deposited_eV` (sensor
hits) and `escaped_or_killed_eV` (cut-offs, bounce limits, bulk losses). The engines
conserve weighted energy: injected = deposited + escaped/killed. The test suite checks
this to 1e-9 for phonon runs and 1e-12 for film cascades.

## Threading

The engine, the lattice and the detector are immutable during transport. To run in
parallel, give each thread its own `Rng`, track stack, hit vector and `SimStats`, and
merge them afterwards. `SimStats` has `operator+=` for this. `examples/phonon_benchmark.cpp`
shows the pattern with `std::thread`. The FEM mesh field keeps a `thread_local`
cell-location cache, so concurrent lookups are safe.

## Performance design

- **Kinematics lookup.** The default mode maps k̂ to (v̂_g, |v_g|) through an octahedral
  256 × 256 FP32 table, with no trigonometric functions on the hot path. The exact
  per-call Christoffel solve is available as `KinematicsMode::Exact`.
- **Samplers.** Anharmonic-decay energy fractions use 1024-bin inverse-CDF tables.
  The Kaplan cascade samplers use variable substitutions that remove the 1/√(x² − Δ²)
  endpoint singularities, so no rejection loop is needed. The Luke emission angle is
  drawn by an analytic inverse CDF.
- **Geometry.** Analytic ray–solid intersection returns the exit distance, the inward
  normal and the surface id in one call. Reflections update the track in place.
- **Extension point.** `transport_custom_excitation<Model>` takes any type that
  satisfies the `ExcitationModel` concept. The model is a template parameter, so its
  methods are inlined into the loop.

## Guarding unvalidated code

`core/validation.hpp` provides `require_validated(allow, module, reason)`. Every entry
point without a reference comparison calls it and throws `UnvalidatedModuleError`
unless the caller opts in, through `DetectorModel::allow_unvalidated` or the explicit
argument of a free function. `test_guards` checks that every guarded path throws by
default and runs when opted in. It also checks that validated paths never throw.
