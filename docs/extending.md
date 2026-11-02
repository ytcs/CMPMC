# Extending CMPMC

The physics lives in small, pure functions and plain data structs, so you can change
or add a model without touching the transport loop, the lookup tables or the
geometry code. Every physics change needs a reference and a test (see
[CONTRIBUTING.md](../CONTRIBUTING.md)).

## Add a material

Return a filled `CrystalConfig`. Start from the closest preset and override the fields
that differ, always with explicit units:

```cpp
#include <cmpmc/materials/presets.hpp>

inline cmpmc::CrystalConfig Diamond() {
    using namespace cmpmc;
    CrystalConfig c = materials::Silicon();   // same crystal system and valley layout
    c.name      = "C";
    c.atomic_Z  = 6.0;
    c.atomic_A  = 12.011;
    c.density   = 3515.0;                     // kg/m^3
    c.C11       = 1079.0 * units::GPa;
    c.C12       =  124.0 * units::GPa;
    c.C44       =  578.0 * units::GPa;
    // ... B_iso, A_anh, DOS fractions, sound speeds, band parameters ...
    return c;
}

const cmpmc::CrystalLattice lattice(Diamond());
```

Constructing `CrystalLattice` builds the kinematics tables and samplers from the
config. Nothing else needs to change. For a non-cubic crystal, set `group` and the
additional constants (`C13`, `C33`, `C66`, …), or call `set_stiffness(i, j, value)` with
1-based Voigt indices. You can also import an existing G4CMP `config.txt` (see
[materials-and-inputs.md](materials-and-inputs.md)).

Add a check to `tests/test_presets.cpp` for every constant you take from a
publication.

## Change a phonon rate or decay distribution

`phonon/phonon_kernels.hpp` holds the rates and distributions as free functions:

- `isotope_scattering_rate(energy_eV, B_iso)`
- `anharmonic_decay_rate(mode, energy_eV, A_anh)`
- `tamura_lt_pdf(x, cfg)` and `tamura_tt_pdf(x, cfg)`: daughter energy-fraction
  densities
- `decay_angle_*`: daughter angles from momentum conservation

`CrystalLattice::build_decay_cdf_tables()` integrates the two PDFs numerically when the
lattice is built and inverts them into 1024-bin tables. A modified PDF therefore costs
nothing at transport time.

## Change a surface or sensor model

`SurfaceProperty` (`geometry/surface_property.hpp`) is a plain struct per surface. It
holds the absorption probability, the reflection model, an optional `KaplanFilm`
sensor, an optional `IsotropicAMMInterface`, and an optional `ElectrodePatternMask`.
Assign one per entry of `DetectorModel::surfaces`.

## Define a new excitation species

Write a type that satisfies `cmpmc::ExcitationModel` (`transport/excitation_model.hpp`)
and pass it to `SimulationEngine::transport_custom_excitation`:

```cpp
struct Exciton {
    double lifetime_s = 15e-9;
    double mfp_m      = 50e-9;

    double min_energy_eV() const noexcept { return 1.5; }

    double bulk_rate(const cmpmc::GenericExcitationTrack& t) const noexcept {
        return t.velocity / mfp_m + 1.0 / lifetime_s;            // 1/s
    }

    void perform_bulk_interaction(cmpmc::GenericExcitationTrack& t,
                                  std::vector<cmpmc::GenericExcitationTrack>& secondaries,
                                  cmpmc::Rng& rng) const noexcept {
        if (rng.uniform() * bulk_rate(t) < 1.0 / lifetime_s) { t.energy = 0.0; return; }  // decay
        t.mom_dir = rng.isotropic_dir();                                                     // scatter
    }
};
static_assert(cmpmc::ExcitationModel<Exciton>);
```

The engine handles free flight, boundary intersection, sensor absorption and energy
accounting. Push new tracks onto `secondaries` to create daughters. Because the model is
a template parameter, its methods are inlined into the loop. For a complete program,
see `examples/custom_excitation.cpp`.

## Promote an experimental module

To move a module out of `experimental/`:

1. Write the reference comparison described in [limitations.md](limitations.md) as a
   test.
2. Remove its `require_validated` call and update `test_guards`.
3. Move the header, add it to `cmpmc/cmpmc.hpp`, and record the change in
   `docs/validation.md` and `CHANGELOG.md`.
