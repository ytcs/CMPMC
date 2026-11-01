# Physicist's Guide to Reading & Extending CMPMC

CMPMC is structured so that **condensed-matter and detector physicists can read, verify, and modify physical models in isolation**, without needing to touch lookup-table (LUT) indexing, multi-threading, or ray-tracing boilerplate.

---

## 1. Directory Layout at a Glance

```
include/cmpmc/
├── core/
│   ├── units.hpp          # Physical units (eV, meV, ueV, THz, GPa, nm, ps) & constants
│   ├── math.hpp           # Vec3, Mat3, SmallVec
│   ├── rng.hpp            # Xoshiro256** PRNG + isotropic & Lambertian angular samplers
│   └── track.hpp          # ExcitationType, PhononTrack, ChargeTrack, GenericExcitationTrack
│
├── materials/             # <--- 1. ADD NEW CRYSTALS OR FILMS HERE
│   ├── crystal_config.hpp # Datasheet structs (CrystalConfig, SuperconductorFilmConfig)
│   └── presets.hpp        # Germanium(), Silicon(), GalliumArsenide(), AluminumFilm(), TungstenTES()
│
├── physics/               # <--- 2. EDIT RATES, PDFs, OR ADD NEW EXCITATIONS HERE
│   ├── phonon_physics.hpp # Pure equations: Christoffel wave solver, Tamura isotope & anharmonic PDFs
│   └── custom_excitation_example.hpp # C++20 Concept + BogoliubovQP & Magnon model templates
│
├── kaplan.hpp             # Superconducting thin-film Cooper-pair breaking & QP cascade
├── charge.hpp             # Herring-Vogt valley drift, Luke Cherenkov emission, IV scattering
├── partition.hpp          # Lindhard nuclear recoil ionization yield & Fano e-/h+ pair sampling
├── geometry.hpp           # Analytic CylinderSolid, BoxSolid & frequency-dependent SurfaceProperty
└── engine.hpp             # High-speed event-based transport loop (Phonons, Charges, Custom Excitations)
```

---

## 2. Common Physics Tasks (Cookbook)

### A. How to Add a New Crystal Material
Open [`include/cmpmc/materials/presets.hpp`](include/cmpmc/materials/presets.hpp) and return a filled `CrystalConfig` struct using explicit units:

```cpp
inline CrystalConfig Sapphire() {
    CrystalConfig c = Germanium();
    c.name       = "Al2O3";
    c.density    = 3980.0; // kg/m^3
    c.C11        = 497.0 * units::GPa;
    c.C12        = 163.0 * units::GPa;
    c.C44        = 147.0 * units::GPa;
    c.B_iso      = 1.25e-43;
    c.A_anh      = 1.10e-56;
    c.debye_freq = 13.0 * units::THz;
    return c;
}
```
When you pass `CrystalLattice(materials::Sapphire())` at startup, the engine **automatically** solves the Christoffel wave equation across all solid angles and builds the 2D group-velocity LUT and 1D anharmonic decay inverse-CDF tables.

---

### B. How to Modify Phonon Scattering Rates or Decay Kinematics
Open [`include/cmpmc/physics/phonon_physics.hpp`](include/cmpmc/physics/phonon_physics.hpp). All equations are pure functions with zero engine boilerplate:
- **Isotope scattering rate:** `isotope_scattering_rate(double energy_eV, double B_iso)`
- **Anharmonic decay rate:** `anharmonic_decay_rate(PhononMode mode, double energy_eV, double A_anh)`
- **Anharmonic energy partition PDFs:** `tamura_lt_pdf(double x, const CrystalConfig& cfg)` and `tamura_tt_pdf(double x, const CrystalConfig& cfg)`

If you modify `tamura_lt_pdf` or `tamura_tt_pdf`, `CrystalLattice::build_decay_cdf_tables()` automatically integrates your new PDF at initialization and compiles it into an O(1) inverse-CDF lookup table—so you can write clear mathematical expressions (including `std::sin`, `std::cos`, `std::pow`) in the PDF without slowing down the runtime Monte Carlo loop!

---

### C. How to Add a Brand New Excitation Species (e.g. Excitons, Magnons, Bogoliubov QPs)
Open [`include/cmpmc/physics/custom_excitation_example.hpp`](include/cmpmc/physics/custom_excitation_example.hpp).
Any struct that defines three simple methods satisfies `ExcitationProcessConcept` and can be simulated at full speed via `engine.transport_custom_excitation(my_model, stack, hits, stats, rng)`:

```cpp
struct SingletExcitonModel {
    double min_energy       = 1.5 * units::eV;
    double radiative_tau_s  = 15.0 * units::ns;
    double diffusion_mfp_m  = 50.0 * units::nm;

    constexpr double min_energy_eV() const noexcept { return min_energy; }

    // 1. Total interaction rate (1/s)
    inline double bulk_rate(const GenericExcitationTrack& trk) const noexcept {
        return (trk.velocity / diffusion_mfp_m) + (1.0 / radiative_tau_s);
    }

    // 2. State update on bulk interaction
    inline void perform_bulk_interaction(GenericExcitationTrack& trk,
                                         std::vector<GenericExcitationTrack>& secondaries,
                                         Rng& rng) const noexcept {
        const double r_rad = 1.0 / radiative_tau_s;
        if (rng.uniform() * bulk_rate(trk) < r_rad) {
            // Radiative recombination -> terminate exciton
            trk.energy = 0.0;
            return;
        }
        // Diffusive scatter
        trk.mom_dir = rng.isotropic_dir();
    }
};
```
Because `transport_custom_excitation` uses C++20 `concept` templates instead of `virtual` base classes, the compiler inlines your physics methods directly into the ray-tracing loop with **zero runtime overhead**.
