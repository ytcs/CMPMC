# CMPMC

**Condensed-Matter Physics Monte Carlo**: a header-only C++20 library that simulates
how phonons, charge carriers and superconducting quasiparticles move through
cryogenic crystals. Typical uses are phonon-mediated particle detectors, TES and KID
calorimeters, and superconducting qubit substrates.

CMPMC takes an energy deposit, splits it into electron-hole pairs and phonons, and
transports every excitation through an anisotropic crystal until it is absorbed. The
result is the sensor hits (position, time, energy, channel) plus a full energy
ledger. The physics follows G4CMP and the published literature, but CMPMC runs
without Geant4, is **15–94× faster than G4CMP** on the same workloads using all cores
of the same machine, and scales almost linearly across threads with bit-reproducible
results ([benchmark](#performance)).

- **No dependencies.** The library needs only a C++20 compiler and the standard
  library. It has no Geant4, ROOT or CLHEP requirement.
- **Anisotropic phonon transport.** It solves the Christoffel equation for any of
  the seven crystal systems or an amorphous solid and follows group velocity, so
  phonon focusing comes out on its own. It models isotope scattering and
  anharmonic decay (Tamura; Leman, arXiv:1109.1193), diffuse, specular and
  polynomial surface models, and acoustic- and diffuse-mismatch interface
  transmission.
- **Charge transport.** Electrons move in Herring–Vogt anisotropic valleys and holes
  with an effective mass. The library models Neganov–Trofimov–Luke phonon emission,
  intervalley scattering, trapping and impact ionisation, and supports uniform fields
  or fields interpolated from a 2D or 3D finite-element mesh.
- **Superconducting sensors.** Phonons absorbed in a thin film run through the
  Kaplan quasiparticle cascade, with pair-breaking efficiency, film-loss models and
  an O(1) surrogate table.
- **Energy partition.** Five nuclear-recoil ionisation-yield models (Lindhard,
  Lewin–Smith, Sarkis, IMPACT@TUNL, an empirical k(E)) and a sampler that
  reproduces the Fano factor exactly.
- **Validated.** Each default model has a test that checks it against a closed form,
  a published reference or a conservation law. See
  [docs/validation.md](docs/validation.md).

## Performance

Measured against G4CMP V10-06-00 on Geant4 11.2.2. Both programs run the detector of
G4CMP's `examples/phonon`: a Ge cylinder of radius 3.81 cm and half-height 1.27 cm with
Al Kaplan sensor films on both faces. The machine was a 4-core Intel Xeon at 2.10 GHz.

| Workload | Single thread: G4CMP | Single thread: CMPMC | Speed-up | All 4 cores: G4CMP | All 4 cores: CMPMC | Speed-up |
|---|---|---|---|---|---|---|
| 7.5 meV phonon, tracked to absorption | 9.85 ms/event | 0.58 ms/event | **17×** | 409 events/s | 6,190 events/s | **15×** |
| e⁻/h⁺ pair drifted at 4 V, no Luke-phonon tracking | 402 ms/event | 5.2 ms/event | **78×** | 9.8 events/s | 915 events/s | **94×** |
| e⁻/h⁺ pair at 4 V with all ~3,800 Luke phonons tracked | 4.0 s/event | 0.14 s/event | **28×** | 0.98 events/s | 29.1 events/s | **30×** |

On all cores, CMPMC runs one process with 4 threads (`cmpmc::run_events`) and
G4CMP runs 4 independent processes, which is G4CMP's ideal-scaling best case. CMPMC
scales 3.9–4.05× from 1 to 4 threads.

The two codes agree on the physics they compute. Deposited energy matches to 0.4 %
and the Luke phonon count per pair to 1.8 %. Most of the gain comes from cheaper
steps: CMPMC has no Geant4 stepping or navigation layer and samples from lookup
tables rather than with rejection loops. Its charge step control also takes about 3×
fewer steps per carrier.

These are CPU figures. CMPMC does not yet have a GPU backend; a CUDA backend is
planned in [issue #1](https://github.com/ytcs/CMPMC/issues/1). The setup, the
alignment of the two configurations, and the reproduction steps are in
[benchmarks/README.md](benchmarks/README.md).

## Quick start

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build          # validation suite
./build/examples/quickstart
```

```cpp
#include <cmpmc/cmpmc.hpp>

int main() {
    using namespace cmpmc;
    const CrystalLattice lattice(materials::Germanium());
    DetectorModel detector;                       // Ge cylinder, Al sensor films top and bottom
    const SimulationEngine engine(lattice, detector);

    Rng rng(1);
    std::vector<PhononTrack> phonons;
    for (int i = 0; i < 200; ++i)
        phonons.push_back(engine.make_track({0, 0, 0}, rng.isotropic_dir(),
                                            PhononMode::L, lattice.debye_energy()));

    std::vector<SensorHit> hits;
    SimStats stats;
    engine.transport(phonons, hits, stats, rng);
}
```

### Running many events on all cores

`cmpmc::run_events` distributes events over threads. Results are bit-identical for
any thread count, because event *i* always draws from the random stream
`Rng::for_stream(seed, i)`.

```cpp
RunOptions options;
options.threads = 0;          // all hardware threads
options.seed    = 42;
const RunResult run = run_events(engine, 100000, [&](EventInput& ev) {
    ev.phonons.push_back(engine.make_track({0, 0, 0}, ev.rng.isotropic_dir(),
                                           PhononMode::L, lattice.debye_energy()));
    ev.injected_eV += lattice.debye_energy();
}, options);
// run.stats: summed statistics; run.hits with run.hit_offsets: the hits of each event
```

### Using CMPMC from another CMake project

```cmake
find_package(cmpmc 0.1 REQUIRED)        # after `cmake --install build --prefix <prefix>`
target_link_libraries(my_app PRIVATE cmpmc::cmpmc)
```

or, vendored:

```cmake
add_subdirectory(third_party/cmpmc)
target_link_libraries(my_app PRIVATE cmpmc::cmpmc)
```

### Build options

| Option | Default | Effect |
|---|---|---|
| `CMPMC_BUILD_TESTS` | ON when top-level | Build the test suite (`ctest`) |
| `CMPMC_BUILD_EXAMPLES` | ON when top-level | Build `examples/` |
| `CMPMC_INSTALL` | ON when top-level | Generate install rules and the `cmpmc` package config |
| `CMPMC_WARNINGS_AS_ERRORS` | OFF | `-Werror` / `/WX` for tests and examples |
| `CMPMC_NATIVE_ARCH` | OFF | `-march=native` for tests and examples |
| `CMPMC_BUILD_BENCHMARKS` | OFF | Build `cmpmc_bench`, the CMPMC side of the G4CMP comparison |
| `CMPMC_G4CMP_CRYSTALMAPS_DIR` | empty | Optional: a G4CMP `CrystalMaps` directory used to cross-check the presets |

Requirements: CMake ≥ 3.21 and a C++20 compiler. Continuous integration builds with GCC and Clang on Linux.

## Documentation

| Document | Contents |
|---|---|
| [docs/architecture.md](docs/architecture.md) | Modules, data flow, track types, threading, the experimental-module guard |
| [docs/physics-models.md](docs/physics-models.md) | Models, equations and references |
| [docs/validation.md](docs/validation.md) | What each model is checked against, and to what tolerance |
| [docs/materials-and-inputs.md](docs/materials-and-inputs.md) | Material presets, `CrystalConfig`, G4CMP `config.txt` import, mesh files |
| [docs/extending.md](docs/extending.md) | Adding materials, changing rates, defining new excitation species |
| [docs/experimental.md](docs/experimental.md) | Modules that are not validated, and how to opt in to them |
| [docs/limitations.md](docs/limitations.md) | Known limitations and planned work |
| [benchmarks/README.md](benchmarks/README.md) | Speed comparison with G4CMP: method, results, how to reproduce |
| [CHANGELOG.md](CHANGELOG.md) | Release notes |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Development workflow and the rules for physics changes |

## Status

CMPMC is pre-1.0 (`0.x`). The API may change between minor versions, and each change
is recorded in the changelog. Modules under `cmpmc/experimental/`, plus the few engine
features listed in [docs/experimental.md](docs/experimental.md), throw
`cmpmc::UnvalidatedModuleError` unless you opt in with
`DetectorModel::allow_unvalidated`.

## Acknowledgements

Several models, default parameters and the material data format come from
[G4CMP](https://github.com/kelseymh/G4CMP), the condensed-matter extension to
Geant4. CMPMC can import G4CMP `config.txt` lattice files directly
(`cmpmc::io::import_g4cmp_config`). The physics references for each model are listed in
[docs/physics-models.md](docs/physics-models.md).

## License

CMPMC is licensed under the **GNU General Public License v3.0 or later**
(`GPL-3.0-or-later`); see [LICENSE](LICENSE). Parts of CMPMC are derived from
G4CMP, which is distributed under the same license; see [NOTICE](NOTICE).

The license allows commercial use and redistribution. If you distribute software
that includes or links CMPMC, you must release it under GPL-compatible terms
and provide its source code.
