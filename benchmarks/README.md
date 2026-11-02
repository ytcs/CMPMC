# CMPMC vs. G4CMP speed comparison

This directory holds two small programs that run the same workloads, one on CMPMC and
one on G4CMP. Neither program writes to disk, and both time only the transport
(initialisation is reported separately, at about 0.08 s for both). Results are given
for a single thread and for all cores.

| File | Contents |
|---|---|
| `cmpmc_bench.cpp` | CMPMC side. Built with `-DCMPMC_BUILD_BENCHMARKS=ON`. |
| `g4cmp/g4cmp_bench.cc`, `g4cmp/CMakeLists.txt` | G4CMP side. A standalone CMake project that needs Geant4 and G4CMP. |
| `run_comparison.py` | Runs both programs (single thread with several seeds, then all cores) and prints the results tables. |

## Workloads

The detector is the one in G4CMP's `examples/phonon`:
- a germanium cylinder of radius 3.81 cm and half-height 1.27 cm, with [100] along z;
- aluminium Kaplan films on both faces, with the example's film parameters;
- the example's polynomial phonon-reflection model on every surface;
- the G4CMP default limit of 100 surface bounces per phonon.

| Mode | One event | Field |
|---|---|---|
| `phonon` | One 7.5 meV phonon from the centre, isotropic, with its mode drawn from the Ge density of states. It is tracked through decay, scattering and reflection until it is absorbed in a film or killed at the bounce limit. | none |
| `charge` | One e⁻/h⁺ pair at the centre, drifted to the electrodes. Luke phonons are not produced. | 4 V across 2.54 cm, uniform |
| `chargefull` | Like `charge`, but every Luke phonon is produced and tracked to absorption. | 4 V across 2.54 cm, uniform |

Several settings were aligned so that both codes do the same work:

- **Sensor hits.** A phonon that reaches a sensor face enters the Kaplan film with
  probability 0.2 (`filmAbsorption`) and otherwise reflects. The stock G4CMP example
  also has a separate absorption without any film physics (`pAbsProb = 0.3`), which
  CMPMC's sensor model does not have. The benchmark therefore sets it to 0.
- **Charges at surfaces.** Charges are absorbed at every surface. This is the default
  in both codes; in G4CMP it follows from a charge-bounce limit of 0.
- **Recombination phonons.** Charge absorption creates none. This is CMPMC's default
  (`ChargeRelease::RecordOnly`); G4CMP is configured with `SetGenPhonons(0)`.
- **Kinematics.** G4CMP uses its (θ, φ) kinematics lookup table, which is its default.
  CMPMC is measured both with its default octahedral table and with its bilinear
  (θ, φ) table, which is the closer analogue.
- **Material constants.** G4CMP reads Ge from `CrystalMaps`. CMPMC's Ge preset uses the
  same phonon constants; `test_g4cmp_import` checks them against `CrystalMaps`.

## Results

Measured on 2026-11-02 on a 4-core Intel Xeon at 2.10 GHz, with GCC 13.3 at `-O3` for
both codes. The versions were Geant4 11.2.2, G4CMP `09b04bc` (release V10-06-00) and
CMPMC 0.1.0.

### Single thread (mean ± standard deviation over 3 seeds)

| Workload | Code | Time per event | Speed-up | Deposited energy per event | Luke phonons per pair |
|---|---|---|---|---|---|
| phonon (2000 events) | G4CMP | 9.85 ± 0.5 ms | 1× | 6.7 meV | – |
| | CMPMC, octahedral table | 0.581 ± 0.019 ms | **17×** | 6.7 meV | – |
| | CMPMC, bilinear table | 0.842 ± 0.014 ms | **12×** | 6.7 meV | – |
| charge (50 events) | G4CMP | 402 ± 28 ms | 1× | see note | – |
| | CMPMC | 5.15 ± 0.52 ms | **78×** | see note | 3769 (emission rate only) |
| chargefull (5 events) | G4CMP | 4.0 ± 0.25 s | 1× | 4.304 ± 0.006 eV | 3838 |
| | CMPMC | 0.142 ± 0.016 s | **28×** | 4.285 ± 0.002 eV | 3773 |

### All 4 cores (throughput)

- **CMPMC** runs one process with 4 threads through `cmpmc::run_events`.
- **G4CMP** runs 4 independent processes at once, each simulating a quarter of the
  events with its own seed. This is G4CMP's best case: it has no threading overhead and
  scales ideally (measured 4.0×), so it is the conservative comparison for CMPMC.

| Workload | G4CMP (4 processes) | CMPMC (4 threads) | Speed-up | CMPMC scaling vs 1 thread |
|---|---|---|---|---|
| phonon | 409 events/s | 6,190 events/s | **15×** | 3.89× |
| charge | 9.76 events/s | 915 events/s | **94×** | 3.89× |
| chargefull | 0.98 events/s | 29.1 events/s | **30×** | 4.05× |

Single-thread speed-ups vary by a few percent between sessions (17× here, 16× in an
earlier run of the same harness), which reflects normal timing noise on a shared machine.

**Physics agreement.** Both codes deposit the same energy in the phonon workload, to
within the printed precision. In `chargefull`, deposited energy agrees to 0.4 % and the
number of Luke phonons per pair to 1.8 %. That residual comes from the codes'
different charge time-stepping schemes.

**Note on `charge`.** When Luke phonons are not produced, G4CMP deposits their energy
locally (4.74 eV per pair, which is the 4 eV of field work plus the 0.74 eV band gap).
CMPMC records only the release energy at collection (0.76 eV per pair). This is a
bookkeeping difference in a mode that discards phonons. It does not affect the
timing, and both codes agree when the phonons are tracked (`chargefull`).

**Where the speed-up comes from.**
- *Phonons.* CMPMC takes about 1.4× fewer steps per phonon (about 3.6k vs. 5.0k);
  Geant4 records additional steps at boundaries. Each step is also about 11× cheaper
  (0.16 µs vs. 1.8 µs): there is no Geant4 stepping or navigation layer, and the
  samplers use lookup tables instead of rejection loops.
- *Charges.* CMPMC's Δk-capped step control takes about 2.8× fewer steps per carrier
  (about 20k vs. 56k). Each step is about 28× cheaper (0.12 µs vs. 3.3 µs).
- *Overall.* Tracking Luke phonons dominates `chargefull`, so its speed-up falls between
  the phonon and charge figures.

**Multi-core.** Both codes parallelise over events and scale almost linearly on this
machine, so the single-thread ratios carry over to all cores. CMPMC's `run_events`
gives bit-identical results for any thread count; each event draws from its own
random stream.

## Reproducing

```bash
# 1. Geant4 (no data sets needed beyond G4ENSDFSTATE) and G4CMP
cmake -S geant4 -B g4-build -DCMAKE_BUILD_TYPE=Release -DGEANT4_BUILD_MULTITHREADED=OFF \
      -DGEANT4_INSTALL_DATA=OFF -DCMAKE_INSTALL_PREFIX=$PWD/g4 && cmake --build g4-build -j && cmake --install g4-build
cmake -S G4CMP -B g4cmp-build -DCMAKE_BUILD_TYPE=Release -DGeant4_DIR=$PWD/g4/lib/cmake/Geant4 \
      -DCMAKE_INSTALL_PREFIX=$PWD/g4cmp && cmake --build g4cmp-build -j && cmake --install g4cmp-build

# 2. The two benchmark programs
cmake -S benchmarks/g4cmp -B g4bench-build -DCMAKE_BUILD_TYPE=Release \
      -DGeant4_DIR=$PWD/g4/lib/cmake/Geant4 -DG4CMP_DIR=$PWD/g4cmp && cmake --build g4bench-build
cmake -B cmpmc-build -DCMPMC_BUILD_BENCHMARKS=ON && cmake --build cmpmc-build --target cmpmc_bench

# 3. Run (single thread, then all cores; --cores defaults to the machine's core count)
export G4LATTICEDATA=$PWD/G4CMP/CrystalMaps G4ENSDFSTATEDATA=<path>/G4ENSDFSTATE2.3
export LD_LIBRARY_PATH=$PWD/g4/lib:$PWD/g4cmp/lib
python3 benchmarks/run_comparison.py --cmpmc cmpmc-build/cmpmc_bench --g4cmp g4bench-build/g4cmp_bench
```
