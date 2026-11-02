# CMPMC vs. G4CMP speed comparison

This directory holds two small programs that run the same workloads, one on CMPMC and
one on G4CMP. Both programs are single-threaded, write nothing to disk, and time only
the transport (initialisation is reported separately, at about 0.08 s for both).

| File | Contents |
|---|---|
| `cmpmc_bench.cpp` | CMPMC side. Built with `-DCMPMC_BUILD_BENCHMARKS=ON`. |
| `g4cmp/g4cmp_bench.cc`, `g4cmp/CMakeLists.txt` | G4CMP side. A standalone CMake project that needs Geant4 and G4CMP. |
| `run_comparison.py` | Runs both programs with several seeds and prints the results table. |

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

Measured on 2026-11-02 on an Intel Xeon at 2.10 GHz, single thread, with GCC 13.3 at
`-O3` for both codes. The versions were Geant4 11.2.2, G4CMP `09b04bc` (release
V10-06-00) and CMPMC 0.1.0. Each entry is the mean ± standard deviation over 3 seeds.

| Workload | Code | Time per event | Speed-up | Deposited energy per event | Luke phonons per pair |
|---|---|---|---|---|---|
| phonon (2000 events) | G4CMP | 9.03 ± 0.6 ms | 1× | 6.7 meV | – |
| | CMPMC, octahedral table | 0.581 ± 0.012 ms | **16×** | 6.7 meV | – |
| | CMPMC, bilinear table | 0.819 ± 0.016 ms | **11×** | 6.7 meV | – |
| charge (50 events) | G4CMP | 367 ± 5 ms | 1× | see note | – |
| | CMPMC | 4.78 ± 0.2 ms | **77×** | see note | 3777 (emission rate only) |
| chargefull (5 events) | G4CMP | 3.75 ± 0.19 s | 1× | 4.304 ± 0.006 eV | 3838 |
| | CMPMC | 0.164 ± 0.032 s | **23×** | 4.287 ± 0.003 eV | 3770 |

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

These are single-thread numbers. Both codes run in parallel over events (G4CMP
through Geant4's multithreading, CMPMC through one engine per thread), so the
ratio should carry over to a multi-core CPU.

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

# 3. Run
export G4LATTICEDATA=$PWD/G4CMP/CrystalMaps G4ENSDFSTATEDATA=<path>/G4ENSDFSTATE2.3
export LD_LIBRARY_PATH=$PWD/g4/lib:$PWD/g4cmp/lib
python3 benchmarks/run_comparison.py --cmpmc cmpmc-build/cmpmc_bench --g4cmp g4bench-build/g4cmp_bench
```
