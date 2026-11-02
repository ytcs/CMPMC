#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 TC Yu and the CMPMC contributors
"""Run the CMPMC-vs-G4CMP comparison and print a Markdown summary.

Usage:
    run_comparison.py --cmpmc build/cmpmc_bench --g4cmp build-g4cmp-bench/g4cmp_bench \
                      [--seeds 3] [--quick]

The G4CMP binary needs G4LATTICEDATA (G4CMP CrystalMaps) and G4ENSDFSTATEDATA in the
environment, plus LD_LIBRARY_PATH for the Geant4 and G4CMP libraries.
"""
import argparse
import platform
import statistics
import subprocess

# (mode, events) per workload; G4CMP sets the scale (about 15-20 s per run).
WORKLOADS = [("phonon", 2000), ("charge", 50), ("chargefull", 5)]
QUICK = [("phonon", 200), ("charge", 5), ("chargefull", 1)]


def run(cmd):
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    line = next(l for l in out.splitlines() if l.startswith("code="))
    return dict(kv.split("=", 1) for kv in line.split())


def cpu_model():
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def mean_sd(xs):
    return statistics.mean(xs), (statistics.stdev(xs) if len(xs) > 1 else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmpmc", required=True)
    ap.add_argument("--g4cmp", required=True)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--voltage", default="4")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    configs = [
        ("G4CMP", lambda m, n, s: [args.g4cmp, m, str(n), str(s), args.voltage]),
        ("CMPMC (octahedral table)", lambda m, n, s: [args.cmpmc, m, str(n), str(s), args.voltage, "octahedral"]),
        ("CMPMC (bilinear table)", lambda m, n, s: [args.cmpmc, m, str(n), str(s), args.voltage, "bilinear"]),
    ]

    print(f"CPU: {cpu_model()}; single thread; {args.seeds} seeds per entry; field workloads at {args.voltage} V\n")
    print("| Workload | Code | Time per event | Speed-up vs G4CMP | Deposited energy per event (eV) | Luke phonons per event |")
    print("|---|---|---|---|---|---|")
    for mode, events in (QUICK if args.quick else WORKLOADS):
        rows = {}
        for name, cmd in configs:
            if mode != "phonon" and "bilinear" in name:
                continue  # charge transport does not use the phonon table choice beyond Luke phonons
            res = [run(cmd(mode, events, 1000 + s)) for s in range(args.seeds)]
            t = [float(r["us_per_event"]) for r in res]
            e = [float(r["edep_eV"]) / events for r in res]
            luke = [int(r.get("luke_phonons", 0)) / events for r in res]
            rows[name] = (mean_sd(t), mean_sd(e), statistics.mean(luke))
        t_ref = rows["G4CMP"][0][0]
        for name, ((t, tsd), (e, esd), luke) in rows.items():
            unit, scale = ("ms", 1e-3) if t >= 1000 else ("µs", 1.0)
            speed = "1×" if name == "G4CMP" else f"{t_ref / t:.0f}×"
            luke_s = f"{luke:.0f}" if mode != "phonon" else "–"
            print(f"| {mode} ({events} events) | {name} | {t * scale:.3g} ± {tsd * scale:.2g} {unit} | {speed} "
                  f"| {e:.4f} ± {esd:.4f} | {luke_s} |")


if __name__ == "__main__":
    main()
