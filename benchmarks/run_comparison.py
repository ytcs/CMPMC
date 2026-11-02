#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 TC Yu and the CMPMC contributors
"""Run the CMPMC-vs-G4CMP comparison and print Markdown tables.

Usage:
    run_comparison.py --cmpmc build/cmpmc_bench --g4cmp build-g4cmp-bench/g4cmp_bench \
                      [--seeds 3] [--cores N] [--quick]

Two comparisons are made:
  1. single thread: one process per code, time per event (mean over seeds);
  2. all cores: CMPMC with N threads (cmpmc::run_events) against N concurrent
     G4CMP processes, each simulating 1/N of the events with its own seed. The
     G4CMP figure is the ideal-scaling throughput of independent jobs, which is
     an upper bound for G4CMP multithreading on the same machine.

The G4CMP binary needs G4LATTICEDATA (G4CMP CrystalMaps) and G4ENSDFSTATEDATA in the
environment, plus LD_LIBRARY_PATH for the Geant4 and G4CMP libraries.
"""
import argparse
import os
import platform
import statistics
import subprocess
import time

# (mode, events): G4CMP sets the scale, about 15-20 s per single-thread run.
WORKLOADS = [("phonon", 2000), ("charge", 50), ("chargefull", 5)]
QUICK = [("phonon", 200), ("charge", 4), ("chargefull", 1)]


def parse(out):
    line = next(l for l in out.splitlines() if l.startswith("code="))
    return dict(kv.split("=", 1) for kv in line.split())


def run(cmd):
    return parse(subprocess.run(cmd, check=True, capture_output=True, text=True).stdout)


def run_concurrent(cmds):
    """Start all commands at once; return (wall seconds, parsed outputs)."""
    t0 = time.perf_counter()
    procs = [subprocess.Popen(c, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True) for c in cmds]
    outs = [p.communicate()[0] for p in procs]
    wall = time.perf_counter() - t0
    if any(p.returncode for p in procs):
        raise RuntimeError("a concurrent run failed")
    return wall, [parse(o) for o in outs]


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


def fmt_time(us, sd):
    if us >= 1e6:
        return f"{us / 1e6:.3g} ± {sd / 1e6:.2g} s"
    if us >= 1e3:
        return f"{us / 1e3:.3g} ± {sd / 1e3:.2g} ms"
    return f"{us:.3g} ± {sd:.2g} µs"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmpmc", required=True)
    ap.add_argument("--g4cmp", required=True)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--cores", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--voltage", default="4")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    workloads = QUICK if args.quick else WORKLOADS

    def g4(mode, n, seed):
        return [args.g4cmp, mode, str(n), str(seed), args.voltage]

    def cm(mode, n, seed, kin="octahedral", threads=1):
        return [args.cmpmc, mode, str(n), str(seed), args.voltage, kin, str(threads)]

    print(f"CPU: {cpu_model()}, {args.cores} cores; field workloads at {args.voltage} V\n")

    # ---- 1. single thread ------------------------------------------------------
    print(f"### Single thread ({args.seeds} seeds per entry)\n")
    print("| Workload | Code | Time per event | Speed-up vs G4CMP | Deposited energy per event (eV) | Luke phonons per pair |")
    print("|---|---|---|---|---|---|")
    configs = [("G4CMP", g4),
               ("CMPMC (octahedral table)", lambda m, n, s: cm(m, n, s, "octahedral")),
               ("CMPMC (bilinear table)", lambda m, n, s: cm(m, n, s, "bilinear"))]
    for mode, events in workloads:
        rows = {}
        for name, cmd in configs:
            if mode != "phonon" and "bilinear" in name:
                continue
            res = [run(cmd(mode, events, 1000 + s)) for s in range(args.seeds)]
            rows[name] = (mean_sd([float(r["us_per_event"]) for r in res]),
                          mean_sd([float(r["edep_eV"]) / events for r in res]),
                          statistics.mean(int(r.get("luke_phonons", 0)) / events for r in res))
        t_ref = rows["G4CMP"][0][0]
        for name, ((t, tsd), (e, esd), luke) in rows.items():
            speed = "1×" if name == "G4CMP" else f"{t_ref / t:.0f}×"
            luke_s = f"{luke:.0f}" if mode != "phonon" else "–"
            print(f"| {mode} ({events} events) | {name} | {fmt_time(t, tsd)} | {speed} | {e:.4f} ± {esd:.4f} | {luke_s} |")

    # ---- 2. all cores ----------------------------------------------------------
    n = args.cores
    print(f"\n### All {n} cores (throughput)\n")
    print("| Workload | Code | Events/s | Speed-up vs G4CMP | CMPMC scaling (vs 1 thread) |")
    print("|---|---|---|---|---|")
    for mode, events in workloads:
        g_events = max(n, events * n // 2)              # enough work per process
        per = g_events // n
        wall, _ = run_concurrent([g4(mode, per, 2000 + i) for i in range(n)])
        g4_rate = per * n / wall
        c_events = g_events * 8                          # CMPMC is faster; keep runs > 1 s
        c1 = float(run(cm(mode, c_events // 4 or 1, 3000))["events_per_s"])
        cn = float(run(cm(mode, c_events, 3000, threads=n))["events_per_s"])
        print(f"| {mode} | G4CMP ({n} processes) | {g4_rate:.3g} | 1× | – |")
        print(f"| {mode} | CMPMC ({n} threads) | {cn:.3g} | {cn / g4_rate:.0f}× | {cn / c1:.2f}× |")


if __name__ == "__main__":
    main()
