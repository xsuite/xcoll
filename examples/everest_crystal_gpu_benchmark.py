# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""CPU-vs-GPU performance benchmark for the Everest *crystal* element.

This is the first time ``EverestCrystal`` can run on the CUDA context at all,
so this script measures the resulting speed-up: it tracks the canonical
``Drift + EverestCrystal + Drift`` line on ``xobjects.ContextCpu`` (single
core) and ``xobjects.ContextCupy`` across a range of bunch sizes and reports
throughput (particles tracked per second) and the GPU/CPU speed-up.

Methodology (fair-timing):
  * kernels are warmed up once per context before timing (the one-off NVRTC /
    cffi compile is excluded);
  * each measurement re-uses a freshly built bunch (build cost excluded from
    the timer) and is repeated; the *median* wall time is reported;
  * GPU timing brackets ``cp.cuda.runtime.deviceSynchronize()`` so it measures
    real device execution, not async launch latency.

Caveats: the CPU baseline is a single core (``ContextCpu`` is serial; OpenMP
would narrow the gap); the Everest physics is branchy + stochastic so the
crystal kernel diverges across GPU threads (less ideal than a plain drift); and
absolute numbers are hardware-specific (here: one consumer GPU).

Run with cupy's CUDA libs on the loader path, e.g.:
    LDP=$(ls -d $(python -c 'import sys;print(sys.prefix)')/lib/python*/site-packages/nvidia/*/lib | tr '\n' ':')
    LD_LIBRARY_PATH="$LDP" python everest_crystal_gpu_benchmark.py
"""

import json
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

from xcoll.scattering_routines.everest import set_crystal_stack_limit

# Canonical crystal config (matches the parity example + the test).
P0C_EV = 400e9
CRYSTAL_LENGTH = 0.002
CRYSTAL_BENDING_ANGLE = 149e-6
CRYSTAL_WIDTH = 0.002
CRYSTAL_HEIGHT = 0.05
CRYSTAL_SIDE = "+"
CRYSTAL_MISCUT = 0.0
CRYSTAL_LATTICE = "strip"
CRYSTAL_JAW = 0.001
DRIFT_LENGTH = 1.0

BEAM_X_MEAN = 1.5e-3
BEAM_X_SPREAD = 75e-6
BEAM_PX_SPREAD = 30e-6
BEAM_Y_SPREAD = 1e-3
BEAM_PY_SPREAD = 5e-6
NUMPY_SEED = 20260603
PARTICLES_SEED = 13579

# Bunch sizes to scan, and timed repeats per size (fewer for the big, slow CPU runs).
SIZES = [1_000, 10_000, 50_000, 100_000, 200_000, 500_000]
REPEATS = {1_000: 7, 10_000: 5, 50_000: 5, 100_000: 3, 200_000: 3, 500_000: 3}


def build_line(context):
    crystal = xc.EverestCrystal(
        length=CRYSTAL_LENGTH, material=xc.materials.Silicon,
        bending_angle=CRYSTAL_BENDING_ANGLE, width=CRYSTAL_WIDTH,
        height=CRYSTAL_HEIGHT, side=CRYSTAL_SIDE, miscut=CRYSTAL_MISCUT,
        lattice=CRYSTAL_LATTICE, jaw=CRYSTAL_JAW, _context=context)
    line = xt.Line(
        elements=[xt.Drift(length=DRIFT_LENGTH), crystal,
                  xt.Drift(length=DRIFT_LENGTH)],
        element_names=["d_upstream", "crystal", "d_downstream"])
    line.particle_ref = xt.Particles(p0c=P0C_EV, mass0=xp.PROTON_MASS_EV)
    line.build_tracker(_context=context)
    return line


def build_bunch(line, context, n_part):
    rng = np.random.default_rng(NUMPY_SEED)
    particles = line.build_particles(
        x=rng.normal(BEAM_X_MEAN, BEAM_X_SPREAD, n_part),
        px=rng.normal(0.0, BEAM_PX_SPREAD, n_part),
        y=rng.normal(0.0, BEAM_Y_SPREAD, n_part),
        py=rng.normal(0.0, BEAM_PY_SPREAD, n_part))
    seeds = PARTICLES_SEED + np.arange(n_part, dtype=np.uint32)
    particles._init_random_number_generator(seeds=seeds)
    return particles


def is_gpu(context):
    return isinstance(context, getattr(xo, "ContextCupy", ()))


def sync(context):
    if is_gpu(context):
        import cupy as cp
        cp.cuda.runtime.deviceSynchronize()


def time_track(line, context, n_part, repeats):
    """Median wall time (s) of one turn through the line for n_part particles."""
    # warm up: triggers kernel build + first-launch costs (excluded from timing)
    warm = build_bunch(line, context, min(n_part, 2000))
    line.track(warm, num_turns=1)
    sync(context)
    samples = []
    for _ in range(repeats):
        particles = build_bunch(line, context, n_part)  # build excluded from timer
        sync(context)
        t0 = time.perf_counter()
        line.track(particles, num_turns=1)
        sync(context)
        samples.append(time.perf_counter() - t0)
    return float(np.median(samples))


def main():
    cpu_ctx = xo.ContextCpu()
    gpu_ctx = xo.ContextCupy()
    eff_stack = set_crystal_stack_limit(gpu_ctx)
    print(f"[setup] GPU per-thread CUDA stack limit = {eff_stack} bytes")

    cpu_line = build_line(cpu_ctx)
    gpu_line = build_line(gpu_ctx)

    rows = []
    print(f"{'N':>9} {'CPU [ms]':>11} {'GPU [ms]':>11} {'speedup':>9} "
          f"{'CPU Mp/s':>9} {'GPU Mp/s':>9}")
    for n_part in SIZES:
        reps = REPEATS[n_part]
        t_cpu = time_track(cpu_line, cpu_ctx, n_part, reps)
        t_gpu = time_track(gpu_line, gpu_ctx, n_part, reps)
        speedup = t_cpu / t_gpu
        cpu_mps = n_part / t_cpu / 1e6
        gpu_mps = n_part / t_gpu / 1e6
        rows.append({"n_part": n_part, "cpu_s": t_cpu, "gpu_s": t_gpu,
                     "speedup": speedup, "cpu_mpps": cpu_mps,
                     "gpu_mpps": gpu_mps})
        print(f"{n_part:>9d} {t_cpu*1e3:>11.2f} {t_gpu*1e3:>11.2f} "
              f"{speedup:>8.2f}x {cpu_mps:>9.3f} {gpu_mps:>9.3f}")

    ns = np.array([r["n_part"] for r in rows])
    sp = np.array([r["speedup"] for r in rows])
    cpu_mpps = np.array([r["cpu_mpps"] for r in rows])
    gpu_mpps = np.array([r["gpu_mpps"] for r in rows])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    ax1.plot(ns, sp, "o-", color="C2")
    ax1.axhline(1.0, color="0.5", ls="--", lw=1, label="break-even (1x)")
    ax1.set_xscale("log")
    ax1.set_xlabel("number of particles")
    ax1.set_ylabel("GPU speed-up  (CPU time / GPU time)")
    ax1.set_title("EverestCrystal track: GPU speed-up vs bunch size")
    for x, y in zip(ns, sp):
        ax1.annotate(f"{y:.1f}x", (x, y), textcoords="offset points",
                     xytext=(0, 7), ha="center", fontsize=8)
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    ax2.plot(ns, cpu_mpps, "o-", color="C0", label="CPU (1 core)")
    ax2.plot(ns, gpu_mpps, "s-", color="C3", label="GPU (Cupy)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("number of particles")
    ax2.set_ylabel("throughput [million particles / s]")
    ax2.set_title("EverestCrystal track throughput")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    out_png = "everest_crystal_gpu_benchmark.png"
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_png}")

    with open("everest_crystal_gpu_benchmark.json", "w") as handle:
        json.dump({"stack_limit_bytes": int(eff_stack), "rows": rows},
                  handle, indent=2)
    best = max(rows, key=lambda r: r["speedup"])
    print(f"peak speed-up {best['speedup']:.1f}x at N={best['n_part']}")


if __name__ == "__main__":
    main()
