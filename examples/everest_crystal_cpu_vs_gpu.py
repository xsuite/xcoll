# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""Visual CPU-vs-GPU comparison of the Everest *crystal* element.

Tracks the SAME fixed-seed bunch through the canonical
``Drift + EverestCrystal + Drift`` line on ``xobjects.ContextCpu`` and
``xobjects.ContextCupy`` and overlays the resulting distributions. This is the
documentation figure for enabling ``EverestCrystal`` on the CUDA context:
agreement is *distributional* (parallel RNG streams + floating-point ordering
differ across contexts, so a per-particle bitwise match is neither expected nor
required -- mirrors
``xtrack/examples/random_number_generator/002_rutherford_cpu_vs_gpu.py``).

It produces:
  1. ``everest_crystal_cpu_vs_gpu_coords.png`` -- per-coordinate (x, px, y, py,
     delta) PDF overlay, empirical-CDF/KS panel, and a GPU/CPU bin-ratio panel
     with a +/-1 sigma counting-noise band.
  2. ``everest_crystal_cpu_vs_gpu_angle_efficiency.png`` -- the channeling
     efficiency vs incoming angle curve on CPU and GPU, plus the residual.

CUDA stack-limit note
---------------------
The Everest crystal kernels build several sizeable per-track structs on the
per-thread CUDA stack which together exceed the 1024 B CUDA default, so a GPU
track traps with cudaErrorIllegalAddress unless the limit is raised first.
``xcoll`` raises it automatically at ``EverestCrystal`` construction via
``xcoll.scattering_routines.everest.set_crystal_stack_limit``; this script
re-asserts it idempotently so it is self-contained.

Run with cupy's CUDA libs on the loader path, e.g.:
    LDP=$(ls -d $(python -c 'import sys;print(sys.prefix)')/lib/python*/site-packages/nvidia/*/lib | tr '\n' ':')
    LD_LIBRARY_PATH="$LDP" python everest_crystal_cpu_vs_gpu.py
"""

import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

from xcoll.scattering_routines.everest import set_crystal_stack_limit


# ---------------------------------------------------------------------------
# Canonical Drift + EverestCrystal + Drift configuration (matches the test in
# tests/test_everest_crystal_gpu.py). Bent silicon strip crystal, channeling
# regime; values match examples/everest_crystal.py.
# ---------------------------------------------------------------------------
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

N_PART = 20000
BEAM_X_MEAN = 1.5e-3
BEAM_X_SPREAD = 75e-6
BEAM_PX_MEAN = 0.0
BEAM_PX_SPREAD = 30e-6
BEAM_Y_SPREAD = 1e-3
BEAM_PY_SPREAD = 5e-6

NUMPY_SEED = 20260603
PARTICLES_SEED = 13579
CHANNELING_KICK_THRESHOLD = 0.5 * CRYSTAL_BENDING_ANGLE

COORDS = ("x", "px", "y", "py", "delta")
COORD_UNITS = {"x": "m", "px": "rad", "y": "m", "py": "rad", "delta": "1"}


def _build_crystal(context):
    return xc.EverestCrystal(
        length=CRYSTAL_LENGTH, material=xc.materials.Silicon,
        bending_angle=CRYSTAL_BENDING_ANGLE, width=CRYSTAL_WIDTH,
        height=CRYSTAL_HEIGHT, side=CRYSTAL_SIDE, miscut=CRYSTAL_MISCUT,
        lattice=CRYSTAL_LATTICE, jaw=CRYSTAL_JAW, _context=context)


def _sample_bunch_coords(n_part, px_mean=BEAM_PX_MEAN, px_spread=BEAM_PX_SPREAD):
    rng = np.random.default_rng(NUMPY_SEED)
    initial_x = rng.normal(loc=BEAM_X_MEAN, scale=BEAM_X_SPREAD, size=n_part)
    initial_px = rng.normal(loc=px_mean, scale=px_spread, size=n_part)
    initial_y = rng.normal(loc=0.0, scale=BEAM_Y_SPREAD, size=n_part)
    initial_py = rng.normal(loc=0.0, scale=BEAM_PY_SPREAD, size=n_part)
    return initial_x, initial_px, initial_y, initial_py


def _build_line(context, initial_px=None):
    crystal = _build_crystal(context)
    line = xt.Line(
        elements=[xt.Drift(length=DRIFT_LENGTH), crystal,
                  xt.Drift(length=DRIFT_LENGTH)],
        element_names=["d_upstream", "crystal", "d_downstream"])
    line.particle_ref = xt.Particles(p0c=P0C_EV, mass0=xp.PROTON_MASS_EV)
    line.build_tracker(_context=context)
    sampled_x, sampled_px, sampled_y, sampled_py = _sample_bunch_coords(N_PART)
    if initial_px is not None:  # angle-scan pencil: override px
        sampled_px = np.full(N_PART, initial_px, dtype=np.float64)
    particles = line.build_particles(x=sampled_x, px=sampled_px,
                                     y=sampled_y, py=sampled_py)
    host_seeds = PARTICLES_SEED + np.arange(N_PART, dtype=np.uint32)
    particles._init_random_number_generator(seeds=host_seeds)
    return line, particles


# ---------------------------------------------------------------------------
# distributional metrics (mirror the Rutherford CPU-vs-GPU example)
# ---------------------------------------------------------------------------
def binned_kl(cpu_values, gpu_values, n_bins=120):
    """Binned KL divergence KL(CPU || GPU) over a shared range (nats)."""
    lo = float(min(cpu_values.min(), gpu_values.min()))
    hi = float(max(cpu_values.max(), gpu_values.max()))
    if hi <= lo:
        return 0.0
    edges = np.linspace(lo, hi, n_bins + 1)
    h_cpu, _ = np.histogram(cpu_values, bins=edges, density=True)
    h_gpu, _ = np.histogram(gpu_values, bins=edges, density=True)
    widths = np.diff(edges)
    p = h_cpu * widths + 1e-12
    q = h_gpu * widths + 1e-12
    p, q = p / p.sum(), q / q.sum()
    return float(np.sum(p * np.log(p / q)))


def empirical_cdf(values, grid):
    return np.searchsorted(np.sort(values), grid) / values.size


def track_canonical(context):
    """Track the canonical fixed-seed bunch one turn; return a coord dict."""
    set_crystal_stack_limit(context)  # no-op on CPU
    line, particles = _build_line(context)
    px_in = np.asarray(context.nparray_from_context_array(particles.px),
                       dtype=np.float64).copy()
    line.track(particles, num_turns=1)
    out = {}
    for field in COORDS + ("state", "particle_id"):
        out[field] = np.asarray(
            context.nparray_from_context_array(getattr(particles, field)))
    out["state"] = out["state"].astype(np.int64)
    out["_px_in_raw"] = px_in
    return out


def fractions(coord_dict):
    state = coord_dict["state"]
    px_out = coord_dict["px"]
    pid = coord_dict["particle_id"].astype(np.int64)
    px_in = coord_dict["_px_in_raw"][pid]
    alive = state > 0
    kick = px_out - px_in
    n_total = int(state.size)
    n_alive = int(np.sum(alive))
    n_chan = int(np.sum(alive & (kick > CHANNELING_KICK_THRESHOLD)))
    return {
        "n_total": n_total, "n_alive": n_alive, "n_lost": n_total - n_alive,
        "n_channeled": n_chan,
        "channeling_fraction_of_alive": (n_chan / n_alive) if n_alive else 0.0,
    }


def channeling_efficiency(context, angle):
    """Fraction of survivors channeled for a pencil at incoming `angle`."""
    set_crystal_stack_limit(context)
    line, particles = _build_line(context, initial_px=angle)
    line.track(particles, num_turns=1)
    state = np.asarray(context.nparray_from_context_array(particles.state))
    px_out = np.asarray(context.nparray_from_context_array(particles.px),
                        dtype=np.float64)
    alive = state > 0
    kick = px_out - angle  # px_in == angle for the whole pencil
    n_alive = int(np.sum(alive))
    n_chan = int(np.sum(alive & (kick > CHANNELING_KICK_THRESHOLD)))
    return (n_chan / n_alive) if n_alive else 0.0


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
def plot_coord_panels(cpu, gpu, ks_stats, kl_stats, path):
    n = len(COORDS)
    fig, axes = plt.subplots(3, n, figsize=(4.0 * n, 11.5))
    for j, field in enumerate(COORDS):
        cpu_alive = cpu[field][cpu["state"] > 0]
        gpu_alive = gpu[field][gpu["state"] > 0]
        lo = float(min(cpu_alive.min(), gpu_alive.min()))
        hi = float(max(cpu_alive.max(), gpu_alive.max()))
        edges = np.linspace(lo, hi, 80)
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax = axes[0, j]
        ax.hist(cpu_alive, bins=edges, density=True, histtype="stepfilled",
                alpha=0.35, color="C0", label="CPU")
        ax.hist(gpu_alive, bins=edges, density=True, histtype="step", lw=1.6,
                color="C3", label="GPU (Cupy)")
        ax.set_title(f"{field} [{COORD_UNITS[field]}]\nKS={ks_stats[field]:.4f}"
                     f"  KL={kl_stats[field]:.2e}", fontsize=10)
        if j == 0:
            ax.set_ylabel("density")
            ax.legend(fontsize=8)

        ax = axes[1, j]
        grid = np.linspace(lo, hi, 1500)
        cdf_cpu = empirical_cdf(cpu_alive, grid)
        cdf_gpu = empirical_cdf(gpu_alive, grid)
        ax.plot(grid, cdf_cpu, color="C0", lw=1.6, label="CPU")
        ax.plot(grid, cdf_gpu, color="C3", lw=1.1, ls="--", label="GPU")
        imax = int(np.argmax(np.abs(cdf_cpu - cdf_gpu)))
        ax.annotate(f"max|dCDF|={ks_stats[field]:.4f}",
                    xy=(grid[imax], cdf_cpu[imax]), xytext=(0.30, 0.20),
                    textcoords="axes fraction", fontsize=8,
                    arrowprops=dict(arrowstyle="->", color="0.3"))
        if j == 0:
            ax.set_ylabel("empirical CDF")

        ax = axes[2, j]
        c_cpu, _ = np.histogram(cpu_alive, bins=edges)
        c_gpu, _ = np.histogram(gpu_alive, bins=edges)
        mask = c_cpu > 0
        ratio = np.full_like(centers, np.nan)
        ratio[mask] = c_gpu[mask] / c_cpu[mask]
        noise = np.full_like(centers, np.nan)
        noise[mask] = np.sqrt(1.0 / c_cpu[mask]
                              + 1.0 / np.maximum(c_gpu[mask], 1))
        ax.axhline(1.0, color="k", lw=0.9)
        ax.fill_between(centers, 1 - noise, 1 + noise, color="0.8",
                        label="+/-1 sigma count noise")
        ax.plot(centers, ratio, ".", ms=3, color="C2", label="GPU/CPU")
        ax.set_ylim(0.0, 2.0)
        ax.set_xlabel(f"{field} [{COORD_UNITS[field]}]")
        if j == 0:
            ax.set_ylabel("GPU/CPU bin ratio")
            ax.legend(fontsize=8)
    fig.suptitle("EverestCrystal CPU vs GPU (Cupy) -- per-coordinate "
                 "distributional agreement (surviving particles)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_angle_curve(angles, eff_cpu, eff_gpu, path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    angles_urad = angles * 1e6
    ax = axes[0]
    ax.plot(angles_urad, eff_cpu, "o-", color="C0", ms=4, label="CPU")
    ax.plot(angles_urad, eff_gpu, "s--", color="C3", ms=4, label="GPU (Cupy)")
    ax.set_xlabel("incoming angle px [urad]")
    ax.set_ylabel("channeling efficiency (channeled / alive)")
    ax.set_title("Channeling-efficiency curve CPU vs GPU")
    ax.legend()
    ax = axes[1]
    ax.plot(angles_urad, eff_gpu - eff_cpu, "d-", color="C2", ms=4)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("incoming angle px [urad]")
    ax.set_ylabel("GPU - CPU efficiency")
    ax.set_title("Curve residual (RNG-stream noise)")
    fig.suptitle("Angle scan: channeling efficiency vs incoming angle",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main():
    cpu_ctx = xo.ContextCpu()
    gpu_ctx = xo.ContextCupy()
    eff_stack = set_crystal_stack_limit(gpu_ctx)
    print(f"[setup] GPU per-thread CUDA stack limit = {eff_stack} bytes")

    print("[1] canonical Drift+EverestCrystal+Drift, identical seeds ...")
    cpu = track_canonical(cpu_ctx)
    gpu = track_canonical(gpu_ctx)
    ks_stats, kl_stats = {}, {}
    for field in COORDS:
        ca = cpu[field][cpu["state"] > 0]
        ga = gpu[field][gpu["state"] > 0]
        ks_stats[field] = float(ks_2samp(ca, ga).statistic)
        kl_stats[field] = binned_kl(ca, ga)
    frac_cpu, frac_gpu = fractions(cpu), fractions(gpu)
    chan_delta = abs(frac_cpu["channeling_fraction_of_alive"]
                     - frac_gpu["channeling_fraction_of_alive"])
    print("  per-coordinate KS / KL (surviving particles):")
    for field in COORDS:
        print(f"    {field:6s} KS={ks_stats[field]:.6f}  "
              f"KL={kl_stats[field]:.3e}")
    print(f"  channeling fraction CPU={frac_cpu['channeling_fraction_of_alive']:.5f}"
          f"  GPU={frac_gpu['channeling_fraction_of_alive']:.5f}"
          f"  |delta|={chan_delta:.5f}")

    print("[2] angle scan (channeling-efficiency curve) ...")
    angles = np.linspace(-15e-6, 15e-6, 13)
    eff_cpu = np.array([channeling_efficiency(cpu_ctx, a) for a in angles])
    eff_gpu = np.array([channeling_efficiency(gpu_ctx, a) for a in angles])
    angle_resid = float(np.max(np.abs(eff_gpu - eff_cpu)))
    print(f"  angle-curve max|delta| = {angle_resid:.4f}")

    print("[3] saving plots ...")
    coord_png = "everest_crystal_cpu_vs_gpu_coords.png"
    angle_png = "everest_crystal_cpu_vs_gpu_angle_efficiency.png"
    plot_coord_panels(cpu, gpu, ks_stats, kl_stats, coord_png)
    plot_angle_curve(angles, eff_cpu, eff_gpu, angle_png)
    print(f"  wrote {coord_png}")
    print(f"  wrote {angle_png}")

    summary = {
        "ks_per_coord": ks_stats, "kl_per_coord": kl_stats,
        "max_coord_ks": max(ks_stats.values()),
        "channeling_fraction_cpu": frac_cpu["channeling_fraction_of_alive"],
        "channeling_fraction_gpu": frac_gpu["channeling_fraction_of_alive"],
        "channeling_fraction_delta": chan_delta,
        "angle_curve_max_abs_delta": angle_resid,
        "stack_limit_bytes": int(eff_stack),
    }
    with open("everest_crystal_cpu_vs_gpu_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(f"max coord KS = {summary['max_coord_ks']:.6f}  "
          f"channeling |delta| = {chan_delta:.5f}  "
          f"angle-curve max|delta| = {angle_resid:.4f}")


if __name__ == "__main__":
    main()
