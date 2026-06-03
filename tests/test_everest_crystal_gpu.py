# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""CPU-vs-GPU distributional agreement for the Everest *crystal* element.

Historically every Everest test in ``test_everest.py`` excluded the GPU
contexts with the note ``# Rutherford RNG not on GPU``. With the upstream
xtrack ``RandomRutherford`` GPU port (PR1) and this xcoll GPU port (PR2) the
``EverestCrystal`` now compiles and tracks on ``xobjects.ContextCupy``, so this
test exercises exactly the path that was previously excluded.

The acceptance is **distributional**, not per-particle bitwise. The parallel
RNG streams and floating-point association ordering differ between
``ContextCpu`` and ``ContextCupy`` (see
``xtrack/tests/test_random_gen_ruth.py::test_cpu_gpu_distribution_match`` and
``xtrack/examples/random_number_generator/002_rutherford_cpu_vs_gpu.py``), so a
bitwise per-particle match across contexts is neither expected nor required.
We instead assert:

  * a small two-sample Kolmogorov-Smirnov statistic per surviving coordinate
    (interpret by magnitude; threshold ``KS < 0.02``), and
  * agreement of the physics observable -- the channeling fraction and the
    surviving fraction.

CUDA stack-limit requirement
----------------------------
The Everest crystal kernels build several sizeable per-track structs on the
per-thread CUDA stack (the geometry descriptor, the collimation data and the
Everest scattering state) which together exceed the CUDA default per-thread
stack frame of 1024 B. A GPU track therefore traps with
``cudaErrorIllegalAddress`` unless the device stack limit is raised first.
``xcoll`` raises it automatically at ``EverestCrystal`` construction time via
``xcoll.scattering_routines.everest.set_crystal_stack_limit`` (16 KiB by
default, ~8x the measured worst case); this test re-asserts it idempotently so
it is self-contained.
"""

import numpy as np
import pytest

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

from xcoll.scattering_routines.everest import set_crystal_stack_limit


# ---------------------------------------------------------------------------
# Canonical Drift + EverestCrystal + Drift line (single source of truth here).
# Bent silicon strip crystal in the channeling regime; values match the clone's
# examples/everest_crystal.py (a physically sensible SPS-style config).
# ---------------------------------------------------------------------------
P0C_EV = 400e9                  # 400 GeV protons (SPS-like)

CRYSTAL_LENGTH = 0.002          # m
CRYSTAL_BENDING_ANGLE = 149e-6  # rad (full bend over the crystal length)
CRYSTAL_WIDTH = 0.002           # m
CRYSTAL_HEIGHT = 0.05           # m
CRYSTAL_SIDE = "+"
CRYSTAL_MISCUT = 0.0
CRYSTAL_LATTICE = "strip"
CRYSTAL_JAW = 0.001             # m

DRIFT_LENGTH = 1.0              # m

N_PART = 20000
BEAM_X_MEAN = 1.5e-3            # m  (just inside the jaw at +1 mm)
BEAM_X_SPREAD = 75e-6           # m
BEAM_PX_MEAN = 0.0             # rad (aimed along the crystal entrance face)
BEAM_PX_SPREAD = 30e-6          # rad (narrow -> near-pencil)
BEAM_Y_SPREAD = 1e-3            # m
BEAM_PY_SPREAD = 5e-6           # rad

NUMPY_SEED = 20260603
PARTICLES_SEED = 13579

# A particle counts as "channeled" if its angular kick px_out - px_in exceeds
# this fraction of the full bending angle (the channeled population picks up
# ~the full bend; the amorphous tail does not).
CHANNELING_KICK_THRESHOLD = 0.5 * CRYSTAL_BENDING_ANGLE

# Distributional acceptance thresholds.
KS_THRESHOLD = 0.02
CHANNELING_FRACTION_TOLERANCE = 0.02

COORDS = ("x", "px", "y", "py", "delta")


def _sample_bunch_coords(n_part):
    """Reproducible numpy sampling of the initial bunch coordinates."""
    rng = np.random.default_rng(NUMPY_SEED)
    initial_x = rng.normal(loc=BEAM_X_MEAN, scale=BEAM_X_SPREAD, size=n_part)
    initial_px = rng.normal(loc=BEAM_PX_MEAN, scale=BEAM_PX_SPREAD, size=n_part)
    initial_y = rng.normal(loc=0.0, scale=BEAM_Y_SPREAD, size=n_part)
    initial_py = rng.normal(loc=0.0, scale=BEAM_PY_SPREAD, size=n_part)
    return initial_x, initial_px, initial_y, initial_py


def _build_line(context):
    """Build the canonical (line, particles) on `context`."""
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

    initial_x, initial_px, initial_y, initial_py = _sample_bunch_coords(N_PART)
    particles = line.build_particles(x=initial_x, px=initial_px,
                                     y=initial_y, py=initial_py)
    # Seed the per-particle RNG on the host: particle_id is the build-order
    # arange, so a host arange reproduces identical seeds on CPU and GPU and
    # keeps the Everest stochastic physics reproducible within a context.
    host_seeds = PARTICLES_SEED + np.arange(N_PART, dtype=np.uint32)
    particles._init_random_number_generator(seeds=host_seeds)
    return line, particles


def _track_canonical(context):
    """Track the canonical fixed-seed bunch one turn; return a coord dict."""
    # Re-assert the per-thread CUDA stack limit (no-op on CPU). EverestCrystal
    # already raises it at construction; doing it here keeps the test
    # self-contained.
    set_crystal_stack_limit(context)
    line, particles = _build_line(context)
    px_in = np.asarray(
        context.nparray_from_context_array(particles.px),
        dtype=np.float64).copy()  # build order == particle_id order

    line.track(particles, num_turns=1)

    out = {}
    for field in COORDS + ("state", "particle_id"):
        out[field] = np.asarray(
            context.nparray_from_context_array(getattr(particles, field)))
    out["state"] = out["state"].astype(np.int64)
    out["_px_in_raw"] = px_in
    return out


def _fractions(coord_dict):
    """Channeled / surviving fractions from a tracked coord dict."""
    state = coord_dict["state"]
    px_out = coord_dict["px"]
    # Realign px_in (sampled in build/particle_id order) to the (possibly
    # unsorted on GPU) output by indexing with each output particle's id.
    pid = coord_dict["particle_id"].astype(np.int64)
    px_in = coord_dict["_px_in_raw"][pid]
    alive = state > 0
    kick = px_out - px_in
    n_total = int(state.size)
    n_alive = int(np.sum(alive))
    n_chan = int(np.sum(alive & (kick > CHANNELING_KICK_THRESHOLD)))
    return {
        "n_total": n_total,
        "n_alive": n_alive,
        "n_lost": n_total - n_alive,
        "n_channeled": n_chan,
        "survival_fraction": n_alive / n_total,
        "channeling_fraction_of_alive": (n_chan / n_alive) if n_alive else 0.0,
    }


def test_everest_crystal_cpu_gpu_distribution_match():
    """EverestCrystal must produce the same distribution on CPU and GPU.

    The previously-excluded GPU path is now exercised: build the canonical
    Drift+EverestCrystal+Drift line, set the CUDA stack limit, track the same
    fixed-seed bunch on ContextCpu and ContextCupy, and assert per-coordinate
    KS < threshold plus channeling-fraction agreement.
    """
    pytest.importorskip("cupy")  # GPU contexts only; skip cleanly without cupy
    from scipy.stats import ks_2samp

    cpu = _track_canonical(xo.ContextCpu())
    gpu = _track_canonical(xo.ContextCupy())

    # The deterministic part (loss bookkeeping) must agree exactly.
    frac_cpu = _fractions(cpu)
    frac_gpu = _fractions(gpu)

    # Per-coordinate KS over the surviving particles.
    ks_stats = {}
    for field in COORDS:
        cpu_alive = cpu[field][cpu["state"] > 0]
        gpu_alive = gpu[field][gpu["state"] > 0]
        ks = ks_2samp(cpu_alive, gpu_alive)
        ks_stats[field] = float(ks.statistic)
        assert ks.statistic < KS_THRESHOLD, (
            f"CPU-vs-GPU EverestCrystal distributions disagree on '{field}': "
            f"KS={ks.statistic:.4f} (p={ks.pvalue:.3g}) "
            f">= threshold {KS_THRESHOLD}")

    # Physics observable: channeling fraction (of survivors) must agree.
    chan_delta = abs(frac_cpu["channeling_fraction_of_alive"]
                     - frac_gpu["channeling_fraction_of_alive"])
    assert chan_delta < CHANNELING_FRACTION_TOLERANCE, (
        f"CPU-vs-GPU channeling fraction disagrees: "
        f"CPU={frac_cpu['channeling_fraction_of_alive']:.4f} "
        f"GPU={frac_gpu['channeling_fraction_of_alive']:.4f} "
        f"|delta|={chan_delta:.4f} >= {CHANNELING_FRACTION_TOLERANCE}")

    # Surviving fraction must agree too (the loss process is the same physics).
    surv_delta = abs(frac_cpu["survival_fraction"]
                     - frac_gpu["survival_fraction"])
    assert surv_delta < CHANNELING_FRACTION_TOLERANCE, (
        f"CPU-vs-GPU survival fraction disagrees: "
        f"CPU={frac_cpu['survival_fraction']:.4f} "
        f"GPU={frac_gpu['survival_fraction']:.4f} "
        f"|delta|={surv_delta:.4f} >= {CHANNELING_FRACTION_TOLERANCE}")

    print(f"\nmax coord KS = {max(ks_stats.values()):.6f}  "
          f"channeling |delta| = {chan_delta:.5f}  "
          f"survival |delta| = {surv_delta:.5f}")


def test_crystal_stack_limit_is_raised_on_gpu():
    """The CUDA per-thread stack must be raised above the 1024 B default.

    A GPU EverestCrystal track traps with cudaErrorIllegalAddress at the CUDA
    default per-thread stack frame (1024 B); xcoll raises it to >= 16 KiB.
    """
    pytest.importorskip("cupy")
    import cupy as cp

    context = xo.ContextCupy()
    effective = set_crystal_stack_limit(context)
    assert effective is not None
    # Must clear the CUDA default and reach at least the documented minimum.
    assert effective >= 16 * 1024, (
        f"crystal stack limit {effective} B is below the documented "
        f"16 KiB minimum")
    # It is genuinely raised above the 1024 B CUDA default.
    assert effective > 1024
    # No-op on CPU (CPU numerics untouched).
    assert set_crystal_stack_limit(xo.ContextCpu()) is None
    # Idempotent / never-lower.
    again = set_crystal_stack_limit(context, nbytes=1024)
    assert again == cp.cuda.runtime.deviceGetLimit(0x00)
    assert again >= effective
