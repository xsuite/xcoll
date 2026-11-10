"""RF frequency sweep setup and tracking in a real LHC Run 3 lattice."""

from pathlib import Path

import numpy as np
import pytest
import xtrack as xt
import xobjects as xo
import xcoll as xc

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"
SEQUENCE = DATA / "sequence_lhc_run3_b1.json"


def _line():
    if not SEQUENCE.exists():
        pytest.skip(f"Missing LHC sequence: {SEQUENCE}")
    env = xt.load(SEQUENCE)
    line = env["lhcb1"]
    line.build_tracker(_context=xo.ContextCpu())
    return line


@pytest.mark.full
@pytest.mark.benchmark(group="RF sweep | preparation")
def bench_rf_sweep_prepare_lhc(benchmark):
    # Loading and tracker construction happen in setup, not in the measured
    # RFSweep preparation, which includes optics and expression installation.
    def setup():
        return (_line(),), {}

    def prepare(line):
        sweep = xc.RFSweep(line)
        sweep.prepare(sweep_per_turn=-0.05)
        return sweep

    result = benchmark.pedantic(prepare, setup=setup, rounds=3, warmup_rounds=1)
    assert result._prepared
    benchmark.extra_info.update(
        protocol=2, physics_case="rf-sweep-prepare-lhc",
        sweep_per_turn_hz=-0.05, beam=1,
    )


@pytest.mark.full
@pytest.mark.benchmark(group="RF sweep | tracking")
def bench_rf_sweep_tracking_lhc(benchmark):
    line = _line()
    sweep = xc.RFSweep(line)
    sweep.prepare(sweep_per_turn=-0.05)
    particles0 = line.build_particles(
        delta=np.linspace(-1e-4, 1e-4, 64),
        x_norm=0, px_norm=0, y_norm=0, py_norm=0,
    )
    line.track(particles0.copy(), num_turns=5)  # warm up

    def setup():
        sweep.reset()  # reset time expression outside measurement
        return (particles0.copy(),), {}

    def track(particles):
        line.track(particles, num_turns=5)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    assert np.all(np.asarray(result.state) > 0)
    benchmark.extra_info.update(
        protocol=2, physics_case="rf-sweep-track-lhc",
        particles=64, turns=5, sweep_per_turn_hz=-0.05,
    )
