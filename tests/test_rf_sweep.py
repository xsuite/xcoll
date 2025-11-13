# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


path = Path(__file__).parent / 'data'


@for_all_test_contexts
@pytest.mark.parametrize("sweep, beam", [[-300, 1], [300, 2], [3500, 3]],
                         ids=["DP pos LHC", "DP neg LHC", "DP neg SPS"])
def test_rf_sweep(sweep, beam, test_context):
    num_turns = 6000
    num_particles = 5
    if beam == 3:
        bh = 3.5e-3
        sweep *= 25
        line = xt.load(path / f'sequence_sps_q20_inj.json')['sps']
        line.insert('cavity_mid', xt.Marker(), at=2822) # Inside thick cavity as to slice it
        tt_c = line.get_table().rows['ac.*']
        assert 'ThickSliceCavity' in tt_c.element_type
    else:
        bh = 3.e-4
        line = xt.load(path / f'sequence_lhc_run3_b{beam}.json')

    line.build_tracker(_context=test_context)

    part = line.build_particles(delta=np.linspace(-0.7*bh, 0.7*bh, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0, _context=test_context)

    rf_sweep = xc.RFSweep(line)
    rf_sweep.prepare(sweep_per_turn=sweep/num_turns)
    rf_sweep.info()
    line.track(particles=part, num_turns=num_turns)

    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    # negative sweep => positive off-momentum etc
    part.move(_context=xo.ContextCpu())
    if sweep < 0:
        assert np.all(part.delta > 3*bh)
    else:
        assert np.all(part.delta < -3*bh)


def test_rf_sweep_old_style():
    num_turns = 6000
    num_particles = 5
    line = xt.load(path / f'sequence_lhc_run3_b2.json')

    line.build_tracker()

    part = line.build_particles(delta=np.linspace(-2e-4, 2e-4, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0)

    sweep = -300
    rf_sweep = xc.RFSweep(line)
    with pytest.warns(DeprecationWarning, match="The `sweep` argument is deprecated."):
        with pytest.warns(DeprecationWarning, match="The `num_turns` argument is deprecated."):
            rf_sweep.info(sweep=sweep, num_turns=num_turns)
    with pytest.warns(DeprecationWarning):
        rf_sweep.track(sweep=sweep, num_turns=num_turns, particles=part)

    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    # negative sweep => positive off-momentum etc
    part.move(_context=xo.ContextCpu())
    assert np.all(part.delta > 1.5e-3)
