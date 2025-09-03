# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import numpy as np
import pytest

import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


path = xc._pkg_root.parent / 'tests' / 'data'


@for_all_test_contexts
@pytest.mark.parametrize("sweep, beam", [[-300, 1], [300, 2], [3500, 3]],
                         ids=["DP pos LHC", "DP neg LHC", "DP neg SPS"])
def test_rf_sweep(sweep, beam, test_context):
    num_turns = 6000
    num_particles = 5
    if beam == 3:
        line = xt.load(path / f'sequence_sps_q20_inj.json')['sps']
        line.insert('markkk', xt.Marker(), at=2822)
        tt_c = line.get_table().rows['ac.*']
        assert 'ThickSliceCavity' in tt_c.element_type
    else:
        line = xt.load(path / f'sequence_lhc_run3_b{beam}.json')

    line.build_tracker()

    part = line.build_particles(delta=np.linspace(-2e-4, 2e-4, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0)

    rf_sweep = xc.RFSweep(line)
    rf_sweep.prepare(sweep_per_turn=sweep/num_turns)
    rf_sweep.info()
    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    line.track(particles=part, num_turns=num_turns)

    # negative sweep => positive off-momentum etc
    if sweep < 0:
        assert np.all(part.delta > 1.5e-3)
    else:
        assert np.all(part.delta < -1.5e-3)

