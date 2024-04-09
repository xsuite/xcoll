# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import numpy as np
import pytest

import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'

@for_all_test_contexts
@pytest.mark.parametrize("sweep, beam", [[-300, 1], [300, 2]],
                         ids=["DP pos", "DP neg"])
def test_rf_sweep(sweep, beam, test_context):
    num_turns = 6000
    num_particles = 5
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

    line.build_tracker()

    part = line.build_particles(delta=np.linspace(-2e-4, 2e-4, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0)

    rf_sweep = xc.RFSweep(line)
    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    rf_sweep.track(sweep=sweep, num_turns=num_turns, particles=part)

    # negative sweep => positive off-momentum etc
    if sweep < 0:
        assert np.all(part.delta > 1.5e-3)
    else:
        assert np.all(part.delta < -1.5e-3)

