from pathlib import Path
import numpy as np
import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'

@for_all_test_contexts
def test_dp_pos(test_context):
    sweep = -300   # negative sweep => positive off-momentum
    num_turns = 6000
    num_particles = 5
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b1.json')

    line.build_tracker()

    part = line.build_particles(delta=np.linspace(-2e-4, 2e-4, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0)

    rf_sweep = xc.RFSweep(line)
    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    rf_sweep.track(sweep=sweep, num_turns=num_turns, particles=part)
    assert np.all(part.delta > 1.5e-3)


@for_all_test_contexts
def test_dp_neg(test_context):
    sweep = 300   # positive sweep => negative off-momentum
    num_turns = 6000
    num_particles = 5
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b2.json')

    line.build_tracker()

    part = line.build_particles(delta=np.linspace(-2e-4, 2e-4, num_particles),
                                x_norm=0, px_norm=0, y_norm=0, py_norm=0)

    rf_sweep = xc.RFSweep(line)
    # This sweep is 3.5 buckets, so check that all particles are at least 3 buckets away
    rf_sweep.track(sweep=sweep, num_turns=num_turns, particles=part)
    assert np.all(part.delta < -1.5e-3)
