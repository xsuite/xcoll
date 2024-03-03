# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import pytest
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'

# https://github.com/xsuite/xtrack/blob/18b1ac33d6a9d87a156e87bfb71cb2c8011085f6/tests/test_radiation.py#LL138C5-L138C29
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@retry()
@pytest.mark.parametrize("beam, plane, npart, interpolation, ignore_crystals", [
                            [1, 'H', 25000, 0.2, True],
                            [2, 'V', 25000, 0.3, True],
                            [1, 'V', 35000, 0.1, False],
                            [2, 'H', 30000, 0.15, False]
                        ], ids=["B1", "B2V", "B1V_crystals", "B2H_crystals"])
def test_run_lossmap(beam, plane, npart, interpolation, ignore_crystals, test_context):

    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

    coll_manager = xc.CollimatorManager.from_yaml(path / f'colldb_lhc_run3_ir7.yaml', line=line, beam=beam,
                                                  ignore_crystals=ignore_crystals, _context=test_context)

    coll_manager.install_everest_collimators()
    coll_manager.build_tracker()
    coll_manager.set_openings()

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=npart)

    coll_manager.enable_scattering()
    line.track(part, num_turns=2)
    coll_manager.disable_scattering()

    with flaky_assertions():
        summ = coll_manager.summary(part, show_zeros=False)
        assert list(summ.columns) == ['collname', 'nabs', 'dE', 'length', 's', 'type']
        assert len(summ[summ.type=='EverestCollimator']) == 10
        # We want at least 5% absorption on the primary
        assert summ.loc[summ.collname==tcp,'nabs'].values[0] > 0.05*npart

        lm = coll_manager.lossmap(part, interpolation=interpolation)
        assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation', 'reversed']
        assert list(lm['collimator'].keys()) == ['s', 'name', 'length', 'n', 'dE']
        assert len(lm['collimator']['s']) == len(summ)
        assert len(lm['collimator']['name']) == len(summ)
        assert len(lm['collimator']['length']) == len(summ)
        assert len(lm['collimator']['n']) == len(summ)
        assert np.all(lm['collimator']['s'] == summ.s)
        assert np.all(lm['collimator']['name'] == summ.collname)
        assert np.all(lm['collimator']['length'] == summ.length)
        assert np.all(lm['collimator']['n'] == summ.nabs)
        assert np.all([nn[:3] in ['tcp', 'tcs'] for nn in lm['collimator']['name']])
        assert np.all([s < lm['machine_length'] for s in lm['collimator']['s']])
        assert list(lm['aperture'].keys()) == ['s', 'name', 'n']
        assert len(lm['aperture']['s']) > 0
        assert len(lm['aperture']['s']) == len(lm['aperture']['name'])
        assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
        assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
        assert lm['interpolation'] == interpolation
        line_is_reversed = True if beam==2 else False
        assert lm['reversed'] == line_is_reversed

