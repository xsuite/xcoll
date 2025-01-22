# copyright ############################### #
# This file is part of the Xcoll package.   #
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


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@retry()
@pytest.mark.parametrize("beam, plane, npart, interpolation, ignore_crystals", [
                            [1, 'H', 2500, 0.2, True],
                            [2, 'V', 500, 0.3, True],
                            [1, 'V', 3500, 0.1, False],
                            [2, 'H', 30000, 0.15, False]
                        ], ids=["B1H", "B2V", "B1V_crystals", "B2H_crystals"])
def test_run_lossmap(beam, plane, npart, interpolation, ignore_crystals, test_context):

    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

    colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3_ir7.yaml',
                                    beam=beam, ignore_crystals=ignore_crystals)

    colldb.install_everest_collimators(line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    line.build_tracker(_context=test_context)
    line.collimators.assign_optics()

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = line[tcp].generate_pencil(npart)

    line.track(part, num_turns=2)
    assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals, 'EverestCollimator', 'EverestCrystal')


def assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals, coll_cls, cry_cls):
    with flaky_assertions():
        line_is_reversed = True if beam == 2 else False
        ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part,
                         interpolation=interpolation)
        print(ThisLM.summary)

        ThisLM.to_json("lossmap.json")
        assert Path("lossmap.json").exists()
        with Path("lossmap.json").open('r') as fid:
            dct = json.load(fid)
            assert xt.line._dicts_equal(dct, ThisLM.lossmap)
        Path("lossmap.json").unlink()
        ThisLM.save_summary("coll_summary.txt")
        assert Path("coll_summary.txt").exists()
        Path("coll_summary.txt").unlink()

        # TODO: check the lossmap quantitaively: rough amount of losses at given positions
        summ = ThisLM.summary
        assert list(summ.columns) == ['collname', 'nabs', 'length', 's', 'type']
        assert len(summ[summ.type==coll_cls]) == 10
        if not ignore_crystals:
            assert len(summ[summ.type==cry_cls]) == 2

        # We want at least 5% absorption on the primary
        assert summ.loc[summ.collname==tcp,'nabs'].values[0] > 0.05*npart

        lm = ThisLM.lossmap
        summ = summ[summ.nabs > 0]
        assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', \
                                   'interpolation', 'reversed']
        assert list(lm['collimator'].keys()) == ['s', 'name', 'length', 'n']
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
        if npart > 5000:
            assert len(lm['aperture']['s']) > 0
            assert len(lm['aperture']['s']) == len(lm['aperture']['name'])
            assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
            assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
        assert lm['interpolation'] == interpolation
        assert lm['reversed'] == line_is_reversed
