# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import numpy as np
import pandas as pd
from pathlib import Path
import xtrack as xt
import xcoll as xc
import pytest
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

try:
    import matplotlib.pyplot as plt
except (ImportError, ModuleNotFoundError):
    plt = None


path = Path(__file__).parent / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("do_plot", [True, False], ids=["with_plot", "without_plot"])
@pytest.mark.parametrize("beam, plane, npart, interpolation, ignore_crystals", [
                            [1, 'H', 2500, 0.2, True],
                            [2, 'V', 500, 0.3, True],
                            [1, 'V', 3500, False, False],
                            [2, 'H', 30000, 0.15, False]
                        ], ids=["B1H", "B2V", "B1V_crystals", "B2H_crystals"])
@retry()
def test_lossmap_everest(beam, plane, npart, interpolation, ignore_crystals, do_plot, test_context):
    if do_plot and plt is None:
        pytest.skip("matplotlib not installed")
    if not do_plot and plt is not None:
        pytest.skip("matplotlib installed")

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
    line.scattering.enable()
    line.track(part, num_turns=2)
    line.scattering.disable()
    this_id = f"B{beam}{plane}-{npart}-{interpolation}-{ignore_crystals}-{test_context}"
    ThisLM = _assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals, 'EverestCollimator', 'EverestCrystal', this_id)
    if do_plot:
        ThisLM.plot(show=False, savefig=f"test-{this_id}.jpg")
        assert Path(f"test-{this_id}.jpg").exists()
        Path(f"test-{this_id}.jpg").unlink()


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("do_plot", [True, False], ids=["with_plot", "without_plot"])
@retry()
@pytest.mark.skipif(not xc.geant4.environment.compiled, reason="BDSIM+Geant4 installation not found")
def test_lossmap_geant4(do_plot, test_context):
    # If a previous test failed, stop the server manually
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop(clean=True)

    npart = 5000
    beam = 2
    plane = 'H'
    interpolation = 0.1
    ignore_crystals = True
    do_plot = True

    if do_plot and plt is None:
        pytest.skip("matplotlib not installed")
    if not do_plot and plt is not None:
        pytest.skip("matplotlib installed")

    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3_ir7.yaml', beam=beam)
    colldb.install_geant4_collimators(line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)
    line.build_tracker(_context=test_context)
    line.collimators.assign_optics()

    xc.geant4.engine.start(line=line, bdsim_config_file=path / 'geant4_protons.gmad')

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = line[tcp].generate_pencil(npart)

    line.scattering.enable()
    line.track(part, num_turns=2)
    line.scattering.disable()
    xc.geant4.engine.stop(clean=True)
    this_id = f"B{beam}{plane}-{npart}-{interpolation}-{ignore_crystals}-{test_context}-geant4"
    ThisLM = _assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals, 'Geant4Collimator', 'Geant4Crystal', this_id)
    if do_plot:
        ThisLM.plot(show=False, savefig=f"test-{this_id}.jpg")
        assert Path(f"test-{this_id}.jpg").exists()
        Path(f"test-{this_id}.jpg").unlink()


def _assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals, coll_cls, cry_cls, this_id):
    with flaky_assertions():
        line_is_reversed = True if beam == 2 else False
        ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part,
                            interpolation=interpolation)
        print(ThisLM.summary)

        ThisLM.to_json(f"lossmap-{this_id}.json")
        assert Path(f"lossmap-{this_id}.json").exists()
        clean_lm_dct = {kk: {kkk: vvv.tolist() for kkk, vvv in vv.items()}
                        if isinstance(vv, dict) else vv
                        for kk, vv in ThisLM.lossmap.items()}
        with Path(f"lossmap-{this_id}.json").open('r') as fid:
            dct = json.load(fid)
            assert dct.pop('xcoll', None) == xc.__version__
            date = dct.pop('date', None)
            assert date is not None
            assert pd.Timestamp.now() - pd.Timestamp(date) < pd.Timedelta('1 minute')
            assert xt.line._dicts_equal(dct, clean_lm_dct)
        Path(f"lossmap-{this_id}.json").unlink()
        ThisLM.save_summary(f"coll_summary-{this_id}.txt")
        assert Path(f"coll_summary-{this_id}.txt").exists()
        Path(f"coll_summary-{this_id}.txt").unlink()

        # TODO: check the lossmap quantitaively: rough amount of losses at given positions
        summ = ThisLM.summary
        assert list(summ.columns) == ['name', 'n', 'e', 'length', 's', 'type']
        assert len(summ[summ.type==coll_cls]) == 10
        if not ignore_crystals:
            assert len(summ[summ.type==cry_cls]) == 2

        # We want at least 5% absorption on the primary
        assert summ.loc[summ.name==tcp,'n'].values[0] > 0.05*npart

        lm = ThisLM.lossmap
        summ = summ[summ.n > 0]
        assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation',
                                'reversed', 'momentum']
        assert list(lm['collimator'].keys()) == ['name', 'n', 'e', 'length', 's', 'type']
        assert len(lm['collimator']['s']) == len(summ)
        assert len(lm['collimator']['name']) == len(summ)
        assert len(lm['collimator']['length']) == len(summ)
        assert len(lm['collimator']['n']) == len(summ)
        assert len(lm['collimator']['e']) == len(summ)
        assert np.all(lm['collimator']['s'] == summ.s)
        assert np.all(lm['collimator']['name'] == summ.name)
        assert np.all(lm['collimator']['length'] == summ.length)
        assert np.all(lm['collimator']['n'] == summ.n)
        assert np.all(lm['collimator']['e'] == summ.e)
        assert np.all([nn[:3] in ['tcp', 'tcs'] for nn in lm['collimator']['name']])
        assert np.all([s < lm['machine_length'] for s in lm['collimator']['s']])
        if interpolation:
            assert list(lm['aperture'].keys()) == ['idx_bins', 's_bins', 'n_bins', 'e_bins', 'length_bins']
            if npart > 5000:
                assert len(lm['aperture']['s_bins']) > 0
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['idx_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['n_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['e_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['length_bins'])
                assert np.all([s < lm['machine_length'] for s in lm['aperture']['s_bins']])
        else:
            assert list(lm['aperture'].keys()) == ['s', 'n', 'e']
            if npart > 5000:
                assert len(lm['aperture']['s']) > 0
                assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
                assert len(lm['aperture']['s']) == len(lm['aperture']['e'])
                assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
        assert lm['interpolation'] == interpolation
        assert lm['reversed'] == line_is_reversed
    return ThisLM
