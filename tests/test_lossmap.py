# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import time
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

import xtrack as xt
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

import xcoll as xc
from xcoll.compare import deep_equal

from _common_api import all_engine_params

try:
    import matplotlib.pyplot as plt
except (ImportError, ModuleNotFoundError):
    plt = None


path = Path(__file__).parent / 'data'


@pytest.mark.parametrize("engine", all_engine_params)
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
def test_lossmap(engine, beam, plane, npart, interpolation, ignore_crystals, do_plot, test_context):
    if do_plot and plt is None:
        pytest.skip("matplotlib not installed")
    if not do_plot and plt is not None:
        pytest.skip("matplotlib installed")

    if engine == "fluka":
        if xc.fluka.engine.is_running():
            xc.fluka.engine.stop(clean=True)
    elif engine == "geant4":
        if xc.geant4.engine.is_running():
            xc.geant4.engine.stop(clean=True)
        if not ignore_crystals:
            pytest.skip("Geant4 crystals not implemented yet")

    env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
    line = env[f'lhcb{beam}']
    colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3_ir7.yaml',
                                    beam=beam, ignore_crystals=ignore_crystals)
    if engine == "everest":
        colldb.install_everest_collimators(line=line)
        assert np.all([isinstance(line[nn], xc.EverestCrystal) for nn in colldb.collimator_families['cry7']])
        assert np.all([isinstance(line[nn], xc.EverestCollimator) for nn in colldb.collimator_names
                    if nn not in colldb.collimator_families['cry7']])
    elif engine == "fluka":
        colldb.install_fluka_collimators(line=line)
        assert np.all([isinstance(line[nn], xc.FlukaCrystal) for nn in colldb.collimator_families['cry7']])
        assert np.all([isinstance(line[nn], xc.FlukaCollimator) for nn in colldb.collimator_names
                    if nn not in colldb.collimator_families['cry7']])
    elif engine == "geant4":
        colldb.install_geant4_collimators(line=line)
        assert np.all([isinstance(line[nn], xc.Geant4Crystal) for nn in colldb.collimator_families['cry7']])
        assert np.all([isinstance(line[nn], xc.Geant4Collimator) for nn in colldb.collimator_names
                    if nn not in colldb.collimator_families['cry7']])

    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)
    line.build_tracker(_context=test_context)
    line.collimators.assign_optics()
    if not ignore_crystals:
        line.collimators.align_to_beam_divergence()

    if engine == "fluka":
        xc.fluka.engine.start(line=line, capacity=2*npart, verbose=True)
    elif engine == "geant4":
        xc.geant4.engine.start(line=line, verbose=True)

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = line[tcp].generate_pencil(npart)
    if engine == "fluka":
        # Make sure all primaries hit the collimator (not trivial because of assembly)
        if plane == 'H':
            part.x[(part.x > 0) & (part.state > -99999)] += 1e-5
            part.x[(part.x < 0) & (part.state > -99999)] -= 1e-5
        elif plane == 'V':
            part.y[(part.y > 0) & (part.state > -99999)] += 1e-5
            part.y[(part.y < 0) & (part.state > -99999)] -= 1e-5
    line.scattering.enable()
    line.track(part, num_turns=2)
    line.scattering.disable()

    if engine == "everest":
        coll_cls = ['EverestCollimator']
        cry_cls  = ['EverestCrystal']
    elif engine == "fluka":
        coll_cls = ['FlukaCollimator']
        cry_cls  = ['FlukaCrystal']
    elif engine == "geant4":
        coll_cls = ['Geant4Collimator', 'Geant4CollimatorTip']
        cry_cls  = ['Geant4Crystal']
    this_id = f"B{beam}{plane}-{npart}-{engine}-{interpolation}-{ignore_crystals}-{test_context}"
    ThisLM = _assert_lossmap(beam, npart, line, part, tcp, interpolation, ignore_crystals,
                             coll_cls, cry_cls, this_id)
    if do_plot:
        ThisLM.plot(show=False, savefig=f"test-{this_id}.jpg")
        assert Path(f"test-{this_id}.jpg").exists()
        Path(f"test-{this_id}.jpg").unlink()

    if engine == "fluka":
        xc.fluka.engine.stop(clean=True)
    elif engine == "geant4":
        xc.geant4.engine.stop(clean=True)


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
            assert dct.pop('xcoll', None) == [xc.__version__]
            date = dct.pop('date', None)
            assert date is not None
            assert pd.Timestamp.now() - pd.Timestamp(date[-1]) < pd.Timedelta('1 minute')
            momentum = dct.pop('momentum', None)
            assert momentum is not None
            assert np.isclose(momentum, line.particle_ref.p0c[0])
            assert dct.pop('beam_type', None) == 2212
            assert dct.pop('num_initial', None) == npart
            assert np.isclose(dct.pop('tot_energy_initial', None), npart*line.particle_ref.energy0[0])
            cold_regions = dct.pop('cold_regions', None)
            warm_regions = dct.pop('warm_regions', None)
            s_range = dct.pop('s_range', None)
            assert deep_equal(cold_regions, ThisLM.cold_regions)
            assert deep_equal(warm_regions, ThisLM.warm_regions)
            assert deep_equal(s_range, ThisLM.s_range)
            assert deep_equal(dct, clean_lm_dct)
        ThisLM2 = xc.LossMap.from_json(f"lossmap-{this_id}.json")
        assert ThisLM == ThisLM2
        Path(f"lossmap-{this_id}.json").unlink()
        ThisLM.save_summary(f"coll_summary-{this_id}.txt")
        assert Path(f"coll_summary-{this_id}.txt").exists()
        Path(f"coll_summary-{this_id}.txt").unlink()

        assert np.isclose(ThisLM.momentum, line.particle_ref.p0c[0])
        assert ThisLM.beam_type == 'proton'
        assert ThisLM.num_initial == npart
        assert np.isclose(ThisLM.tot_energy_initial, npart*line.particle_ref.energy0[0])

        # TODO: check the lossmap quantitaively: rough amount of losses at given positions
        summ = ThisLM.summary
        assert list(summ.columns) == ['name', 'n', 'e', 'length', 's', 'type']
        assert len(summ[[tt in coll_cls for tt in summ.type]]) == 10
        if not ignore_crystals:
            assert len(summ[[tt in cry_cls for tt in summ.type]]) == 2

        # We want at least 5% absorption on the primary
        assert summ.loc[summ.name==tcp,'n'].values[0] > 0.05*npart

        lm = ThisLM.lossmap
        summ = summ[summ.n > 0]
        assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation',
                                   'reversed']
        assert lm['interpolation'] == interpolation
        assert lm['reversed'] == line_is_reversed
        assert np.isclose(lm['machine_length'], line.get_length())
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
            assert list(lm['aperture'].keys()) == ['idx_bins', 's_bins', 'length_bins', 'n_bins', 'e_bins']
            if npart > 5000:
                assert len(lm['aperture']['s_bins']) > 0
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['idx_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['n_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['e_bins'])
                assert len(lm['aperture']['s_bins']) == len(lm['aperture']['length_bins'])
                assert np.all([s < lm['machine_length'] for s in lm['aperture']['s_bins']])
        else:
            assert list(lm['aperture'].keys()) == ['name', 'n', 'e', 'length', 's', 'type']
            if npart > 5000:
                assert len(lm['aperture']['s']) > 0
                assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
                assert len(lm['aperture']['s']) == len(lm['aperture']['e'])
                assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
    return ThisLM
