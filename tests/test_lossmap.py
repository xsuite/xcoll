import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc


path = Path.cwd() / 'data'


def _run_lossmap(beam, plane, npart, interpolation, ignore_crystals=True):

    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

    coll_manager = xc.CollimatorManager.from_yaml(path / f'colldb_lhc_run3_ir7.yaml', line=line, beam=beam, ignore_crystals=ignore_crystals)

    coll_manager.install_everest_collimators()
    coll_manager.build_tracker()
    coll_manager.set_openings()

    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    part = coll_manager.generate_pencil_on_collimator(tcp, num_particles=npart)

    coll_manager.enable_scattering()
    line.track(part, num_turns=2)
    coll_manager.disable_scattering()

    summ = coll_manager.summary(part, show_zeros=False)
    assert list(summ.columns) == ['collname', 'nabs', 'length', 's', 'type']
    assert len(summ) == 10
    # We want at least 5% absorption on the primary
    assert summ.loc[summ.collname==tcp,'nabs'].values[0] > 0.05*npart

    lm = coll_manager.lossmap(part, interpolation=interpolation)
    assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation', 'reversed']
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
    assert len(lm['aperture']['s']) > 0
    assert len(lm['aperture']['s']) == len(lm['aperture']['name'])
    assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
    assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
    assert lm['interpolation'] == interpolation
    line_is_reversed = True if beam==2 else False
    assert lm['reversed'] == line_is_reversed


def test_lossmap_B1H():
    _run_lossmap(1, 'H', 25000, 0.2)

def test_lossmap_B2V():
    _run_lossmap(2, 'V', 25000, 0.3)

def test_lossmap_crystals_B1V():
    _run_lossmap(1, 'V', 35000, 0.1, ignore_crystals=False)

def test_lossmap_crystals_B2H():
    _run_lossmap(2, 'H', 30000, 0.15, ignore_crystals=False)

