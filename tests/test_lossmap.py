import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry

path = Path.cwd() / 'data'

# https://github.com/xsuite/xtrack/blob/18b1ac33d6a9d87a156e87bfb71cb2c8011085f6/tests/test_radiation.py#LL138C5-L138C29
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


@retry
def test_lossmap_B1H():
    with flaky_assertions():
        _run_lossmap(1, 'H', 25000, 0.2)

@retry
def test_lossmap_B2V():
    with flaky_assertions():
        _run_lossmap(2, 'V', 25000, 0.3)

@retry
def test_lossmap_crystals_B1V():
    with flaky_assertions():
        _run_lossmap(1, 'V', 35000, 0.1, ignore_crystals=False)

@retry
def test_lossmap_crystals_B2H():
    with flaky_assertions():
        _run_lossmap(2, 'H', 30000, 0.15, ignore_crystals=False)



# def test_off_momentum_lossmap():
#     beam           = 1
#     npart          = 200
#     sweep          = 300
#     sweep          = -abs(sweep) if plane == 'DPpos' else abs(sweep)
#     pretrack_turns = 50
#     num_turns      = int(10*abs(sweep))
#     at_element     = 'ip3'

#     line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
#     coll_manager = xc.CollimatorManager.from_yaml(path / f'colldb_lhc_run3.yaml', line=line, beam=beam)

#     coll_manager.install_everest_collimators()
#     coll_manager.build_tracker()
#     coll_manager.set_openings()

#     x_norm, px_norm, _, _ = xp.generate_2D_uniform_circular_sector(npart, r_range=(4, 6))
#     y_norm, py_norm, _, _ = xp.generate_2D_uniform_circular_sector(npart, r_range=(4, 6))
#     zeta, delta = xp.generate_longitudinal_coordinates(
#                         num_particles=num_particles, distribution='gaussian', sigma_z=7.55e-2, line=line
#                 )
#     theta = np.random.uniform(0, np.pi/2, npart)
#     x_norm  *= np.cos(theta)
#     px_norm *= np.cos(theta)
#     y_norm  *= np.sin(theta)
#     py_norm *= np.sin(theta)
#     part = xp.build_particles(
#             x_norm=x_norm, px_norm=px_norm, y_norm=y_norm, py_norm=py_norm, zeta=zeta, delta=delta,
#             nemitt_x=3.5e-6, nemitt_y=3.5e-6, line=line, at_element=at_element
#     )

#     line.optimize_for_tracking(keep_markers=[at_element])
#     idx = line.element_names.index(at_element)
#     part.at_element = idx
#     part.start_tracking_at_element = idx

#     coll_manager.enable_scattering()
#     line.track(part, num_turns=pretrack_turns)
#     part = part.filter(part.state == 1)
#     coll_manager.rf_sweep(sweep=sweep, num_turns=num_turns, particles=part)
#     coll_manager.disable_scattering()

#     tcp3  = f"tcp.6{'l' if beam==1 else 'r'}3.b{beam}"
#     summ = coll_manager.summary(part, show_zeros=False)
#     assert list(summ.columns) == ['collname', 'nabs', 'length', 's', 'type']
#     assert len(summ) == 10
#     # We want at least 5% absorption on the primary
#     assert summ.loc[summ.collname==tcp3,'nabs'].values[0] > 0.05*npart

#     lm = coll_manager.lossmap(part, interpolation=interpolation)
#     assert list(lm.keys()) == ['collimator', 'aperture', 'machine_length', 'interpolation', 'reversed']
#     assert list(lm['collimator'].keys()) == ['s', 'name', 'length', 'n']
#     assert len(lm['collimator']['s']) == len(summ)
#     assert len(lm['collimator']['name']) == len(summ)
#     assert len(lm['collimator']['length']) == len(summ)
#     assert len(lm['collimator']['n']) == len(summ)
#     assert np.all(lm['collimator']['s'] == summ.s)
#     assert np.all(lm['collimator']['name'] == summ.collname)
#     assert np.all(lm['collimator']['length'] == summ.length)
#     assert np.all(lm['collimator']['n'] == summ.nabs)
#     assert np.all([nn[:3] in ['tcp', 'tcs'] for nn in lm['collimator']['name']])
#     assert np.all([s < lm['machine_length'] for s in lm['collimator']['s']])
#     assert list(lm['aperture'].keys()) == ['s', 'name', 'n']
#     assert len(lm['aperture']['s']) > 0
#     assert len(lm['aperture']['s']) == len(lm['aperture']['name'])
#     assert len(lm['aperture']['s']) == len(lm['aperture']['n'])
#     assert np.all([s < lm['machine_length'] for s in lm['aperture']['s']])
#     assert lm['interpolation'] == interpolation
#     line_is_reversed = True if beam==2 else False
#     assert lm['reversed'] == line_is_reversed
