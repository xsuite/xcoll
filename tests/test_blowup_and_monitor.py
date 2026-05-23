# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xpart as xp
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts
import xobjects as xo

num_turns = 25
num_part = 5000
nemitt_x = 3.5e-6
nemitt_y = 2.5e-6

path = Path(__file__).parent / 'data'


# TODO: need to test frev and sampling_frequency with coasting beams!


def _assert_monitor(mon, dct={}):
    expected = {
        'particle_id_start': 0,
        'particle_id_stop': -1,
        'start_at_turn': 0,
        'stop_at_turn': 1,
        'horizontal': True,
        'vertical': True,
        'longitudinal': True,
        'sampling_frequency': 1,
        'frev': 1,
        'suppress_warnings': False,
    }
    expected.update(dct)
    for key in expected:
        assert getattr(mon, key) == expected[key]
    size = expected['stop_at_turn'] - expected['start_at_turn']
    size *= expected['sampling_frequency'] / expected['frev']
    for col in ['count', 'cached', 'cached_modes']:
        assert len(mon.data[col]) == size
    for col in ['x_sum1', 'px_sum1', 'x_x_sum2', 'x_px_sum2', 'px_px_sum2']:
        if expected['horizontal']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1
    for col in ['y_sum1', 'py_sum1', 'y_y_sum2', 'y_py_sum2', 'py_py_sum2']:
        if expected['vertical']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1
    for col in ['zeta_sum1', 'pzeta_sum1', 'zeta_zeta_sum2', 'zeta_pzeta_sum2', 'pzeta_pzeta_sum2']:
        if expected['longitudinal']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1
    for col in ['x_y_sum2', 'x_py_sum2', 'px_y_sum2', 'px_py_sum2']:
        if expected['horizontal'] and expected['vertical']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1
    for col in ['x_zeta_sum2', 'x_pzeta_sum2', 'px_zeta_sum2', 'px_pzeta_sum2']:
        if expected['horizontal'] and expected['longitudinal']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1
    for col in ['y_zeta_sum2', 'y_pzeta_sum2', 'py_zeta_sum2', 'py_pzeta_sum2']:
        if expected['vertical'] and expected['longitudinal']:
            assert len(mon.data[col]) == size
        else:
            assert len(mon.data[col]) == 1


@pytest.mark.parametrize("cls", [xc.EmittanceMonitor], ids=["EmittanceMonitor"])
def test_monitor_instance(cls):
    mon = cls()
    _assert_monitor(mon)
    mon = cls(suppress_warnings=True)
    _assert_monitor(mon, {'suppress_warnings': True})
    mon = cls(particle_id_start=7)
    _assert_monitor(mon, {'particle_id_start': 7})
    mon = cls(particle_id_stop=5)
    _assert_monitor(mon, {'particle_id_stop': 5})
    mon = cls(particle_id_start=11, particle_id_stop=21)
    _assert_monitor(mon, {'particle_id_start': 11, 'particle_id_stop': 21})
    mon = cls(particle_id_range=[4, 19])
    _assert_monitor(mon, {'particle_id_start': 4, 'particle_id_stop': 19})
    mon = cls(num_particles=33)
    _assert_monitor(mon, {'particle_id_stop': 33})
    mon = cls(num_particles=33, particle_id_start=5)
    _assert_monitor(mon, {'particle_id_start': 5, 'particle_id_stop': 38})
    mon = cls(num_particles=33, particle_id_stop=35)
    _assert_monitor(mon, {'particle_id_start': 2, 'particle_id_stop': 35})
    with pytest.raises(ValueError):
        mon = cls(particle_id_start=11, particle_id_stop=3)
    with pytest.raises(ValueError):
        mon = cls(num_particles=33, particle_id_stop=5)
    with pytest.raises(ValueError):
        mon = cls(num_particles=33, particle_id_range=[66,999]) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=66) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=[66]) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=[888, 66]) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=[66, 888, 777]) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=[66,999], particle_id_start=5) # this should raise
    with pytest.raises(ValueError):
        mon = cls(particle_id_range=[66,999], particle_id_stop=5) # this should raise

    mon = cls(start_at_turn=5)
    _assert_monitor(mon, {'start_at_turn': 5, 'stop_at_turn': 6})
    mon = cls(stop_at_turn=10)
    _assert_monitor(mon, {'stop_at_turn': 10})
    mon = cls(start_at_turn=5, stop_at_turn=10)
    _assert_monitor(mon, {'start_at_turn': 5, 'stop_at_turn': 10})
    with pytest.raises(ValueError):
        mon = cls(start_at_turn=5, stop_at_turn=4)

    mon = cls(horizontal=True)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': False})
    mon = cls(horizontal=False)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': True})
    mon = cls(vertical=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': False})
    mon = cls(vertical=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': True})
    mon = cls(longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': False, 'longitudinal': True})
    mon = cls(longitudinal=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': True, 'longitudinal': False})
    mon = cls(horizontal=True, vertical=True)
    _assert_monitor(mon, {'horizontal': True, 'vertical': True, 'longitudinal': False})
    mon = cls(horizontal=True, vertical=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': False})
    mon = cls(horizontal=False, vertical=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': False})
    mon = cls(horizontal=False, vertical=False)
    _assert_monitor(mon, {'horizontal': False, 'vertical': False, 'longitudinal': True})
    mon = cls(horizontal=True, longitudinal=True)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': True})
    mon = cls(horizontal=True, longitudinal=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': False})
    mon = cls(horizontal=False, longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': False, 'longitudinal': True})
    mon = cls(horizontal=False, longitudinal=False)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': False})
    mon = cls(vertical=True, longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': True})
    mon = cls(vertical=True, longitudinal=False)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': False})
    mon = cls(vertical=False, longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': False, 'longitudinal': True})
    mon = cls(vertical=False, longitudinal=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': False})
    mon = cls(horizontal=True, vertical=True,longitudinal=True)
    _assert_monitor(mon, {'horizontal': True, 'vertical': True, 'longitudinal': True})
    mon = cls(horizontal=True, vertical=True,longitudinal=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': True, 'longitudinal': False})
    mon = cls(horizontal=True, vertical=False,longitudinal=True)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': True})
    mon = cls(horizontal=True, vertical=False,longitudinal=False)
    _assert_monitor(mon, {'horizontal': True, 'vertical': False, 'longitudinal': False})
    mon = cls(horizontal=False, vertical=True,longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': True})
    mon = cls(horizontal=False, vertical=True,longitudinal=False)
    _assert_monitor(mon, {'horizontal': False, 'vertical': True, 'longitudinal': False})
    mon = cls(horizontal=False, vertical=False,longitudinal=True)
    _assert_monitor(mon, {'horizontal': False, 'vertical': False, 'longitudinal': True})
    with pytest.raises(ValueError):
        mon = cls(horizontal=False, vertical=False,longitudinal=False)


@retry()
@for_all_test_contexts
@pytest.mark.parametrize("aper", [None, "auto", "single", "both"],
                         ids=["without_aper", "auto_aper", "single_aper", "both_aper"])
@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']],
                         ids=["B1H", "B1V", "B2H", "B2V"])
def test_blowup_install(beam, plane, aper, test_context):
    aperture = None
    if aper == 'auto':
        env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
    else:
        env = xt.load(path / f'sequence_lhc_run3_b{beam}_no_aper.json')
        if aper == 'single':
            aperture = xt.LimitEllipse(a=0.01, b=0.01)
        elif aper == 'both':
            aperture = [xt.LimitEllipse(a=0.01, b=0.01), xt.LimitEllipse(a=0.02, b=0.02)]
    need_apertures = aper is not None
    line = env[f'lhcb{beam}']
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    xc.BlowUp.install(line, name=f'{name}_blowup', at=adt_pos, need_apertures=need_apertures,
                      aperture=aperture, plane=plane, stop_at_turn=num_turns,
                      use_individual_kicks=True, _context=test_context)


@retry()
@for_all_test_contexts
@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']],
                         ids=["B1H", "B1V", "B2H", "B2V"])
def test_blowup(beam, plane, test_context):
    env = xt.load(path / f'sequence_lhc_run3_b{beam}_no_aper.json')
    line = env[f'lhcb{beam}']
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt = xc.BlowUp.install(line, name=f'{name}_blowup', at=adt_pos, need_apertures=False, plane=plane,
                            stop_at_turn=num_turns, use_individual_kicks=True, _context=test_context)
    mon = xc.EmittanceMonitor.install(line, name="monitor", at=adt_pos, stop_at_turn=num_turns, _context=test_context)

    line.build_tracker(_context=test_context)
    if plane == 'H':
        adt.calibrate_by_emittance(nemitt=nemitt_x)
    else:
        adt.calibrate_by_emittance(nemitt=nemitt_y)

    part = xp.generate_matched_gaussian_bunch(num_particles=num_part, total_intensity_particles=1.6e11,
                                              nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)

    adt.activate()
    line.track(part, num_turns=num_turns, with_progress=1)

    # Verify emittances
    with flaky_assertions():
        if plane == 'H':
            nemitt_expected = 6.43e-6 if beam == 1 else 6.46e-6
            xo.assert_allclose(mon.nemitt_x[-1], nemitt_expected, rtol=7.5e-2)
            xo.assert_allclose(mon.nemitt_I[-1], nemitt_expected, rtol=7.5e-2)
            assert mon.nemitt_x.argmin() < 0.05*num_turns
            assert mon.nemitt_x.argmax() > 0.95*num_turns
            # y should not have changed (max 10%):
            assert all([abs(nn-nemitt_y)/nemitt_y < 1.e-1 for nn in mon.nemitt_y])
        else:
            nemitt_expected = 4.32e-6 if beam == 1 else 4.49e-6
            xo.assert_allclose(mon.nemitt_y[-1], nemitt_expected, rtol=7.5e-2)
            xo.assert_allclose(mon.nemitt_II[-1], nemitt_expected, rtol=7.5e-2)
            assert mon.nemitt_y.argmin() < 0.05*num_turns
            assert mon.nemitt_y.argmax() > 0.95*num_turns
            # x should not have changed (max 10%):
            assert all([abs(nn-nemitt_x)/nemitt_x < 1.e-1 for nn in mon.nemitt_x])


@for_all_test_contexts
def test_monitor_reset(test_context):
    beam = 1
    plane = 'V'
    env = xt.load(path / f'sequence_lhc_run3_b{beam}_no_aper.json')
    line = env[f'lhcb{beam}']
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt = xc.BlowUp.install(line, name=f'{name}_blowup', at=adt_pos, need_apertures=False, plane=plane,
                            stop_at_turn=num_turns, use_individual_kicks=True)
    mon = xc.EmittanceMonitor.install(line, name="monitor", at=adt_pos, stop_at_turn=num_turns)

    line.build_tracker(_context=test_context)
    if plane == 'H':
        adt.calibrate_by_emittance(nemitt=nemitt_x)
    else:
        adt.calibrate_by_emittance(nemitt=nemitt_y)

    part_init = xp.generate_matched_gaussian_bunch(num_particles=int(num_part/5), total_intensity_particles=1.6e11,
                                                   nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)

    adt.activate()
    for _ in range(2):
        part = part_init.copy()
        line.track(part, num_turns=num_turns, with_progress=1)
        part.move() # Ensure particles are back on CPU context
        mon2 = mon.copy()
        mon2.move()  # Ensure monitor is back on CPU context

        assert np.unique(part.state) == np.array([1])  # Sanity check for test
        assert np.all(mon2.count == int(num_part/5) * np.ones(num_turns))  # Ensure no new particles added

        # Reset monitor for next tracking
        mon.reset()
