# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
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
path = xc._pkg_root.parent / 'examples' / 'machines'


@retry()
@for_all_test_contexts
@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']],
                         ids=["B1H", "B1V", "B2H", "B2V"])
def test_blow_up(beam, plane, test_context):
    line = xt.Line.from_json(path / f'lhc_run3_b{beam}_no_aper.json')
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt = xc.BlowUp.install(line, name=f'{name}_blowup', at_s=adt_pos, need_apertures=False, plane=plane,
                            stop_at_turn=num_turns, use_individual_kicks=True, _context=test_context)
    mon = xc.EmittanceMonitor.install(line, name="monitor", at_s=adt_pos, stop_at_turn=num_turns, _context=test_context)

    line.build_tracker(_context=test_context)
    if plane == 'H':
        adt.calibrate_by_emittance(nemitt=nemitt_x)
    else:
        adt.calibrate_by_emittance(nemitt=nemitt_y)

    part = xp.generate_matched_gaussian_bunch(num_particles=num_part, total_intensity_particles=1.6e11,
                                              nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line,
                                              _context=test_context)

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

    # Quick test to check storing
    file = Path(f'monitor_{beam}{plane}{test_context}.json')
    assert not file.exists()
    mon.to_json(file)
    assert file.exists()
    mon2 = xc.EmittanceMonitor.from_json(file)
    dct1 = mon.to_dict()
    dct2 = mon2.to_dict()
    # TODO: this compares ArrNFloat64 ids, however, should do to_nparray()
    # assert xt.line._dicts_equal(dct1, dct2)
    assert set(dct1.keys()) == set(dct2.keys())
    assert set(dct1['data'].keys()) == set(dct2['data'].keys())
    for kk in dct1['data'].keys():
        assert np.allclose(dct1['data'][kk].to_nparray(), dct2['data'][kk].to_nparray(),
                           rtol=1e-15, atol=1e-15)
    dct1.pop('data')
    dct2.pop('data')
    assert xt.line._dicts_equal(dct1, dct2)

    file.unlink()
