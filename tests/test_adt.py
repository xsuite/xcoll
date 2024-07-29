# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest

import xtrack as xt
import xpart as xp
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

num_turns = 25
num_part = 5000
nemitt_x = 3.5e-6
nemitt_y = 2.5e-6


@retry()
@for_all_test_contexts
@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']],
                         ids=["B1H", "B1V", "B2H", "B2V"])
def test_blow_up(beam, plane, test_context):
    line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / f'lhc_run3_b{beam}_no_aper.json')
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt = xc.BlowUp.install(line, name=name, at_s=adt_pos, need_apertures=False, plane=plane,
                            stop_at_turn=num_turns, use_individual_kicks=True)
    mon = xc.EmittanceMonitor.install(line, name="monitor", at_s=adt_pos, stop_at_turn=num_turns)

    line.build_tracker(_context=test_context)
    tw = line.twiss()
    if plane == 'H':
        adt.calibrate_by_emittance(nemitt=nemitt_x, twiss=tw)
    else:
        adt.calibrate_by_emittance(nemitt=nemitt_y, twiss=tw)
    mon.set_closed_orbit(twiss=tw)

    part = xp.generate_matched_gaussian_bunch(num_particles=num_part, total_intensity_particles=1.6e11,
                                              nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)

    adt.activate()
    line.track(part, num_turns=num_turns, with_progress=1)

    # Verify emittances
    with flaky_assertions():
        if plane == 'H':
            nemitt_expected = 9.37e-6 if beam == 1 else 9.70e-6
            assert abs(mon.nemitt_x[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert abs(mon.nemitt_I[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert mon.nemitt_x.argmin() < 0.05*num_turns
            assert mon.nemitt_x.argmax() > 0.95*num_turns
            # y should not have changed (max 10%):
            assert all([abs(nn-nemitt_y)/nemitt_y < 1.e-1 for nn in mon.nemitt_y])
        else:
            nemitt_expected = 6.06e-6 if beam == 1 else 6.50e-6
            assert abs(mon.nemitt_y[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert abs(mon.nemitt_II[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert mon.nemitt_y.argmin() < 0.05*num_turns
            assert mon.nemitt_y.argmax() > 0.95*num_turns
            # x should not have changed (max 10%):
            assert all([abs(nn-nemitt_x)/nemitt_x < 1.e-1 for nn in mon.nemitt_x])
