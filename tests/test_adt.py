# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

from xpart.test_helpers import flaky_assertions, retry

n_turns = 500
num_part = 5000
nemitt = 3.5e-6

# Helper function to calculate the emittance
def _calculate_nemitt(part):
    cov_x = np.cov(part.x, part.px)
    cov_y = np.cov(part.y, part.py)
    nemitt_x = part.gamma0[0]*np.sqrt(cov_x[0,0]*cov_x[1,1]-cov_x[1,0]*cov_x[0,1])
    nemitt_y = part.gamma0[0]*np.sqrt(cov_y[0,0]*cov_y[1,1]-cov_y[1,0]*cov_y[0,1])
    return nemitt_x, nemitt_y

@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']], ids=["B1H", "B1V", "B2H", "B2V"])
@retry()
def test_blow_up(beam, plane):
    line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / f'lhc_run3_b{beam}_no_aper.json')

    adt = xc.ADT(plane=plane, amplitude=1)
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt.install(line, name=name, at_s=adt_pos, need_apertures=False)

    line.build_tracker()
    line.optimize_for_tracking()
    part = xp.generate_matched_gaussian_bunch(num_particles=num_part, total_intensity_particles=1.6e11,
                                      nemitt_x=nemitt, nemitt_y=nemitt, sigma_z=7.55e-2, line=line)

    ex, ey = _calculate_nemitt(part)
    nemitt_x = [ex]
    nemitt_y = [ey]
    assert abs(nemitt_x[0]-nemitt)/nemitt < 5.e-2
    assert abs(nemitt_y[0]-nemitt)/nemitt < 5.e-2

    tw = line.twiss()
    part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt, nemitt_y=nemitt)
    x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
    x_norm_mean = [np.mean(x_norm)]
    x_norm_std  = [np.std(x_norm)]
    y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
    y_norm_mean = [np.mean(y_norm)]
    y_norm_std  = [np.std(y_norm)]
    mean_norm_original = 1.2
    assert abs(x_norm_mean[0]-mean_norm_original)/mean_norm_original < 1.e-1
    assert abs(y_norm_mean[0]-mean_norm_original)/mean_norm_original < 1.e-1

    adt.calibrate(twiss=tw, nemitt_x=nemitt, nemitt_y=nemitt)
    adt.activate()

    line.discard_tracker()
    line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))

    for _ in range(n_turns):
        line.track(part)
        ex, ey = _calculate_nemitt(part)
        nemitt_x.append(ex)
        nemitt_y.append(ey)
        part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt, nemitt_y=nemitt)
        x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
        x_norm_mean.append(np.mean(x_norm))
        x_norm_std.append(np.std(x_norm))
        y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
        y_norm_mean.append(np.mean(y_norm))
        y_norm_std.append(np.std(y_norm))

    with flaky_assertions():
        nemitt_expected = 45.e-6
        if beam == 1 and plane == 'V':
            nemitt_expected = 42.e-6
        if plane == 'H':
            assert abs(nemitt_x[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert np.array(nemitt_x).argmin() < 0.05*n_turns
            assert np.array(nemitt_x).argmax() > 0.95*n_turns
            # y should not have changed (max 10%):
            assert all([abs(nn-nemitt)/nemitt < 1.e-1 for nn in nemitt_y])
        else:
            assert abs(nemitt_y[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
            assert np.array(nemitt_y).argmin() < 0.05*n_turns
            assert np.array(nemitt_y).argmax() > 0.95*n_turns
            # x should not have changed (max 10%):
            assert all([abs(nn-nemitt)/nemitt < 1.e-1 for nn in nemitt_x])

        mean_norm_expected = 4.5
        if beam == 1 and plane == 'V':
            mean_norm_expected = 4.2
        if plane == 'H':
            assert abs(x_norm_mean[-1]-mean_norm_expected)/mean_norm_expected < 1.e-1
            assert np.array(x_norm_mean).argmin() < 0.05*n_turns
            assert np.array(x_norm_mean).argmax() > 0.95*n_turns
            # y should not have changed (max 10%):
            assert all([abs(nn-mean_norm_original)/mean_norm_original < 1.e-1 for nn in y_norm_mean])
        else:
            assert abs(y_norm_mean[-1]-mean_norm_expected)/mean_norm_expected < 1.e-1
            assert np.array(y_norm_mean).argmin() < 0.05*n_turns
            assert np.array(y_norm_mean).argmax() > 0.95*n_turns
            # x should not have changed (max 10%):
            assert all([abs(nn-mean_norm_original)/mean_norm_original < 1.e-1 for nn in x_norm_mean])
