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

n_turns = 50
num_part = int(1e5)
nemitt_x = 3.5e-6
nemitt_y = 2.5e-6


# Helper function to calculate the emittance
def calculate_nemitt(part):
    cov_x = np.cov(part.x, part.px)
    cov_y = np.cov(part.y, part.py)
    nemitt_x = part.beta0[0]*part.gamma0[0]*np.sqrt(cov_x[0,0]*cov_x[1,1]-cov_x[1,0]*cov_x[0,1])
    nemitt_y = part.beta0[0]*part.gamma0[0]*np.sqrt(cov_y[0,0]*cov_y[1,1]-cov_y[1,0]*cov_y[0,1])
    return nemitt_x, nemitt_y


@pytest.mark.parametrize("beam, plane", [[1,'H'], [1,'V'], [2,'H'], [2,'V']], ids=["B1H", "B1V", "B2H", "B2V"])
@retry()
def test_blow_up(beam, plane):
    line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / f'lhc_run3_b{beam}_no_aper.json')
    adt = xc.BlowUp(plane=plane, amplitude=0.25)
    pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
    pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
    name = f'adtk{plane.lower()}.{pos}.b{beam}'
    tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
    tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
    adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
    adt.install(line, name=name, at_s=adt_pos, need_apertures=False)

    line.build_tracker()
    line.optimize_for_tracking()
    tw = line.twiss()
    if plane == 'H':
        adt.calibrate_by_emittance(nemitt=nemitt_x, twiss=tw)
    else:
        adt.calibrate_by_emittance(nemitt=nemitt_y, twiss=tw)

    part = xp.generate_matched_gaussian_bunch(num_particles=num_part, total_intensity_particles=1.6e11,
                                              nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)

    ex, ey = calculate_nemitt(part)
    ex = [ex]
    ey = [ey]
    assert abs(ex[0]-nemitt_x)/nemitt_x < 5.e-2
    assert abs(ey[0]-nemitt_y)/nemitt_y < 5.e-2
    part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt_x, nemitt_y=nemitt_y)
    x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
    x_norm_mean = [np.mean(x_norm)]
    x_norm_std  = [np.std(x_norm)]
    y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
    y_norm_mean = [np.mean(y_norm)]
    y_norm_std  = [np.std(y_norm)]
    mean_norm_original = 1.2
    assert abs(x_norm_mean[0]-mean_norm_original)/mean_norm_original < 1.e-1
    assert abs(y_norm_mean[0]-mean_norm_original)/mean_norm_original < 1.e-1

    line.discard_tracker()
    line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))

    adt.activate()

    for i in range(n_turns):
        line.track(part)
        this_ex, this_ey = calculate_nemitt(part)
        ex.append(this_ex)
        ey.append(this_ey)
        part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt_x, nemitt_y=nemitt_y)
        x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
        x_norm_mean.append(np.mean(x_norm))
        y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
        y_norm_mean.append(np.mean(y_norm))
        if i % 10 == 9:
            print(f"Tracked turn {i+1}/{n_turns}.")

    # with flaky_assertions():
    if plane == 'H':
        nemitt_expected = 46.6e-6 if beam == 1 else 47.2e-6
        assert abs(ex[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
        assert np.array(ex).argmin() < 0.05*n_turns
        assert np.array(ex).argmax() > 0.95*n_turns
        # y should not have changed (max 10%):
        assert all([abs(nn-nemitt_y)/nemitt_y < 1.e-1 for nn in ey])
    else:
        nemitt_expected = 29.6e-6 if beam == 1 else 32.3e-6
        assert abs(ey[-1]-nemitt_expected)/nemitt_expected < 7.5e-2
        assert np.array(ey).argmin() < 0.05*n_turns
        assert np.array(ey).argmax() > 0.95*n_turns
        # x should not have changed (max 10%):
        assert all([abs(nn-nemitt_x)/nemitt_x < 1.e-1 for nn in ex])

    if plane == 'H':
        mean_norm_expected = 4.57 if beam == 1 else 4.60
        assert abs(x_norm_mean[-1]-mean_norm_expected)/mean_norm_expected < 1.e-1
        assert np.array(x_norm_mean).argmin() < 0.05*n_turns
        assert np.array(x_norm_mean).argmax() > 0.95*n_turns
        # y should not have changed (max 10%):
        assert all([abs(nn-mean_norm_original)/mean_norm_original < 1.e-1 for nn in y_norm_mean])
    else:
        mean_norm_expected = 4.31 if beam == 1 else 4.50
        assert abs(y_norm_mean[-1]-mean_norm_expected)/mean_norm_expected < 1.e-1
        assert np.array(y_norm_mean).argmin() < 0.05*n_turns
        assert np.array(y_norm_mean).argmax() > 0.95*n_turns
        # x should not have changed (max 10%):
        assert all([abs(nn-mean_norm_original)/mean_norm_original < 1.e-1 for nn in x_norm_mean])
