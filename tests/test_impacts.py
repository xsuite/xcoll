# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts
import xcoll.constants as xcc


num_part = 10000
num_turns = 3
path = Path(__file__).parent / 'data'


@pytest.mark.everest
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam, plane", [
                            [1, 'H'],
                            [2, 'V'],
                            [1, 'V'],
                            [2, 'H']], ids=["B1H", "B2V", "B1V", "B2H"])
def test_impacts_from_line(beam, plane, test_context):
    env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
    line = env[f'lhcb{beam}']
    colldb = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=beam)
    colldb.install_everest_collimators(verbose=True, line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    impacts = xc.InteractionRecord.start(line=line, record_impacts=True, record_exits=True)
    line.build_tracker(_context=test_context)

    line.xcoll.collimators.assign_optics()
    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    tw = line.twiss()
    part = line[tcp].generate_pencil(num_part, twiss=tw)

    line.xcoll.scattering.enable()
    line.track(part, num_turns=num_turns, time=True, with_progress=1)
    line.xcoll.scattering.disable()

    _assert_impacts(impacts, lengths=line.xcoll.collimators.length)


@pytest.mark.everest
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_impacts_single_collimator(test_context):
    coll = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.MolybdenumGraphite,
                                emittance=3.5e-6, _context=test_context)

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part     = xt.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    impacts = xc.InteractionRecord.start(elements=[coll], names='TCP', record_impacts=True, record_exits=True)
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    _assert_impacts(impacts, lengths=coll.length)


@pytest.mark.everest
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("R, side", [
                            [1, '+'],
                            [-1, '+'],
                            [1, '-'],
                            [-1, '-']], ids=["R>0 side=+ ", "R<0 side=+ ", "R>0 side=- ", "R<0 side=- "])
def test_impacts_single_crystal(R, side, test_context):
    coll = xc.EverestCrystal(length=0.002, material=xc.materials.Silicon, bending_angle=R*149e-6,
                        width=0.002, height=0.05, side=side, lattice='strip', jaw=0.001, _context=test_context)

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part     = xt.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    impacts = xc.InteractionRecord.start(elements=[coll], names='TCPCH', record_impacts=True, record_exits=True)
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    _assert_impacts(impacts, expected_types=['Enter Jaw L', 'Exit Jaw'])


def test_impacts_to_pandas_default_frame():
    impacts, _, _ = _make_impacts_for_frame_tests()

    assert impacts._record_all_columns == 1
    assert impacts.to_pandas().equals(impacts.to_pandas(frame='jaw'))
    with pytest.raises(ValueError, match="Invalid frame"):
        impacts.to_pandas(frame='wrong')


def test_impacts_selected_columns_to_pandas():
    impacts, _, data = _make_impacts_for_frame_tests(
        columns=['particle_id_before', 'x_before'])

    assert impacts._recorded_columns == (
        '_index', 'at_turn', 'at_element', 'shape_id', '_inter',
        'particle_id_before', 'x_before')
    assert impacts._record_all_columns == 0
    assert len(impacts.at_turn) == 2
    assert len(impacts.shape_id) == 2
    assert len(impacts.particle_id_before) == 2
    assert len(impacts.x_before) == 2
    assert len(impacts.s_before) == 0
    assert len(impacts.particle_id_after) == 0

    df = impacts.to_pandas()
    assert list(df.columns) == [
        'turn', 'collimator', 'interaction_type', 'particle_id_before', 'x_before']
    assert np.all(df.particle_id_before.values == data['particle_id_before'])
    assert np.all(df.x_before.values == data['x_before'])

    with pytest.raises(ValueError, match="columns .* were not recorded"):
        impacts.to_pandas(frame='collimator')


def test_impacts_selected_columns_tracking():
    coll = xc.TransparentCollimator(length=0.6, jaw=0.001, name='TCP')
    part = xt.Particles(
        p0c=4e11,
        x=[1.1e-3, 1.2e-3, 1.3e-3],
        px=[0, 0, 0],
        y=[0, 0, 0],
        py=[0, 0, 0],
        delta=[0, 0, 0])

    impacts = xc.InteractionRecord.start(
        elements=[coll], columns=['particle_id_before', 'x_before'],
        record_impacts=True, record_exits=True, capacity=10)
    coll.track(part)

    n_rows = impacts._index.num_recorded
    assert n_rows > 0
    assert impacts.capacity == 10
    assert impacts.io_buffer_capacity >= 10
    assert len(impacts.particle_id_before) == 10
    assert len(impacts.x_before) == 10
    assert len(impacts.particle_id_after) == 0
    assert len(impacts.s_before) == 0

    df = impacts.to_pandas()
    assert len(df) == n_rows
    assert 'particle_id_before' in df.columns
    assert 'x_before' in df.columns
    assert 'particle_id_after' not in df.columns
    assert 's_before' not in df.columns


def test_impacts_selected_columns_unknown():
    coll = xc.TransparentCollimator(length=0.6, jaw=0.001, name='TCP')

    with pytest.raises(ValueError, match="Unknown InteractionRecord columns"):
        xc.InteractionRecord.start(elements=[coll], columns=['not_a_column'])


def test_impacts_to_pandas_collimator_frame():
    impacts, coll, data = _make_impacts_for_frame_tests()
    df = impacts.to_pandas(frame='collimator')
    expected = _expected_frame(data, coll, frame='collimator')

    assert np.all(df.interaction_type.values == ['Enter Jaw L', 'Enter Jaw R'])
    assert np.all(df.collimator.values == 'TCP')
    for coord in ['s', 'x', 'px', 'y', 'py']:
        assert np.allclose(df[f'{coord}_before'], expected[f'{coord}_before'])
        assert np.allclose(df[f'{coord}_after'],  expected[f'{coord}_after'])


def test_impacts_to_pandas_lattice_frame():
    impacts, coll, data = _make_impacts_for_frame_tests()
    df = impacts.to_pandas(frame='lattice')
    expected = _expected_frame(data, coll, frame='lattice')

    for coord in ['s', 'x', 'px', 'y', 'py']:
        assert np.allclose(df[f'{coord}_before'], expected[f'{coord}_before'])
        assert np.allclose(df[f'{coord}_after'],  expected[f'{coord}_after'])


def _assert_impacts(impacts, expected_types=['Enter Jaw L', 'Enter Jaw R', 'Exit Jaw'], lengths=None):
    df = impacts.to_pandas()
    types = np.unique(df.interaction_type)
    assert np.all([type in expected_types for type in types])

    for this_type in ['Enter Jaw L', 'Enter Jaw R']:
        if this_type in expected_types:
            mask = df.interaction_type == this_type
            assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                        np.isclose(df.x_before[mask], 0.0, atol=1e-12))
            mask_all = mask & np.isclose(df.s_before, 0.0, atol=1e-12) & \
                       np.isclose(df.x_before, 0.0, atol=1e-12)
            assert mask_all.sum() < 3 # Allow maximally two particles at the corner
        else:
            assert this_type not in types

    if lengths:
        mask = df.interaction_type == 'Exit Jaw'
        if not isinstance(lengths, dict):
            lengths = {coll: lengths for coll in np.unique(df.collimator[mask])}
        assert np.all(np.isclose(df.s_before[mask], [lengths[coll] for coll in df.collimator[mask]], atol=1e-12) |
                    np.isclose(df.x_before[mask], 0.0, atol=1e-12))


def _make_impacts_for_frame_tests(columns=None):
    coll = xc.TransparentCollimator(length=4., jaw=[[0.2, 1.0], [-0.9, -0.3]],
                                    angle=[30, -45], name='TCP')
    impacts = xc.InteractionRecord.start(elements=[coll], record_impacts=True,
                                         record_exits=True, capacity=2,
                                         columns=columns)
    data = {
        'at_turn':      np.array([0, 1]),
        'at_element':   np.array([0, 0]),
        'shape_id':     np.array([1, -1]),
        '_inter':       np.array([xcc.ENTER_JAW_L,
                                  xcc.ENTER_JAW_R]),
        'particle_id_before':    np.array([10, 11]),
        's_before':     np.array([0.30, 0.45]),
        'x_before':     np.array([0.04, 0.06]),
        'px_before':    np.array([0.010, 0.020]),
        'y_before':     np.array([0.07, -0.08]),
        'py_before':    np.array([0.005, -0.006]),
        'delta_before': np.array([0.10, -0.20]),
        'particle_id_after':     np.array([-1, 11]),
        's_after':      np.array([-1., 1.20]),
        'x_after':      np.array([-1., 0.10]),
        'px_after':     np.array([-1., -0.015]),
        'y_after':      np.array([-1., 0.11]),
        'py_after':     np.array([-1., 0.012]),
        'delta_after':  np.array([-1., 0.25]),
    }
    impacts._index.num_recorded = len(data['shape_id'])
    for field, values in data.items():
        if impacts._column_is_recorded(field):
            for ii, val in enumerate(values):
                getattr(impacts, field)[ii] = val
    return impacts, coll, data


def _expected_frame(data, coll, frame):
    expected = {}
    for at in ['before', 'after']:
        s, x, px, y, py = _expected_collimator_frame(data, coll, at)
        if frame == 'lattice':
            x, px, y, py = _expected_lattice_frame(data, coll, x, px, y, py)
        expected[f's_{at}']  = s
        expected[f'x_{at}']  = x
        expected[f'px_{at}'] = px
        expected[f'y_{at}']  = y
        expected[f'py_{at}'] = py
    return expected


def _expected_collimator_frame(data, coll, at):
    s     = data[f's_{at}'].copy()
    x     = data[f'x_{at}'].copy()
    px    = data[f'px_{at}'].copy()
    y     = data[f'y_{at}'].copy()
    py    = data[f'py_{at}'].copy()
    delta = data[f'delta_{at}']

    for ii, shape_id in enumerate(data['shape_id']):
        if shape_id >= 0:
            sin_y, cos_y, tilt, jaw = coll._sin_yL, coll._cos_yL, coll.tilt_L, coll.jaw_LU
        else:
            sin_y, cos_y, tilt, jaw = coll._sin_yR, coll._cos_yR, coll.tilt_R, coll.jaw_RU
            if x[ii] != -1:
                x[ii] = -x[ii]
            if px[ii] != -1:
                px[ii] = -px[ii]
        if s[ii] != -1 and x[ii] != -1:
            old_s = s[ii]
            old_x = x[ii]
            s[ii] = old_s*cos_y - old_x*sin_y - coll.length/2*(1 - cos_y)
            x[ii] = old_s*sin_y + old_x*cos_y
        if x[ii] != -1:
            x[ii] += jaw
        if px[ii] != -1:
            px[ii] += tilt*(1 + delta[ii])
    return s, x, px, y, py


def _expected_lattice_frame(data, coll, x, px, y, py):
    x  = x.copy()
    px = px.copy()
    y  = y.copy()
    py = py.copy()
    for ii, shape_id in enumerate(data['shape_id']):
        if shape_id >= 0:
            sin_z, cos_z = coll._sin_zL, coll._cos_zL
        else:
            sin_z, cos_z = coll._sin_zR, coll._cos_zR
        if x[ii] != -1 and y[ii] != -1:
            old_x = x[ii]
            old_y = y[ii]
            x[ii] = old_x*cos_z - old_y*sin_z
            y[ii] = old_x*sin_z + old_y*cos_z
        if px[ii] != -1 and py[ii] != -1:
            old_px = px[ii]
            old_py = py[ii]
            px[ii] = old_px*cos_z - old_py*sin_z
            py[ii] = old_px*sin_z + old_py*cos_z
    return x, px, y, py
