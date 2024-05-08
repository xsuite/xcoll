# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xcoll as xc
import xtrack as xt
from xtrack.line import _dicts_equal
from xobjects.test_helpers import for_all_test_contexts

import numpy as np
import pytest


# Example values to test
# ======================
#
# This testing script might look a bit unorthodox, as we try to capture all cases
# where we read/write fields of the collimators.
#
# The dictonary 'base_fields contains all xofields of BaseCollimator, with some
# values to set each field to (which will then be verified if correctly written).
#
# The dictionay 'base_dict_fields' contains all fields in the to_dict that are
# not base xofields. Often there are multiple allowed formats for such a to_dict
# field (like is the case for jaw=0.015, jaw=[0.015], jaw=[0.015,-0.015] which all
# represent the same). In 'base_dict_fields', the value for each field is a list
# of allowed formats. Each of those is a dictionary, with 'val' the value to set
# the field to, and 'expected' a dictionary of the xofields that will be updated
# and their expected values.
#
# The dictionary 'base_user_fields' works the same as 'base_dict_fields', but is
# for fields that are not saved in the to_dict but only directly accessible in
# python (like 'angle_L').
#
# The list 'base_user_fields_read_only' is for fields that should not be written
# to: it is verified that an exception is raised.


# BaseCollimator

base_fields = {
    'length':  1.3,
    '_jaw_LU': 0.01544,
    '_jaw_LD': 0.0152,
    '_jaw_RU': -0.0152,
    '_jaw_RD': -0.0152,
    '_sin_zL': np.sin(137.8*np.pi/180),
    '_cos_zL': np.cos(137.8*np.pi/180),
    '_sin_zR': np.sin(125.3*np.pi/180),
    '_cos_zR': np.cos(125.3*np.pi/180),
    '_sin_zDiff': np.sin(-12.5*np.pi/180),
    '_cos_zDiff': np.cos(-12.5*np.pi/180),
    '_sin_yL': -0.00022/1.3,
    '_cos_yL': np.sqrt(1. - 0.00022**2/1.3**2),
    '_tan_yL': -0.00022/1.3/np.sqrt(1. - 0.00022**2/1.3**2),
    '_sin_yR': 0.,
    '_cos_yR': 1.,
    '_tan_yR': 0.,
    '_side':   0,
    'active':  1,
    'record_touches':      0,
    'record_scatterings': 0
}
base_dict_fields = {
    'angle': [
        {'val':  40,        'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180),
                                         '_sin_zR': np.sin(40*np.pi/180), '_cos_zR': np.cos(40*np.pi/180)}},
        {'val':  [40],      'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180),
                                         '_sin_zR': np.sin(40*np.pi/180), '_cos_zR': np.cos(40*np.pi/180)}},
        {'val':  [40, 47],  'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180),
                                         '_sin_zR': np.sin(47*np.pi/180), '_cos_zR': np.cos(47*np.pi/180)}}],
    'jaw': [
        {'val':  0.015,            'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.015, 'jaw_LD': 0.01488, 'jaw_RD': -0.015,
                                                'jaw_L':  0.015,   'jaw_R':  -0.015, 'jaw': [[0.01512, -0.015], [0.01488, -0.015]],
                                                'tilt_L': -0.0001692308, 'tilt_R': 0}},  # Existing tilt from base_fields
        {'val':  [0.015],          'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.015, 'jaw_LD': 0.01488, 'jaw_RD': -0.015,
                                                'jaw_L':  0.015,   'jaw_R':  -0.015, 'jaw': [[0.01512, -0.015], [0.01488, -0.015]],
                                                'tilt_L': -0.0001692308, 'tilt_R': 0}},
        {'val':  [0.015, -0.014],  'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.014, 'jaw_LD': 0.01488, 'jaw_RD': -0.014,
                                                'jaw_L':  0.015,   'jaw_R':  -0.014, 'jaw': [[0.01512, -0.014], [0.01488, -0.014]],
                                                'tilt_L': -0.0001692308, 'tilt_R': 0}},
        {'val':  [[0.015, -0.014],[0.0152, -0.0143]], 
                                   'expected': {'jaw_LU': 0.015,  'jaw_RU': -0.014,   'jaw_LD': 0.0152, 'jaw_RD': -0.0143,
                                                'jaw_L':  0.0151, 'jaw_R':  -0.01415, 'jaw': [[0.015, -0.014], [0.0152, -0.0143]],
                                                'tilt_L': 0.0001538462, 'tilt_R': -0.0002307692}}],
    'tilt': [
        {'val':  5.2e-3,           'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3),
                                                '_sin_yR': np.sin(5.2e-3), '_cos_yR': np.cos(5.2e-3),
                                                '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(5.2e-3),
                                                'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01752998476748726, 
                                                'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.010770015232512739,
                                                'jaw_L':  0.0151, 'jaw_R':  -0.01415}},
        {'val':  [5.2e-3],         'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3),
                                                '_sin_yR': np.sin(5.2e-3), '_cos_yR': np.cos(5.2e-3),
                                                '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(5.2e-3),
                                                'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01752998476748726, 
                                                'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.010770015232512739,
                                                'jaw_L':  0.0151, 'jaw_R':  -0.01415}},
        {'val':  [5.2e-3, 3.7e-3], 'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3),
                                                '_sin_yR': np.sin(3.7e-3), '_cos_yR': np.cos(3.7e-3),
                                                '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(3.7e-3),
                                                'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01655499451259542,
                                                'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.011745005487404576,
                                                'jaw_L':  0.0151, 'jaw_R':  -0.01415}}],
    'gap': [
        {'val':  5,       'expected': {'gap': 5,       'gap_L': 5, 'gap_R': -5}},
        {'val':  [5],     'expected': {'gap': 5,       'gap_L': 5, 'gap_R': -5}},
        {'val':  [5, -3], 'expected': {'gap': [5, -3], 'gap_L': 5, 'gap_R': -3}}],
    'side': [
        {'val':  'left',  'expected': {'_side': 1,  'jaw_LU': 0.01172001523251274, 'jaw_RU': None, '_jaw_RU': -3.0024049945125952, 'gap_L': 5,
                                                    'jaw_LD': 0.018479984767487263,'jaw_RD': None, '_jaw_RD': -2.9975950054874047, 'gap_R': None}},
        {'val':  'both',  'expected': {'_side': 0,  'jaw_LU': 0.01172001523251274, 'jaw_RU': None, '_jaw_RU': -3.0024049945125952, 'gap_L': 5,
                                                    'jaw_LD': 0.018479984767487263,'jaw_RD': None, '_jaw_RD': -2.9975950054874047, 'gap_R': None}},
        {'val':  'right', 'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  'L',     'expected': {'_side': 1,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  'R',     'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  '+',     'expected': {'_side': 1,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  '-',     'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  '+-',    'expected': {'_side': 0,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
        {'val':  '-+',    'expected': {'_side': 0,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, 'gap_L': None,
                                                    '_jaw_LD': 3.0033799847674874, '_jaw_RD': -2.9975950054874047, 'gap_R': None,
                                                    'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}}]
}
base_user_fields = {
    'angle_L': [{'val':  40, 'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180)}}],
    'angle_R': [{'val':  43, 'expected': {'_sin_zR': np.sin(43*np.pi/180), '_cos_zR': np.cos(43*np.pi/180)}}],
    'tilt_L':  [{'val':  0.5*np.pi/180, 'expected': {
                '_sin_yL': np.sin(0.5*np.pi/180), '_cos_yL': np.cos(0.5*np.pi/180), '_tan_yL': np.tan(0.5*np.pi/180)}}],
    'tilt_R':  [{'val':  0.7*np.pi/180, 'expected': {
                '_sin_yR': np.sin(0.7*np.pi/180), '_cos_yR': np.cos(0.7*np.pi/180), '_tan_yR': np.tan(0.7*np.pi/180)}}],
    'jaw_L':   [{'val':  0.013,  'expected': {'jaw_LU': 0.007327751926056947, 'jaw_RU': None, '_jaw_RU': -3.0079410505429107,
                                              'jaw_LD': 0.018672248073943076, 'jaw_RD': None, '_jaw_RD': -2.9920589494570894,
                                              'tilt':  [0.0087266463, 0.0122173048]}}],
    'jaw_R':   [{'val':  -0.011, 'expected': {'jaw_LU': 0.007327751926056947, 'jaw_RU': -0.01894105054291073,
                                              'jaw_LD': 0.018672248073943076, 'jaw_RD': -0.0030589494570893994,
                                              'tilt':  [0.0087266463, 0.0122173048]}}],
    'gap_L':   [{'val':  5,      'expected': {'gap': [5, None],  'gap_L': 5, 'gap_R': None}}],
    'gap_R':   [{'val':  -3,     'expected': {'gap': [5, -3],    'gap_L': 5, 'gap_R': -3}}]
}
base_user_fields_read_only = []


# BlackAbsorber

absorber_fields = base_fields
absorber_dict_fields = base_dict_fields
absorber_user_fields = base_user_fields
absorber_user_fields_read_only = base_user_fields_read_only


# EverestCollimator

everest_fields = {**base_fields,
    '_material':         xc.materials.Copper,
    'rutherford_rng':    xt.RandomRutherford(lower_val=0.002, upper_val=0.98, A=34, B=0.1, Newton_iterations=20),
    '_tracking':         False
}
everest_dict_fields = {**base_dict_fields,
    'material':          [{'val': xc.materials.Copper, 'expected': {'_material': xc.materials.Copper}}]
}
everest_user_fields = base_user_fields
everest_user_fields_read_only = base_user_fields_read_only


# EverestCrystal

everest_crystal_fields = {**everest_fields,
    'align_angle':       3.7e-5,
    '_bending_radius':   61.54,
    '_bending_angle':    50.e-6,
    'xdim':              2.0e-3,
    'ydim':              50.0e-3,
    'thick':             1e-6,
    'miscut':            1.2e-6,
    '_orient':           2
}
everest_crystal_fields['_material'] = xc.materials.TungstenCrystal  # needs to be a CrystalMaterial, else it will fail
everest_crystal_dict_fields = {**everest_dict_fields,
    'lattice': [
        {'val': 'strip',             'expected': {'_orient': 1}},
        {'val': '110',               'expected': {'_orient': 1}},
        {'val': 110,                 'expected': {'_orient': 1}},
        {'val': 'quasi-mosaic',      'expected': {'_orient': 2}},
        {'val': '111',               'expected': {'_orient': 2}},
        {'val': 111,                 'expected': {'_orient': 2}}],
    'bending_radius': [
        {'val': 61.54,               'expected': {'_bending_radius': 61.54, '_bending_angle': 0.02112604331283213}}],
    'bending_angle': [
        {'val': 0.02112604331283213, 'expected': {'_bending_radius': 61.54, '_bending_angle': 0.02112604331283213}}]
}
everest_crystal_dict_fields['material'] = [
    {'val': xc.materials.TungstenCrystal, 'expected': {'_material': xc.materials.TungstenCrystal}}]
everest_crystal_user_fields = everest_user_fields
# everest_crystal_user_fields = {kk: vv for kk, vv in everest_user_fields.items() if '_L' not in kk}
# everest_crystal_user_fields = {kk: [{'val': vv[0]['val'], 'expected': {kkk: vvv for kkk, vvv in vv[0]['expected'].items() if '_L' not in kkk}}]
#                                for kk, vv in everest_user_fields.items() if '_L' not in kk}
everest_crystal_user_fields_read_only = everest_user_fields_read_only



def assert_all_close(expected, setval):
    if hasattr(expected, 'to_dict'):
        assert _dicts_equal(expected.to_dict(), setval.to_dict())
    elif hasattr(expected, '__iter__') and not isinstance(expected, str):
        for exp, stv in zip(expected, setval):
            assert_all_close(exp, stv)
    else:
        if expected is None:
            assert setval is None
        elif isinstance(expected, str):
            assert setval == expected
        else:
            assert np.isclose(expected, setval, atol=1e-12, rtol=0)  


# Tests
# =====
@for_all_test_contexts
def test_black_absorber(test_context):
    # Test instantiation
    elem = xc.BlackAbsorber(length=1, _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in absorber_fields])
    elem.jaw = 0.8 # Need to give non-default value to jaw and gap, otherwise it is None and not in to_dict()
    elem.gap = 5
    assert np.all([key in elem.to_dict() for key in absorber_dict_fields])
    elem.jaw = 1
    elem.gap = None

    # Test reading fields
    for field in list(absorber_fields.keys()) + list(absorber_dict_fields.keys()) \
               + list(absorber_user_fields.keys()) + absorber_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in absorber_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        assert_all_close(val, setval)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for field, vals in {**absorber_dict_fields, **absorber_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            print(val['val'])
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                setval = getattr(elem, basefield)
                assert_all_close(expected, setval)

    # Writing to a read-only field should fail
    for field in absorber_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest_block(test_context):
    # Test instantiation
    elem = xc.EverestBlock(length=1.3, material=xc.materials.Carbon, _context=test_context)
    assert np.isclose(elem.length, 1.3)
    assert elem._tracking == True
    assert xt.line._dicts_equal(elem.material.to_dict(), xc.materials.Carbon.to_dict())
    assert np.isclose(elem.rutherford_rng.lower_val, 0.0009982)
    assert np.isclose(elem.rutherford_rng.upper_val, 0.02)
    assert np.isclose(elem.rutherford_rng.A, 0.0012280392539122623)
    assert np.isclose(elem.rutherford_rng.B, 53.50625)
    assert elem.rutherford_rng.Newton_iterations == 7
    elem.material = xc.materials.Tungsten
    assert np.isclose(elem.rutherford_rng.lower_val, 0.0009982)
    assert np.isclose(elem.rutherford_rng.upper_val, 0.01)
    assert np.isclose(elem.rutherford_rng.A, 0.0018637950841805943)
    assert np.isclose(elem.rutherford_rng.B, 231.48944000000003)
    assert elem.rutherford_rng.Newton_iterations == 7


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest(test_context):
    # Test instantiation
    elem = xc.EverestCollimator(length=1, material=xc.materials.Carbon, _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in everest_fields])
    elem.jaw = 0.8 # Need to give non-default value to jaw and gap, otherwise it is None and not in to_dict()
    elem.gap = 5
    assert np.all([key in elem.to_dict() for key in everest_dict_fields])
    elem.jaw = 1
    elem.gap = None

    # Test reading fields
    for field in list(everest_fields.keys()) + list(everest_dict_fields.keys()) \
               + list(everest_user_fields.keys()) + everest_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in everest_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        assert_all_close(val, setval)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for field, vals in {**everest_dict_fields, **everest_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            print(val['val'])
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                setval = getattr(elem, basefield)
                assert_all_close(expected, setval)

    # Writing to a read-only field should fail
    for field in everest_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)

            
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest_crystal(test_context):
    # Test instantiation
    elem = xc.EverestCrystal(length=1, jaw=0.99, material=xc.materials.SiliconCrystal, side='-', _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in everest_crystal_fields])
    elem.jaw = 0.8 # Need to give non-default value to jaw and gap, otherwise it is None and not in to_dict()
    elem.gap = 5
    assert np.all([key in elem.to_dict() for key in everest_crystal_dict_fields])
    elem.jaw = 1
    elem.gap = None

    # Test reading fields
    for field in list(everest_crystal_fields.keys()) + list(everest_crystal_dict_fields.keys()) \
               + list(everest_crystal_user_fields.keys()) + everest_crystal_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in everest_crystal_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        assert_all_close(val, setval)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    side_set = False
    for field, vals in {**everest_crystal_dict_fields, **everest_crystal_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            if field == 'side' and (val['val'] == '+-' or val['val'] == '-+' or val['val'] == 'both'):
                side_set = True
                continue
            print(val['val'])
            if side_set and field.endswith('_L'):
                continue
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                if not side_set or not field.endswith('_L'):
                    setval = getattr(elem, basefield)
                    assert_all_close(expected, setval)

    # Writing to a read-only field should fail
    for field in everest_crystal_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)

# TODO:
# def test_jaw_func():
