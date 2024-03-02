# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
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
    'inactive_front':    0.7,
    'active_length':     1.3,
    'inactive_back':     0.2,
    'jaw_L':             0.01544,
    'jaw_R':             -0.0152,
    'ref_x':             4.78e-5,
    'ref_y':             -3.2e-5,
    'sin_zL':            np.sin(137.8*np.pi/180),
    'cos_zL':            np.cos(137.8*np.pi/180),
    'sin_zR':            np.sin(42.5*np.pi/180),
    'cos_zR':            np.cos(42.5*np.pi/180),
    'sin_yL':            np.sin(2.78*np.pi/180),
    'cos_yL':            np.cos(2.78*np.pi/180),
    'tan_yL':            np.tan(2.78*np.pi/180),
    'sin_yR':            np.sin(3.5*np.pi/180),
    'cos_yR':            np.cos(3.5*np.pi/180),
    'tan_yR':            np.tan(3.5*np.pi/180),
    '_side':             2,
    'active':            0
}
base_dict_fields = {
    'jaw': [
        {'val':  0.015,            'expected': {'jaw_L': 0.015, 'jaw_R': -0.015}},
        {'val':  [0.015],          'expected': {'jaw_L': 0.015, 'jaw_R': -0.015}},
        {'val':  [0.015, -0.014],  'expected': {'jaw_L': 0.015, 'jaw_R': -0.014}}],
    'angle': [
        {'val':  40,        'expected': {'sin_zL': np.sin(40*np.pi/180), 'cos_zL': np.cos(40*np.pi/180),
                                         'sin_zR': np.sin(40*np.pi/180), 'cos_zR': np.cos(40*np.pi/180)}},
        {'val':  [40],      'expected': {'sin_zL': np.sin(40*np.pi/180), 'cos_zL': np.cos(40*np.pi/180),
                                         'sin_zR': np.sin(40*np.pi/180), 'cos_zR': np.cos(40*np.pi/180)}},
        {'val':  [40, 47],  'expected': {'sin_zL': np.sin(40*np.pi/180), 'cos_zL': np.cos(40*np.pi/180),
                                         'sin_zR': np.sin(47*np.pi/180), 'cos_zR': np.cos(47*np.pi/180)}}],
    'tilt': [
        {'val':  5.2e-2,           'expected': {'sin_yL': np.sin(5.2e-2), 'cos_yL': np.cos(5.2e-2),
                                                'sin_yR': np.sin(5.2e-2), 'cos_yR': np.cos(5.2e-2),
                                                'tan_yL': np.tan(5.2e-2), 'tan_yR': np.tan(5.2e-2)}},
        {'val':  [5.2e-2],         'expected': {'sin_yL': np.sin(5.2e-2), 'cos_yL': np.cos(5.2e-2),
                                                'sin_yR': np.sin(5.2e-2), 'cos_yR': np.cos(5.2e-2),
                                                'tan_yL': np.tan(5.2e-2), 'tan_yR': np.tan(5.2e-2)}},
        {'val':  [5.2e-2, 3.7e-2], 'expected': {'sin_yL': np.sin(5.2e-2), 'cos_yL': np.cos(5.2e-2),
                                                'sin_yR': np.sin(3.7e-2), 'cos_yR': np.cos(3.7e-2),
                                                'tan_yL': np.tan(5.2e-2), 'tan_yR': np.tan(3.7e-2)}}],
    'reference_center': [
        {'val':  0,                  'expected': {'ref_x': 0,       'ref_y': 0}},
        {'val':  [4.78e-5, -3.2e-5], 'expected': {'ref_x': 4.78e-5, 'ref_y': -3.2e-5}}],
    'side': [
        {'val':  'left',  'expected': {'_side': 1}},
        {'val':  'both',  'expected': {'_side': 0}},
        {'val':  'right', 'expected': {'_side': 2}},
        {'val':  'L',     'expected': {'_side': 1}},
        {'val':  'R',     'expected': {'_side': 2}},
        {'val':  '+',     'expected': {'_side': 1}},
        {'val':  '-',     'expected': {'_side': 2}},
        {'val':  '+-',    'expected': {'_side': 0}},
        {'val':  '-+',    'expected': {'_side': 0}}]
}
base_user_fields = {
    'angle_L': [{'val':  40, 'expected': {'sin_zL': np.sin(40*np.pi/180), 'cos_zL': np.cos(40*np.pi/180)}}],
    'angle_R': [{'val':  40, 'expected': {'sin_zR': np.sin(40*np.pi/180), 'cos_zR': np.cos(40*np.pi/180)}}],
    'tilt_L':  [{'val':  0.5*np.pi/180, 'expected': {
                'sin_yL': np.sin(0.5*np.pi/180), 'cos_yL': np.cos(0.5*np.pi/180), 'tan_yL': np.tan(0.5*np.pi/180)}}],
    'tilt_R':  [{'val':  0.7*np.pi/180, 'expected': {
                'sin_yR': np.sin(0.7*np.pi/180), 'cos_yR': np.cos(0.7*np.pi/180), 'tan_yR': np.tan(0.7*np.pi/180)}}],
    # Assuming jaw is previously set to [0.015, -0.014] and tilt to [0.5,0.7] and active_length to 1.3
    # In other words, the current setting is: jaw_LU = 0.009327751926056942
    #                                         jaw_LD = 0.020672248073943057
    #                                         jaw_RU = -0.02194105054291066
    #                                         jaw_RD = -0.00605894945708934
    'jaw_LU':  [{'val':  0.013,  'expected': {'jaw_L':  0.016836124036971527, 'sin_yL': 0.005901729287648506,
                                              'cos_yL': 0.9999825846440603,   'tan_yL': 0.005901832070154704}}],
    'jaw_LD':  [{'val':  0.011,  'expected': {'jaw_L':  0.012,                'sin_yL': -0.0015384615384615385,
                                              'cos_yL': 0.9999988165673471,   'tan_yL': -0.001538463359129313}}],
    'jaw_RU':  [{'val':  -0.015, 'expected': {'jaw_R':  -0.01052947472854467, 'sin_yR': 0.006877731186854354,
                                              'cos_yR': 0.9999763481271551,   'tan_yR': 0.006877893861925417}}],
    'jaw_RD':  [{'val':  -0.019, 'expected': {'jaw_R':  -0.017,               'sin_yR': -0.003076923076923077,
                                              'cos_yR': 0.9999952662609852,   'tan_yR': -0.0030769376423428405}}]
}
base_user_fields_read_only = ['length']


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
everest_crystal_user_fields_read_only = everest_user_fields_read_only



# Tests
# =====
@for_all_test_contexts
def test_black_absorber(test_context):
    # Test instantiation
    elem = xc.BlackAbsorber(length=1, _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in absorber_fields])
    assert np.all([key in elem.to_dict() for key in absorber_dict_fields])

    # Test reading fields
    for field in list(absorber_fields.keys()) + list(absorber_dict_fields.keys()) \
               + list(absorber_user_fields.keys()) + absorber_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in absorber_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        if hasattr(val, 'to_dict'):
            assert _dicts_equal(val.to_dict(), setval.to_dict())
        else:
            assert np.allclose(val, setval, atol=1e-12, rtol=0)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for field, vals in {**absorber_dict_fields, **absorber_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            print(val['val'])
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                setval = getattr(elem, basefield)
                if hasattr(expected, 'to_dict'):
                    assert _dicts_equal(expected.to_dict(), setval.to_dict())
                else:
                    assert np.allclose(expected, setval, atol=1e-12, rtol=0)

    # Writing to a read-only field should fail
    for field in absorber_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)

            
@for_all_test_contexts
def test_black_absorber_length_instantiation(test_context):
    # Test instantiation with different ways of specifying the length
    elem = xc.BlackAbsorber(length=1.1, _context=test_context)
    assert np.allclose(elem.inactive_front, 0, rtol=0)
    assert np.allclose(elem.inactive_back, 0, rtol=0)
    assert np.allclose(elem.active_length, 1.1, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, active_length=0.5, _context=test_context)
    assert np.allclose(elem.inactive_front, 0.3, rtol=0)
    assert np.allclose(elem.inactive_back, 0.3, rtol=0)
    assert np.allclose(elem.active_length, 0.5, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, active_length=0.5, inactive_front=0.4, _context=test_context)
    assert np.allclose(elem.inactive_front, 0.4, rtol=0)
    assert np.allclose(elem.inactive_back, 0.2, rtol=0)
    assert np.allclose(elem.active_length, 0.5, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, active_length=0.5, inactive_back=0.4, _context=test_context)
    assert np.allclose(elem.inactive_front, 0.2, rtol=0)
    assert np.allclose(elem.inactive_back, 0.4, rtol=0)
    assert np.allclose(elem.active_length, 0.5, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, inactive_front=0.4, _context=test_context)
    assert np.allclose(elem.inactive_front, 0.4, rtol=0)
    assert np.allclose(elem.inactive_back, 0., rtol=0)
    assert np.allclose(elem.active_length, 0.7, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, inactive_back=0.4, _context=test_context)
    assert np.allclose(elem.inactive_front, 0., rtol=0)
    assert np.allclose(elem.inactive_back, 0.4, rtol=0)
    assert np.allclose(elem.active_length, 0.7, rtol=0)
    elem = xc.BlackAbsorber(length=1.1, inactive_front=0.4, inactive_back=0.3, _context=test_context)
    assert np.allclose(elem.inactive_front, 0.4, rtol=0)
    assert np.allclose(elem.inactive_back, 0.3, rtol=0)
    assert np.allclose(elem.active_length, 0.4, rtol=0)
    with pytest.raises(Exception) as e_info:
        elem = xc.BlackAbsorber(length=1.1, active_length=0.5, inactive_front=0.4, inactive_back=0.2, _context=test_context)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest(test_context):
    # Test instantiation
    elem = xc.EverestCollimator(length=1, material=xc.materials.Carbon, _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in everest_fields])
    assert np.all([key in elem.to_dict() for key in everest_dict_fields])

    # Test reading fields
    for field in list(everest_fields.keys()) + list(everest_dict_fields.keys()) \
               + list(everest_user_fields.keys()) + everest_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in everest_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        if hasattr(val, 'to_dict'):
            assert _dicts_equal(val.to_dict(), setval.to_dict())
        else:
            assert np.allclose(val, setval, atol=1e-12, rtol=0)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for field, vals in {**everest_dict_fields, **everest_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            print(val['val'])
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                setval = getattr(elem, basefield)
                if hasattr(expected, 'to_dict'):
                    assert _dicts_equal(expected.to_dict(), setval.to_dict())
                else:
                    assert np.allclose(expected, setval, atol=1e-12, rtol=0)

    # Writing to a read-only field should fail
    for field in everest_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)

            
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest_crystal(test_context):
    # Test instantiation
    elem = xc.EverestCrystal(length=1, material=xc.materials.SiliconCrystal, _context=test_context)

    # Test existence of fields
    assert np.all([key in dir(elem) for key in everest_crystal_fields])
    assert np.all([key in elem.to_dict() for key in everest_crystal_dict_fields])

    # Test reading fields
    for field in list(everest_crystal_fields.keys()) + list(everest_crystal_dict_fields.keys()) \
               + list(everest_crystal_user_fields.keys()) + everest_crystal_user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in everest_crystal_fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        if hasattr(val, 'to_dict'):
            assert _dicts_equal(val.to_dict(), setval.to_dict())
        else:
            assert np.allclose(val, setval, atol=1e-12, rtol=0)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for field, vals in {**everest_crystal_dict_fields, **everest_crystal_user_fields}.items():
        print(f"Writing field {field}...")
        for val in vals:
            print(val['val'])
            setattr(elem, field, val['val'])
            for basefield, expected in val['expected'].items():
                setval = getattr(elem, basefield)
                if hasattr(expected, 'to_dict'):
                    assert _dicts_equal(expected.to_dict(), setval.to_dict())
                else:
                    assert np.allclose(expected, setval, atol=1e-12, rtol=0)

    # Writing to a read-only field should fail
    for field in everest_crystal_user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)

# def test_jaw_func():
