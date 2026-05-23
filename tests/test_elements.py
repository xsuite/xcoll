# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt
from xobjects.test_helpers import for_all_test_contexts

import xcoll as xc
from xcoll.compare import deep_equal

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

# BaseBlock
base_fields = {
    'length':  1.3,
    'active':  1,
    '_record_interactions': 0
}
base_dict_fields = [
    {'field': 'record_impacts',     'val':  True,  'expected': {'record_impacts': True,  'record_exits': False, 'record_scatterings': False, '_record_interactions': 1}},
    {'field': 'record_impacts',     'val':  False, 'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': False, '_record_interactions': 0}},
    {'field': 'record_exits',       'val':  True,  'expected': {'record_impacts': False, 'record_exits': True,  'record_scatterings': False, '_record_interactions': 2}},
    {'field': 'record_exits',       'val':  False, 'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': False, '_record_interactions': 0}},
    {'field': 'record_scatterings', 'val':  True,  'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': True,  '_record_interactions': 4}},
    {'field': 'record_scatterings', 'val':  False, 'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': False, '_record_interactions': 0}},
    {'field': 'record_scatterings', 'val':  True,  'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': True,  '_record_interactions': 4}},
    {'field': 'record_impacts',     'val':  True,  'expected': {'record_impacts': True,  'record_exits': False, 'record_scatterings': True,  '_record_interactions': 5}},
    {'field': 'record_impacts',     'val':  False, 'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': True,  '_record_interactions': 4}},
    {'field': 'record_exits',       'val':  True,  'expected': {'record_impacts': False, 'record_exits': True,  'record_scatterings': True,  '_record_interactions': 6}},
    {'field': 'record_exits',       'val':  False, 'expected': {'record_impacts': False, 'record_exits': False, 'record_scatterings': True,  '_record_interactions': 4}},
    {'field': 'record_exits',       'val':  True,  'expected': {'record_impacts': False, 'record_exits': True,  'record_scatterings': True,  '_record_interactions': 6}},
    {'field': 'record_scatterings', 'val':  False, 'expected': {'record_impacts': False, 'record_exits': True,  'record_scatterings': False, '_record_interactions': 2}},
    {'field': 'record_impacts',     'val':  True,  'expected': {'record_impacts': True,  'record_exits': True,  'record_scatterings': False, '_record_interactions': 3}},
    {'field': 'record_impacts',     'val':  False, 'expected': {'record_impacts': False, 'record_exits': True,  'record_scatterings': False, '_record_interactions': 2}},
    {'field': 'record_impacts',     'val':  True,  'expected': {'record_impacts': True,  'record_exits': True,  'record_scatterings': False, '_record_interactions': 3}},
    {'field': 'record_exits',       'val':  False, 'expected': {'record_impacts': True,  'record_exits': False, 'record_scatterings': False, '_record_interactions': 1}},
    {'field': 'record_scatterings', 'val':  True,  'expected': {'record_impacts': True,  'record_exits': False, 'record_scatterings': True,  '_record_interactions': 5}},
    {'field': 'record_exits',       'val':  True,  'expected': {'record_impacts': True,  'record_exits': True,  'record_scatterings': True,  '_record_interactions': 7}}
]
base_user_fields = {}
base_user_fields_read_only = ['gemitt_x', 'gemitt_y']

# BaseCollimator
base_coll_fields = {**base_fields,
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
    '_jaws_parallel': False,
    '_sin_yL': -0.00022/1.3,
    '_cos_yL': np.sqrt(1. - 0.00022**2/1.3**2),
    '_tan_yL': -0.00022/1.3/np.sqrt(1. - 0.00022**2/1.3**2),
    '_sin_yR': 0.,
    '_cos_yR': 1.,
    '_tan_yR': 0.,
    '_side':   0,
    '_align':    0,
    '_gap_L':    3,
    '_gap_R':    4,
    '_nemitt_x': 3.5e-6,
    '_nemitt_y': 2.5e-6
}
base_coll_dict_fields = [*base_dict_fields,
    {'field': 'angle', 'val':  40,       'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180), '_sin_zR': np.sin(40*np.pi/180), '_cos_zR': np.cos(40*np.pi/180)}},
    {'field': 'angle', 'val':  [40],     'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180), '_sin_zR': np.sin(40*np.pi/180), '_cos_zR': np.cos(40*np.pi/180)}},
    {'field': 'angle', 'val':  [40, 47], 'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180), '_sin_zR': np.sin(47*np.pi/180), '_cos_zR': np.cos(47*np.pi/180)}},
    {'field': 'jaw',   'val':  0.015,                                'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.015, 'jaw_LD': 0.01488, 'jaw_RD': -0.015,  'jaw_L':  0.015,  'jaw_R':  -0.015,   'jaw': [[0.01512, 0.01488], [-0.015, -0.015]], 'tilt_L': -0.0001692308, 'tilt_R': 0}},  # Existing tilt from base_fields
    {'field': 'jaw',   'val':  [0.015],                              'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.015, 'jaw_LD': 0.01488, 'jaw_RD': -0.015,  'jaw_L':  0.015,  'jaw_R':  -0.015,   'jaw': [[0.01512, 0.01488], [-0.015, -0.015]], 'tilt_L': -0.0001692308, 'tilt_R': 0}},
    {'field': 'jaw',   'val':  [0.015, -0.014],                      'expected': {'jaw_LU': 0.01512, 'jaw_RU': -0.014, 'jaw_LD': 0.01488, 'jaw_RD': -0.014,  'jaw_L':  0.015,  'jaw_R':  -0.014,   'jaw': [[0.01512, 0.01488], [-0.014, -0.014]], 'tilt_L': -0.0001692308, 'tilt_R': 0}},
    {'field': 'jaw',   'val':  [[0.015, 0.0152], [-0.014, -0.0143]], 'expected': {'jaw_LU': 0.015,   'jaw_RU': -0.014, 'jaw_LD': 0.0152,  'jaw_RD': -0.0143, 'jaw_L':  0.0151, 'jaw_R':  -0.01415, 'jaw': [[0.015, 0.0152], [-0.014, -0.0143]],   'tilt_L': 0.0001538462,  'tilt_R': -0.0002307692}},
    {'field': 'tilt',  'val':  5.2e-3,           'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3), '_sin_yR': np.sin(5.2e-3), '_cos_yR': np.cos(5.2e-3), '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(5.2e-3), 'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01752998476748726, 'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.010770015232512739, 'jaw_L':  0.0151, 'jaw_R':  -0.01415}},
    {'field': 'tilt',  'val':  [5.2e-3],         'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3), '_sin_yR': np.sin(5.2e-3), '_cos_yR': np.cos(5.2e-3), '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(5.2e-3), 'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01752998476748726, 'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.010770015232512739, 'jaw_L':  0.0151, 'jaw_R':  -0.01415}},
    {'field': 'tilt',  'val':  [5.2e-3, 3.7e-3], 'expected': {'_sin_yL': np.sin(5.2e-3), '_cos_yL': np.cos(5.2e-3), '_sin_yR': np.sin(3.7e-3), '_cos_yR': np.cos(3.7e-3), '_tan_yL': np.tan(5.2e-3), '_tan_yR': np.tan(3.7e-3), 'jaw_LU': 0.01172001523251274,  'jaw_RU': -0.01655499451259542, 'jaw_LD': 0.018479984767487263, 'jaw_RD': -0.011745005487404576, 'jaw_L':  0.0151, 'jaw_R':  -0.01415}},
    {'field': 'gap',   'val':  5,       'expected': {'gap': 5,       'gap_L': 5, 'gap_R': -5}},
    {'field': 'gap',   'val':  [5],     'expected': {'gap': 5,       'gap_L': 5, 'gap_R': -5}},
    {'field': 'gap',   'val':  [5, -3], 'expected': {'gap': [5, -3], 'gap_L': 5, 'gap_R': -3}},
    {'field': 'side',  'val':  'left',  'expected': {'_side': 1,  'jaw_LU': 0.01172001523251274, '_jaw_RU': -3.0024049945125952, 'jaw_LD': 0.018479984767487263, '_jaw_RD': -2.9975950054874047, 'gap_L': 5,    'gap_R': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  'both',  'expected': {'_side': 0,  'jaw_LU': 0.01172001523251274, '_jaw_RU': -3.0024049945125952, 'jaw_LD': 0.018479984767487263, '_jaw_RD': -2.9975950054874047, 'gap_L': 5,    'gap_R': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  'right', 'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  'L',     'expected': {'_side': 1,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  'R',     'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  '+',     'expected': {'_side': 1,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  '-',     'expected': {'_side': -1, '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  '+-',    'expected': {'_side': 0,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'side',  'val':  '-+',    'expected': {'_side': 0,  '_jaw_LU': 2.9966200152325128, '_jaw_RU': -3.0024049945125952, '_jaw_LD': 3.0033799847674874,  '_jaw_RD': -2.9975950054874047, 'gap_L': None, 'gap_R': None, 'jaw_LU': None, 'jaw_LD': None, 'jaw_RU': None, 'jaw_RD': None}}
]
base_coll_user_fields = [*base_user_fields,
    {'field': 'angle_L',   'val':  40,               'expected': {'_sin_zL': np.sin(40*np.pi/180), '_cos_zL': np.cos(40*np.pi/180)}},
    {'field': 'angle_R',   'val':  43,               'expected': {'_sin_zR': np.sin(43*np.pi/180), '_cos_zR': np.cos(43*np.pi/180)}},
    {'field': 'tilt_L',    'val':  0.5*np.pi/180,    'expected': {'_sin_yL': np.sin(0.5*np.pi/180), '_cos_yL': np.cos(0.5*np.pi/180), '_tan_yL': np.tan(0.5*np.pi/180)}},
    {'field': 'tilt_R',    'val':  0.7*np.pi/180,    'expected': {'_sin_yR': np.sin(0.7*np.pi/180), '_cos_yR': np.cos(0.7*np.pi/180), '_tan_yR': np.tan(0.7*np.pi/180)}},
    {'field': 'jaw_L',     'val':  0.013,            'expected': {'jaw_LU': 0.007327751926056947, '_jaw_RU': -3.0079410505429107, 'jaw_LD': 0.018672248073943076, '_jaw_RD': -2.9920589494570894,   'tilt':  [0.0087266463, 0.0122173048], 'jaw_RU': None, 'jaw_RD': None}},
    {'field': 'jaw_R',     'val':  -0.011,           'expected': {'jaw_LU': 0.007327751926056947, 'jaw_RU': -0.01894105054291073, 'jaw_LD': 0.018672248073943076, 'jaw_RD': -0.0030589494570893994, 'tilt':  [0.0087266463, 0.0122173048]}},
    {'field': 'gap_L',     'val':  5,                'expected': {'gap': [5, None],  'gap_L': 5, 'gap_R': None}},
    {'field': 'gap_R',     'val':  -3,               'expected': {'gap': [5, -3],    'gap_L': 5, 'gap_R': -3}},
    {'field': 'align',     'val':  'upstream',       'expected': {'align': 'upstream'}},
    {'field': 'align',     'val':  'downstream',     'expected': {'align': 'downstream'}},
    {'field': 'emittance', 'val':  3.5e-6,           'expected': {'emittance': 3.5e-6,  'nemitt_x': 3.5e-6,  'nemitt_y': 3.5e-6}},
    {'field': 'emittance', 'val':  [3.25e-6],        'expected': {'emittance': 3.25e-6, 'nemitt_x': 3.25e-6, 'nemitt_y': 3.25e-6}},
    {'field': 'emittance', 'val':  [3.0e-6, 3.0e-6], 'expected': {'emittance': 3.0e-6,  'nemitt_x': 3.0e-6,  'nemitt_y': 3.0e-6}},
    {'field': 'emittance', 'val':  [2.0e-6, 2.5e-6], 'expected': {'emittance': [2.0e-6, 2.5e-6], 'nemitt_x': 2.0e-6, 'nemitt_y': 2.5e-6}},
    {'field': 'nemitt_x',  'val':  1.5e-6,           'expected': {'emittance': [1.5e-6, 2.5e-6], 'nemitt_x': 1.5e-6, 'nemitt_y': 2.5e-6}},
    {'field': 'nemitt_y',  'val':  1.5e-6,           'expected': {'emittance': 1.5e-6, 'nemitt_x': 1.5e-6, 'nemitt_y': 1.5e-6}}
]
base_coll_user_fields_read_only = base_user_fields_read_only

# BaseCrystal
base_crystal_fields = {**base_fields,
    '_jaw_U': 0.01544,
    '_sin_z': np.sin(137.8*np.pi/180),
    '_cos_z': np.cos(137.8*np.pi/180),
    '_sin_y': -0.00022/1.3,
    '_cos_y': np.sqrt(1. - 0.00022**2/1.3**2),
    '_tan_y': -0.00022/1.3/np.sqrt(1. - 0.00022**2/1.3**2),
    '_side':  1,
    '_align': 0,
    '_gap':   3,
    '_nemitt_x': 3.5e-6,
    '_nemitt_y': 2.5e-6,
    '_bending_radius':  12,
    '_bending_angle':   0.10854636232780158,
    'width':            0.05,
    'height':           0.02
}
base_crystal_dict_fields = [*base_dict_fields,
    {'field': 'angle',          'val':  40,      'expected': {'_sin_z': np.sin(40*np.pi/180), '_cos_z': np.cos(40*np.pi/180)}},
    {'field': 'jaw',            'val':  0.015,   'expected': {'jaw': 0.015, 'jaw_U': 0.015, 'jaw_D': 0.08540449144429954}},
    {'field': 'tilt',           'val':  5.2e-3,  'expected': {'_sin_y': np.sin(5.2e-3), '_cos_y': np.cos(5.2e-3), '_tan_y': np.tan(5.2e-3), 'jaw_U': 0.015,  'jaw_D': 0.09238350714959694}},
    {'field': 'gap',            'val':  5,       'expected': {'gap': 5, '_gap': 5}},
    {'field': 'side',           'val':  'left',  'expected': {'_side': 1,  'jaw_D': 0.09238350714959694}},
    {'field': 'side',           'val':  'right', 'expected': {'_side': -1, 'jaw_D': 0.09206107586980695}},
    {'field': 'side',           'val':  'L',     'expected': {'_side': 1,  'jaw_D': 0.09238350714959694}},
    {'field': 'side',           'val':  'R',     'expected': {'_side': -1, 'jaw_D': 0.09206107586980695}},
    {'field': 'side',           'val':  '+',     'expected': {'_side': 1,  'jaw_D': 0.09238350714959694}},
    {'field': 'side',           'val':  '-',     'expected': {'_side': -1, 'jaw_D': 0.09206107586980695}},
    {'field': 'bending_radius', 'val':  47,      'expected': {'_bending_radius': 47.0,              '_bending_angle': 0.027663102518570612, 'jaw_U': 0.015, 'jaw_D': 0.039715568642416266}},
    {'field': 'bending_angle',  'val':  0.07,    'expected': {'_bending_radius': 18.58660391285491, '_bending_angle': 0.07,                 'jaw_U': 0.015, 'jaw_D': 0.06713730900985289}}
]
base_crystal_user_fields = [*base_user_fields,
    {'field': 'jaw_U',     'val':  0.013,            'expected': {'jaw': 0.013, 'jaw_U': 0.013, 'jaw_D': 0.06513730900985289}},
    {'field': 'align',     'val':  'upstream',       'expected': {'align': 'upstream'}},
    {'field': 'emittance', 'val':  3.5e-6,           'expected': {'emittance': 3.5e-6,           'nemitt_x': 3.5e-6,  'nemitt_y': 3.5e-6}},
    {'field': 'emittance', 'val':  [3.25e-6],        'expected': {'emittance': 3.25e-6,          'nemitt_x': 3.25e-6, 'nemitt_y': 3.25e-6}},
    {'field': 'emittance', 'val':  [3.0e-6, 3.0e-6], 'expected': {'emittance': 3.0e-6,           'nemitt_x': 3.0e-6,  'nemitt_y': 3.0e-6}},
    {'field': 'emittance', 'val':  [2.0e-6, 2.5e-6], 'expected': {'emittance': [2.0e-6, 2.5e-6], 'nemitt_x': 2.0e-6,  'nemitt_y': 2.5e-6}},
    {'field': 'nemitt_x',  'val':  1.5e-6,           'expected': {'emittance': [1.5e-6, 2.5e-6], 'nemitt_x': 1.5e-6,  'nemitt_y': 2.5e-6}},
    {'field': 'nemitt_y',  'val':  1.5e-6,           'expected': {'emittance': 1.5e-6,           'nemitt_x': 1.5e-6,  'nemitt_y': 1.5e-6}}
]
base_crystal_user_fields_read_only = base_user_fields_read_only


# BlackAbsorber
absorber_fields = base_coll_fields
absorber_dict_fields = base_coll_dict_fields
absorber_user_fields = base_coll_user_fields
absorber_user_fields_read_only = base_coll_user_fields_read_only

# BlackCrystal
black_crystal_fields = base_crystal_fields
black_crystal_dict_fields = base_crystal_dict_fields
black_crystal_user_fields = base_crystal_user_fields
black_crystal_user_fields_read_only = base_crystal_user_fields_read_only

# EverestCollimator
everest_fields = {**base_coll_fields,
    '_material':         xc.materials.Copper,
    'rutherford_rng':    xt.RandomRutherford(lower_val=0.002, upper_val=0.98, A=34, B=0.1, Newton_iterations=20),
    '_tracking':         False
}
everest_dict_fields = [*base_coll_dict_fields,
    {'field': 'material', 'val': xc.materials.Copper, 'expected': {'_material': xc.materials.Copper}}
]
everest_user_fields = base_coll_user_fields
everest_user_fields_read_only = base_coll_user_fields_read_only

# EverestCrystal
everest_crystal_fields = {**base_crystal_fields,
    'miscut':            1.2e-6,
    '_orient':           2,
    '_critical_angle':   1.3e-6,
    '_material':         xc.materials.Tungsten,
    'rutherford_rng':    xt.RandomRutherford(lower_val=0.002, upper_val=0.98, A=34, B=0.1, Newton_iterations=20),
    '_tracking':         False
}
everest_crystal_dict_fields = [*base_crystal_dict_fields,
    {'field': 'material', 'val': xc.materials.Silicon, 'expected': {'_material': xc.materials.Silicon}},
    {'field': 'lattice', 'val': 'strip',        'expected': {'_orient': 1}},
    {'field': 'lattice', 'val': '110',          'expected': {'_orient': 1}},
    {'field': 'lattice', 'val': 110,            'expected': {'_orient': 1}},
    {'field': 'lattice', 'val': 'quasi-mosaic', 'expected': {'_orient': 2}},
    {'field': 'lattice', 'val': '111',          'expected': {'_orient': 2}},
    {'field': 'lattice', 'val': 111,            'expected': {'_orient': 2}}
]
everest_crystal_user_fields = base_crystal_user_fields
everest_crystal_user_fields_read_only = [*base_crystal_user_fields_read_only, 'critical_radius', 'critical_angle']

# FlukaCollimator
fluka_fields = {**base_coll_fields,
    'fluka_id':              137,
    'length_front':          0.4125,
    'length_back':           0.738,
    '_tracking':             1,
    '_acc_ionisation_loss':  4.349234923e18
}
fluka_dict_fields = [*base_coll_dict_fields[:-9],
    {'field': 'assembly', 'val': 'hilumi_tcppm',
        'expected': {'material': None, 'height': None, 'width': None, 'side': None}}
]
fluka_generic_dict_fields = [*base_coll_dict_fields,
    {'field': 'material', 'val': xc.materials.Yttrium, 'expected': {'material': xc.materials.Yttrium}},
    {'field': 'height', 'val': 0.03, 'expected': {'height': 0.03}},
    {'field': 'width',  'val': 0.02, 'expected': {'width': 0.02}},
    {'field': 'side',   'val': 'right', 'expected': {'side': 'right'}},
]

# FlukaCrystal
fluka_crystal_fields = {**base_crystal_fields,
    'fluka_id':              137,
    'length_front':          0.4125,
    'length_back':           0.738,
    '_tracking':             1,
    '_acc_ionisation_loss':  4.349234923e18
}
fluka_crystal_dict_fields = [*base_crystal_dict_fields,
    {'field': 'material', 'val': xc.materials.Yttrium, 'expected': {'material': xc.materials.Yttrium}},
    {'field': 'height', 'val': 0.03, 'expected': {'height': 0.03}},
    {'field': 'width',  'val': 0.02, 'expected': {'width': 0.02}},
    {'field': 'side',   'val': 'right', 'expected': {'side': 'right'}},
]


# Tests
# =====
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # BlackAbsorber not on GPU
)
def test_black_absorber(test_context):
    # Test instantiation
    elem = xc.BlackAbsorber(length=1, _context=test_context)
    _check_all_elements(elem, absorber_fields, absorber_dict_fields, absorber_user_fields, \
                        absorber_user_fields_read_only)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # not on GPU
)
def test_black_crystal(test_context):
    # Test instantiation
    elem = xc.BlackCrystal(length=1, jaw=0.99, side='-', _context=test_context)
    _check_all_elements(elem, black_crystal_fields, black_crystal_dict_fields, \
                        black_crystal_user_fields, black_crystal_user_fields_read_only)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest_block(test_context):
    # Test instantiation
    elem = xc.EverestBlock(length=1.3, material=xc.materials.CarbonFibreCarbon, _context=test_context)
    assert np.isclose(elem.length, 1.3)
    assert elem._tracking == True
    assert deep_equal(elem.material.to_dict(), xc.materials.CarbonFibreCarbon.to_dict())
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
    elem = xc.EverestCollimator(length=1, material=xc.materials.CarbonFibreCarbon, _context=test_context)
    _check_all_elements(elem, everest_fields, everest_dict_fields, everest_user_fields, \
                        everest_user_fields_read_only)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_everest_crystal(test_context):
    # Test instantiation
    elem = xc.EverestCrystal(length=1, jaw=0.99, material=xc.materials.Silicon, side='-', _context=test_context)
    _check_all_elements(elem, everest_crystal_fields, everest_crystal_dict_fields, \
                        everest_crystal_user_fields, everest_crystal_user_fields_read_only)

@pytest.mark.fluka
def test_fluka():
    # Test instantiation
    elem = xc.FlukaCollimator(length=1, assembly='lhc_tcp')
    _check_all_elements(elem, fluka_fields, fluka_dict_fields, [], [])

@pytest.mark.fluka
def test_fluka_generic():
    # Test instantiation
    elem = xc.FlukaCollimator(length=1, material=xc.materials.CarbonFibreCarbon)
    _check_all_elements(elem, fluka_fields, fluka_generic_dict_fields, [], [])

# @pytest.mark.fluka
# def test_fluka_crystal():
#     # Test instantiation
#     elem = xc.FlukaCrystal(length=1, jaw=0.99, material=xc.materials.Silicon, bending_radius=20, side='+')
#     _check_all_elements(elem, fluka_crystal_fields, fluka_crystal_dict_fields, [], [])


def _assert_all_close(expected, setval):
    if hasattr(expected, 'to_dict'):
        assert deep_equal(expected.to_dict(), setval.to_dict())
    elif hasattr(expected, '__iter__') and not isinstance(expected, str):
        for exp, stv in zip(expected, setval):
            _assert_all_close(exp, stv)
    else:
        if expected is None:
            assert setval is None
        elif isinstance(expected, str):
            assert setval == expected
        else:
            assert np.isclose(expected, setval, atol=1e-12, rtol=0)


def _check_all_elements(elem, fields, dict_fields, user_fields, user_fields_read_only):
    # Test existence of fields
    assert np.all([key in dir(elem) for key in fields])
    elem.jaw = 0.8 # Need to give non-default value to jaw and gap, otherwise it is None and not in to_dict()
    elem.gap = 5
    assert np.all([key['field'] in elem.to_dict() for key in dict_fields])
    elem.jaw = 1
    elem.gap = None

    # Test reading fields
    for field in list(fields.keys()) + [vv['field'] for vv in dict_fields] \
               + [vv['field'] for vv in user_fields] + user_fields_read_only:
        print(f"Reading field {field}: {getattr(elem, field)}")

    # Test writing xofields
    for field, val in fields.items():
        print(f"Writing field {field}")
        setattr(elem, field, val)
        setval = getattr(elem, field)
        _assert_all_close(val, setval)

    # Test writing the to_dict and user-friendly fields (can be multiple options per field)
    for data in [*dict_fields, *user_fields]:
        field = data['field']
        val = data['val']
        setattr(elem, field, val)
        for basefield, expected in data['expected'].items():
            setval = getattr(elem, basefield)
            print(f"Writing field {field} with {val}:  {basefield} expected {expected}, got {setval}")
            _assert_all_close(expected, setval)

    # Writing to a read-only field should fail
    for field in user_fields_read_only:
        with pytest.raises(Exception) as e_info:
            setattr(elem, field, 0.3)


# TODO:
# def test_jaw_func():
