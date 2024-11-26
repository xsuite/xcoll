# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path
import xtrack as xt
import xcoll as xc
import xpart as xp
import pytest

from xcoll.scattering_routines.fluka.fluka_input import _coll_dict, create_fluka_input

_HEADER_START = "*  XCOLL START  **"
_HEADER_END  = "*  XCOLL END  **"
EXPECTED_HEADER = {
    "tctpv.4l2.b1": {
        "fluka_id": 1,
        "length": 1.482,
        "jaw": [0.00648527, -0.00997401],
    }
}
EXPECTED_ROT_ANGLE = "90.0"
EXPECTED_CENTERED_JAW = [-0.822964, 0.822964]
EXPECTED_SHIFT = 0.17443658

def get_collimators_from_input_file(input_file):
    with open(input_file, 'r') as fp:     
        data = fp.read() 
    if _HEADER_START not in data or _HEADER_END not in data:
        raise ValueError("No XCOLL header found in input file. Regenerate input file!")
    commented_dict = data.split(_HEADER_START)[1].split(_HEADER_END)[0].split('\n')[1:-1]
    cleaned_dict = "".join([val[3:] for val in commented_dict])        
    return json.loads(cleaned_dict) 

def validate_header(header):
    # Compare the values of the dictionary in the header with the expected values
    assert header['tctpv.4l2.b1']['fluka_id'] == EXPECTED_HEADER['tctpv.4l2.b1']['fluka_id']
    assert header['tctpv.4l2.b1']['length'] == EXPECTED_HEADER['tctpv.4l2.b1']['length']
    # Check jaw with a tolerance of 1e-8
    for jaw, expected_jaw in zip(header['tctpv.4l2.b1']['jaw'], EXPECTED_HEADER['tctpv.4l2.b1']['jaw']):
        assert round(jaw, 8) == expected_jaw

    # assert round(header['tctpv.4l2.b1']['jaw'][0], 8) == 0.00648527
    # assert round(header['tctpv.4l2.b1']['jaw'][1], 8) == -0.00997401

def validate_fluka_cards(input_file):
    # Check that the input file contains the expected cards
    with open(input_file, 'r') as fp:
        data = fp.read()

    # Rot-defi lines for jaw container
    rot_defi_lines_container = [line for line in data.split('\n') if 'ROT-DEFI' in line and 'TCTPV____1' in line]

    # Check angle for container
    assert rot_defi_lines_container[1].split()[3] == EXPECTED_ROT_ANGLE

    # Rot-defi lines for jaw left
    rot_defi_lines_jaw_l = [line for line in data.split('\n') if 'ROT-DEFI' in line and 'TCTPV____2' in line]

    # Check angle, jaw opening and shift for jaw left
    assert rot_defi_lines_jaw_l[1].split()[3] == EXPECTED_ROT_ANGLE
    assert round(float(rot_defi_lines_jaw_l[2].split()[4]), 6) == EXPECTED_CENTERED_JAW[0] # jaw opening
    assert round(float(rot_defi_lines_jaw_l[3].split()[4]), 8) == EXPECTED_SHIFT # jaw shift

    # Rot-defi lines for jaw right
    rot_defi_lines_jaw_r = [line for line in data.split('\n') if 'ROT-DEFI' in line and 'TCTPV____3' in line]

    # Check angle, jaw opening and shift for jaw right
    assert rot_defi_lines_jaw_r[1].split()[3] == EXPECTED_ROT_ANGLE
    assert round(float(rot_defi_lines_jaw_r[2].split()[4]), 6) == EXPECTED_CENTERED_JAW[1] # jaw opening
    assert round(float(rot_defi_lines_jaw_r[3].split()[4]), 8) == EXPECTED_SHIFT # jaw shift
                                          
def test_fluka_input():
    beam = 1
    path = xc._pkg_root.parent / 'examples'

    # Load from json
    line = xt.Line.from_json(path / 'machines' / f'lhc_run3_b{beam}.json')

    # Load collimators
    colldb = xc.CollimatorDatabase.from_yaml(path / 'colldb' / f'lhc_run3_fluka.yaml', beam=beam)
    colldb.install_fluka_collimators(line=line, verbose=True)
    line.build_tracker()
    line.collimators.assign_optics()

    prototypes_file = xc._pkg_root.parent / 'xcoll' / 'scattering_routines' / 'fluka' / 'data' / 'prototypes.lbp'
    #prototypes_file = "/afs/cern.ch/work/a/adonadon/public/fellow/templates/fordevelopment/prototypes_b1.lbp"
    include_files = [
        xc._pkg_root.parent / 'xcoll' / 'scattering_routines' / 'fluka' / 'data' / 'include_settings_beam.inp',
        xc._pkg_root.parent / 'xcoll' / 'scattering_routines' / 'fluka' / 'data' / 'include_settings_physics.inp',
        xc._pkg_root.parent / 'xcoll' / 'scattering_routines' / 'fluka' / 'data' / 'include_custom_scoring.inp'
    ]

    from xcoll.beam_elements import FlukaCollimator

    elements, names = line.get_elements_of_type(FlukaCollimator)

    coll = elements[4]
    coll_name = names[4] # tctpv.4l2.b1

    input_file = create_fluka_input(elements=coll, names=coll_name, prototypes_file=prototypes_file, include_files=include_files)
    header = get_collimators_from_input_file(input_file)

    # Test header
    validate_header(header)
    # Test fluka cards
    validate_fluka_cards(input_file)