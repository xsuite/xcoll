# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
import sys
import os
from pathlib import Path
import shutil

from ...beam_elements import FlukaCollimator
from ...beam_elements.base import OPEN_GAP, OPEN_JAW
from .paths import fluka_builder, fedb, linebuilder


_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"


# TODO check that prototype is valid and its sides
def _fluka_builder(elements, names, prototypes_file):
    # Save system state
    old_sys_path = sys.path.copy()
    old_os_env = os.environ.copy()

    os.environ['FEDB_PATH'] = fedb.as_posix()
    os.environ['LB_PATH'] = linebuilder.as_posix()

    sys.path.append(fluka_builder.as_posix())
    file_path = fluka_builder / "FLUKA_builder.py"
    if file_path.exists():
        try:
            import FLUKA_builder as fb
        except ImportError as e:
            raise EnvironmentError(f"Cannot import FLUKA_builder: {e}")
    else:
        raise EnvironmentError("FLUKA_builder.py not found at:", file_path)

    collimator_dict = {}
    for ee, name in zip(elements, names):
        if ee.side == 'left':
            if ee.jaw_L is None:
                nsig = OPEN_GAP
                half_gap = OPEN_JAW
            else:
                nsig = ee.gap_L
                half_gap = ee.jaw_L
            offset = 0
            tilt_1 = ee.tilt_L
            tilt_2 = 0
        elif ee.side == 'right':
            if ee.jaw_R is None:
                nsig = OPEN_GAP
                half_gap = OPEN_JAW
            else:
                nsig = ee.gap_R
                half_gap = -ee.jaw_R   #  TODO: is the sign correct?
            offset = 0
            tilt_1 = 0
            tilt_2 = ee.tilt_R
        else:
            if ee.jaw_L is None and ee.jaw_R is None:
                nsig = OPEN_GAP
                half_gap = OPEN_JAW
                offset = 0
            else:
                nsig = ee.gap_L if ee.gap_L is not None else ee.gap_R
                half_gap = (ee._jaw_LU + ee._jaw_LD - ee._jaw_RU - ee._jaw_RD) / 4
                offset   = (ee._jaw_LU + ee._jaw_LD + ee._jaw_RU + ee._jaw_RD) / 4
            tilt_1 = ee.tilt_L
            tilt_2 = ee.tilt_R

        collimator_dict[name] = {
            'name': name,
            'betx': 1,
            'bety': 1,
            'material': ee.material,
            'length': ee.length,
            'angle': np.deg2rad(ee.angle),
            'sigma_x': 1,
            'sigma_y': 1,
            'offset': offset,
            'tilt_1': tilt_1,
            'tilt_2': tilt_2,
            'nsig': nsig,
            'half_gap': half_gap
        }

    collimatorList = fb.CollimatorList()
    collimatorList.acquireCollxsuite(collimator_dict)

    args_fb = fb.args_fluka_builder()
    args_fb.collimatorList = collimatorList
    args_fb.geometrical_emittance = 1
    args_fb.prototype_file = 'prototypes.lbp'
    args_fb.output_name = 'fluka_input'

    shutil.copy(prototypes_file, Path.cwd() / 'prototypes.lbp')
    input_file, coll_dict = fb.fluka_builder(args_fb)

    # Restore system state
    sys.path = old_sys_path
    os.environ = old_os_env

    return input_file, coll_dict


def _write_xcoll_header_to_fluka_input(input_file, collimator_dict):
    header = ["*  DO NOT CHANGE THIS HEADER", _header_start, "*  {"]
    for kk, vv in collimator_dict.items():
        header.append(f'*  "{kk}": ' + json.dumps(vv).replace('" jaw"', '\n*          "jaw"') + ',')
    header[-1] = header[-1][:-1]  # remove last comma
    header.append("*  }")
    header.append(_header_stop)

    with open(input_file, 'r') as fp:
        data = fp.read()
    with open(input_file, 'w') as fp:
        fp.write("\n".join(header) + "\n*\n" + data)


def create_fluka_input(line, prototypes_file, include_files=[], cwd=None):
    elements, names = line.get_elements_of_type(FlukaCollimator)
    if len(elements) == 0:
        raise ValueError('No FlukaCollimator elements found in line')
    prototypes_file = Path(prototypes_file).resolve()
    if not prototypes_file.exists():
        raise FileNotFoundError(f"Prototypes file not found: {prototypes_file}")
    prev_cwd = Path.cwd()
    tempdir = Path.cwd() / 'temp_fluka_input'
    tempdir.mkdir(exist_ok=True)
    os.chdir(tempdir)
    input_file, collimator_dict = _fluka_builder(elements, names, prototypes_file)
    input_file = Path(input_file).resolve()
    # copy include files
    # fedb / 'expand.sh
    _write_xcoll_header_to_fluka_input(input_file, collimator_dict)
    os.chdir(prev_cwd)
    return input_file


def get_collimators_from_input_file(input_file):
    with open(input_file, 'r') as fp:
        data = fp.read()
    if _header_start not in data or _header_stop not in data:
        raise ValueError("No XCOLL header found in input file. Regenerate input file!")
    commented_dict = data.split(_header_start)[1].split(_header_stop)[0].split('\n')[1:-1]
    cleaned_dict = "".join([val[3:] for val in commented_dict])
    return json.loads(cleaned_dict)


def verify_insertion_file(insertion_file, collimator_dict):
    all_fluka_ids = []
    with open(insertion_file, 'r') as fid:
        for line in fid.readlines():
            all_fluka_ids.append(int(line.split()[0]))
    for name, val in collimator_dict.items():
        if val['fluka_id'] not in all_fluka_ids:
            raise ValueError(f'FlukaCollimator {name} not found in insertion file!')

