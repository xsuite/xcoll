# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
import sys
import os
import shutil
from subprocess import run, PIPE

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except ImportError:
    from ...xaux import FsPath

from ...beam_elements import FlukaCollimator
from ...beam_elements.base import OPEN_GAP, OPEN_JAW
from ...general import _pkg_root
from .paths import fedb, linebuilder, flukafile_resolve


_lb_souce_files = [
    fedb / "tools" / "expand.sh",
    fedb / "structure.py",
    *list((linebuilder / "src").glob('*.py')),
    *list((linebuilder / "lib").glob('*.py')),
    *list((linebuilder / "db").glob('*.py')),
]

_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"


# TODO check that prototype is valid and its sides
def _element_dict_to_fluka(element_dict, dump=False):
    collimator_dict = {}
    for name, ee in element_dict.items():
        # nsig = OPEN_GAP
        nsig = 1
        if ee.side == 'left':
            if ee.jaw_L is None:
                half_gap = OPEN_JAW
            else:
                nsig = ee.gap_L
                half_gap = ee.jaw_L
            offset = 0
            tilt_1 = ee.tilt_L
            tilt_2 = 0
        elif ee.side == 'right':
            if ee.jaw_R is None:
                half_gap = OPEN_JAW
            else:
                nsig = ee.gap_R
                half_gap = -ee.jaw_R   #  TODO: is the sign correct?
            offset = 0
            tilt_1 = 0
            tilt_2 = ee.tilt_R
        else:
            if ee.jaw_L is None and ee.jaw_R is None:
                half_gap = OPEN_JAW
                offset = 0
            else:
                if ee.gap_L is not None:
                    nsig = ee.gap_L
                elif ee.gap_R is not None:
                    nsig = ee.gap_R
                half_gap = (ee._jaw_LU + ee._jaw_LD - ee._jaw_RU - ee._jaw_RD) / 4
                offset   = (ee._jaw_LU + ee._jaw_LD + ee._jaw_RU + ee._jaw_RD) / 4
            tilt_1 = ee.tilt_L
            tilt_2 = ee.tilt_R

        collimator_dict[name] = {
            'name': name,
            'betx': 1,
            'bety': 1,
            'material': 'stub',
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
    if dump:
        # dump coll_dictionary in json format
        with open('collimator_dict.json', 'w') as fp:
            json.dump(collimator_dict, fp, indent=4)
    return collimator_dict


def _fluka_builder(element_dict):
    # Save system state
    old_sys_path = sys.path.copy()
    old_os_env = os.environ.copy()

    this_fedb = flukafile_resolve(fedb)
    if this_fedb is not None:
        os.environ['FEDB_PATH'] = this_fedb.as_posix()
    else:
        raise ValueError(f"Could not find fedb folder {fedb}!")
    this_linebuilder = flukafile_resolve(linebuilder)
    if this_linebuilder is not None:
        os.environ['LB_PATH'] = this_linebuilder.as_posix()
    else:
        raise ValueError(f"Could not find linebuilder folder {linebuilder}!")

    for f in _lb_souce_files:
        # This forces syncing of the files
        if not f.exists():
            raise FileNotFoundError(f"Linebuilder source file not found: {f}.")

    sys.path.append((linebuilder / "src").as_posix())
    file_path = linebuilder / "src" / "FLUKA_builder.py"
    if file_path.exists():
        try:
            import FLUKA_builder as fb
        except ImportError as e:
            raise EnvironmentError(f"Cannot import FLUKA_builder: {e}")
    else:
        raise EnvironmentError("FLUKA_builder.py not found at:", file_path)

    collimator_dict = _element_dict_to_fluka(element_dict)
    collimatorList = fb.CollimatorList()
    collimatorList.acquireCollxsuite(collimator_dict)

    args_fb = fb.args_fluka_builder()
    args_fb.collimatorList = collimatorList
    args_fb.geometrical_emittance = None
    args_fb.prototype_file = 'prototypes.lbp'
    args_fb.output_name = 'fluka_input'

    input_file, coll_dict = fb.fluka_builder(args_fb, auto_accept=True)

    # Restore system state
    sys.path = old_sys_path
    os.environ = old_os_env

    return input_file, coll_dict


def _write_xcoll_header_to_fluka_input(input_file, collimator_dict):
    header = ["*  DO NOT CHANGE THIS HEADER", _header_start, "*  {"]
    for kk, vv in collimator_dict.items():
        header.append(f'*  "{kk}": ' + json.dumps(vv).replace('"jaw"', '\n*          "jaw"') + ',')
    header[-1] = header[-1][:-1]  # remove last comma
    header.append("*  }")
    header.append(_header_stop)

    with open(input_file, 'r') as fp:
        data = fp.read()
    with open(input_file, 'w') as fp:
        fp.write("\n".join(header) + "\n*\n" + data)


def create_fluka_input(element_dict, prototypes_file, include_files, *, filename=None,
                       cwd=None, verbose=False):
    prototypes_file = FsPath(prototypes_file).resolve()
    if not prototypes_file.exists():
        raise FileNotFoundError(f"Prototypes file not found: {prototypes_file}.")
    include_files = [FsPath(ff).resolve() for ff in include_files]
    for ff in include_files:
        if not ff.exists():
            raise FileNotFoundError(f"Include file not found: {ff}.")
    required_includes = ['include_settings_beam.inp', 'include_settings_physics.inp',
                         'include_custom_scoring.inp']
    for ff in required_includes:
        if ff not in [file.name for file in include_files]:
            raise ValueError(f"Missing include file {ff}.")
    for ff in (_pkg_root / 'scattering_routines' / 'fluka' / 'data').glob('include_*'):
        if ff.name not in [file.name for file in include_files]:
            include_files.append(ff)

    # Change to the directory of the input file
    if filename is not None:
        filename = FsPath(filename).expanduser().resolve().with_suffix('.inp')
        if cwd is None:
            cwd = filename.parent
    if cwd is None:
        cwd = FsPath.cwd()
    cwd.mkdir(parents=True, exist_ok=True)
    prev_cwd = FsPath.cwd()
    os.chdir(cwd)
    shutil.copy(prototypes_file, FsPath.cwd() / 'prototypes.lbp')
    for ff in include_files:
        shutil.copy(ff, FsPath.cwd() / ff.name)
    # Call FLUKA_builder
    input_file, fluka_dict = _fluka_builder(element_dict)
    input_file = FsPath(input_file).resolve()
    insertion_file = (input_file.parent / 'insertion.txt').resolve()
    assert input_file.exists()
    assert insertion_file.exists()
    # Expand using include files
    if verbose:
        print(f"Expanding {input_file} using {include_files}.")
    cmd = run([(fedb / 'tools' / 'expand.sh').as_posix(), input_file.name],
              cwd=FsPath.cwd(), stdout=PIPE, stderr=PIPE)
    if cmd.returncode == 0:
        if verbose:
            print("Expanded include files.")
    else:
        stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
        raise RuntimeError(f"Could not expand include files!\nError given is:\n{stderr}")
    for name, ee in element_dict.items():
        if name not in fluka_dict:
            continue
        fluka_dict[name]['length'] /= 100
        fluka_dict[name]['angle'] = ee.angle
        fluka_dict[name]['tilt'] = [ee.tilt_L, ee.tilt_R]
        fluka_dict[name]['jaw'] = [ee.jaw_L, ee.jaw_R]
    # Delete prototypes and include files to avoid confusion
    if prototypes_file.parent != FsPath.cwd():
        (FsPath.cwd() / 'prototypes.lbp').unlink()
    for ff in include_files:
        if ff.parent != FsPath.cwd():
            (FsPath.cwd() / ff.name).unlink()
    # Return to previous directory
    os.chdir(prev_cwd)
    if filename is None:
        filename = input_file.name
    input_file     = input_file.parent / f'{input_file.stem}_exp.inp'
    input_file     = input_file.rename(cwd / filename)
    insertion_file = insertion_file.rename(cwd / 'insertion.txt')
    _write_xcoll_header_to_fluka_input(input_file, fluka_dict)
    assert input_file.exists()
    assert insertion_file.exists()
    if verbose:
        print(f"Created FLUKA input file {input_file}.")
    return input_file, insertion_file


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
        if val.fluka_id not in all_fluka_ids:
            raise ValueError(f'FlukaCollimator {name} not found in insertion file!')

