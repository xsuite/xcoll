# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import numpy as np
import sys
import os
from pathlib import Path
import shutil
from subprocess import run, PIPE

from ...beam_elements import FlukaCollimator
from ...beam_elements.base import OPEN_GAP, OPEN_JAW
from ...general import _pkg_root
from .paths import fedb, linebuilder, flukafile_resolve
from .prototypes import FlukaAssembly


_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"


# TODO check that prototype is valid and its sides
def _coll_dict(elements, names, dump=False):
    collimator_dict = {}
    for ee, name in zip(elements, names):
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
            tilt_1 = np.round(ee.tilt_L, 9)
            tilt_2 = np.round(ee.tilt_R, 9)
        if abs(tilt_1) > 1.e-12 or abs(tilt_2) > 1.e-12:
            raise NotImplementedError(f"Collimator {name}: Tilts are not (yet) supported in FLUKA-Xcoll!")

        if nsig is None:
            nsig = 1

        collimator_dict[name] = {
            'name': name,
            'betx': 1,
            'bety': 1,
            'material': 'stub',
            'length': ee.length,
            'angle': np.round(np.deg2rad(ee.angle) - ee.assembly.angle, 9),
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


def _brute_force_environment():
    # Add the paths to the environment
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
    # Brute-force the system paths
    sys.path.append(fedb.as_posix())
    sys.path.append((fedb / "tools").as_posix())
    sys.path.append((fedb / "tools" / "materials" / "cables").as_posix())
    sys.path.append((linebuilder / "src").as_posix())
    sys.path.append((linebuilder / "lib").as_posix())
    # Force-resolve all dependent files
    flukafile_resolve(fedb / "structure.py")
    flukafile_resolve(fedb / "tools")
    for ff in (fedb / "tools").glob("*.py"):
        flukafile_resolve(ff)
    flukafile_resolve(fedb / "tools" / "materials" / "cables")
    for ff in (fedb / "tools" / "materials" / "cables").glob("*.py"):
        flukafile_resolve(ff)
    flukafile_resolve(linebuilder / "lib")
    for ff in (linebuilder / "lib").glob("*.py"):
        flukafile_resolve(ff)
    flukafile_resolve(linebuilder / "additionals" / "generic_frame.fluka")
    flukafile_resolve(linebuilder / "src")
    for ff in (linebuilder / "src").glob("*.py"):
        flukafile_resolve(linebuilder / "src" / ff)


def _fluka_builder(elements, names):
    # Save system state
    old_sys_path = sys.path.copy()
    old_os_env = os.environ.copy()
    _brute_force_environment()
    file_path = linebuilder / "src" / "FLUKA_builder.py"
    if file_path.exists():
        try:
            import FLUKA_builder as fb
        except ImportError as e:
            raise EnvironmentError(f"Cannot import FLUKA_builder: {e}")
    else:
        raise EnvironmentError(f"FLUKA_builder.py not found at: {file_path.as_posix()}")

    collimator_dict = _coll_dict(elements, names)
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


def _select_active_collimators(elements, names):
    this_names = []
    this_elements = []
    for name, ee in zip(names, elements):
        if (ee.jaw_L is None and ee.jaw_R is None):
            print(f"Collimator {name} is fully open. Deactivated.")
            ee.active = False
            ee.assembly.remove_element(name)
        elif not ee.active:
            ee.assembly.remove_element(name)
        else:
            this_names.append(name)
            this_elements.append(ee)
    return this_elements, this_names


def create_fluka_input(include_files, *, line=None, elements=None, names=None,
                       prototypes_file=None, filename=None, cwd=None):
    # Checks
    if elements is None or names is None:
        if line is None:
            raise ValueError("Need to provide either `line` or `elements` and `names`.")
        elements, names = line.get_elements_of_type(FlukaCollimator)
    elif line is not None:
        print("Warning: `line`, `elements`, and `names` are provided. Only the "
            + "latter two will be used.")
    if not hasattr(elements, '__iter__') or isinstance(elements, str):
        elements = [elements]
    if not hasattr(names, '__iter__') or isinstance(names, str):
        names = [names]
    assert len(elements) == len(names)
    elements, names = _select_active_collimators(elements, names)
    if len(elements) == 0:
        raise ValueError('No active FlukaCollimator elements found in line!')
    # Includes
    include_files = [Path(ff).resolve() for ff in include_files]
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
        filename = Path(filename).expanduser().resolve().with_suffix('.inp')
        if cwd is None:
            cwd = filename.parent
    if cwd is None:
        cwd = Path.cwd()
    cwd.mkdir(parents=True, exist_ok=True)
    prev_cwd = Path.cwd()
    os.chdir(cwd)
    # Prototypes
    if prototypes_file is None:
        for name, el in zip(names, elements):
            el.assembly.add_element(name)
        FlukaAssembly.make_prototypes()
    else:
        prototypes_file = Path(prototypes_file).resolve()
        shutil.copy(prototypes_file, Path.cwd() / 'prototypes.lbp')
    for ff in include_files:
        shutil.copy(ff, Path.cwd() / ff.name)
    for name, el in zip(names, elements):
        if not el.assembly.in_file(prototypes_file):
            raise ValueError(f"Prototype {el.assembly.name} for {name} not found "
                           + f"in prototypes file!")
    # Call FLUKA_builder
    input_file, collimator_dict = _fluka_builder(elements, names)
    input_file = Path(input_file).resolve()
    insertion_file = (input_file.parent / 'insertion.txt').resolve()
    assert input_file.exists()
    assert insertion_file.exists()
    # Expand using include files
    print(f"Expanding {input_file} using {include_files}.")
    cmd = run([(fedb / 'tools' / 'expand.sh').as_posix(), input_file.name],
              cwd=Path.cwd(), stdout=PIPE, stderr=PIPE)
    if cmd.returncode == 0:
        print("Expanded include files.")
    else:
        stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
        raise RuntimeError(f"Could not expand include files!\nError given is:\n{stderr}")
    for name,ee in zip(names, elements):
        if name not in collimator_dict:
            continue
        collimator_dict[name]['length'] /= 100
        collimator_dict[name]['angle'] = ee.angle
        collimator_dict[name]['tilt'] = [ee.tilt_L, ee.tilt_R]
        collimator_dict[name]['jaw'] = [ee.jaw_L, ee.jaw_R]
    # If a prototypes file was provided (instead of generated), delete it
    if prototypes_file is not None and prototypes_file.parent != Path.cwd():
        (Path.cwd() / 'prototypes.lbp').unlink()
    # Delete include files
    for ff in include_files:
        if ff.parent != Path.cwd():
            (Path.cwd() / ff.name).unlink()
    # Return to previous directory
    os.chdir(prev_cwd)
    if filename is None:
        filename = input_file.name
    input_file     = input_file.parent / f'{input_file.stem}_exp.inp'
    input_file     = input_file.rename(cwd / filename)
    insertion_file = insertion_file.rename(cwd / 'insertion.txt')
    _write_xcoll_header_to_fluka_input(input_file, collimator_dict)
    assert input_file.exists()
    assert insertion_file.exists()
    print(f"Created FLUKA input file {input_file}.")
    return input_file, elements, names


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

