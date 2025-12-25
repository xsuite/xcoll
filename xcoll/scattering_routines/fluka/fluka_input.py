# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import numpy as np
from subprocess import run, PIPE
from contextlib import redirect_stdout

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ...beam_elements.base import OPEN_GAP, OPEN_JAW
from ...general import _pkg_root
from .prototype import FlukaPrototype, FlukaAssembly
from .includes import get_include_files

_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"


def create_fluka_input(element_dict, particle_ref, prototypes_file=None,
                       verbose=True, **kwargs):
    import xcoll as xc
    _create_prototypes_file(element_dict, prototypes_file)
    for prototype in {**FlukaPrototype._assigned_registry, **FlukaAssembly._assigned_registry}.values():
        if prototype.is_crystal:
            kwargs.setdefault('use_crystals', True)
    include_files, kwargs = get_include_files(particle_ref, verbose=verbose,
                                              **kwargs)
    # Call FLUKA_builder
    collimator_dict = _element_dict_to_fluka(element_dict)
    input_file, fluka_dict = _fluka_builder(collimator_dict)
    input_file = FsPath(input_file).resolve()
    insertion_file = (input_file.parent / 'insertion.txt').resolve()
    assert input_file.exists()
    assert insertion_file.exists()
    # Expand using include files
    fedb = xc.fluka.environment.fedb
    cmd = run([(fedb / 'tools' / 'expand.sh').as_posix(), input_file.name],
              cwd=FsPath.cwd(), stdout=PIPE, stderr=PIPE)
    if cmd.returncode == 0:
        if verbose:
            print("Expanded include files.")
    else:
        stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
        raise RuntimeError(f"Could not expand include files!\nError given is:\n{stderr}")
    new_input_file = input_file.parent / f'{input_file.stem}_exp.inp'
    assert new_input_file.exists()
    input_file.rename(input_file.parent / f'{input_file.stem}_orig.inp')
    new_input_file.rename(input_file)
    # Check that all collimators were treated and write header
    for name, ee in element_dict.items():
        if name not in fluka_dict:
            if verbose:
                print(f"Warning: Collimator {name} was requested but not treated by "
                    + f"the LineBuilder.")
        else:
            fluka_dict[name]['length'] /= 100
            fluka_dict[name]['angle'] = ee.angle
            if ee.assembly.is_crystal:
                if ee.side == 'left':
                    fluka_dict[name]['tilt'] = [ee.tilt, None]
                    fluka_dict[name]['jaw'] = [ee.jaw, None]
                elif ee.side == 'right':
                    fluka_dict[name]['tilt'] = [None, ee.tilt]
                    fluka_dict[name]['jaw'] = [None, ee.jaw]
            else:
                fluka_dict[name]['tilt'] = [ee.tilt_L, ee.tilt_R]
                fluka_dict[name]['jaw'] = [ee.jaw_L, ee.jaw_R]
    _write_xcoll_header_to_fluka_input(input_file, fluka_dict)
    if verbose:
        print(f"Created FLUKA input file {input_file}.")
    return [input_file, insertion_file], kwargs


def get_collimators_from_input_file(input_file):
    with open(input_file, 'r') as fp:
        data = fp.read()
    if _header_start not in data or _header_stop not in data:
        raise ValueError("No Xcoll header found in input file. Regenerate input file!")
    commented_dict = data.split(_header_start)[1].split(_header_stop)[0].split('\n')[1:-1]
    cleaned_dict = "".join([val[3:] for val in commented_dict])
    return json.loads(cleaned_dict)


def verify_insertion_file(insertion_file, element_dict):
    all_fluka_ids = []
    with open(insertion_file, 'r') as fid:
        for line in fid.readlines():
            all_fluka_ids.append(int(line.split()[0]))
    for name, val in element_dict.items():
        if val.fluka_id not in all_fluka_ids:
            raise ValueError(f'FlukaCollimator {name} not found in insertion file!')


def _create_prototypes_file(element_dict, prototypes_file=None):
    for name, ee in element_dict.items():
        if ee.assembly is None:
            raise ValueError(f"Collimator {name} has no assembly!")
    if prototypes_file is None:
        FlukaPrototype.make_prototypes()
    else:
        prototypes_file = FsPath(prototypes_file).resolve()
        prototypes_file.copy_to(FsPath.cwd() / 'prototypes.lbp', method='mount')
    FlukaPrototype.inspect_prototypes_file(prototypes_file)


# TODO check that prototype is valid and its sides
def _element_dict_to_fluka(element_dict, dump=False):
    collimator_dict = {}
    for name, ee in element_dict.items():
        if ee.length < 1.e-12:
            raise ValueError(f"Collimator {name} has zero length!")

        nsig = 1 # TODO can remove?
        if ee.assembly.is_crystal:
            if ee.jaw is None:
                half_gap = OPEN_JAW
            else:
                nsig = ee.gap
                half_gap = -ee.jaw  # Mistake in LineBuilder?
            offset = 0
            tilt_1 = ee.tilt
            tilt_2 = 0
        else:
            tilt_1 = ee.tilt_L
            tilt_2 = ee.tilt_R
            if ee.side == 'left':
                if ee.jaw_L is None:
                    half_gap = OPEN_JAW
                else:
                    nsig = ee.gap_L
                    half_gap = ee.jaw_L
                offset = 0
                tilt_2 = 0
            elif ee.side == 'right':
                if ee.jaw_R is None:
                    half_gap = OPEN_JAW
                else:
                    nsig = ee.gap_R
                    half_gap = -ee.jaw_R   #  TODO: is the sign correct?
                offset = 0
                tilt_1 = 0
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
        tilt_1 = round(tilt_1, 9)
        tilt_2 = round(tilt_2, 9)

        if nsig is None:
            nsig = 1

        collimator_dict[name] = {
            'name': name,
            'betx': 1,
            'bety': 1,
            'material': 'stub',
            'length': ee.length,
            'angle': round(np.deg2rad(ee.angle) - ee.assembly.angle, 9),
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


def _fluka_builder(collimator_dict):
    import xcoll as xc
    # Save system state
    xc.fluka.environment.set_fedb_environment()
    file_path = xc.fluka.environment.linebuilder / "src" / "FLUKA_builder.py"
    if file_path.exists():
        try:
            import FLUKA_builder as fb
        except ImportError as e:
            raise EnvironmentError(f"Cannot import FLUKA_builder: {e}")
    else:
        raise EnvironmentError(f"FLUKA_builder.py not found at: {file_path.as_posix()}")
    collimatorList = fb.CollimatorList()
    collimatorList.acquireCollxsuite(collimator_dict)

    args_fb = fb.args_fluka_builder()
    args_fb.collimatorList = collimatorList
    args_fb.geometrical_emittance = None
    args_fb.prototype_file = 'prototypes.lbp'
    args_fb.output_name = 'fluka_input'

    with open('linebuilder.log', 'w') as f:
        with redirect_stdout(f):
            input_file, coll_dict = fb.fluka_builder(args_fb, auto_accept=True)

    # Restore system state
    xc.fluka.environment.restore_environment()

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
