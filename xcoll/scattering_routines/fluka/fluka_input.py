# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import sys
import json
import numpy as np
import pathlib as Path
from subprocess import run, PIPE
from contextlib import redirect_stdout

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ...beam_elements.base import OPEN_GAP, OPEN_JAW
from ...general import _pkg_root
from .environment import FlukaEnvironment
from .prototypes import FlukaAssembly
from .includes import get_include_files


_header_start = "*  XCOLL START  **"
_header_stop  = "*  XCOLL END  **"


def create_fluka_input(element_dict, particle_ref, prototypes_file=None, include_files=[],
                       verbose=True, **kwargs):
    # _create_generic_colls(element_dict)
    include_files = get_include_files(particle_ref, include_files, verbose=verbose, **kwargs)
    # Call FLUKA_builder
    collimator_dict = _element_dict_to_fluka(element_dict)
    input_file, fluka_dict = _fluka_builder(collimator_dict)
    input_file = FsPath(input_file).resolve()
    insertion_file = (input_file.parent / 'insertion.txt').resolve()
    assert input_file.exists()
    assert insertion_file.exists()
    # Expand using include files
    cmd = run([(FlukaEnvironment.fedb_base / 'tools' / 'expand.sh').as_posix(), input_file.name],
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
            fluka_dict[name]['tilt'] = [ee.tilt_L, ee.tilt_R]
            fluka_dict[name]['jaw'] = [ee.jaw_L, ee.jaw_R]
    _rite_xcoll_header_to_fluka_input(input_file, fluka_dict)
    if verbose:
        print(f"Created FLUKA input file {input_file}.")
    return input_file, insertion_file

def xcoll_to_fluka_material(material):
    # XXX EXPAND THIS DICT

    material_dict = {
        "iner": "INERM180",
        "c":    "  CARBON",
        "Si":   " SILICON",
        "cu":   "  COPPER",
        "mogr": "AC150GPH"
    }

    if material not in material_dict.keys():
        raise ValueError(f"Material {material} not found in material dictionary!")
    else:
        return material_dict[material]

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

def _create_generic_assemblies(_generic_colls):
    for gen in _generic_colls.keys():
        template_assembly = f"""\
# --------------------------------------------------------------------------------------------------------------
PROTOTYPE       {gen}_T
# --------------------------------------------------------------------------------------------------------------
# prototype of the tank:
FEDB_TAG        {gen}_T                # tag in the fedb (for the name of the file)
FEDB_SERIES     test                    # series of the fedb (for the name of the file)
# this rot-defi is only for local use:
ROT-DEFI         0.0       0.0       0.0     100.0   -3000.0    1000.0 proto
#
# --------------------------------------------------------------------------------------------------------------
PROTOTYPE       {gen}_B
# --------------------------------------------------------------------------------------------------------------
# prototype of the jaw:
FEDB_TAG        {gen}_B                # tag in the fedb (for the name of the file)
FEDB_SERIES     test                    # series of the fedb (for the name of the file)
# this rot-defi is only for local use:
ROT-DEFI         0.0       0.0       0.0       0.0   -3000.0    1000.0 proto
#
#
#
# --------------------------------------------------------------------------------------------------------------
ASSEMBLY        {gen}
# --------------------------------------------------------------------------------------------------------------
#
# needed bodies:
#               rename          rotname         name_in_file   fedb_series     fedb_tag       index
BODY            CONTAINO        CONTAINO        {gen}_T        test            {gen}_T        1
BODY            CONTAINI        CONTAINO        {gen}_TIN      test            {gen}_T        1
* jaw on positive x:
BODY            JAW_POS         JAW_POS         {gen}_B        test            {gen}_B        1
* jaw on negative x:
BODY            JAW_NEG         JAW_NEG         {gen}_B        test            {gen}_B        1
#
# define regions:
#               rename          material        rotbody         defition
REGION          *               EXTERNAL        *               -CONTAINO
* tank:
REGION          TANK            LATTICE         CONTAINO        +CONTAINO -CONTAINI
* vacuum between jaws:
REGION          INNERVAC        VACUUM          *               +CONTAINO -JAW_POS  -JAW_NEG
* jaw on positive x:
REGION          JAW_POS         LATTICE         JAW_POS         +JAW_POS
* jaw on negative x:
REGION          JAW_NEG         LATTICE         JAW_NEG         +JAW_NEG
#
ROT-DEFI         0.0       0.0       0.0       0.0       0.0       0.0  CONTAINO
#
ROT-DEFI         0.0       0.0       0.0       0.0       0.0       0.0  JAW_POS
* rotate by 180 deg the negative jaw:
ROT-DEFI       300.0       0.0     180.0       0.0       0.0       0.0  JAW_NEG
"""
        filepath = FlukaEnvironment.fedb_base / "assemblies" / f"test_{gen}.lbp"
        with filepath.open('w') as fp:
            fp.write(template_assembly)

def _create_generic_bodies(_generic_colls):
    for gen in _generic_colls.keys():
        template_body = f"""\
RPP {gen}_B   0.0 10.0 -5 5 -{_generic_colls[gen].length*100/2} {_generic_colls[gen].length*100/2}
"""
        filepath = FlukaEnvironment.fedb_base / "bodies" / f"test_{gen}_B.bodies"
        with filepath.open('w') as fp:
            # dump file in path FLUKA_environment.linebuilder / src
            fp.write(template_body)

        template_tank = f"""\
RPP {gen}_T   -15.0 15.0 -10 10 -{_generic_colls[gen].length*100/2} {_generic_colls[gen].length*100/2}
RPP {gen}_TIN -15.0 15.0 -10 10 -{_generic_colls[gen].length*100/2} {_generic_colls[gen].length*100/2}
"""
        filepath = FlukaEnvironment.fedb_base / "bodies" / f"test_{gen}_T.bodies"
        with filepath.open('w') as fp:
            # dump file in path FLUKA_environment.linebuilder / src
            fp.write(template_tank)
def _create_generic_regions(_generic_colls):
    for gen in _generic_colls.keys():
        template_body_reg = f"""\
{gen}_B     25 +{gen}_B
"""
        filepath = FlukaEnvironment.fedb_base / "regions" / f"test_{gen}_B.regions"
        with filepath.open('w') as fp:
            fp.write(template_body_reg)

        template_tank_reg = f"""\
{gen}_T     25 +{gen}_T -{gen}_TIN
{gen}_TIN   25 +{gen}_TIN
"""
        filepath = FlukaEnvironment.fedb_base / "regions" / f"test_{gen}_T.regions"
        with filepath.open('w') as fp:
            fp.write(template_tank_reg)

def _create_generic_materials(_generic_colls):
    for gen in _generic_colls.keys():
        mat = xcoll_to_fluka_material(_generic_colls[gen].material)
        template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat}    {gen}_B
"""
        filepath = FlukaEnvironment.fedb_base / "materials" / f"test_{gen}_B.assignmat"
        with filepath.open('w') as fp:
            fp.write(template_body_mat)
        template_tank_mat = f"""\
ASSIGNMA      VACUUM    {gen}_T
ASSIGNMA      VACUUM  {gen}_TIN
"""
        filepath = FlukaEnvironment.fedb_base / "materials" / f"test_{gen}_T.assignmat"
        with filepath.open('w') as fp:
            fp.write(template_tank_mat)

def _create_gen_files(_generic_colls):
    _create_generic_assemblies(_generic_colls)
    _create_generic_bodies(_generic_colls)
    _create_generic_regions(_generic_colls)
    _create_generic_materials(_generic_colls)

def _create_crystal_files(_cristals):
    _create_crystal_assemblies(_cristals)


def _create_generic_colls(element_dict):

    class GenericCollimator:
        def __init__(self, name, length, material):
            self.name = name
            self.length = length
            self.material = material
            self.collimator_list = []

    # Check collimator with material and go through cases with same length

    _generic_colls = {}
    for name, ee in element_dict.items():
        # if ee.generic:
        if ee.material:
            matched = False
            for generic in _generic_colls.values():
                if generic.length == ee.length and generic.material == ee.material:
                    generic.collimator_list.append(name)
                    matched = True
                    break
            if not matched:
                new_name = f"GEN{len(_generic_colls) + 2}"
                new_generic = GenericCollimator(new_name, ee.length, ee.material)
                new_generic.collimator_list.append(name)
                _generic_colls[new_name] = new_generic

    _create_gen_files(_generic_colls)


def _create_prototypes_file(element_dict, prototypes_file=None):
    for name, ee in element_dict.items():
        if ee.assembly is None:
            raise ValueError(f"Collimator {name} has no assembly!")
    if prototypes_file is None:
        FlukaAssembly.make_prototypes()
    else:
        prototypes_file = FsPath(prototypes_file).resolve()
        prototypes_file.copy_to(FsPath.cwd() / 'prototypes.lbp')
    FlukaAssembly.inspect_prototypes_file(prototypes_file)


# TODO check that prototype is valid and its sides
def _element_dict_to_fluka(element_dict, dump=False):
    collimator_dict = {}
    for name, ee in element_dict.items():
        if ee.length < 1.e-12:
            raise ValueError(f"Collimator {name} has zero length!")

        nsig = 1 # TODO can remove?
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
            tilt_1 = round(ee.tilt_L, 9)
            tilt_2 = round(ee.tilt_R, 9)
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
    # Save system state
    FlukaEnvironment.set_fedb_environment()
    file_path = FlukaEnvironment.linebuilder / "src" / "FLUKA_builder.py"
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
    FlukaEnvironment.unset_fedb_environment()

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
