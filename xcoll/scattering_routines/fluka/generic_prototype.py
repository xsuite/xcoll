# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import atexit
import numpy as np

from .prototype import FlukaPrototype, FlukaAssembly
from ...compare import deep_equal


# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
#
# DOES NOT WORK WITH PARALLEL JOBS: Exit handler of one job will delete prototype
# as soon as finished, even if other job still running...
#
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# def exit_handler():
#     """Remove generic assemblies on exit."""
#     for assembly in FlukaPrototype._registry:
#         if isinstance(assembly, FlukaAssembly) \
#         and assembly.fedb_series == 'generic':
#             if len(assembly.elements) == 0:
#                 print(f"Deleting generic FlukaAssembly {assembly.fedb_tag}...")
#                 assembly.delete()
#     # If there are broken assemblies, there might be prototypes left
#     for prototype in FlukaPrototype._registry:
#         if isinstance(prototype, FlukaPrototype) \
#         and prototype.fedb_series == 'generic':
#             if len(prototype.elements) == 0:
#                 print(f"Deleting generic FlukaPrototype {prototype.fedb_tag}...")
#                 prototype.delete()
# atexit.register(exit_handler)


_generic_required_fields = ['material', 'length']
_generic_optional_fields = {'side': 'both', 'width': 0.2, 'height': 0.2}
_generic_crystal_required_fields = ['material', 'length', 'bending_radius']
_generic_crystal_optional_fields = {'side': 'left', 'width': 0.02, 'height': 0.05, 'tip_l': None, 'tip_mat': None}


def create_generic_assembly(**kwargs):
    _validate_kwargs(kwargs)
    # Check if the assembly already exists
    for prototype in FlukaPrototype._registry:
        if not isinstance(prototype, FlukaAssembly):
            continue
        if prototype.fedb_series == 'generic':
            found = True
            for field in _generic_required_fields + list(_generic_optional_fields.keys()):
                if not(deep_equal(kwargs[field], getattr(prototype, field))):
                    found = False
                    break
            if kwargs['is_crystal']:
                for field in _generic_crystal_required_fields + list(_generic_crystal_optional_fields.keys()):
                    if not(deep_equal(kwargs[field], getattr(prototype, field))):
                        found = False
                        break
            if found and prototype.exists():
                return prototype
    # Get an ID
    existing_ids = [int(p.fedb_tag[3:]) for p in FlukaPrototype._registry
                    if p.fedb_series == 'generic' and isinstance(p, FlukaAssembly)]
    max_id = max(existing_ids) if existing_ids else 0
    new_id = min(set(range(max_id+2)) - set(existing_ids))
    fedb_series = 'generic'
    fedb_tag = f'GEN{new_id:03d}'
    # Create the objects
    tank = FlukaPrototype(fedb_series, f'{fedb_tag}_T', _allow_generic=True)
    kwargs_body = kwargs.copy()
    kwargs_body.pop('side')
    body = FlukaPrototype(fedb_series, f'{fedb_tag}_B', _allow_generic=True, **kwargs_body)
    assm = FlukaAssembly(fedb_series, fedb_tag, _allow_generic=True, **kwargs)
    assm._prototypes = [tank, body]
    return assm


def _validate_kwargs(kwargs):
    kwargs.setdefault('is_crystal', False)
    from ...materials import _resolve_material
    kwargs['material'] = _resolve_material(kwargs['material'], ref='fluka', allow_none=False)
    if kwargs['is_crystal']:
        for field in _generic_crystal_required_fields:
            if field not in kwargs or kwargs[field] is None:
                raise ValueError(f"Need to provide {field}!")
        for field, opt_value in _generic_crystal_optional_fields.items():
            if field not in kwargs or kwargs[field] is None:
                kwargs[field] = opt_value
    else:
        for field in _generic_required_fields:
            if field not in kwargs or kwargs[field] is None:
                raise ValueError(f"Need to provide {field}!")
        for field, opt_value in _generic_optional_fields.items():
            if field not in kwargs or kwargs[field] is None:
                kwargs[field] = opt_value
    if kwargs.get('side') not in ['both', 'left', 'right']:
        raise ValueError(f"Side must be 'both', 'left' or 'right', but got {kwargs.get('side')}!")
    if kwargs['width'] > 0.25:
        kwargs['width'] = 0.25
    if kwargs['height'] > 0.25:
        kwargs['height'] = 0.25
    kwargs.pop('fedb_series', None) # fedb_series is always 'generic'
    kwargs.pop('fedb_tag', None)    # fedb_tag is always 'GENnnn'
    kwargs['angle'] = 0             # Only horizontal assembly (angle will be provided by LineBuilder)


def _assembly_file(fedb, fedb_tag, side, **kwargs):
    template_assembly = f"""\
# --------------------------------------------------------------------------------------------------------------
PROTOTYPE       {fedb_tag}_T
# --------------------------------------------------------------------------------------------------------------
# prototype of the tank:
FEDB_TAG        {fedb_tag}_T                # tag in the fedb (for the name of the file)
FEDB_SERIES     generic                # series of the fedb (for the name of the file)
# this rot-defi is only for local use:
ROT-DEFI         0.0       0.0       0.0     100.0   -3000.0    1000.0 proto
#
# --------------------------------------------------------------------------------------------------------------
PROTOTYPE       {fedb_tag}_B
# --------------------------------------------------------------------------------------------------------------
# prototype of the jaw:
FEDB_TAG        {fedb_tag}_B                # tag in the fedb (for the name of the file)
FEDB_SERIES     generic                # series of the fedb (for the name of the file)
# this rot-defi is only for local use:
ROT-DEFI         0.0       0.0       0.0       0.0   -3000.0    1000.0 proto
#
#
#
# --------------------------------------------------------------------------------------------------------------
ASSEMBLY        {fedb_tag}
# --------------------------------------------------------------------------------------------------------------
#
# needed bodies:
#           rename      rotname     name_in_file fedb_series fedb_tag    index
BODY        CONTAINO    CONTAINO    {fedb_tag}_T     generic     {fedb_tag}_T    1
BODY        CONTAINI    CONTAINO    {fedb_tag}_I     generic     {fedb_tag}_T    1"""
    if side in [None, 'both', 'left']:
        template_assembly += f"""
* jaw on positive x:
BODY        JAW_POS     JAW_POS     {fedb_tag}_B     generic     {fedb_tag}_B    1"""
    if side in [None, 'both', 'right']:
        template_assembly += f"""
* jaw on negative x:
BODY        JAW_NEG     JAW_NEG     {fedb_tag}_B     generic     {fedb_tag}_B    1"""
    template_assembly += f"""
#
# define regions:
#           rename      material    rotbody     defition
REGION      *           EXTERNAL    *           -CONTAINO
* tank:
REGION      TANK        LATTICE     CONTAINO    +CONTAINO -CONTAINI
* vacuum between jaws:"""
    subtract_jaws = ''
    if side in [None, 'both', 'left']:
        subtract_jaws += ' -JAW_POS'
    if side in [None, 'both', 'right']:
        subtract_jaws += ' -JAW_NEG'
    template_assembly += f"""
REGION      INNERVAC    VACUUM      *           +CONTAINI{subtract_jaws}"""
    if side in [None, 'both', 'left']:
        template_assembly += f"""
* jaw on positive x:
REGION      JAW_POS     LATTICE     JAW_POS     +JAW_POS"""
    if side in [None, 'both', 'right']:
        template_assembly += f"""
* jaw on negative x:
REGION      JAW_NEG     LATTICE     JAW_NEG     +JAW_NEG"""
    template_assembly += f"""
#
ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 CONTAINO"""
    if side in [None, 'both', 'left']:
        template_assembly += f"""
ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 JAW_POS"""
    if side in [None, 'both', 'right']:
        template_assembly += f"""
* rotate by 180 deg the negative jaw:
ROT-DEFI           300.0         0.0       180.0         0.0         0.0         0.0 JAW_NEG
"""
    return _write_file(fedb, "assemblies", f"generic_{fedb_tag}.lbp",
                       template_assembly)


def _inp_prot_file(fedb, fedb_tag, length, material, width, height, **kwargs):

    bodies_start =f"""\
TITLE
Test element
GLOBAL                                         1.0       1.0
DEFAULTS                                                              NEW-DEFA
BEAM
BEAMPOS
GEOBEGIN                                                              COMBNAME
    0    0
* Black body
SPH blkbody    0.0 0.0 0.0 10000000.0
* Void sphere
SPH void       0.0 0.0 0.0 100000.0
* Bodies from the include file
* START_CUT_BODIES
"""
    bodies_end =f"""\
* END_CUT_BODIES
END
* Black hole
BLKBODY      5 +blkbody -void
* Void around
VOID         5 +void -TCP_BODY
"""
    region_start = f"""\
END
* Black hole
BLKBODY      5 +blkbody -void
* Void around
VOID         5 +void -TCP_BODY
* Region from the include file
* START_CUT_REGIONS
"""
    region_end = f"""\
* END_CUT_REGIONS
"""
    mat_start = f"""\
ASSIGNMA    BLCKHOLE   BLKBODY
ASSIGNMA      VACUUM      VOID
#include ../materials/materials.inp
* START_CUT_MATERIALS
"""
    mat_end = f"""\
* END_CUT_MATERIALS
"""
    end_template = f"""\
RANDOMIZ         1.0
START
STOP
"""

    if kwargs["is_crystal"]:
        body_file, tank_file = _crystal_body_file(fedb, fedb_tag, length, 
                                kwargs['bending_radius'],
                                width, height,                  
                                bodies_start = bodies_start,  bodies_end = bodies_end)
        body_region_file, tank_region_file = _crystal_region_file(fedb, fedb_tag,
                                                region_start = region_start, region_end = region_end)
        body_mat_file, tank_mat_file = _crystal_material_file(fedb, fedb_tag, material, 
                                          mat_start = mat_start, mat_end = mat_end)
    else:
        body_file, tank_file = _body_file(fedb, fedb_tag, length, width, height, 
                bodies_start = bodies_start,  bodies_end = bodies_end, tip_l = kwargs['tip_l'])
        body_region_file, tank_region_file = _region_file(fedb, fedb_tag, 
                region_start = region_start, region_end = region_end, tip_l = kwargs['tip_l'])
        body_mat_file, tank_mat_file = _material_file(fedb, fedb_tag, material,
                mat_start = mat_start, mat_end = mat_end, tip_mat = kwargs['tip_mat'])

    inp_body_file = body_file + body_region_file + body_mat_file
    inp_tank_file = tank_file + tank_region_file + tank_mat_file

    # XXX end needed to have a complete fluka inp (for flair)
    body_file = inp_body_file + end_template
    tank_file = inp_tank_file + end_template

    body_file = _write_file(fedb, "prototypes", f"generic_{fedb_tag}_B.inp",
                                body_file)
    tank_file = _write_file(fedb, "prototypes", f"generic_{fedb_tag}_T.inp",
                                tank_file)

    return body_file, tank_file

def _body_file(fedb, fedb_tag, length, width, height, **kwargs):

    tip_l = kwargs['tip_l']
    template_body = f"""\
*
RPP {fedb_tag}_B   0.0 {100*width} -{100*height/2} {100*height/2} -{length*100/2} {length*100/2}
*
"""
    if tip_l:
        template_body += f"""\
YZP   tip      {tip_l*100}
"""

    # Tank body should fit in blackhole (0.8m x 0.8m) for any angle, so maximally 0.8*sqrt(2)/2 = 0.565 for each side
    template_tank = f"""\
RPP {fedb_tag}_T  -28 28 -28 28 -{length*100/2 + 5} {length*100/2 + 5}
RPP {fedb_tag}_I  -28 28 -28 28 -{length*100/2 + 5} {length*100/2 + 5}
* For BB, need less margin!
*RPP {fedb_tag}_T  -28 28 -28 28 -{length*100/2 + 1e-12} {length*100/2 + 1e-12}
*RPP {fedb_tag}_I  -28 28 -28 28 -{length*100/2 + 1e-12} {length*100/2 + 1e-12}
"""
    body_file = kwargs["bodies_start"] + template_body + kwargs["bodies_end"]
    tank_file = kwargs["bodies_start"] + template_tank + kwargs["bodies_end"]

    return body_file, tank_file


def _region_file(fedb, fedb_tag, **kwargs):

    if kwargs['tip_l']:
        template_body_reg = f"""\
{fedb_tag}_B     5 +{fedb_tag}_B -tip
{fedb_tag}_C     5 +{fedb_tag}_B +tip
"""
    else:
        template_body_reg = f"""\
{fedb_tag}_B     5 +{fedb_tag}_B
"""
#     template_body_reg = f"""\
# {fedb_tag}_B     5 +{fedb_tag}_B -tip
# {fedb_tag}_C     5 +{fedb_tag}_B +tip
# """
    template_tank_reg = f"""\
{fedb_tag}_T     5 +{fedb_tag}_T -{fedb_tag}_I
{fedb_tag}_I     5 +{fedb_tag}_I
"""
    body_file = kwargs["region_start"] + template_body_reg + kwargs["region_end"]
    tank_file = kwargs["region_start"] + template_tank_reg + kwargs["region_end"]

    return body_file, tank_file


def _material_file(fedb, fedb_tag, material, **kwargs):
    mat = material.fluka_name
    tip_mat = kwargs['tip_mat'].fluka_name

    template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat:>8}  {fedb_tag:>6}_B
"""
    if tip_mat:
        template_body_mat += f"""\
ASSIGNMA    {tip_mat:>8}  {fedb_tag:>6}_C
"""

    template_tank_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA      VACUUM  {fedb_tag:>6}_T
ASSIGNMA      VACUUM  {fedb_tag:>6}_I
"""
    body_file = kwargs["mat_start"] + template_body_mat + kwargs["mat_end"]
    tank_file = kwargs["mat_start"] + template_tank_mat + kwargs["mat_end"]

    return body_file, tank_file


def _crystal_body_file(fedb, fedb_tag, length, bending_radius, width, height, **kwargs):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {width*(100+10)} -{height*(100+10)/2} {height*(100+10)/2} -{length*(100+20)} {length*(100+20)}
YCC {fedb_tag}Z1  0.0 {bending_radius*100} {bending_radius*100}
YCC {fedb_tag}Z2  0.0 {bending_radius*100} {bending_radius*100-width*100}
PLA {fedb_tag}P1  1.0 0.0 {np.cos(length/bending_radius)/np.sin(length/bending_radius)} {bending_radius*100} 0.0 0.0
XYP {fedb_tag}P2  0.0
"""
    template_tank = f"""\
RPP {fedb_tag}_T  -28 28 -28 28 -{length*(100+20)/2 + 5} {length*(100+20)/2 + 5}
RPP {fedb_tag}_I  -28 28 -28 28 -{length*(100+20)/2 + 5} {length*(100+20)/2 + 5}
"""
    body_file = kwargs["bodies_start"] + template_body + kwargs["bodies_end"]
    tank_file = kwargs["bodies_start"] + template_tank + kwargs["bodies_end"]

    return body_file, tank_file


def _crystal_region_file(fedb, fedb_tag, **kwargs):
    template_body_reg = f"""\
{fedb_tag}_B     5 | +{fedb_tag}_B +{fedb_tag}Z1 -{fedb_tag}Z2 +{fedb_tag}P1 - {fedb_tag}P2
{fedb_tag}B2     5 | +{fedb_tag}_B +{fedb_tag}Z2
                   | +{fedb_tag}_B -{fedb_tag}Z1
                   | +{fedb_tag}_B -{fedb_tag}P1
                   | +{fedb_tag}_B +{fedb_tag}P2 -{fedb_tag}Z2
"""
    template_tank_reg = f"""\
{fedb_tag}_T     5 | +{fedb_tag}_T -{fedb_tag}_I
{fedb_tag}_I     5 | +{fedb_tag}_I
"""

    body_file = kwargs["region_start"] + template_body_reg + kwargs["region_end"]
    tank_file = kwargs["region_start"] + template_tank_reg + kwargs["region_end"]

    return body_file, tank_file


def _crystal_material_file(fedb, fedb_tag, material, **kwargs):
    mat = material.fluka_name
    template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat:>8}  {fedb_tag:>6}_B
ASSIGNMA      VACUUM  {fedb_tag:>6}B2
"""
    
    template_tank_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA      VACUUM  {fedb_tag:>6}_T
ASSIGNMA      VACUUM  {fedb_tag:>6}_I
"""

    body_file = kwargs["mat_start"] + template_body_mat + kwargs["mat_end"]
    tank_file = kwargs["mat_start"] + template_tank_mat + kwargs["mat_end"]

    return body_file, tank_file


def _write_file(fedb, directory, filename, content):
    path = fedb / directory / filename
    with path.open('w') as fid:
        fid.write(content)
    return path
