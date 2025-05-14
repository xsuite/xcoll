# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import atexit
import numpy as np

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .prototype import FlukaPrototype, FlukaAssembly


def xcoll_to_fluka_material(material):
    # XXX EXPAND THIS DICT

    material_dict = {
        "iner": "INERM180",
        "c":    "  CARBON",
        "si":   " SILICON",
        "cu":   "  COPPER",
        "mogr": "AC150GPH"
    }

    if material not in material_dict.keys():
        raise ValueError(f"Material {material} not found in material dictionary!")
    else:
        return material_dict[material]


def exit_handler():
    """Remove generic assemblies on exit."""
    for assembly in FlukaPrototype._registry:
        if isinstance(assembly, FlukaAssembly) \
        and assembly.fedb_series == 'generic':
            for el in assembly.elements:
                assembly.remove_element(el)
            assembly.delete()
atexit.register(exit_handler)


_generic_required_fields = ['material', 'length']
_generic_optional_fields = {'side': 'both', 'width': 0.2, 'height': 0.2}
_generic_crystal_required_fields = ['material', 'length', 'bending_radius']
_generic_crystal_optional_fields = {'side': 'left', 'width': 0.02, 'height': 0.05}


def create_generic_assembly(**kwargs):
    _validate_kwargs(kwargs)
    # Check if the assembly already exists
    for prototype in FlukaPrototype._registry:
        if prototype.fedb_series == 'generic':
            found = True
            for field in _generic_required_fields:
                if kwargs[field] != getattr(prototype, field):
                    found = False
                    break
            for field, opt_value in _generic_optional_fields.items():
                if kwargs.get(field, opt_value) != getattr(prototype, field):
                    found = False
                    break
            if found:
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
    ass = FlukaAssembly(fedb_series, fedb_tag, _allow_generic=True, **kwargs)
    # Create and assign the files
    if kwargs['is_crystal']:
            body_file, tank_file = _crystal_body_file(fedb_tag, **kwargs)
            _crystal_region_file(fedb_tag, **kwargs)
            _crystal_material_file(fedb_tag, **kwargs)
    else:
        body_file, tank_file = _body_file(fedb_tag, **kwargs)
        _region_file(fedb_tag, **kwargs)
        _material_file(fedb_tag, **kwargs)
    body.body_file = body_file
    tank.body_file = tank_file
    ass.assembly_file = _assembly_file(fedb_tag, **kwargs)
    return ass


def _validate_kwargs(kwargs):
    kwargs.setdefault('is_crystal', False)
    if kwargs['is_crystal']:
        for field in _generic_crystal_required_fields:
            if field not in kwargs:
                raise ValueError(f"Need to provide {field}!")
        for field, opt_value in _generic_crystal_optional_fields.items():
            kwargs.setdefault(field, opt_value)
    else:
        for field in _generic_required_fields:
            if field not in kwargs:
                raise ValueError(f"Need to provide {field}!")
        for field, opt_value in _generic_optional_fields.items():
            kwargs.setdefault(field, opt_value)
    if kwargs.get('side') not in ['both', 'left', 'right']:
        raise ValueError("Side must be 'both', 'left' or 'right'!")
    if kwargs['width'] > 0.25:
        kwargs['width'] = 0.25
    if kwargs['height'] > 0.25:
        kwargs['height'] = 0.25
    kwargs.pop('fedb_series', None) # fedb_series is always 'generic'
    kwargs.pop('fedb_tag', None)    # fedb_tag is always 'GENnnn'
    kwargs['angle'] = 0             # Only horizontal assembly (angle will be provided by LineBuilder)


def _assembly_file(fedb_tag, side, **kwargs):
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
    return _write_file("assemblies", f"generic_{fedb_tag}.lbp",
                       template_assembly)


def _body_file(fedb_tag, length, width, height, **kwargs):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {100*width} -{100*height/2} {100*height/2} -{length*100/2} {length*100/2}
"""
    body_file = _write_file("bodies", f"generic_{fedb_tag}_B.bodies",
                        template_body)

    # Tank body should fit in blackhole (0.8m x 0.8m) for any angle, so maximally 0.8*sqrt(2)/2 = 0.565 for each side
    template_tank = f"""\
RPP {fedb_tag}_T  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
RPP {fedb_tag}_I  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
"""
    tank_file = _write_file("bodies", f"generic_{fedb_tag}_T.bodies",
                        template_tank)
    return body_file, tank_file


def _region_file(fedb_tag, **kwargs):
    template_body_reg = f"""\
{fedb_tag}_B     5 +{fedb_tag}_B
"""
    _write_file("regions", f"generic_{fedb_tag}_B.regions",
                template_body_reg)

    template_tank_reg = f"""\
{fedb_tag}_T     5 +{fedb_tag}_T -{fedb_tag}_I
{fedb_tag}_I     5 +{fedb_tag}_I
"""
    _write_file("regions", f"generic_{fedb_tag}_T.regions",
                template_tank_reg)


def _material_file(fedb_tag, material, **kwargs):
    mat = xcoll_to_fluka_material(material)
    template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat:>8}  {fedb_tag:>6}_B
"""
    _write_file("materials", f"generic_{fedb_tag}_B.assignmat",
                template_body_mat)

    template_tank_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA      VACUUM  {fedb_tag:>6}_T
ASSIGNMA      VACUUM  {fedb_tag:>6}_I
"""
    _write_file("materials", f"generic_{fedb_tag}_T.assignmat",
                template_tank_mat)


def _crystal_body_file(fedb_tag, length, bending_radius, width, height, **kwargs):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {width*(100+10)} -{height*(100+10)/2} {height*(100+10)/2} 0.00001 {length*(100+20)}
ZCC {fedb_tag}Z1  {bending_radius*100} 0.0 {bending_radius*100}
ZCC {fedb_tag}Z2  {bending_radius*100} 0.0 {bending_radius*100-width*100}
PLA {fedb_tag}P1  1.0 0.0 {np.cos(length/bending_radius)/np.sin(length/bending_radius)} {bending_radius*100} 0.0 0.0
"""
    body_file = _write_file("bodies", f"generic_{fedb_tag}_B.bodies",
                            template_body)
    template_tank = f"""\
RPP {fedb_tag}_T  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
RPP {fedb_tag}_I  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
"""
    tank_file = _write_file("bodies", f"generic_{fedb_tag}_T.bodies",
                            template_tank)
    return body_file, tank_file


def _crystal_region_file(fedb_tag, **kwargs):
    template_body_reg = f"""\
{fedb_tag}_B     5 | +{fedb_tag}_B +{fedb_tag}Z1 -{fedb_tag}Z2 +{fedb_tag}P1
{fedb_tag}B2     5 | +{fedb_tag}_B +{fedb_tag}Z2
                   | +{fedb_tag}_B -{fedb_tag}Z1
                   | +{fedb_tag}_B -{fedb_tag}P1
"""
    _write_file("regions", f"generic_{fedb_tag}_B.regions",
                template_body_reg)
    template_body_tank = f"""\
{fedb_tag}_T     5 | +{fedb_tag}_T -{fedb_tag}_I
{fedb_tag}_I     5 | +{fedb_tag}_I
"""
    _write_file("regions", f"generic_{fedb_tag}_T.regions",
                template_body_tank)


def _crystal_material_file(fedb_tag, material, **kwargs):
    mat = xcoll_to_fluka_material(material)
    template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat:>8}  {fedb_tag:>6}_B
ASSIGNMA     VACUUM  {fedb_tag:>6}B2
"""
    _write_file("materials", f"generic_{fedb_tag}_B.assignmat",
                template_body_mat)
    template_tank_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA      VACUUM  {fedb_tag:>6}_T
ASSIGNMA      VACUUM  {fedb_tag:>6}_I
"""
    _write_file("materials", f"generic_{fedb_tag}_T.assignmat",
                template_tank_mat)


def _write_file(directory, filename, content):
    from xcoll import fluka
    path = fluka.environment.fedb / directory / filename
    with path.open('w') as fid:
        fid.write(content)
    return path
