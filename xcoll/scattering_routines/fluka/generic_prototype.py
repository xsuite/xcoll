# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .prototype import FlukaPrototype, FlukaAssembly

import numpy as np

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


class FlukaGenericAssembly(FlukaAssembly):
    _generic_required_fields = ['material', 'length']
    _generic_optional_fields = {'side': 'both', 'x_dim': 0.2, 'y_dim': 0.2}

    def __new__(cls, **kwargs):
        for field in cls._generic_required_fields:
            if field not in kwargs:
                raise ValueError(f"Need to provide {field}!")
        if kwargs.get('side') not in [None, 'both', 'left', 'right']:
            raise ValueError("Side must be 'both', 'left' or 'right'!")
        for prototype in FlukaPrototype._registry:
            if isinstance(prototype, FlukaGenericAssembly):
                found = True
                for field in cls._generic_required_fields:
                    if kwargs[field] != getattr(prototype, field):
                        found = False
                        break
                for field, opt_value in cls._generic_optional_fields.items():
                    if kwargs.get(field, opt_value) != getattr(prototype, field):
                        found = False
                        break
                if found:
                    return prototype
        # Create new generic assembly
        return super().__new__(cls, fedb_series='NOPE', fedb_tag='NOPENOPE', **kwargs)

    def __init__(self, **kwargs):
        if not hasattr(self, 'fedb_tag') and not hasattr(self, 'fedb_series'):
            self._init(kwargs)
            self._x_dim = kwargs['x_dim']
            self._y_dim = kwargs['y_dim']
            if self._x_dim > 0.25:
                self._x_dim = 0.25
            if self._y_dim > 0.25:
                self._y_dim = 0.25
            create_assembly_file(self.fedb_tag, self.side)
            create_body_file(self.fedb_tag, self.length, self.x_dim, self.y_dim)
            create_region_file(self.fedb_tag)
            create_material_file(self.fedb_tag, self.material)

    def _init(self, kwargs):
        num_generic_assemblies = len([p for p in FlukaPrototype._registry
                                    if isinstance(p, FlukaGenericAssembly)])
        kwargs['fedb_series'] = 'generic'
        kwargs['fedb_tag'] = f'GEN{num_generic_assemblies:03d}'
        for field, opt_value in self._generic_optional_fields.items():
            kwargs.setdefault(field, opt_value)
        super().__init__(**kwargs)

    @property
    def x_dim(self):
        return self._x_dim

    @property
    def y_dim(self):
        return self._y_dim


class FlukaGenericCrystalAssembly(FlukaGenericAssembly):
    _generic_required_fields = ['material', 'length', 'bending_radius', 'side']
    _generic_optional_fields = {'x_dim': 0.2, 'y_dim': 0.2}

    def __new__(cls,  *, bending_radius=None, **kwargs):
        if bending_radius is None:
            raise ValueError("Need to provide bending_radius!")
        kwargs['bending_radius'] = bending_radius
        return super().__new__(cls, **kwargs)

    def __init__(self, **kwargs):
        if not hasattr(self, 'fedb_tag') and not hasattr(self, 'fedb_series'):
            self._init(kwargs)
            self.bending_radius = kwargs['bending_radius']
            # Set dimensions
            self._x_dim = kwargs['x_dim']
            self._y_dim = kwargs['y_dim']
            # create files
            create_assembly_file(self.fedb_tag, self.side)
            create_crystal_body_file(self.fedb_tag, self.length, self.bending_radius,  self.x_dim, self.y_dim)
            create_crystal_region_file(self.fedb_tag)
            create_crystal_material_file(self.fedb_tag, self.material)


def create_assembly_file(fedb_tag, side):
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
    _write_file("assemblies", f"generic_{fedb_tag}.lbp",
                template_assembly)


def create_body_file(fedb_tag, length, x_dim, y_dim):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {100*x_dim} -{100*y_dim/2} {100*y_dim/2} -{length*100/2} {length*100/2}
"""
    _write_file("bodies", f"generic_{fedb_tag}_B.bodies",
                template_body)

    # Tank body should fit in blackhole (0.8m x 0.8m) for any angle, so maximally 0.8*sqrt(2)/2 = 0.565 for each side
    template_tank = f"""\
RPP {fedb_tag}_T  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
RPP {fedb_tag}_I  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
"""
    _write_file("bodies", f"generic_{fedb_tag}_T.bodies",
                template_tank)


def create_region_file(fedb_tag):
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


def create_material_file(fedb_tag, material):
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


def create_crystal_body_file(fedb_tag, length, bending_radius, x_dim, y_dim):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {x_dim*(100+10)} -{y_dim*(100+10)/2} {y_dim*(100+10)/2} 0.00001 {length*(100+20)}
ZCC {fedb_tag}Z1  {bending_radius*100} 0.0 {bending_radius*100}
ZCC {fedb_tag}Z2  {bending_radius*100} 0.0 {bending_radius*100-x_dim*100}
PLA {fedb_tag}P1  1.0 0.0 {np.cos(length/bending_radius)/np.sin(length/bending_radius)} {bending_radius*100} 0.0 0.0
"""
    _write_file("bodies", f"generic_{fedb_tag}_B.bodies",
                template_body)
    template_tank = f"""\
RPP {fedb_tag}_T  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
RPP {fedb_tag}_I  -{28} {28} -{28} {28} -{length*100 + 5} {length*100 + 5}
"""
    _write_file("bodies", f"generic_{fedb_tag}_T.bodies",
                template_tank)


def create_crystal_region_file(fedb_tag):
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


def create_crystal_material_file(fedb_tag, material):
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
    fluka.environment._add_to_index(directory, filename)
