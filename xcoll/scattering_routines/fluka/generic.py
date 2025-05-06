# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .environment import FlukaEnvironment
from .prototypes import FlukaAssembly


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


class FlukaGenericAssembly(xc.FlukaAssembly):
    def __new__(cls,  *, material=None, length=None, side='both', **kwargs):
        if material is None:
            raise ValueError("Need to provide material!")
        if length is None:
            raise ValueError("Need to provide length!")
        if side not in [None, 'both', 'left', 'right']:
            raise ValueError("Side must be 'both', 'left' or 'right'!")
        for prototype in FlukaPrototype._registry:
            if isinstance(prototype, FlukaGenericAssembly):
                if prototype.material == material and prototype.length == length \
                and if prototype.side == side:
                    return prototype
        # Create new generic assembly
        return super().__new__(cls, fedb_series='NOPE', fedb_tag='NOPENOPE', **kwargs)


    def __init__(self, **kwargs):
        num_generic_assemblies = len([p for p in FlukaPrototype._registry
                                      if isinstance(p, FlukaGenericAssembly)])
        kwargs['fedb_series'] = 'generic'
        kwargs['fedb_tag'] = f'GEN{num_generic_assemblies:03d}'
        super().__init__(**kwargs)
        self._x_dim = kwargs.get('x_dim', 1)
        self._y_dim = kwargs.get('y_dim', 1)
        create_fedb_files(fedb_tag, self.material, self.length, self.side, self.x_dim, self.y_dim)

    @property
    def x_dim(self):
        return self._x_dim

    @property
    def y_dim(self):
        return self._y_dim


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
REGION      INNERVAC    VACUUM      *           +CONTAINO{subtract_jaws}"""
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
    filepath = FlukaEnvironment.fedb_base / "assemblies" / f"generic_{fedb_tag}.lbp"
    with filepath.open('w') as fp:
        fp.write(template_assembly)


def create_body_file(fedb_tag, length, x_dim, y_dim):
    template_body = f"""\
RPP {fedb_tag}_B   0.0 {100*x_dim} -{100*y_dim/2} {100*y_dim/2} -{length*100/2} {length*100/2}
"""
    filepath = FlukaEnvironment.fedb_base / "bodies" / f"generic_{fedb_tag}_B.bodies"
    with filepath.open('w') as fp:
        fp.write(template_body)

    template_tank = f"""\
RPP {fedb_tag}_T  -{200*x_dim} {200*x_dim} -{100*y_dim} {100*y_dim} -{length*120} {length*120}
RPP {fedb_tag}_I  -{200*x_dim} {200*x_dim} -{100*y_dim} {100*y_dim} -{length*120} {length*120}
"""
    filepath = FlukaEnvironment.fedb_base / "bodies" / f"generic_{fedb_tag}_T.bodies"
    with filepath.open('w') as fp:
        fp.write(template_tank)


def create_region_file(fedb_tag):
    template_body_reg = f"""\
{fedb_tag}_B     5 +{fedb_tag}_B
"""
    filepath = FlukaEnvironment.fedb_base / "regions" / f"generic_{fedb_tag}_B.regions"
    with filepath.open('w') as fp:
        fp.write(template_body_reg)

    template_tank_reg = f"""\
{fedb_tag}_T     5 +{fedb_tag}_T -{fedb_tag}_I
{fedb_tag}_I     5 +{fedb_tag}_I
"""
    filepath = FlukaEnvironment.fedb_base / "regions" / f"generic_{fedb_tag}_T.regions"
    with filepath.open('w') as fp:
        fp.write(template_tank_reg)


def create_material_file(fedb_tag, material):
    mat = xcoll_to_fluka_material(material)
    template_body_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA    {mat:8}   {fedb_tag}_B
"""
    filepath = FlukaEnvironment.fedb_base / "materials" / f"generic_{fedb_tag}_B.assignmat"
    with filepath.open('w') as fp:
        fp.write(template_body_mat)

    template_tank_mat = f"""\
* ..+....1....+....2....+....3....+....4....+....5....+....6....+....7..
ASSIGNMA      VACUUM  {fedb_tag}_T
ASSIGNMA      VACUUM  {fedb_tag}_I
"""
    filepath = FlukaEnvironment.fedb_base / "materials" / f"generic_{fedb_tag}_T.assignmat"
    with filepath.open('w') as fp:
        fp.write(template_tank_mat)
