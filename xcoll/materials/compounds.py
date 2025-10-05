# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import CompoundMaterial
from .atoms import *
from .database import db


CarbonDioxide = CompoundMaterial(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.951e-3, state='gas', temperature=273.15, pressure=1)
DryIce        = CompoundMaterial(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.564, state='solid')
Water         = CompoundMaterial(components=[Hydrogen, Oxygen], n_atoms=[2,1], density=1.0, state='liquid')

# From FLUKA
BoronNitride  = CompoundMaterial(components=[Boron, Nitrogen],              n_atoms=[1,1], density=2.25, state='solid')
Silica        = CompoundMaterial(components=[Silicon, Oxygen],              n_atoms=[1,2], density=1.92, state='solid')
BoronTrioxide = CompoundMaterial(components=[Boron, Oxygen],                n_atoms=[2,3], density=2.46, state='solid')
AluminumOxide = CompoundMaterial(components=[Aluminium, Oxygen],            n_atoms=[2,3], density=3.98, state='solid')
NiZnFerrite   = CompoundMaterial(components=[Nickel, Zinc, Iron, Oxygen],   n_atoms=[0.45,0.55,2,4], density=5.2, state='solid')
CopperDiamond = CompoundMaterial(components=[Copper, Carbon, Boron],        n_atoms=[0.2359, 0.7548, 0.0093], density=5.4, state='solid')
MoGr6400      = CompoundMaterial(components=[Molybdenum, Carbon],           n_atoms=[0.0146, 0.9854], density=2.48, state='solid')
MoGr          = CompoundMaterial(components=[Molybdenum, Carbon, Titanium], n_atoms=[0.0184, 0.9809, 0.07], density=2.55, state='solid')
TiZrMo        = CompoundMaterial(components=[Molybdenum, Carbon, Titanium, Zirconium], density=10.22, state='solid',
                                    n_atoms=[0.98762, 0.00158, 0.00996, 0.00084])

# Give name and store in database
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, CompoundMaterial):
        obj.name = name

db['Sand'] = Silica

Water.fluka_name = 'WATER'
BoronNitride.fluka_name = 'hBNpure'
Silica.fluka_name = 'SiO2'
BoronTrioxide.fluka_name = 'B2O3'
AluminumOxide.fluka_name = 'Al2O3'
NiZnFerrite.fluka_name = 'NIZNFER2'
CopperDiamond.fluka_name = 'CUDIAM75'
MoGr.fluka_name = 'MG6403Fc'
MoGr6400.fluka_name = 'MoGRMG64'
TiZrMo.fluka_name = 'TZM'
