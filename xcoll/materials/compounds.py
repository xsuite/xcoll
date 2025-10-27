# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material
from .atoms import Carbon, Oxygen, Hydrogen, Boron, Nitrogen, Silicon, Aluminium, Iron, Nickel, Zinc
from .database import db


# General compounds
# =================

CarbonDioxide = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.951e-3, state='gas', temperature=273.15, pressure=1)
DryIce        = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.564,    state='solid')
Water         = Material(components=[Hydrogen, Oxygen], n_atoms=[2,1], density=1.0,      state='liquid')


# Compound definitions from FLUKA
# ===============================

BoronNitride  = Material(components=[Boron, Nitrogen],            n_atoms=[1, 1],             density=2.25, state='solid')
Silica        = Material(components=[Silicon, Oxygen],            n_atoms=[1, 2],             density=1.92, state='solid')
BoronTrioxide = Material(components=[Boron, Oxygen],              n_atoms=[2, 3],             density=2.46, state='solid')
AluminumOxide = Material(components=[Aluminium, Oxygen],          n_atoms=[2, 3],             density=3.98, state='solid')
NiZnFerrite   = Material(components=[Nickel, Zinc, Iron, Oxygen], n_atoms=[0.45, 0.55, 2, 4], density=5.2,  state='solid')


# Metadata for database
# =====================

del Carbon, Oxygen, Hydrogen, Boron, Nitrogen, Silicon, Aluminium, Iron, Nickel, Zinc
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and obj.name is None:
        obj.name = name

db['Sand'] = Silica

Water.fluka_name = 'WATER'
BoronNitride.fluka_name = 'hBNpure'
Silica.fluka_name = 'SiO2'
BoronTrioxide.fluka_name = 'B2O3'
AluminumOxide.fluka_name = 'Al2O3'
NiZnFerrite.fluka_name = 'NIZNFER2'


# Clean up namespace
del name, obj
del Material
del db
