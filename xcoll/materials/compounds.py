# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material
from .atoms import Carbon, Oxygen, Hydrogen, Boron, Nitrogen, Silicon, Aluminium, Iron, Nickel, Zinc
from .database import db, _manually_add_material_to_db


# General compounds
# =================

CarbonDioxide = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.951e-3, state='gas', temperature=273.15, pressure=1)
DryIce        = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.564,    state='solid')
Water         = Material(components=[Hydrogen, Oxygen], n_atoms=[2,1], density=1.0,      state='liquid')


# Compound definitions from FLUKA
# ===============================

BoronNitride   = Material(components=[Boron, Nitrogen],            n_atoms=[1, 1],             density=2.25, state='solid')
Silica         = Material(components=[Silicon, Oxygen],            n_atoms=[1, 2],             density=1.92, state='solid')
BoronTrioxide  = Material(components=[Boron, Oxygen],              n_atoms=[2, 3],             density=2.46, state='solid')
AluminiumOxide = Material(components=[Aluminium, Oxygen],          n_atoms=[2, 3],             density=3.98, state='solid')
NiZnFerrite    = Material(components=[Nickel, Zinc, Iron, Oxygen], n_atoms=[0.45, 0.55, 2, 4], density=5.2,  state='solid')


# Metadata for database
# =====================

_manually_add_material_to_db(CarbonDioxide,  'CarbonDioxide',  short_name='CO2')
_manually_add_material_to_db(DryIce,         'DryIce')
_manually_add_material_to_db(Water,          'Water',          short_name='H2O',  fluka_name='WATER')
_manually_add_material_to_db(BoronNitride,   'BoronNitride',                      fluka_name='hBNpure')
_manually_add_material_to_db(Silica,         'Silica',         short_name='Sand', fluka_name='SiO2')
_manually_add_material_to_db(BoronTrioxide,  'BoronTrioxide',                     fluka_name='B2O3')
_manually_add_material_to_db(AluminiumOxide, 'AluminiumOxide',                    fluka_name='Al2O3')
_manually_add_material_to_db(NiZnFerrite,    'NiZnFerrite',                       fluka_name='NIZNFER2')


# Clean up namespace
del Carbon, Oxygen, Hydrogen, Boron, Nitrogen, Silicon, Aluminium, Iron, Nickel, Zinc
del Material
del db
