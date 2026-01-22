# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .atoms import Helium, Nitrogen, Oxygen, Argon, Carbon
from .material import Material
from .database import db, _manually_add_material_to_db


# Liquid gases
LiquidHelium   = Helium.adapt(state='liquid',   density=0.146,  temperature=1.9,   pressure=None)
LiquidNitrogen = Nitrogen.adapt(state='liquid', density=0.808,  temperature=77.36, pressure=None)
LiquidOxygen   = Oxygen.adapt(state='liquid',   density=1.141,  temperature=90.20, pressure=None)
LiquidArgon    = Argon.adapt(state='liquid',    density=1.3954, temperature=87.3,  pressure=None)


# Carbon variants
Carbon140 = Carbon.adapt(density=1.40)
Carbon175 = Carbon.adapt(density=1.75)
Carbon180 = Carbon.adapt(density=1.80)
Carbon185 = Carbon.adapt(density=1.85)
Diamond   = Carbon.adapt(density=3.52, excitation_energy=88.5)
GraphiteR4550 = Carbon.adapt(density=1.83, excitation_energy=78.0)
CarbonFibreCarbon = GraphiteR4550.adapt(density=1.67, nuclear_radius=0.25, nuclear_elastic_slope=70.0,
                                        cross_section=[0.337, 0.232, 0, 0, 0, 0.0076e-2], hcut=0.02,
                                        crystal_plane_distance=0.63e-7, eta=0.9,
                                        crystal_potential=21.0, nuclear_collision_length=1.e-12)

# Metadata for database
# =====================

_manually_add_material_to_db(LiquidHelium,      'LiquidHelium')
_manually_add_material_to_db(LiquidNitrogen,    'LiquidNitrogen')
_manually_add_material_to_db(LiquidOxygen,      'LiquidOxygen')
_manually_add_material_to_db(LiquidArgon,       'LiquidArgon')
_manually_add_material_to_db(Carbon140,         'Carbon140',                           fluka_name='CC_1_40')
_manually_add_material_to_db(Carbon175,         'Carbon175',                           fluka_name='CC_1_75')
_manually_add_material_to_db(Carbon180,         'Carbon180',                           fluka_name='CC_1_80')
_manually_add_material_to_db(Carbon185,         'Carbon185',                           fluka_name='CC_1_85')
_manually_add_material_to_db(Diamond,           'Diamond')
_manually_add_material_to_db(GraphiteR4550,     'GraphiteR4550',                       fluka_name='GRAR4550')
_manually_add_material_to_db(CarbonFibreCarbon, 'CarbonFibreCarbon', short_name='CFC', fluka_name='AC150GPH')


Carbon140.info         = 'Carbon-based material. Used in LHC dump.'
Carbon175.info         = 'Carbon-based material. Used in LHC dump.'
Carbon180.info         = 'Carbon-based material. Used in LHC dump.'
Carbon185.info         = 'Carbon-based material. Used in LHC dump.'
GraphiteR4550.info     = 'Graphite material. Was a candidate for LHC collimators.'
CarbonFibreCarbon.info = 'Carbon-fibre composite. Used in LHC collimators.'


# Clean up namespace
del Helium, Nitrogen, Oxygen, Argon, Carbon
del Material
del db
