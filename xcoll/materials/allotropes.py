# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .atoms import Helium, Nitrogen, Oxygen, Argon, Carbon
from .material import Material
from .database import db

# Liquid gases
LiquidHelium   = Helium.adapt(state='liquid',   density=0.146,  temperature=1.9,   pressure=None)
LiquidNitrogen = Nitrogen.adapt(state='liquid', density=0.808,  temperature=77.36, pressure=None)
LiquidOxygen   = Oxygen.adapt(state='liquid',   density=1.141,  temperature=90.20, pressure=None)
LiquidArgon    = Argon.adapt(state='liquid',    density=1.3954, temperature=87.3,  pressure=None)


# Carbon variants (need to provide name explicitly to avoid fluka_name being written as main name)
Carbon140 = Carbon.adapt(name='Carbon140', density=1.40, fluka_name='CC_1_40')
Carbon175 = Carbon.adapt(name='Carbon175', density=1.75, fluka_name='CC_1_75')
Carbon180 = Carbon.adapt(name='Carbon180', density=1.80, fluka_name='CC_1_80')
Carbon185 = Carbon.adapt(name='Carbon185', density=1.85, fluka_name='CC_1_85')
Diamond   = Carbon.adapt(name='Diamond', density=3.52, excitation_energy=88.5)
GraphiteR4550 = Carbon.adapt(name='GraphiteR4550', density=1.83, excitation_energy=78.0, fluka_name='GRAR4550')
CarbonFibreCarbon = GraphiteR4550.adapt(name='CarbonFibreCarbon', density=1.67, fluka_name='AC150GPH')
db['CFC'] = CarbonFibreCarbon


# Metadata for database
# =====================

del Helium, Nitrogen, Oxygen, Argon, Carbon
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and obj.name is None:
        obj.name = name


# Clean up namespace
del name, obj
del Material
del db
