# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .database import db
from .atoms import Helium, Nitrogen, Oxygen, Argon, Carbon


# Liquid gases
LiquidHelium   = Helium.adapt(name='LiquidHelium',     state='liquid', density=0.146,  temperature=1.9,   pressure=None)
LiquidNitrogen = Nitrogen.adapt(name='LiquidNitrogen', state='liquid', density=0.808,  temperature=77.36, pressure=None)
LiquidOxygen   = Oxygen.adapt(name='LiquidOxygen',     state='liquid', density=1.141,  temperature=90.20, pressure=None)
LiquidArgon    = Argon.adapt(name='LiquidArgon',       state='liquid', density=1.3954, temperature=87.3,  pressure=None)


# Carbon variants
Carbon140 = Carbon.adapt(name='Carbon140', density=1.40, fluka_name='CC_1_40')
Carbon175 = Carbon.adapt(name='Carbon175', density=1.75, fluka_name='CC_1_75')
Carbon180 = Carbon.adapt(name='Carbon180', density=1.80, fluka_name='CC_1_80')
Carbon185 = Carbon.adapt(name='Carbon185', density=1.85, fluka_name='CC_1_85')
Diamond   = Carbon.adapt(name='Diamond', density=3.52)
GraphiteR4550 = Carbon.adapt(name='GraphiteR4550', density=1.83, excitation_energy=78.0e-9,
                             nuclear_radius=0.25, nuclear_elastic_slope=70.0,
                             cross_section=[0.337, 0.232, 0, 0, 0, 0.0076e-2], fluka_name='GRAR4550')
CarbonFibreCarbon = GraphiteR4550.adapt(name='CarbonFibreCarbon', density=1.67, fluka_name='AC150GPH')
db['CFC'] = CarbonFibreCarbon

# Was in SixTrack, not sure what it is
Carbon2 = GraphiteR4550.adapt(name='Carbon2', density=4.52)
db['C2'] = Carbon2
