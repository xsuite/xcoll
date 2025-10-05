# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import MixtureMaterial
from .atoms import *
from .compounds import *


# Air = MixtureMaterial(name='Air', components=[Nitrogen, Oxygen, Argon, CarbonDioxide], state='gas',
#                       mass_fractions=[78.084, 20.946, 0.934, 0.036], density=0.001225, pressure=101.325, temperature=15.0)

# INERM180 : matdef, density=18.0, components=["W", "Ni", "Cu"], componentsFractions={0.95, 0.035, 0.015};
# MG6403Fc : matdef, density=2.55, components=["Mo", "C", "Ti"], componentsFractions={0.130012156609076, 0.8675206104800844, 0.002467232910839576};
# AC150GPH : matdef, Z=6.0, A=12.011, density=1.67;
# GPH      : matdef, Z=6.0, A=12.011, density=1.83;
# CUDIAM75 : matdef, density=5.4, components=["Cu", "C", "B"], componentsFractions={0.6205462018968333, 0.3752917530367507, 0.004162045066415989};



# From FLUKA
BoronNitride5000    = MixtureMaterial(components=[BoronNitride, Silica, Oxygen, BoronTrioxide, Calcium],
                                      state='solid', mass_fractions=[0.929, 0.002, 0.04, 0.004, 0.025], density=1.925)
FeNiCo              = MixtureMaterial(components=[Iron, Nickel, Cobalt], state='solid',
                                      mass_fractions=[0.54, 0.29, 0.17], density=8.35)
GlidCop15           = MixtureMaterial(components=[Copper, AluminumOxide], state='solid',
                                      mass_fractions=[0.9985, 0.0015], density=8.9)
CuNiFeMn            = MixtureMaterial(components=[Copper, Nickel, Iron, Manganese], state='solid',
                                      mass_fractions=[0.875, 0.1, 0.015, 0.01], density=8.87)
Inermet180          = MixtureMaterial(components=[Tungsten, Nickel, Copper], state='solid',
                                      mass_fractions=[0.95, 0.035, 0.015], density=18.0)
StainLessSteel304L  = MixtureMaterial(components=[Iron, Chromium, Nickel, Manganese, Silicon, Phosphorus,
                                       Sulfur, Carbon], state='solid', mass_fractions=[0.67145, 0.185, 0.1125, 0.02, 0.01, 0.00045,
                                       0.0003, 0.0003], density=8.02)
StainLessSteel316L  = MixtureMaterial(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon,
                                       Phosphorus, Sulfur, Carbon], state='solid', mass_fractions=[0.61395, 0.185, 0.14, 0.03, 0.02,
                                       0.01, 0.00045, 0.0003, 0.0003], density=8.0)
StainLessSteel316LN = MixtureMaterial(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon,
                                       Titanium, Nitrogen, Niobium, Copper, Cobalt, Phosphorus, Carbon, Sulfur, Tantalum, Boron],
                                      state='solid', mass_fractions=[0.65093, 0.17, 0.12, 0.025, 0.02, 0.0075, 0.0015, 0.0014,
                                       0.001, 0.001, 0.0005, 0.00045, 0.0003, 0.0003, 0.0001, 0.00002], density=8.03)
TitaniumNitride     = MixtureMaterial(components=[Titanium, Nitrogen], state='solid',
                                      mass_fractions=[0.77, 0.23], density=5.22)
Antico              = MixtureMaterial(components=[Aluminium, Silicon, Iron, Magnesium, Manganese, Chromium,
                                       Titanium, Zinc, Copper], state='solid', density=2.7,
                                      mass_fractions=[0.9725, 0.006, 0.005, 0.004, 0.004, 0.0035, 0.002, 0.002, 0.001])
CopperChromium      = MixtureMaterial(components=[Copper, Chromium], state='solid',
                                      mass_fractions=[0.991, 0.009], density=8.8)
CastIron            = MixtureMaterial(components=[Iron, Carbon, Silicon], state='solid',
                                      mass_fractions=[0.94, 0.04, 0.02], density=7.2)
CuCrZr              = MixtureMaterial(components=[Copper, Chromium, Iron, Silicon, Zirconium],
                                      state='solid', mass_fractions=[0.98805, 0.0085, 0.0008, 0.001, 0.00165],  density=8.9)
TiGr2               = MixtureMaterial(components=[Titanium, Iron, Carbon, Oxygen, Nitrogen, Hydrogen], density=4.51,
                                      state='solid', mass_fractions=[0.99325, 0.003, 0.0008, 0.0025, 0.0003, 0.00015])


# WATER Water   1.0
# POLYSTYR Polystyrene   1.06
# PLASCINT Plastic scintillator   1.032
# PMMA Polymethyl methacrylate, Plexiglas, Lucite, Perspex   1.19
# BONECOMP Compact bone   1.85
# BONECORT Cortical bone   1.85
# MUSCLESK Skeletal muscle   1.04
# MUSCLEST Striated muscle   1.04
# ADTISSUE Adipose tissue   0.92
# KAPTON Kapton polyimide film   1.42
# POLYETHY Polyethylene   0.94
# AIR Dry air at NTP conditions   0.00120479


for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, MixtureMaterial):
        obj.name = name

BoronNitride5000.fluka_name = 'BN5000'
FeNiCo.fluka_name = 'ASTMF-15'
GlidCop15.fluka_name = 'GLIDCP15'
CuNiFeMn.fluka_name = 'CuNiFeMn'
Inermet180.fluka_name = 'INERM180'
StainLessSteel304L.fluka_name = 'SS304L'
StainLessSteel316L.fluka_name = 'SS316L'
StainLessSteel316LN.fluka_name = 'SS316LN'
TitaniumNitride.fluka_name = 'TITNAT'
Antico.fluka_name = 'ANTICO'
CopperChromium.fluka_name = 'CuCr'
CastIron.fluka_name = 'CASTIRON'
CuCrZr.fluka_name = 'CuCrZr'
TiGr2.fluka_name = 'TiGr2'
