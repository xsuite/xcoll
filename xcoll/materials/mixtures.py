# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material
from .database import db
from .atoms import (Carbon, Oxygen, Argon, Hydrogen, Boron, Nitrogen, Silicon,
                    Aluminium, Iron, Nickel, Cobalt, Zinc, Chromium, Titanium,
                    Zirconium, Manganese, Phosphorus, Sulfur, Copper, Molybdenum,
                    Tungsten, Calcium, Niobium, Magnesium, Tantalum)
from .compounds import CarbonDioxide, BoronNitride, Silica, BoronTrioxide, AluminiumOxide, NiZnFerrite


# TODO: add more from Geant4: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html


# Mixture definitions from FLUKA
# ===============================

# By molar fractions (like number of atoms but not a compound)
CopperDiamond          = Material(components=[Copper, Carbon, Boron], molar_fractions=[0.2359, 0.7548, 0.0093], density=5.4, state='solid')
TiZrMo                 = Material(components=[Molybdenum, Carbon, Titanium, Zirconium], molar_fractions=[0.98762, 0.00158, 0.00996, 0.00084], density=10.22, state='solid')
MolybdenumGraphite     = Material(components=[Molybdenum, Carbon, Titanium], molar_fractions=[0.0184, 0.9809, 0.0007], density=2.55, state='solid')
MolybdenumGraphite6400 = Material(components=[Molybdenum, Carbon], molar_fractions=[0.0146, 0.9854], density=2.48, state='solid')

MolybdenumGraphite.info     = 'Molybdenum-Graphite composite used in LHC collimators.'
MolybdenumGraphite6400.info = 'Molybdenum-Graphite composite. Candidate material for LHC collimators, but not used'

# By mass fractions
Air                 = Material(components=[Nitrogen, Oxygen, Argon, CarbonDioxide], state='gas', pressure=1, temperature=288.15,
                               mass_fractions=[78.084, 20.946, 0.934, 0.036], density=0.001225)
BoronNitride5000    = Material(components=[BoronNitride, Silica, Oxygen, BoronTrioxide, Calcium], state='solid',
                               mass_fractions=[0.929, 0.002, 0.04, 0.004, 0.025], density=1.925)
FeNiCo              = Material(components=[Iron, Nickel, Cobalt], state='solid',
                               mass_fractions=[0.54, 0.29, 0.17], density=8.35)
Glidcop15           = Material(components=[Copper, AluminiumOxide], state='solid',
                               mass_fractions=[0.9985, 0.0015], density=8.93)
CuNiFeMn            = Material(components=[Copper, Nickel, Iron, Manganese], state='solid',
                               mass_fractions=[0.875, 0.1, 0.015, 0.01], density=8.87)
Inermet180          = Material(components=[Tungsten, Nickel, Copper], state='solid',
                               mass_fractions=[0.95, 0.035, 0.015], density=18.0)
StainLessSteel304L  = Material(components=[Iron, Chromium, Nickel, Manganese, Silicon, Phosphorus, Sulfur, Carbon], state='solid',
                               mass_fractions=[0.67145, 0.185, 0.1125, 0.02, 0.01, 0.00045, 0.0003, 0.0003], density=8.02)
StainLessSteel316L  = Material(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon, Phosphorus, Sulfur, Carbon], state='solid',
                               mass_fractions=[0.61395, 0.185, 0.14, 0.03, 0.02, 0.01, 0.00045, 0.0003, 0.0003], density=8.0)
StainLessSteel316LN = Material(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon, Titanium, Nitrogen, Niobium, Copper, Cobalt, Phosphorus, Carbon, Sulfur, Tantalum, Boron], state='solid',
                               mass_fractions=[0.65093, 0.17, 0.12, 0.025, 0.02, 0.0075, 0.0015, 0.0014, 0.001, 0.001, 0.0005, 0.00045, 0.0003, 0.0003, 0.0001, 0.00002], density=8.03)
TitaniumNitride     = Material(components=[Titanium, Nitrogen], state='solid',
                               mass_fractions=[0.77, 0.23], density=5.22)
Antico              = Material(components=[Aluminium, Silicon, Iron, Magnesium, Manganese, Chromium, Titanium, Zinc, Copper], state='solid',
                               mass_fractions=[0.9725, 0.006, 0.005, 0.004, 0.004, 0.0035, 0.002, 0.002, 0.001], density=2.7)
CopperChromium      = Material(components=[Copper, Chromium], state='solid',
                               mass_fractions=[0.991, 0.009], density=8.8)
CastIron            = Material(components=[Iron, Carbon, Silicon], state='solid',
                               mass_fractions=[0.94, 0.04, 0.02], density=7.2)
CuCrZr              = Material(components=[Copper, Chromium, Iron, Silicon, Zirconium], state='solid',
                               mass_fractions=[0.98805, 0.0085, 0.0008, 0.001, 0.00165],  density=8.9)
TiGr2               = Material(components=[Titanium, Iron, Carbon, Oxygen, Nitrogen, Hydrogen], state='solid',
                               mass_fractions=[0.99325, 0.003, 0.0008, 0.0025, 0.0003, 0.00015], density=4.51)


# Extra parameters for Everest
# ============================

MolybdenumGraphite.adapt(inplace=True, nuclear_radius=0.25, nuclear_elastic_slope=76.7,
                         cross_section=[0.362, 0.247, 0, 0, 0, 0.0094e-2], hcut=0.02)
CopperDiamond.adapt(inplace=True, nuclear_radius=0.308, nuclear_elastic_slope=115.0,
                    cross_section=[0.572, 0.370, 0, 0, 0, 0.0279e-2], hcut=0.02)
Glidcop15.adapt(inplace=True, nuclear_radius=0.418, nuclear_elastic_slope=208.7,
                cross_section=[1.246, 0.765, 0, 0, 0, 0.1390e-2], hcut=0.02)
Inermet180.adapt(inplace=True, nuclear_radius=0.578, nuclear_elastic_slope=392.1,
                 cross_section=[2.548, 1.473, 0, 0, 0, 0.5740e-2], hcut=0.02)


# Existing FLUKA mixtures
# =======================

# POLYSTYR Polystyrene 1.06
# PLASCINT Plastic scintillator 1.032
# PMMA Polymethyl methacrylate, Plexiglas, Lucite, Perspex 1.19
# BONECOMP Compact bone 1.85
# BONECORT Cortical bone 1.85
# MUSCLESK Skeletal muscle 1.04
# MUSCLEST Striated muscle 1.04
# ADTISSUE Adipose tissue 0.92
# KAPTON Kapton polyimide film 1.42
# POLYETHY Polyethylene 0.94


# Metadata for database
# =====================

del Carbon, Oxygen, Argon, Hydrogen, Boron, Nitrogen, Silicon, Aluminium, Iron, Nickel, Cobalt, Zinc, Chromium, Titanium
del Zirconium, Manganese, Phosphorus, Sulfur, Copper, Molybdenum, Tungsten, Calcium, Niobium, Magnesium, Tantalum
del CarbonDioxide, BoronNitride, Silica, BoronTrioxide, AluminiumOxide, NiZnFerrite
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and obj.name is None:
        obj.name = name

CopperDiamond.fluka_name = 'CUDIAM75'
MolybdenumGraphite.fluka_name = 'MG6403Fc'
MolybdenumGraphite6400.fluka_name = 'MoGRMG64'
TiZrMo.fluka_name = 'TZM'
Air.fluka_name = 'AIR'
BoronNitride5000.fluka_name = 'BN5000'
FeNiCo.fluka_name = 'ASTMF-15'
Glidcop15.fluka_name = 'GLIDCP15'
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

db['MoGR'] = MolybdenumGraphite
db['CuCD'] = CopperDiamond
db['Glid'] = Glidcop15
db['Iner'] = Inermet180


# Clean up namespace
del name, obj
del Material
del db
