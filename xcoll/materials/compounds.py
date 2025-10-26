# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material
from .atoms import *
from .database import db


# TODO: add more from Geant4: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html

# Air = MixtureMaterial(name='Air', components=[Nitrogen, Oxygen, Argon, CarbonDioxide], state='gas',
#                       mass_fractions=[78.084, 20.946, 0.934, 0.036], density=0.001225, pressure=101.325, temperature=15.0)

# INERM180 : matdef, density=18.0, components=["W", "Ni", "Cu"], componentsFractions={0.95, 0.035, 0.015};
# MG6403Fc : matdef, density=2.55, components=["Mo", "C", "Ti"], componentsFractions={0.130012156609076, 0.8675206104800844, 0.002467232910839576};
# AC150GPH : matdef, Z=6.0, A=12.011, density=1.67;
# GPH      : matdef, Z=6.0, A=12.011, density=1.83;
# CUDIAM75 : matdef, density=5.4, components=["Cu", "C", "B"], componentsFractions={0.6205462018968333, 0.3752917530367507, 0.004162045066415989};



# General compounds
# =================

CarbonDioxide = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.951e-3, state='gas', temperature=273.15, pressure=1)
DryIce        = Material(components=[Carbon, Oxygen],   n_atoms=[1,2], density=1.564, state='solid')
Water         = Material(components=[Hydrogen, Oxygen], n_atoms=[2,1], density=1.0, state='liquid')


# Compound definitions from FLUKA
# ===============================

# By number of atoms
BoronNitride  = Material(components=[Boron, Nitrogen],              n_atoms=[1,1], density=2.25, state='solid')
Silica        = Material(components=[Silicon, Oxygen],              n_atoms=[1,2], density=1.92, state='solid')
BoronTrioxide = Material(components=[Boron, Oxygen],                n_atoms=[2,3], density=2.46, state='solid')
AluminumOxide = Material(components=[Aluminium, Oxygen],            n_atoms=[2,3], density=3.98, state='solid')
NiZnFerrite   = Material(components=[Nickel, Zinc, Iron, Oxygen],   n_atoms=[0.45,0.55,2,4], density=5.2, state='solid')
CopperDiamond = Material(components=[Copper, Carbon, Boron],        n_atoms=[0.2359, 0.7548, 0.0093], density=5.4, state='solid')
TiZrMo        = Material(components=[Molybdenum, Carbon, Titanium, Zirconium], density=10.22, state='solid',
                                    n_atoms=[0.98762, 0.00158, 0.00996, 0.00084])
MolybdenumGraphite6400 = Material(components=[Molybdenum, Carbon],           n_atoms=[0.0146, 0.9854], density=2.48, state='solid')
MolybdenumGraphite     = Material(components=[Molybdenum, Carbon, Titanium], n_atoms=[0.0184, 0.9809, 0.0007], density=2.55, state='solid')

# By mass fractions
BoronNitride5000    = Material(components=[BoronNitride, Silica, Oxygen, BoronTrioxide, Calcium],
                                       state='solid', mass_fractions=[0.929, 0.002, 0.04, 0.004, 0.025], density=1.925)
FeNiCo              = Material(components=[Iron, Nickel, Cobalt], state='solid',
                                       mass_fractions=[0.54, 0.29, 0.17], density=8.35)
GlidCop15           = Material(components=[Copper, AluminumOxide], state='solid',
                                       mass_fractions=[0.9985, 0.0015], density=8.93)
CuNiFeMn            = Material(components=[Copper, Nickel, Iron, Manganese], state='solid',
                                       mass_fractions=[0.875, 0.1, 0.015, 0.01], density=8.87)
Inermet180          = Material(components=[Tungsten, Nickel, Copper], state='solid',
                                       mass_fractions=[0.95, 0.035, 0.015], density=18.0)
StainLessSteel304L  = Material(components=[Iron, Chromium, Nickel, Manganese, Silicon, Phosphorus,
                                        Sulfur, Carbon], state='solid', mass_fractions=[0.67145, 0.185, 0.1125, 0.02, 0.01, 0.00045,
                                       0.0003, 0.0003], density=8.02)
StainLessSteel316L  = Material(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon,
                                        Phosphorus, Sulfur, Carbon], state='solid', mass_fractions=[0.61395, 0.185, 0.14, 0.03, 0.02,
                                       0.01, 0.00045, 0.0003, 0.0003], density=8.0)
StainLessSteel316LN = Material(components=[Iron, Chromium, Nickel, Molybdenum, Manganese, Silicon,
                                        Titanium, Nitrogen, Niobium, Copper, Cobalt, Phosphorus, Carbon, Sulfur, Tantalum, Boron],
                                      state='solid', mass_fractions=[0.65093, 0.17, 0.12, 0.025, 0.02, 0.0075, 0.0015, 0.0014,
                                       0.001, 0.001, 0.0005, 0.00045, 0.0003, 0.0003, 0.0001, 0.00002], density=8.03)
TitaniumNitride     = Material(components=[Titanium, Nitrogen], state='solid',
                                       mass_fractions=[0.77, 0.23], density=5.22)
Antico              = Material(components=[Aluminium, Silicon, Iron, Magnesium, Manganese, Chromium,
                                        Titanium, Zinc, Copper], state='solid', density=2.7,
                                      mass_fractions=[0.9725, 0.006, 0.005, 0.004, 0.004, 0.0035, 0.002, 0.002, 0.001])
CopperChromium      = Material(components=[Copper, Chromium], state='solid',
                                       mass_fractions=[0.991, 0.009], density=8.8)
CastIron            = Material(components=[Iron, Carbon, Silicon], state='solid',
                                       mass_fractions=[0.94, 0.04, 0.02], density=7.2)
CuCrZr              = Material(components=[Copper, Chromium, Iron, Silicon, Zirconium],
                                       state='solid', mass_fractions=[0.98805, 0.0085, 0.0008, 0.001, 0.00165],  density=8.9)
TiGr2               = Material(components=[Titanium, Iron, Carbon, Oxygen, Nitrogen, Hydrogen], density=4.51,
                                       state='solid', mass_fractions=[0.99325, 0.003, 0.0008, 0.0025, 0.0003, 0.00015])


# Extra parameters for Everest
print(f"{MolybdenumGraphite6400.Z=}  {MolybdenumGraphite6400.A=}  {MolybdenumGraphite6400.excitation_energy=}  {MolybdenumGraphite6400.radiation_length=}")
print(f"{MolybdenumGraphite.Z=}  {MolybdenumGraphite.A=}  {MolybdenumGraphite.excitation_energy=}  {MolybdenumGraphite.radiation_length=}")
print(f"{CopperDiamond.Z=}  {CopperDiamond.A=}  {CopperDiamond.excitation_energy=}  {CopperDiamond.radiation_length=}")
MolybdenumGraphite.adapt(inplace=True, excitation_energy=87.1, radiation_length=0.1193, nuclear_radius=0.25,
                         nuclear_elastic_slope=76.7, cross_section=[0.362, 0.247, 0, 0, 0, 0.0094e-2], hcut=0.02)
CopperDiamond.adapt(inplace=True, excitation_energy=152.9, radiation_length=0.0316, nuclear_radius=0.308,
                    nuclear_elastic_slope=115.0, cross_section=[0.572, 0.370, 0, 0, 0, 0.0279e-2], hcut=0.02)
print(f"{GlidCop15.Z=}  {GlidCop15.A=}  {GlidCop15.excitation_energy=}  {GlidCop15.radiation_length=}")
print(f"{Inermet180.Z=}  {Inermet180.A=}  {Inermet180.excitation_energy=}  {Inermet180.radiation_length=}")
GlidCop15.adapt(inplace=True, excitation_energy=320.8e-9, radiation_length=0.0144, nuclear_radius=0.418,
                nuclear_elastic_slope=208.7, cross_section=[1.246, 0.765, 0, 0, 0, 0.1390e-2], hcut=0.02)
Inermet180.adapt(inplace=True, excitation_energy=682.2e-9, radiation_length=0.00385, nuclear_radius=0.578,
                 nuclear_elastic_slope=392.1, cross_section=[2.548, 1.473, 0, 0, 0, 0.5740e-2], hcut=0.02)


# Metadata for database
# =====================

for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material):
        if obj.name is None:
            obj.name = name

db['Sand'] = Silica

Water.fluka_name = 'WATER'
BoronNitride.fluka_name = 'hBNpure'
Silica.fluka_name = 'SiO2'
BoronTrioxide.fluka_name = 'B2O3'
AluminumOxide.fluka_name = 'Al2O3'
NiZnFerrite.fluka_name = 'NIZNFER2'
CopperDiamond.fluka_name = 'CUDIAM75'
MolybdenumGraphite.fluka_name = 'MG6403Fc'
MolybdenumGraphite6400.fluka_name = 'MoGRMG64'
TiZrMo.fluka_name = 'TZM'

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
