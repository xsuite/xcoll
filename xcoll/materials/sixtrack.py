# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material, CrystalMaterial
from .atoms import Beryllium, Aluminium, Silicon, Copper, Germanium, Molybdenum, Tungsten, Lead
from .allotropes import CarbonFibreCarbon
from .mixtures import MolybdenumGraphite, CopperDiamond, Glidcop15, Inermet180
from .database import _manually_add_material_to_db


# K2 variants of elements
K2Beryllium  = Beryllium.adapt(A=9.01, radiation_length=0.353)
K2Beryllium._num_nucleons_eff = 3.3668216738035692

K2Aluminium  = Aluminium.adapt(A=26.98, radiation_length=0.089)
K2Aluminium._num_nucleons_eff = 4.852801185429141

K2Silicon    = Silicon.adapt(A=28.08, density=2.33, radiation_length=1)
K2Silicon._num_nucleons_eff = 4.917875746143141

K2Copper     = Copper.adapt(A=63.55, radiation_length=0.0143)
K2Copper._num_nucleons_eff = 6.456795558713786

K2Germanium  = Germanium.adapt(radiation_length=1)
K2Germanium._num_nucleons_eff = 6.750726791090729

K2Molybdenum = Molybdenum.adapt(A=95.96, radiation_length=0.0096, density=10.22)
K2Molybdenum._num_nucleons_eff = 7.4075614639458625

K2Tungsten   = Tungsten.adapt(A=183.85, radiation_length=0.0035)
K2Tungsten._num_nucleons_eff = 9.200252118694163

K2Lead       = Lead.adapt(density=11.35, radiation_length=0.0056)
K2Lead._num_nucleons_eff = 9.57417689618882



# K2 variants of allotropes
K2CarbonFibreCarbon = CarbonFibreCarbon.adapt(A=12.01, radiation_length=0.2557)
K2CarbonFibreCarbon._num_nucleons_eff = 3.705323974123165

K2Carbon2 = K2CarbonFibreCarbon.adapt(density=4.52, A=12.01, radiation_length=0.094)  # Was in SixTrack, not sure what it is
K2Carbon2._num_nucleons_eff = 3.705323974123165


# K2 variants of compounds
K2MolybdenumGraphite = MolybdenumGraphite.adapt(excitation_energy=87.1, radiation_length=0.1193, density=2.5)
K2MolybdenumGraphite._ZA_mean = 0.4915003695491501
K2MolybdenumGraphite._Z2_eff = 44.222500000000004
K2MolybdenumGraphite._num_nucleons_eff = 3.8554740361466595
K2MolybdenumGraphite._atoms_per_volume = 1.1127384996304508e+29

K2CopperDiamond = CopperDiamond.adapt(excitation_energy=152.9, radiation_length=0.0316)
K2CopperDiamond._ZA_mean = 0.47147385103011097
K2CopperDiamond._Z2_eff = 141.61
K2CopperDiamond._num_nucleons_eff = 4.74615190534514
K2CopperDiamond._atoms_per_volume = 1.2884136332805071e+29

K2Glidcop15 = Glidcop15.adapt(excitation_energy=320.8, radiation_length=0.0144)
K2Glidcop15._ZA_mean = 0.45605700712589076
K2Glidcop15._Z2_eff = 829.44
K2Glidcop15._num_nucleons_eff = 6.4432201272142064
K2Glidcop15._atoms_per_volume = 8.515869673285828e+28

K2Inermet180 = Inermet180.adapt(excitation_energy=682.2, radiation_length=0.00385)
K2Inermet180._ZA_mean = 0.40611877624475107
K2Inermet180._Z2_eff = 4583.29
K2Inermet180._num_nucleons_eff = 8.904790721425886
K2Inermet180._atoms_per_volume = 6.5026114985003005e+28


_manually_add_material_to_db(K2Beryllium,          'K2Beryllium')
_manually_add_material_to_db(K2Aluminium,          'K2Aluminium')
_manually_add_material_to_db(K2Silicon,            'K2Silicon')
_manually_add_material_to_db(K2Copper,             'K2Copper')
_manually_add_material_to_db(K2Germanium,          'K2Germanium')
_manually_add_material_to_db(K2Molybdenum,         'K2Molybdenum')
_manually_add_material_to_db(K2Tungsten,           'K2Tungsten')
_manually_add_material_to_db(K2Lead,               'K2Lead')
_manually_add_material_to_db(K2CarbonFibreCarbon,  'K2CarbonFibreCarbon')
_manually_add_material_to_db(K2Carbon2,            'K2Carbon2')
_manually_add_material_to_db(K2MolybdenumGraphite, 'K2MolybdenumGraphite')
_manually_add_material_to_db(K2CopperDiamond,      'K2CopperDiamond')
_manually_add_material_to_db(K2Glidcop15,          'K2Glidcop15')
_manually_add_material_to_db(K2Inermet180,         'K2Inermet180')


# K2 variants of crystal materials
K2CarbonCrystal    = CrystalMaterial.from_material(K2CarbonFibreCarbon,
                        name='K2CarbonCrystal', crystal_plane_distance=0.63e-7,
                        crystal_potential=21.0, nuclear_collision_length=0,
                        radiation_length=0.188, eta=0.9)
K2CarbonCrystal._num_nucleons_eff = 3.705323974123165

K2SiliconCrystal   = CrystalMaterial.from_material(K2Silicon,
                        crystal_plane_distance=0.96e-7, crystal_potential=21.34,
                        nuclear_collision_length=0.3016, radiation_length=0.0937,
                        eta=0.9)
K2SiliconCrystal._num_nucleons_eff = 4.917875746143141

K2GermaniumCrystal = CrystalMaterial.from_material(K2Germanium,
                        crystal_plane_distance=1.0e-7, crystal_potential=40.0,
                        nuclear_collision_length=0.1632, radiation_length=0.02302,
                        eta=0.9)
K2GermaniumCrystal._num_nucleons_eff = 6.750726791090729

K2TungstenCrystal  = CrystalMaterial.from_material(K2Tungsten,
                        crystal_plane_distance=0.56e-7, crystal_potential=21.0,
                        nuclear_collision_length=0, eta=0.9)
K2TungstenCrystal._num_nucleons_eff = 9.200252118694163


# Clean up namespace
del (Beryllium, Aluminium, Silicon, Copper, Germanium, Molybdenum, Tungsten, Lead,
     CarbonFibreCarbon, MolybdenumGraphite, CopperDiamond, Glidcop15, Inermet180)
del Material, CrystalMaterial
