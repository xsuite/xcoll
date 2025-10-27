# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material, CrystalMaterial
from .atoms import Beryllium, Aluminium, Silicon, Copper, Germanium, Molybdenum, Tungsten, Lead
from .allotropes import CarbonFibreCarbon
from .mixtures import MolybdenumGraphite, CopperDiamond, Glidcop15, Inermet180


# for mat in xc.materials.db:
#     if mat.name is not None and mat.name.startswith('K2'):
#         old_dct = mat._xobject._to_dict()
#         new_dct = xc.Material.from_dict(xc.materials.db[mat.name[2:]].to_dict())._xobject._to_dict()
#     for key in set(old_dct.keys()) | set(new_dct.keys()):
#         if key in old_dct and key in new_dct:
#             old_val = old_dct[key]
#             new_val = new_dct[key]
#             if old_val != new_val:
#                 print(f"Material {mat.name}: field {key} differs: {old_val} != {new_val}")


# K2 variants of elements
K2Beryllium  = Beryllium.adapt(nuclear_radius=0.22, nuclear_elastic_slope=74.7,
                    cross_section=[0.271, 0.192, 0, 0, 0, 0.0035e-2], hcut=0.02,
                    A=9.01, radiation_length=0.353)
K2Beryllium._num_nucleons_eff = 3.3668216738035692

K2Aluminium  = Aluminium.adapt(nuclear_radius=0.302, nuclear_elastic_slope=120.3,
                    cross_section=[0.643, 0.418, 0, 0, 0, 0.0340e-2], hcut=0.02,
                    A=26.98, radiation_length=0.089)
K2Aluminium._num_nucleons_eff = 4.852801185429141

K2Silicon    = Silicon.adapt(nuclear_radius=0.441, nuclear_elastic_slope=120.14,
                    cross_section=[0.664, 0.430, 0, 0, 0, 0.0390e-2], hcut=0.02,
                    A=28.08, density=2.33, radiation_length=1)
K2Silicon._num_nucleons_eff = 4.917875746143141

K2Copper     = Copper.adapt(nuclear_radius=0.366, nuclear_elastic_slope=217.8,
                    cross_section=[1.253, 0.769, 0, 0, 0, 0.1530e-2], hcut=0.01,
                    A=63.55, radiation_length=0.0143)
K2Copper._num_nucleons_eff = 6.456795558713786

K2Germanium  = Germanium.adapt(nuclear_radius=0.605, nuclear_elastic_slope=226.35,
                    cross_section=[1.388, 0.844, 0, 0, 0, 0.1860e-2], hcut=0.02,
                    radiation_length=1)
K2Germanium._num_nucleons_eff = 6.750726791090729

K2Molybdenum = Molybdenum.adapt(nuclear_radius=0.481, nuclear_elastic_slope=273.9,
                    cross_section=[1.713, 1.023, 0, 0, 0, 0.2650e-2], hcut=0.02,
                    A=95.96, radiation_length=0.0096, density=10.22)
K2Molybdenum._num_nucleons_eff = 7.4075614639458625

K2Tungsten   = Tungsten.adapt(nuclear_radius=0.520, nuclear_elastic_slope=440.3,
                    cross_section=[2.765, 1.591, 0, 0, 0, 0.7680e-2], hcut=0.01,
                    A=183.85, radiation_length=0.0035)
K2Tungsten._num_nucleons_eff = 9.200252118694163

K2Lead       = Lead.adapt(nuclear_radius=0.542, nuclear_elastic_slope=455.3,
                    cross_section=[3.016, 1.724, 0, 0, 0, 0.9070e-2], hcut=0.01,
                    density=11.35, radiation_length=0.0056)
K2Lead._num_nucleons_eff = 9.57417689618882



# K2 variants of allotropes
K2CarbonFibreCarbon = CarbonFibreCarbon.adapt(excitation_energy=78.0,
                    nuclear_radius=0.25, nuclear_elastic_slope=70.0,
                    cross_section=[0.337, 0.232, 0, 0, 0, 0.0076e-2], hcut=0.02,
                    A=12.01, radiation_length=0.2557)
K2CarbonFibreCarbon._num_nucleons_eff = 3.705323974123165

K2Carbon2 = K2CarbonFibreCarbon.adapt(density=4.52, A=12.01, radiation_length=0.094)  # Was in SixTrack, not sure what it is
K2Carbon2._num_nucleons_eff = 3.705323974123165


# K2 variants of compounds
K2MolybdenumGraphite = MolybdenumGraphite.adapt(excitation_energy=87.1,
                    radiation_length=0.1193, nuclear_radius=0.25,
                    nuclear_elastic_slope=76.7, hcut=0.02, density=2.5,
                    cross_section=[0.362, 0.247, 0, 0, 0, 0.0094e-2])
K2MolybdenumGraphite._ZA_mean = 0.4915003695491501
K2MolybdenumGraphite._Z2_eff = 44.222500000000004
K2MolybdenumGraphite._num_nucleons_eff = 3.8554740361466595
K2MolybdenumGraphite._atoms_per_volume = 1.1127384996304508e+29

K2CopperDiamond = CopperDiamond.adapt(excitation_energy=152.9,
                    radiation_length=0.0316, nuclear_radius=0.308,
                    nuclear_elastic_slope=115.0, hcut=0.02,
                    cross_section=[0.572, 0.370, 0, 0, 0, 0.0279e-2])
K2CopperDiamond._ZA_mean = 0.47147385103011097
K2CopperDiamond._Z2_eff = 141.61
K2CopperDiamond._num_nucleons_eff = 4.74615190534514
K2CopperDiamond._atoms_per_volume = 1.2884136332805071e+29

K2Glidcop15 = Glidcop15.adapt(excitation_energy=320.8e-9,
                    radiation_length=0.0144, nuclear_radius=0.418,
                    nuclear_elastic_slope=208.7, hcut=0.02,
                    cross_section=[1.246, 0.765, 0, 0, 0, 0.1390e-2])
K2Glidcop15._ZA_mean = 0.45605700712589076
K2Glidcop15._Z2_eff = 829.44
K2Glidcop15._num_nucleons_eff = 6.4432201272142064
K2Glidcop15._atoms_per_volume = 8.515869673285828e+28

K2Inermet180 = Inermet180.adapt(excitation_energy=682.2e-9,
                    radiation_length=0.00385, nuclear_radius=0.578,
                    nuclear_elastic_slope=392.1, hcut=0.02,
                    cross_section=[2.548, 1.473, 0, 0, 0, 0.5740e-2])
K2Inermet180._ZA_mean = 0.40611877624475107
K2Inermet180._Z2_eff = 4583.29
K2Inermet180._num_nucleons_eff = 8.904790721425886
K2Inermet180._atoms_per_volume = 6.5026114985003005e+28


del Beryllium, Aluminium, Silicon, Copper, Germanium, Molybdenum, Tungsten, Lead
del CarbonFibreCarbon
del MolybdenumGraphite, CopperDiamond, Glidcop15, Inermet180
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and obj.name is None:
        obj.name = name


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
del name, obj
del Material, CrystalMaterial
