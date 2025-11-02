# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xcoll as xc


# General Usage
# =============

# Xcoll ships with a buit-in database of materials (which cannot be modified at runtime).
xc.materials.show(full=True)

# Materials can be accessed from the module by their full name only:
print(xc.materials.Aluminium)

# or from the database using their name or aliases
mat1 = xc.materials.db['Aluminium']
mat2 = xc.materials.db['Aluminum']
mat3 = xc.materials.db['Al']
print(mat1)
print(mat2)
print(mat3)
print(mat1 == mat2)
print(mat1 == mat3)

# When a material is known to FLUKA or Geant4, it has a fluka_name resp geant4_name attribute
print(f"Molybdenum-Graphite in FLUKA is called: {xc.materials.MolybdenumGraphite.fluka_name}")
print(f"Aluminium in Geant4 is called: {mat1.geant4_name}")

# They can be looked up in a sub-database by these names
print(f"'MG6403Fc' is in the materials database: {'MG6403Fc' in xc.materials.db}")
print(f"'MG6403Fc' is in the FLUKA sub-database: {'MG6403Fc' in xc.materials.db.fluka}")
print(f"'MG6403Fc' is in the Geant4 sub-database: {'MG6403Fc' in xc.materials.db.geant4}")
print(f"{xc.materials.db.fluka['MG6403Fc']=}")
print(f"'G4_Mo' is in the materials database: {'G4_Mo' in xc.materials.db}")
print(f"'G4_Mo' is in the FLUKA sub-database: {'G4_Mo' in xc.materials.db.fluka}")
print(f"'G4_Mo' is in the Geant4 sub-database: {'G4_Mo' in xc.materials.db.geant4}")
print(f"{xc.materials.db.geant4['G4_Mo']=}")

# The material in an Everest element can be specified using a Material instance or a name known to the database
coll1 = xc.EverestBlock(length=1, material=xc.materials.Beryllium)
coll2 = xc.EverestBlock(length=1, material='Beryllium')
coll3 = xc.EverestBlock(length=1, material='Be')
coll4 = xc.EverestBlock(length=1, material='G4_Be')
coll5 = xc.EverestBlock(length=1, material='BERYLLIU')
print(coll1.material == coll2.material)
print(coll1.material == coll3.material)
print(coll1.material == coll4.material)
print(coll1.material == coll5.material)


# Defining New Materials
# ======================

# All elements in the periodic table are predefined in the database.
# When defining a new elemental material (i.e. an allotrope), the fields
# 'Z', 'A', and 'density are required:
WhitePhosphorus = xc.Material(Z=15, A=123.895/4, density=1.823, name='WhitePhosphorus',
                              state='solid', info="P4, but defined as element instead of compound.")

# It is not possible to redefine existing elements in the database, but it
# is possible to adapt them with the `adapt` method.
Ozone = xc.materials.Oxygen.adapt(density=2.144e-3, name='Ozone',
                                  info="O3, but defined as element instead of compound.")
# When adapting a material, unspecified fields are taken from the original material,
# like the state, and temperature and pressure at which the density applies:
print(xc.materials.Oxygen)
print(Ozone)
print(xc.materials.Oxygen.to_dict())
print(Ozone.to_dict())
# But names and info are not copied:
print(f"Original name: {xc.materials.Oxygen.name}, adapted name: {Ozone.name}")
print(f"Original short name: {xc.materials.Oxygen.short_name}, adapted short name: {Ozone.short_name}")
print(f"Original FLUKA name: {xc.materials.Oxygen.fluka_name}, adapted FLUKA name: {Ozone.fluka_name}")
print(f"Original Geant4 name: {xc.materials.Oxygen.geant4_name}, adapted Geant4 name: {Ozone.geant4_name}")
print(f"Original info: `{xc.materials.Oxygen.info}`")
print(f"Adapted info: `{Ozone.info}`")


# Adding Materials to the Database
# ================================

# When providing a name, the material is automatically added to the database:
print(f"WhitePhosphorus in database: {'WhitePhosphorus' in xc.materials.db}")
print(f"Ozone in database: {'Ozone' in xc.materials.db}")
print(f"Ozone in FLUKA sub-database: {'Ozone' in xc.materials.db.fluka}")
Ozone.fluka_name = 'OZONE'
print(f"Ozone in FLUKA sub-database: {'OZONE' in xc.materials.db.fluka}")

# And hence, materials without a name are not in the database:
HeavyGas = xc.Material(Z=10, A=20.1797, density=1.784e-3, state='gas')
print(f"HeavyGas in database: {'HeavyGas' in xc.materials.db}")


# Defining New Compounds and Mixtures
# ===================================

# Compounds are defined by specifying their components as a chemical formula,
# and hence the fields `components` and `n_atoms`:
Ethanol = xc.Material(components=['C', 'H', 'O'], n_atoms=[2, 6, 1], density=0.78945, name='Ethanol',
                      state='liquid', temperature=293.15)
print(Ethanol)
print(Ethanol.composition)
# Any doubled components are automatically combined:
EthanolBis = xc.Material(components=['C', 'H', 'C', 'H', 'O', 'H'], n_atoms=[1, 3, 1, 2, 1, 1], density=0.78945,
                  name='EthanolBis', state='liquid', temperature=293.15)
print(f"{Ethanol == EthanolBis=}")


# Mixtures are not a chemical formula, but a physical mixture of different materials.
# They are defined by specifying their components and either their mass fractions,
# volume fractions, or molar fractions:
Concrete = xc.Material(components=['H', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Fe'],
                       mass_fractions=[0.01, 0.001, 0.529107, 0.016, 0.002, 0.033872, 0.337021, 0.013, 0.044, 0.014],
                       name='Concrete', density=2.35, state='solid')
print(Concrete)
print(Concrete.composition)
print(Concrete.to_dict())

# It is possible to combine any previously defined material as components of a mixture
# (even mixtures themselves). The code will recursively resolve the elemental composition:
Glucose      = xc.Material(components=['C', 'H', 'O'], n_atoms=[6, 12, 6], density=1.54, name='Glucose')
Anethole     = xc.Material(components=['C', 'H', 'O'], n_atoms=[10, 12, 1], density=1.05, name='Anethole',
                           info="Main flavor component (anise) of absinthe.")
components=[Ethanol, xc.materials.Water, Glucose, Anethole]
volume_fractions=[0.7, 0.25, 0.045, 0.005]
FakeAbsinthe = xc.Material(components=components, volume_fractions=volume_fractions,
                           density = sum([vf * el.density for el, vf in zip(components, volume_fractions)]),
                           name='FakeAbsinthe', state='liquid')
print(FakeAbsinthe)
print(FakeAbsinthe.composition)
print(FakeAbsinthe.to_dict())


# Everest Compatibility
# =====================

# Not all materials are 100% compatible with Everest. Multiple Coulomb scattering and ionisation loss are
# always supported, but nuclear interactions are only supported for materials known to Everest.
# This can be checeked with the `full_everest_supported` attribute:
print(f"CFC full Everest compatibility: {xc.materials.CarbonFibreCarbon.full_everest_supported}")
print(f"Ethanol full Everest compatibility: {Ethanol.full_everest_supported}")
print(f"FakeAbsinthe full Everest compatibility: {FakeAbsinthe.full_everest_supported}")
