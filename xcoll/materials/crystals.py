# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import CrystalMaterial
from .atoms import Silicon, Germanium, Tungsten
from .atom_variants import Graphite


CarbonCrystal    = CrystalMaterial.from_material(Graphite, crystal_plane_distance=0.63e-7,
                   crystal_nuclear_length=0.400, crystal_potential=21.0, nuclear_collision_length=0)

SiliconCrystal   = CrystalMaterial.from_material(Silicon, crystal_plane_distance=0.96e-7,
                   crystal_nuclear_length=0.4652, crystal_potential=21.34, nuclear_collision_length=0.3016)

GermaniumCrystal = CrystalMaterial.from_material(Germanium, crystal_plane_distance=1.0e-7,
                   crystal_nuclear_length=0.2686, crystal_potential=40.0, nuclear_collision_length=0.1632)

TungstenCrystal  = CrystalMaterial.from_material(Tungsten, crystal_plane_distance=0.56e-7,
                   crystal_nuclear_length=0.096, crystal_potential=21.0, nuclear_collision_length=0)
