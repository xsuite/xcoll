# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .material import Material, CrystalMaterial, CompoundMaterial, MixtureMaterial, RefMaterial
from .atoms import *
from .atom_variants import *
from .compounds import *
from .mixtures import *
from .crystals import *

# Freeze all material instances (to prevent accidental modification)
for name, obj in list(globals().items()):  # Have to wrap in list to take a snapshot (avoid updating in-place)
    if isinstance(obj, Material) and name != 'obj':
        obj._frozen = True
del name, obj

from .database import db
def show():
    db.show()


_DEFAULT_MATERIAL = Material(A=0, Z=0, density=0, radiation_length=0)
_DEFAULT_CRYSTALMATERIAL = CrystalMaterial(A=0, Z=0, density=0, radiation_length=0)
