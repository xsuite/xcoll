# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .general import _pkg_root, __version__, citation

from .beam_elements import BlackAbsorber, BlackCrystal, EverestBlock, EverestCollimator, EverestCrystal, \
                           BlowUp, EmittanceMonitor, collimator_classes, crystal_classes, element_classes
from .geometry import LineSegment, HalfOpenLineSegment, CircularSegment, BezierSegment, Shape2D, Shape2DV, \
                      create_jaw, create_polygon, create_polygon_jaw, create_crystal
from .scattering_routines.everest import materials, Material, CrystalMaterial
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
from .rf_sweep import RFSweep
from .lossmap import LossMap

# Deprecated
from ._manager import CollimatorManager
from .install import install_elements
from .line_tools import assign_optics_to_collimators, open_collimators, send_to_parking, enable_scattering, disable_scattering
def generate_pencil_on_collimator(line, name, *args, **kwargs):
    from warnings import warn
    warn("`xcoll.generate_pencil_on_collimator()` is deprecated and will be removed. Use "
       + "`line[coll].generate_pencil()` instead.", FutureWarning)
    return line[name].generate_pencil(*args, **kwargs)
def generate_delta_from_dispersion(line, at_element, *args, **kwargs):
    from warnings import warn
    warn("`xcoll.generate_delta_from_dispersion()` is deprecated and will be removed. Use "
       + "`line[at_element].generate_delta()` instead.", FutureWarning)
    return line[name].generate_delta(*args, **kwargs)

# print("If you use Xcoll in your simulations, please cite us :-)")
# print(citation)

