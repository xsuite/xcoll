# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..wrapper import BaseWrapper
from .engine import Geant4Engine
from .environment import Geant4Environment
from .reference_masses import geant4_masses_meta


class Geant4Wrapper(BaseWrapper):
    """Wrapper for all Geant4 and BDSIM functions."""

    _engine_cls = Geant4Engine
    _environment_cls = Geant4Environment
    _particle_masses_meta = geant4_masses_meta
