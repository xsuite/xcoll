# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import xtrack as xt

from ..wrapper import BaseWrapper
from .engine import Geant4Engine
from .environment import Geant4Environment


class Geant4Wrapper(BaseWrapper):
    """Wrapper for all Geant4 and BDSIM functions."""

    _engine_cls = Geant4Engine
    _environment_cls = Geant4Environment
