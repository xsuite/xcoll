# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import os
import json
from pathlib import Path
from subprocess import run

from ..environment import BaseEnvironment
from ...general import _pkg_root


class Geant4Environment(BaseEnvironment):

    def __init__(self):
        super().__init__()

    def geant4_found(self):
        cmd = run(['which', 'geant4-config'], capture_output=True)
        return cmd.returncode == 0

    def bdsim_found(self):
        cmd = run(['which', 'bdsim'], capture_output=True)
        return cmd.returncode == 0

    def collimasim_found(self):
        try:
            from collimasim import XtrackInterface
        except (ModuleNotFoundError, ImportError):
            return False
        else:
            return True

    @property
    def compiled(self):
        return self.geant4_found() and self.bdsim_found() and self.collimasim_found()
