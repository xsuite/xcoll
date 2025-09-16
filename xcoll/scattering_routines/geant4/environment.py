# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from subprocess import run

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from ..environment import BaseEnvironment


class Geant4Environment(BaseEnvironment):
    _read_only_paths = {'bdsim': 0, 'geant4': 0, 'collimasim': 0}

    def __init__(self):
        cmd = run(['which', 'geant4-config'], capture_output=True)
        self._geant4 = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        cmd = run(['which', 'bdsim'], capture_output=True)
        self._bdsim = FsPath(cmd.stdout.decode().strip()) if cmd.returncode == 0 else None
        try:
            from collimasim import XtrackInterface
            import collimasim as cs
        except (ModuleNotFoundError, ImportError):
            self._collimasim = None
        else:
            self._collimasim = FsPath(cs.__path__)
        super().__init__()

    @property
    def geant4(self):
        return self._geant4

    @property
    def bdsim(self):
        return self._bdsim

    @property
    def collimasim(self):
        return self._collimasim

    @property
    def compiled(self):
        return self.geant4_config is not None and self.bdsim is not None and self.collimasim is not None
