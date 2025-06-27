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
    # The paths to be set. The value is the parent depth that needs to be brute-forced (0 = file itself, None = no brute-force)
    _paths = {}
    _read_only_paths = {}

    def __init__(self):
        super().__init__()

    def geant4_found(self):
        cmd = run(['which', 'geant4-config'])
        if cmd.returncode != 0:
            return False
        else:
            return True

    def bdsim_found(self):
        cmd = run(['which', 'bdsim'])
        if cmd.returncode != 0:
            return False
        else:
            return True

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




# _geant4_environ = {
#     'G4NEUTRONHPDATA': 'G4NDL*',
#     'G4LEDATA': 'G4EMLOW*',
#     'G4LEVELGAMMADATA': 'PhotonEvaporation*',
#     'G4RADIOACTIVEDATA': 'RadioactiveDecay*',
#     'G4NEUTRONXSDATA': 'G4NEUTRONXS*',
#     'G4PIIDATA': 'G4PII*',
#     'G4REALSURFACEDATA': 'RealSurface*',
#     'G4SAIDXSDATA': 'G4SAIDDATA*',
#     'G4ABLADATA': 'G4ABLA*',
#     'G4ENSDFSTATEDATA': 'G4ENSDFSTATE*'
# }
# _bdsim_environ = {
#     'ROOT_INCLUDE_PATH': ['include/bdsim',
#                           'include/bdsim/analysis',
#                           'include/bdsim/parser'],
# }


# def check_env():
#     os_environ = dict(os.environ)
#     for kk, vv in _geant4_environ.items():
#         if kk not in os_environ:
#             raise ValueError(f"Geant4 environment variable {kk} not found! "
#                               "Please source the Geant4 environment script.")
#         if not Path(os_environ[kk]).exists():
#             raise ValueError(f"Geant4 data path {os_environ[kk]} (for environment variable {kk}) "
#                            + f"not found! Please source the Geant4 environment script.")
#         if not Path(os_environ[kk]).name.startswith(vv[:-1]):
#             print(f"Warning: Geant4 data path {os_environ[kk]} (for environment variable {kk}) "
#                 + f"does not match the expected pattern {vv[-1]}!")
#     for kk, vv in _bdsim_environ.items():
#         if kk not in os_environ:
#             raise ValueError(f"BDSIM environment variable {kk} not found! "
#                               "Please source the BDSIM environment script.")
#         for path in os_environ[kk].split(':'):
#             if not Path(path).exists():
#                 raise ValueError(f"BDSIM data path {vv} (for environment variable {kk}) "
#                                + f"not found! Please source the BDSIM environment script.")
#             if not any([Path(path).as_posix().endswith(vvv) for vvv in vv]):
#                 print(f"Warning: BDSIM data path {vv} (for environment variable {kk}) "
#                     + f"does not match the expected pattern {vv}!")


# def geant4_installed():
#     cmd = run(['which', 'geant4-config'])
#     if cmd.returncode != 0:
#         return False
#     else:
#         return True


# def geant4_env(geant4_path=None, bdsim_path=None):
#     if not geant4_installed():
#         geant4_path, version = _get_geant4_path(geant4_path)
#         bdsim_path = _get_bdsim_path(bdsim_path, version)
#     os_environ = {}
#     if geant4_path:
#         os_environ.update({kk: list(geant4_path.glob(vv)) for kk, vv in _geant4_environ.items()})
#     if bdsim_path:
#         os_environ.update({kk: [bdsim_path / vvv for vvv in vv] for kk, vv in _bdsim_environ.items()})
#     path = os.environ.get('PATH', None)
#     path = [path] if path is not None else []
#     if geant4_path:
#         path = [geant4_path / 'bin', *path]
#     if bdsim_path:
#         path = [bdsim_path / 'bin', *path]
#     if path:
#         os_environ['PATH'] = path
#     ld_library_path = os.environ.get('LD_LIBRARY_PATH', None)
#     ld_library_path = [ld_library_path] if ld_library_path is not None else []
#     if bdsim_path:
#         os_environ['LD_LIBRARY_PATH'] = [bdsim_path / 'lib', *ld_library_path]
#     for kk, vv in os_environ.items():
#         if kk.endswith('PATH'):
#             for vvv in vv:
#                 os_environ[kk] = ':'.join([v.as_posix() if isinstance(v, Path) else v
#                                        for v in vv])
#         else:
#             os_environ[kk] = vv[0].as_posix()
#     return os_environ


# def set_geant4_env(geant4_path=None, bdsim_path=None):
#     old_os_environ = dict(os.environ)
#     for kk, vv in geant4_env(geant4_path, bdsim_path).items():
#        os.environ[kk] = vv
#     check_env()
#     if not geant4_installed():
#         raise ValueError("Geant4 sourcing failed!")
#     return old_os_environ


# def unset_geant4_env(old_os_environ):
#     os.environ.clear()
#     os.environ.update(old_os_environ)


# def _get_geant4_path(geant4_path):
#     version = None
#     if geant4_path is None:
#         geant4_path = list((_pkg_root / "scattering_routines" / "geant4" ).glob("v*/geant4*"))
#         if not geant4_path:
#             raise ValueError("Geant4 installation not found! Please provide " \
#                            + "'geant4_path' explicitly.")
#         elif len(geant4_path) > 1:
#             versions = [vv.parent.name for vv in geant4_path]
#             vals = [1000_000*vv.split('.')[0] + 1000*vv.split('.')[1] + vv.split('.')[2]
#                     for vv in versions]
#             idx_max = max(range(len(vals)), key=vals.__getitem__)
#             return geant4_path[idx_max], versions[idx_max]
#         else:
#             return geant4_path[0], None
#     return geant4_path, version


# def _get_bdsim_path(bdsim_path, version):
#     version = version or 'v*'
#     if bdsim_path is None:
#         bdsim_path = list((_pkg_root / "scattering_routines" / "geant4" ).glob(f"{version}/bdsim*"))
#         if not bdsim_path:
#             raise ValueError("BDSIM installation not found! Please provide " \
#                            + "'bdsim_path' explicitly.")
#         elif len(bdsim_path) > 1:
#             raise ValueError("Multiple BDSIM installations found! Please provide " \
#                           + f"'bdsim_path' explicitly.\n{bdsim_path}")
#         return bdsim_path[0]
#     return bdsim_path
