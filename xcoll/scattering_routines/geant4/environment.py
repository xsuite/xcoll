# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import os
from pathlib import Path

from ...general import _pkg_root

def set_geant4_env(geant4_path=None, bdsim_path=None):
    if geant4_path is None:
        geant4_path = list((_pkg_root / "scattering_routines" / "geant4").glob("v1*"))
        if not geant4_path:
            raise ValueError("No Geant4 installation found! Please provide 'geant4_path' " \
                           + "and 'bdsim_path' explicitly.")
        elif len(geant4_path) > 1:
            raise ValueError("Multiple Geant4 installations found! Please provide 'geant4_path' " \
                           + "and 'bdsim_path' explicitly.")
        print(geant4_path)
        print(geant4_path[0])
        geant4_path = list(geant4_path[0].glob("geant4*"))[0]
        bdsim_path = list(geant4_path.parent.glob("bdsim*"))[0]
    old_os_environ = dict(os.environ)
    ld_library_path = old_os_environ.get('LD_LIBRARY_PATH')
    ld_library_path = f":LD_LIBRARY_PATH" if ld_library_path else ''
    geant4_env = {
            'G4NEUTRONHPDATA': list(geant4_path.glob('G4NDL*')),
            'G4LEDATA': list(geant4_path.glob('G4EMLOW*')),
            'G4LEVELGAMMADATA': list(geant4_path.glob('PhotonEvaporation*')),
            'G4RADIOACTIVEDATA': list(geant4_path.glob('RadioactiveDecay*')),
            'G4NEUTRONXSDATA': list(geant4_path.glob('G4NEUTRONXS*')),
            'G4PIIDATA': list(geant4_path.glob('G4PII*')),
            'G4REALSURFACEDATA': list(geant4_path.glob('RealSurface*')),
            'G4SAIDXSDATA': list(geant4_path.glob('G4SAIDDATA*')),
            'G4ABLADATA': list(geant4_path.glob('G4ABLA*')),
            'G4ENSDFSTATEDATA': list(geant4_path.glob('G4ENSDFSTATE*')),
            'LD_LIBRARY_PATH': bdsim_path / 'lib',
            'ROOT_INCLUDE_PATH': [bdsim_path / 'include' / 'bdsim',
                                  bdsim_path / 'include' / 'bdsim' / 'analysis',
                                  bdsim_path / 'include' / 'bdsim' / 'parser',]
    }
    for kk, vv in geant4_env.items():
        if kk == 'LD_LIBRARY_PATH':
            if not vv.exists():
                raise ValueError(f"BDSIM path {vv} not found!")
            os.environ[kk] = vv.as_posix() + ld_library_path
        elif kk == 'ROOT_INCLUDE_PATH':
            for vvv in vv:
                if not vvv.exists():
                    raise ValueError(f"BDSIM path {vvv} not found!")
            os.environ[kk] = ':'.join([v.as_posix() for v in vv])
        else:
            if not vv or not vv[0].exists():
                raise ValueError(f"Geant4 environment variable {kk} not found!")
            os.environ[kk] = vv[0].as_posix()
    return old_os_environ


def unset_geant4_env(old_os_environ):
    os.environ.clear()
    os.environ.update(old_os_environ)
