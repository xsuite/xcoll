# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import os
from pathlib import Path

def set_geant4_env(geant4_path):
    old_os_environ = dict(os.environ)
    ld_library_path = old_os_environ.get('LD_LIBRARY_PATH')
    ld_library_path = f":{ld_library_path}" if ld_library_path else ''
    geant4_path = Path(geant4_path).as_posix()
    geant4_env = {
            'G4NEUTRONHPDATA': geant4_path + 'G4NDL4.5',
            'G4LEDATA': geant4_path + 'G4EMLOW7.3',
            'G4LEVELGAMMADATA': geant4_path + 'PhotonEvaporation5.2',
            'G4RADIOACTIVEDATA': geant4_path + 'RadioactiveDecay5.2',
            'G4NEUTRONXSDATA': geant4_path + 'G4NEUTRONXS1.4',
            'G4PIIDATA': geant4_path + 'G4PII1.3',
            'G4REALSURFACEDATA': geant4_path + 'RealSurface2.1.1',
            'G4SAIDXSDATA': geant4_path + 'G4SAIDDATA1.1',
            'G4ABLADATA': geant4_path + 'G4ABLA3.1',
            'G4ENSDFSTATEDATA': geant4_path + 'G4ENSDFSTATE2.2',
            'LD_LIBRARY_PATH': geant4_path + 'bdsim/lib:LD_LIBRARY_PATH',
            'ROOT_INCLUDE_PATH': geant4_path + 'bdsim//include/bdsim/:' + geant4_path + 'bdsim/include/bdsim/analysis/:' + geant4_path + 'bdsim/include/bdsim/parser/'
    }
    for kk, val in geant4_env.items():
        os.environ[kk] = val
    return old_os_environ


def unset_geant4_env(old_os_environ):
    os.environ.clear()
    os.environ.update(old_os_environ)
