# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import os
import json
from pathlib import Path
from subprocess import run

from ...general import _pkg_root


def geant4_env(geant4_path=None, bdsim_path=None):
    geant4_path = _get_geant4_path(geant4_path)
    bdsim_path = _get_bdsim_path(bdsim_path)
    os_environ = {
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
        'ROOT_INCLUDE_PATH': [bdsim_path / 'include' / 'bdsim',
                                bdsim_path / 'include' / 'bdsim' / 'analysis',
                                bdsim_path / 'include' / 'bdsim' / 'parser'],
    }
    path = os.environ.get('PATH', None)
    if path is not None:
        os_environ['PATH'] = [bdsim_path / 'bin', geant4_path / 'bin', path]
    else:
        os_environ['PATH'] = [bdsim_path / 'bin', geant4_path / 'bin']
    ld_library_path = os.environ.get('LD_LIBRARY_PATH', None)
    if ld_library_path is not None:
        os_environ['LD_LIBRARY_PATH'] = [bdsim_path / 'lib', ld_library_path]
    else:
        os_environ['LD_LIBRARY_PATH'] = [bdsim_path / 'lib']
    for kk, vv in os_environ.items():
        if kk.endswith('PATH'):
            for vvv in vv:
                if isinstance(vvv, Path) and not vvv.exists():
                    raise ValueError(f"Geant4/BDSIM path {vvv} not found!")
            os_environ[kk] = ':'.join([v.as_posix() if isinstance(v, Path) else v
                                       for v in vv])
        else:
            if not vv or not vv[0].exists():
                raise ValueError(f"Geant4 environment variable {kk} not found!")
            os_environ[kk] = vv[0].as_posix()
    return os_environ


def set_geant4_env(geant4_path=None, bdsim_path=None):
    old_os_environ = dict(os.environ)
    for kk, vv in geant4_env(geant4_path, bdsim_path).items():
       os.environ[kk] = vv
    return old_os_environ


def unset_geant4_env(old_os_environ):
    os.environ.clear()
    os.environ.update(old_os_environ)


# Use this when using a new version of Geant4 or BDSIM, to check that 'set_geant4_env'
# is doing the same as the sourcing scripts. Then adapt 'set_geant4_env' if necessary.
def compare_environs(geant4_path=None, bdsim_path=None):
    # Get the existing environment, updated with the paths as we know them
    old_os_environ = _get_environ()#extra_dict=geant4_env(geant4_path, bdsim_path))
    old_os_environ.update(_get_environ(geant4_env(geant4_path, bdsim_path)))

    # Source the environment scripts manually
    geant4_path = _get_geant4_path(geant4_path)
    bdsim_path = _get_bdsim_path(bdsim_path)
    if geant4_path is None:
        geant4_path = list((_pkg_root / "scattering_routines" / "geant4" / "lib").glob("geant4*"))
        if not geant4_path:
            raise ValueError("Geant4 installation not found! Please provide " \
                           + "'geant4_path' explicitly.")
        elif len(geant4_path) > 1:
            raise ValueError("Multiple Geant4 installations found! Please provide " \
                           + "'geant4_path' explicitly.")
        geant4_path = geant4_path[0]
    if bdsim_path is None:
        bdsim_path = list((_pkg_root / "scattering_routines" / "geant4" / "lib").glob("bdsim*"))
        if not bdsim_path:
            raise ValueError("BDSIM installation not found! Please provide " \
                           + "'bdsim_path' explicitly.")
        elif len(bdsim_path) > 1:
            raise ValueError("Multiple BDSIM installations found! Please provide " \
                           + "'bdsim_path' explicitly.")
        bdsim_path = bdsim_path[0]
    new_os_environ = Path('new_os_environ.json').resolve()
    run(f"source ./geant4.sh; unset LD_LIBRARY_PATH; source {bdsim_path}/bin/bdsim.sh; "\
      + "python -c 'import os, json; fp = open(\"new_os_environ.json\", \"w\"); "\
      + "json.dump(dict(os.environ), fp); fp.close();' ",
      shell=True, cwd=(geant4_path / 'bin').as_posix())
    file = Path(geant4_path / 'bin' / 'new_os_environ.json')
    with file.resolve().open('r') as f:
        new_os_environ = _get_environ(json.load(f))
    file.unlink()

    # Compare the two environments
    old_os_environ = set(old_os_environ.items())
    new_os_environ = set(new_os_environ.items())
    difference = old_os_environ ^ new_os_environ
    return difference#[[diff[0], diff[1].split(':')] for diff in difference]


def _get_geant4_path(geant4_path=None):
    if geant4_path is None:
        geant4_path = list((_pkg_root / "scattering_routines" / "geant4" / "lib").glob("geant4*"))
        if not geant4_path:
            raise ValueError("Geant4 installation not found! Please provide " \
                           + "'geant4_path' explicitly.")
        elif len(geant4_path) > 1:
            raise ValueError("Multiple Geant4 installations found! Please provide " \
                           + "'geant4_path' explicitly.")
        return geant4_path[0]
    return geant4_path


def _get_bdsim_path(bdsim_path=None):
    if bdsim_path is None:
        bdsim_path = list((_pkg_root / "scattering_routines" / "geant4" / "lib").glob("bdsim*"))
        if not bdsim_path:
            raise ValueError("BDSIM installation not found! Please provide " \
                           + "'bdsim_path' explicitly.")
        elif len(bdsim_path) > 1:
            raise ValueError("Multiple BDSIM installations found! Please provide " \
                           + "'bdsim_path' explicitly.")
        return bdsim_path[0]
    return bdsim_path


def _get_environ(os_environ=None, extra_dict={}):
    if os_environ is None:
        os_environ = dict(os.environ)
    os_environ.pop('_', None)
    os_environ.pop('PWD', None)
    os_environ.pop('OLDPWD', None)
    os_environ.pop('DISPLAY', None)
    os_environ.pop('SSH_TTY', None)
    os_environ.pop('SSH_CLIENT', None)
    os_environ.pop('SSH_CONNECTION', None)
    os_environ.pop('XDG_SESSION_ID', None)
    os_environ.update(extra_dict)
    for kk, vv in os_environ.items():
        if kk.endswith('PATH'):
            os_environ[kk] = tuple(Path(vvv).resolve().as_posix() for vvv in vv.split(':'))
        elif kk.startswith('G4'):
            os_environ[kk] = Path(vv).resolve().as_posix()
    return os_environ
