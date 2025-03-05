# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import sys
import os
from subprocess import run, PIPE

from .paths import flukafile_resolve, fedb, linebuilder
try:
    from xaux import singleton, FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import singleton, FsPath


@singleton
class FlukaEnvironment:
    def __init__(self):
         self._gfortran_installed = False

    def test_gfortran(self, verbose=False):
        try:
            cmd = run(["gfortran", "-dumpversion"], stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            self._gfortran_installed = False
            raise RuntimeError("Could not find gfortran installation! Need gfortran 9 or higher.")
        if cmd.returncode == 0:
            version = cmd.stdout.decode('UTF-8').strip().split('\n')[0]
            if int(version.split('.')[0]) < 9:
                self._gfortran_installed = False
                raise RuntimeError(f"Need gfortran 9 or higher, but found gfortran {version}!")
            self._gfortran_installed = True
            if verbose:
                cmd2 = run(["which", "gfortran"], stdout=PIPE, stderr=PIPE)
                if cmd2.returncode == 0:
                    file = cmd2.stdout.decode('UTF-8').strip().split('\n')[0]
                    print(f"Found gfortran {version} in {file}", flush=True)
                else:
                    stderr = cmd2.stderr.decode('UTF-8').strip().split('\n')
                    raise RuntimeError(f"Could not run `which gfortran`!\nError given is:\n{stderr}")
        else:
            stderr = cmd.stderr.decode('UTF-8').strip().split('\n')
            self._gfortran_installed = False
            raise RuntimeError(f"Could not run gfortran! Verify its installation.\nError given is:\n{stderr}")


    def brute_force_fedb_environment(self):
        # Add the paths to the environment
        if flukafile_resolve(fedb):
            os.environ['FEDB_PATH'] = fedb.as_posix()
        else:
            raise ValueError(f"Could not find fedb folder {fedb}!")
        if flukafile_resolve(linebuilder):
            os.environ['LB_PATH'] = linebuilder.as_posix()
        else:
            raise ValueError(f"Could not find linebuilder folder {linebuilder}!")
        # Brute-force the system paths
        sys.path.append(fedb.as_posix())
        sys.path.append((fedb / "tools").as_posix())
        sys.path.append((fedb / "tools" / "materials" / "cables").as_posix())
        sys.path.append((linebuilder / "src").as_posix())
        sys.path.append((linebuilder / "lib").as_posix())
        # Force-resolve all dependent files
        flukafile_resolve(fedb / "structure.py")
        flukafile_resolve(fedb / "tools")
        for ff in (fedb / "tools").glob("*.py"):
            flukafile_resolve(ff)
        flukafile_resolve(fedb / "tools" / "materials" / "cables")
        for ff in (fedb / "tools" / "materials" / "cables").glob("*.py"):
            flukafile_resolve(ff)
        flukafile_resolve(linebuilder / "lib")
        for ff in (linebuilder / "lib").glob("*.py"):
            flukafile_resolve(ff)
        flukafile_resolve(linebuilder / "additionals" / "generic_frame.fluka")
        flukafile_resolve(linebuilder / "src")
        for ff in (linebuilder / "src").glob("*.py"):
            flukafile_resolve(linebuilder / "src" / ff)
