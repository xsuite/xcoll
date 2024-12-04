# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import time
from pathlib import Path


_fluka_coupling = Path('/eos/project/f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
# linebuilder = (_fluka_coupling / 'linebuilder').resolve()
# TODO if MR on gitlab accepted, update path
linebuilder = Path("/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()


# Trying to find fluka and flukaserver executables
# Wait for 3 minutes if they are not found to allow EOS to sync
# TODO this should be done with xaux.FsPath
def flukafile_resolve(fluka_file, timeout=180):
        start_time = time.time()
        fluka_file = Path(fluka_file).resolve()
        _alt_fluka_file = Path(fluka_file.as_posix().replace("project/f", "project-f"))

        while time.time() - start_time < timeout:
            if fluka_file.exists():
                return fluka_file
            elif _alt_fluka_file.exists():
                return _alt_fluka_file
            time.sleep(1)

        return None
