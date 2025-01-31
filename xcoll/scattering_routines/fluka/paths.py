# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import time
from xaux import FsPath


_fluka_coupling = FsPath('/eos/project/f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
# linebuilder = (_fluka_coupling / 'linebuilder').resolve()
# TODO if MR on gitlab accepted, update path
linebuilder = FsPath("/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()


# Trying to find fluka and flukaserver executables
# Wait for 3 minutes if they are not found to allow EOS to sync
def flukafile_resolve(fluka_file, timeout=180):
        start_time = time.time()
        fluka_file = FsPath(fluka_file).resolve()

        while time.time() - start_time < timeout:
            if hasattr(fluka_file, 'getfid'):
                fluka_file.getfid()
            if fluka_file.exists():
                return fluka_file
            time.sleep(1)

        return None
