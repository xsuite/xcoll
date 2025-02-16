# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
from xaux import FsPath

from ...general import _pkg_root


_fluka_coupling = FsPath('/eos/project/f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
# linebuilder = (_fluka_coupling / 'linebuilder').resolve()
# TODO if MR on gitlab accepted, update path
linebuilder = FsPath("/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
fedb = (_pkg_root / 'scattering_routines' / 'fluka' / 'fedb').resolve()


# Trying to find fluka and flukaserver executables
# Wait for 1 minute if they are not found to allow EOS to sync
def flukafile_resolve(fluka_file, timeout=60):
    start_time = time.time()
    fluka_file = FsPath(fluka_file).resolve()
    while time.time() - start_time < timeout:
        if hasattr(fluka_file, 'getfid'):
            fluka_file.getfid()
        if fluka_file.exists():
            return fluka_file
        time.sleep(1)
    return None
