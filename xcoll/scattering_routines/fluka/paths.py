# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path


_fluka_coupling = Path('/eos/project-f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
# linebuilder = (_fluka_coupling / 'linebuilder').resolve()
# TODO if MR on gitlab accepted, update path
linebuilder = Path("/eos/project-c/collimation-team/software/fluka_coupling_tmp_patch_xsuite").resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()
