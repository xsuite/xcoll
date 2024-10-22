# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path


_fluka_coupling = Path('/eos/project-f/flukafiles/fluka-coupling').resolve()

# TODO if MR on gitlab accepted, update path
# fluka_builder = Path("/eos/project-c/collimation-team/software/fluka_coupling_tmp_patch_xsuite/tools/preprocess").resolve()
#fluka_builder = Path("/afs/cern.ch/user/a/adonadon/conf/dev/fluka_coupling/tools/preprocess").resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
#linebuilder = (_fluka_coupling / 'linebuilder').resolve()
linebuilder = Path("/afs/cern.ch/user/a/adonadon/conf/dev/linebuilder_ads").resolve()
fluka_builder = (linebuilder / "src").resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()
