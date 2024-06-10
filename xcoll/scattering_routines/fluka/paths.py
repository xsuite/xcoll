# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path


_fluka_coupling = Path('/eos/project-f/flukafiles/fluka-coupling').resolve()

# _fluka_coupling_ads = Path("/afs/cern.ch/work/a/adonadon/public/fellow/fluka_builder_development/fluka_coupling/").resolve()

_fluka_coupling_ads = Path("/afs/cern.ch/work/a/adonadon/public/fellow/development-xcoll-fluka-coupling/fluka_coupling_ads/fluka_coupling/").resolve()

fluka = (_fluka_coupling / 'fluka4-3.4' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
linebuilder = (_fluka_coupling / 'linebuilder').resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()
generic_frame = (linebuilder / 'additionals' / 'generic_frame.fluka').resolve()
fluka_builder = (_fluka_coupling_ads/ 'tools' / 'preprocess').resolve()
# Path("tools/preprocess/FLUKA_builder_with_main_ads.py")

# sys.path.append(linebuilder / 'src') 
# sys.path.append(linebuilder / 'lib')
# sys.path.append(linebuilder / 'additionals')
# sys.path.append(_fluka_coupling / 'fluka_coupling' / 'tools')
