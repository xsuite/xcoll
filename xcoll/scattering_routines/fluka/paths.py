# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path


_fluka_coupling = Path('/eos/project-f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-3.4' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
linebuilder = (_fluka_coupling / 'linebuilder').resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()
generic_frame = (linebuilder / 'additionals' / 'generic_frame.fluka').resolve()

# sys.path.append(linebuilder / 'src') 
# sys.path.append(linebuilder / 'lib')
# sys.path.append(linebuilder / 'additionals')
# sys.path.append(_fluka_coupling / 'fluka_coupling' / 'tools')
