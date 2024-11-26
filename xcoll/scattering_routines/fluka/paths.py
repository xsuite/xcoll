# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path


_fluka_coupling = Path('/eos/project-f/flukafiles/fluka-coupling').resolve()

fluka = (_fluka_coupling / 'fluka4-4.1' / 'bin' / 'rfluka').resolve()
flukaserver = (_fluka_coupling / 'fluka_coupling' / 'fluka' / 'flukaserver').resolve()
fedb = (_fluka_coupling / 'fedb_coupling').resolve()
