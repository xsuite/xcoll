# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xcoll import __version__

def test_version():
    assert __version__ == '0.6.2.dev1+geant4'
