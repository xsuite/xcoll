# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from xcoll import __version__

def test_version():
    assert __version__ == '0.8.0rc0'
