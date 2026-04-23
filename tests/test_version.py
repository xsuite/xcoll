# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
from xcoll import __version__


@pytest.mark.xcother
def test_version():
    assert __version__ == '0.10.0rc0'
