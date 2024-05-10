# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import warnings

from xsuite.kernel_definitions import kernel_definitions
from xsuite.prebuild_kernels import regenerate_kernels


def test_init():
    regenerate_kernels()
