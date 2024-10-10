# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from xsuite.prebuild_kernels import regenerate_kernels


def test_init():
    regenerate_kernels()
