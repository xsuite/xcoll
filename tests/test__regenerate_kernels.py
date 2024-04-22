# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import warnings

from xsuite_kernels.kernel_definitions import kernel_definitions
from xsuite_kernels.prebuild_kernels import regenerate_kernels


def test_init():
    xcoll_kernels = [name for name, ker in kernel_definitions if 'xcoll' in name]
    if len(xcoll_kernels) > 0:
        regenerate_kernels(kernels=xcoll_kernels)
    else:
        warnings.warn('No Xcoll kernels found!')
