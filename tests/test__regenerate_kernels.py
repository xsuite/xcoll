# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xtrack as xt

def test_init():
    from xtrack.prebuilt_kernels.kernel_definitions import kernel_definitions
    xcoll_kernels = [ker for ker in kernel_definitions if 'xcoll' in ker]
    if len(xcoll_kernels) > 0:
        xt.prebuild_kernels.regenerate_kernels(kernels=xcoll_kernels)
