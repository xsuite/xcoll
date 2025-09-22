# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xcoll.headers.constants import Constants, constant, group

from .type1 import GR1
from .type4 import TestType4, GR2, THING_OTHER, YAY_NO_META_2, YAY_NO_META_3


class TestOtherType(Constants):
    __category__  = "other_type"
    __reverse__   = "unique"
    __c_prefix__  = "XFREOO"

    VAS1 = constant(2, "Vas 1")
    VAS2 = constant(-3, "Vas 2")
    GRO  = group(GR1, VAS1, TestType4.GR2, TestType4.YAY_NO_META_3, info='Fancyyyy')  # cross-module group
    GROO = group(GR1, VAS1, GR2, THING_OTHER, YAY_NO_META_2, YAY_NO_META_3)  # cross-module group
