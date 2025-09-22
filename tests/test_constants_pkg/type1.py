# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xcoll.xoconstants import Constants, constant, group


class TestType1(Constants):
    __category__ = "type"
    __reverse__  = None
    __c_prefix__ = "XF"

    ENABLED = constant(True, 'Important info', c_name='OLA_ENABLED')
    BIG     = constant(2**64 - 1, "Big int")
    MASS    = constant(0.938, "Mass in GeV/c^2")
    COOL    = group(BIG, MASS, info="Isn't this cool?")


class TestType2(Constants):
    __category__ = "unique_type"
    __reverse__  = "unique"
    __c_prefix__ = "UU"

    VAR1 = constant(2, "Var 1")
    VAR2 = constant(-3, "Var 2")
    VAR3 = constant(9, "Var 3")
    GR1  = group(VAR1, VAR2, TestType1.COOL)


class TestType3(Constants):
    __category__ = "multi_type"
    __reverse__  = "multi"
    __c_prefix__ = "XF"

    VAR4 = constant(4, "Var 4")
    VAR5 = constant(-5, "Var 5")
    VAR6 = constant(4, "Var 6")
