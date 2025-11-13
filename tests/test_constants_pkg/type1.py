# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xcoll.xoconstants import Constants, constant, group


class TestType1(Constants):
    _category_ = "type"
    _reverse_  = None
    _c_prefix_ = "XF"

    ENABLED = constant(True, 'Important info', c_name='OLA_ENABLED')
    BIG     = constant(2**64 - 1, "Big int")
    MASS    = constant(0.938, "Mass in GeV/c^2")
    COOL    = group(BIG, MASS, info="Isn't this cool?")


class TestType2(Constants):
    _category_ = "unique_type"
    _reverse_  = "unique"
    _c_prefix_ = "UU"

    VAR1 = constant(2, "Var 1")
    VAR2 = constant(-3, "Var 2")
    VAR3 = constant(9, "Var 3")
    GR1  = group(VAR1, VAR2, TestType1.COOL)


class TestType3(Constants):
    _category_ = "multi_type"
    _reverse_  = "multi"
    _c_prefix_ = "XF"

    VAR4 = constant(4, "Var 4")
    VAR5 = constant(-5, "Var 5")
    VAR6 = constant(4, "Var 6")
