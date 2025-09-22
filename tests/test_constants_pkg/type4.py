# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from xcoll.xoconstants import Constants, constant, group

from .type1 import VAR4


class TestType4(Constants):
    __category__ = "unique_type"   # auto-plural -> "states"
    __reverse__  = "unique"        # builds particle_state_names
    __c_prefix__ = "UU"

    THING_SOME_VAL = constant(0.734, "Type4 thing with some value.")
    THING_OTHER    = constant(1.23, "Another type4 thing.")
    THING_INT      = constant(12., "Integer type4 thing -> has to be float")
    YAY_NO_META    = 9.81
    YAY_NO_META_2  = 23.
    YAY_NO_META_3  = True
    GR2            = group(THING_SOME_VAL, YAY_NO_META, YAY_NO_META_2)

class TestType5(Constants):
    __category__ = "multi_type"
    __reverse__  = "multi"
    __c_prefix__ = "RR"

    VAR7 = constant(4, "Var 7")
    VAR8 = constant(-5, "Var 8")
    VAR9 = constant(10, "Var 9")


class TestType6(Constants):
    __category__ = "multi_type"
    __reverse__  = "multi"
    __c_prefix__ = "RR"

    VAR10 = constant(4, "Var 10")
    VAR11 = constant(88, "Var 11")
    VAR12 = constant(10, "Var 12")
    GR3   = group(VAR4, TestType5.VAR7, VAR12)
