# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from ..xoconstants import Constants, constant, group

class XcollInteractions(Constants):
    _category_ = "interaction" # auto-plural -> "interactions"
    _reverse_  = "unique"      # builds interaction_names
    _c_prefix_ = "XC"

    UNITIALISED                 = constant(0,   shortcut="NAN",  type="point", info="Internal; do not use.")

    ENTER_JAW_L                 = constant(-1,  shortcut="JI",   type="point", info="Impact on left jaw. If not on the front, location is given by 'ds'.")
    ENTER_JAW_R                 = constant(-2,  shortcut="JI",   type="point", info="Impact on right jaw. If not on the front, location is given by 'ds'.")
    EXIT_JAW                    = constant(-3,  shortcut="JO",   type="point", info="Exit location.")
    ENTER_JAW                   = constant(-4,  shortcut="JI",   type="point", info="Deprecated. Impact on jaw.")
    ABSORBED                    = constant(1,   shortcut="A",    type="point", info="Particle died in material.")
    MULTIPLE_COULOMB_SCATTERING = constant(13,  shortcut="MCS",  type="continuous", info="MCS transport.")
    PN_ELASTIC                  = constant(14,  shortcut="PN",   type="point", info="Hadronic elastic interaction on a nucleon.")
    PP_ELASTIC                  = constant(15,  shortcut="PP",   type="point", info="Hadronic elastic interaction on a proton.")
    SINGLE_DIFFRACTIVE          = constant(16,  shortcut="SD",   type="point", info="Hadronic single diffractive interaction.")
    COULOMB                     = constant(17,  shortcut="C",    type="point", info="Electromagnetic interaction.")
    CHANNELLING                 = constant(100, shortcut="CH",   type="continuous", info="Channeling in crystal.")
    DECHANNELLING               = constant(101, shortcut="DCH",  type="point", info="Dechanneling in crystal.")
    VOLUME_REFLECTION_TRANS_CH  = constant(102, shortcut="VRCH", type="point", info="Transition region in crystal around +-xpcrit.")
    VOLUME_REFLECTION           = constant(103, shortcut="VR",   type="point", info="Volume reflection in crystal.")
    VOLUME_REFLECTION_TRANS_MCS = constant(104, shortcut="VRAM", type="point", info="Transition region in crystal around t_P.")
    MULTIPLE_COULOMB_TRANS_VR   = constant(105, shortcut="AMVR", type="continuous", info="Transition region in crystal around t_P + 2 xpcrit.")
    VOLUME_CAPTURE              = constant(106, shortcut="VC",   type="point", info="Volume capture in crystal.")

    POINT_INTERACTIONS          = group(ENTER_JAW_L, ENTER_JAW_R, EXIT_JAW, ENTER_JAW, ABSORBED, PN_ELASTIC, PP_ELASTIC,
                                        SINGLE_DIFFRACTIVE, COULOMB, DECHANNELLING, VOLUME_REFLECTION_TRANS_CH,
                                        VOLUME_REFLECTION, VOLUME_REFLECTION_TRANS_MCS, VOLUME_CAPTURE,
                                        info="All interactions that happen at a single point.")
    CONTINUOUS_INTERACTIONS     = group(MULTIPLE_COULOMB_SCATTERING, CHANNELLING, MULTIPLE_COULOMB_TRANS_VR,
                                        info="All interactions that happen continuously.")
    CRYSTAL_INTERACTIONS        = group(CHANNELLING, DECHANNELLING, VOLUME_REFLECTION_TRANS_CH, VOLUME_REFLECTION,
                                        VOLUME_REFLECTION_TRANS_MCS, MULTIPLE_COULOMB_TRANS_VR, VOLUME_CAPTURE,
                                        info="All interactions that can happen in a crystal.")
