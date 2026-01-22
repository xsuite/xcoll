// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

from ..xoconstants import Constants, constant, group

class XcollInteractions(Constants):
    _category_ = "interaction" # auto-plural -> "interactions"
    _reverse_  = "unique"      # builds interaction_names
    _c_prefix_ = "XC"

    UNITIALISED                 = constant(0, "NAN  // Do not use")

    ENTER_JAW_L                 = constant(-1, "JI   // point (no children)    Set ds > 0 if entering later")
    ENTER_JAW_R                 = constant(-2, "JI   // point (no children)    Set ds > 0 if entering later")
    EXIT_JAW                    = constant(-3, "JO   // point (no children)")
    ENTER_JAW                   = constant(-4, "JI   // point (no children)    Set ds > 0 if entering later    still here for compatibility")
    ABSORBED                    = constant(1,  "A    // point (no children)    Don't use 0 (is default for unitialised)")
    MULTIPLE_COULOMB_SCATTERING = constant(13, "MCS  // continuous")
    PN_ELASTIC                  = constant(14, "PN   // point (no children)")
    PP_ELASTIC                  = constant(15, "PP   // point (no children)")
    SINGLE_DIFFRACTIVE          = constant(16, "SD   // point (no children)")
    COULOMB                     = constant(17, "C    // point (no children)")
    CHANNELLING                 = constant(100, "CH   // continuous")
    DECHANNELLING               = constant(101, "DCH  // point (no children)")
    VOLUME_REFLECTION_TRANS_CH  = constant(102, "VRCH // point (no children)    Transition region around +-xpcrit")
    VOLUME_REFLECTION           = constant(103, "VR   // point (no children)")
    VOLUME_REFLECTION_TRANS_MCS = constant(104, "VRAM // point (no children)    Transition region around t_P")
    MULTIPLE_COULOMB_TRANS_VR   = constant(105, "AMVR // continuous             Transition region around t_P + 2 xpcrit")
    VOLUME_CAPTURE              = constant(106, "VC   // point (no children)")

    POINT_INTERACTIONS          = group(ENTER_JAW_L, ENTER_JAW_R, EXIT_JAW, ENTER_JAW, ABSORBED, PN_ELASTIC, PP_ELASTIC,
                                        SINGLE_DIFFRACTIVE, COULOMB, DECHANNELLING, VOLUME_REFLECTION_TRANS_CH,
                                        VOLUME_REFLECTION, VOLUME_REFLECTION_TRANS_MCS, VOLUME_CAPTURE,
                                        info="All interactions that happen at a single point.")
    CONTINUOUS_INTERACTIONS     = group(MULTIPLE_COULOMB_SCATTERING, CHANNELLING, MULTIPLE_COULOMB_TRANS_VR,
                                        info="All interactions that happen continuously.")
    CRYSTAL_INTERACTIONS        = group(CHANNELLING, DECHANNELLING, VOLUME_REFLECTION_TRANS_CH, VOLUME_REFLECTION,
                                        VOLUME_REFLECTION_TRANS_MCS, MULTIPLE_COULOMB_TRANS_VR, VOLUME_CAPTURE,
                                        info="All interactions that can happen in a crystal.")

# TODO: should be metadata as with pdg
shortcuts = {
    int(ll[0]): ll[1].split('//')[0].strip()
    for ll in [(ii['value'], ii['info']) for ii in interactions_meta.values()]
}
