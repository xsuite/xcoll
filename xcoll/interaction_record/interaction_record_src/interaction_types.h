// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_INTERACTIONS_H
#define XCOLL_INTERACTIONS_H

#define  XC_UNITIALISED                     0        // NAN  // Do not use

#define  XC_ENTER_JAW_L                    -1        // JI   // point (no children)    Set ds > 0 if entering later
#define  XC_ENTER_JAW_R                    -2        // JI   // point (no children)    Set ds > 0 if entering later
#define  XC_EXIT_JAW                       -3        // JO   // point (no children)
#define  XC_ENTER_JAW                      -4        // JI   // point (no children)    Set ds > 0 if entering later    still here for compatibility
#define  XC_ABSORBED                        1        // A    // point (no children)    Don't use 0 (is default for unitialised)
#define  XC_MULTIPLE_COULOMB_SCATTERING    13        // MCS  // continuous
#define  XC_PN_ELASTIC                     14        // PN   // point (no children)
#define  XC_PP_ELASTIC                     15        // PP   // point (no children)
#define  XC_SINGLE_DIFFRACTIVE             16        // SD   // point (no children)
#define  XC_COULOMB                        17        // C    // point (no children)
#define  XC_CHANNELLING                   100        // CH   // continuous
#define  XC_DECHANNELLING                 101        // DCH  // point (no children)
#define  XC_VOLUME_REFLECTION_TRANS_CH    102        // VRCH // point (no children)    Transition region around +-xpcrit
#define  XC_VOLUME_REFLECTION             103        // VR   // point (no children)
#define  XC_VOLUME_REFLECTION_TRANS_MCS   104        // VRAM // point (no children)    Transition region around t_P
#define  XC_MULTIPLE_COULOMB_TRANS_VR     105        // AMVR // continuous             Transition region around t_P + 2 xpcrit
#define  XC_VOLUME_CAPTURE                106        // VC   // point (no children)

#endif /* XCOLL_INTERACTIONS_H */