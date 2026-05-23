// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CONSTANTS_H
#define XCOLL_EVEREST_CONSTANTS_H

#if !defined( XC_PROTON_MASS )
    #define   XC_PROTON_MASS ( 938.271998 )
#endif

#if !defined( XC_ELECTRON_MASS )
    #define   XC_ELECTRON_MASS ( 0.51099890 )
#endif

// Constant in front bethe-bloch [mev g^-1 cm^2]
#if !defined( XC_BETHE_BLOCH )
    #define   XC_BETHE_BLOCH ( 0.307075 )
#endif

#if !defined( XC_PP_CS_REF )
    #define   XC_PP_CS_REF ( 0.04 )
#endif

// -------------------------
// #### Crystal-specific ###
// -------------------------

// Screening function [m]
#if !defined( XC_SCREENING )
    #define   XC_SCREENING ( 0.194e-10 )
#endif

// Distance between planes (110) [m]
#if !defined( XC_PLANE_DISTANCE )
    #define   XC_PLANE_DISTANCE ( 1.920e-10 )
#endif

// Thermal vibrations amplitude
#if !defined( XC_THERMAL_VIBRATIONS )
    #define   XC_THERMAL_VIBRATIONS ( 0.075e-10 )
#endif 

// Classical radius of electron [m]
#if !defined( XC_CRADE )
    #define   XC_CRADE ( 2.817940285e-15 )
#endif


#endif /* XCOLL_EVEREST_CONSTANTS_H */
