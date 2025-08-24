// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_STATES_H
#define XCOLL_STATES_H

#define  XC_LOST_ON_EVEREST_BLOCK   -330
#define  XC_LOST_ON_EVEREST_COLL    -331
#define  XC_LOST_ON_EVEREST_CRYSTAL -332

#define  XC_LOST_ON_FLUKA_BLOCK     -333
#define  XC_LOST_ON_FLUKA_COLL      -334
#define  XC_LOST_ON_FLUKA_CRYSTAL   -335

#define  XC_LOST_ON_GEANT4_BLOCK    -336
#define  XC_LOST_ON_GEANT4_COLL     -337
#define  XC_LOST_ON_GEANT4_CRYSTAL  -338

#define  XC_LOST_ON_ABSORBER        -340

#define  XC_MASSLESS_OR_NEUTRAL     -341
#define  XC_ACC_IONISATION_LOSS     -342 // Not a real particle, but accumulation of some ionisation losses that are not accounted for (per collimator)
#define  XC_VIRTUAL_ENERGY          -343 // Not a real particle, but energy that is deposited

#define  XC_ERR_INVALID_TRACK       -390
#define  XC_ERR_NOT_IMPLEMENTED     -391
#define  XC_ERR_INVALID_XOFIELD     -392
#define  XC_ERR                     -399

#endif /* XCOLL_STATES_H */
