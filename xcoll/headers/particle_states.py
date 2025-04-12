# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


XC_LOST_ON_EVEREST_BLOCK   = -330
XC_LOST_ON_EVEREST_COLL    = -331
XC_LOST_ON_EVEREST_CRYSTAL = -332

XC_LOST_ON_FLUKA_BLOCK     = -333
XC_LOST_ON_FLUKA_COLL      = -334
XC_LOST_ON_FLUKA_CRYSTAL   = -335

XC_LOST_ON_GEANT4_BLOCK    = -336
XC_LOST_ON_GEANT4_COLL     = -337
XC_LOST_ON_GEANT4_CRYSTAL  = -338

XC_LOST_ON_ABSORBER        = -340

XC_MASSLESS_OR_NEUTRAL     = -341
XC_ACC_IONISATION_LOSS     = -342  # Not a real particle, but accumulation of some ionisation losses that are not accounted for (per collimator)
XC_VIRTUAL_ENERGY          = -343  # Not a real particle, but energy that is deposited

XC_ERR_INVALID_TRACK       = -390
XC_ERR_NOT_IMPLEMENTED     = -391
XC_ERR_INVALID_XOFIELD     = -392
XC_ERR                     = -399

particle_states_src = f"""
#ifndef XCOLL_STATES_H
#define XCOLL_STATES_H
#define  XC_LOST_ON_EVEREST_BLOCK   {XC_LOST_ON_EVEREST_BLOCK}
#define  XC_LOST_ON_EVEREST_COLL    {XC_LOST_ON_EVEREST_COLL}
#define  XC_LOST_ON_EVEREST_CRYSTAL {XC_LOST_ON_EVEREST_CRYSTAL}
#define  XC_LOST_ON_FLUKA_BLOCK     {XC_LOST_ON_FLUKA_BLOCK}
#define  XC_LOST_ON_FLUKA_COLL      {XC_LOST_ON_FLUKA_COLL}
#define  XC_LOST_ON_FLUKA_CRYSTAL   {XC_LOST_ON_FLUKA_CRYSTAL}
#define  XC_LOST_ON_GEANT4_BLOCK    {XC_LOST_ON_GEANT4_BLOCK}
#define  XC_LOST_ON_GEANT4_COLL     {XC_LOST_ON_GEANT4_COLL}
#define  XC_LOST_ON_GEANT4_CRYSTAL  {XC_LOST_ON_GEANT4_CRYSTAL}
#define  XC_LOST_ON_ABSORBER        {XC_LOST_ON_ABSORBER}
#define  XC_MASSLESS_OR_NEUTRAL     {XC_MASSLESS_OR_NEUTRAL}
#define  XC_ACC_IONISATION_LOSS     {XC_ACC_IONISATION_LOSS}
#define  XC_VIRTUAL_ENERGY          {XC_VIRTUAL_ENERGY}
#define  XC_ERR_INVALID_TRACK       {XC_ERR_INVALID_TRACK}
#define  XC_ERR_NOT_IMPLEMENTED     {XC_ERR_NOT_IMPLEMENTED}
#define  XC_ERR_INVALID_XOFIELD     {XC_ERR_INVALID_XOFIELD}
#define  XC_ERR                     {XC_ERR}
#endif /* XCOLL_STATES_H */
"""