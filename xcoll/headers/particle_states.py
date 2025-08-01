# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #


LOST_ON_EVEREST_BLOCK   = -330
LOST_ON_EVEREST_COLL    = -331
LOST_ON_EVEREST_CRYSTAL = -332

LOST_ON_FLUKA_BLOCK     = -333
LOST_ON_FLUKA_COLL      = -334
LOST_ON_FLUKA_CRYSTAL   = -335

LOST_ON_GEANT4_BLOCK    = -336
LOST_ON_GEANT4_COLL     = -337
LOST_ON_GEANT4_CRYSTAL  = -338

LOST_ON_ABSORBER        = -340

MASSLESS_OR_NEUTRAL     = -341
ACC_IONISATION_LOSS     = -342  # Not a real particle, but accumulation of some ionisation losses that are not accounted for (per collimator)
VIRTUAL_ENERGY          = -343  # Not a real particle, but energy that is deposited

ERR_INVALID_TRACK       = -390
ERR_NOT_IMPLEMENTED     = -391
ERR_INVALID_XOFIELD     = -392
ERR                     = -399

HIT_ON_FLUKA_COLL       = 334

particle_states_src = f"""
#ifndef XCOLL_STATES_H
#define XCOLL_STATES_H
#define  XC_LOST_ON_EVEREST_BLOCK   {LOST_ON_EVEREST_BLOCK}
#define  XC_LOST_ON_EVEREST_COLL    {LOST_ON_EVEREST_COLL}
#define  XC_LOST_ON_EVEREST_CRYSTAL {LOST_ON_EVEREST_CRYSTAL}
#define  XC_LOST_ON_FLUKA_BLOCK     {LOST_ON_FLUKA_BLOCK}
#define  XC_LOST_ON_FLUKA_COLL      {LOST_ON_FLUKA_COLL}
#define  XC_HIT_ON_FLUKA_COLL       {HIT_ON_FLUKA_COLL}
#define  XC_LOST_ON_FLUKA_CRYSTAL   {LOST_ON_FLUKA_CRYSTAL}
#define  XC_LOST_ON_GEANT4_BLOCK    {LOST_ON_GEANT4_BLOCK}
#define  XC_LOST_ON_GEANT4_COLL     {LOST_ON_GEANT4_COLL}
#define  XC_LOST_ON_GEANT4_CRYSTAL  {LOST_ON_GEANT4_CRYSTAL}
#define  XC_LOST_ON_ABSORBER        {LOST_ON_ABSORBER}
#define  XC_MASSLESS_OR_NEUTRAL     {MASSLESS_OR_NEUTRAL}
#define  XC_ACC_IONISATION_LOSS     {ACC_IONISATION_LOSS}
#define  XC_VIRTUAL_ENERGY          {VIRTUAL_ENERGY}
#define  XC_ERR_INVALID_TRACK       {ERR_INVALID_TRACK}
#define  XC_ERR_NOT_IMPLEMENTED     {ERR_NOT_IMPLEMENTED}
#define  XC_ERR_INVALID_XOFIELD     {ERR_INVALID_XOFIELD}
#define  XC_ERR                     {ERR}
#endif /* XCOLL_STATES_H */
"""