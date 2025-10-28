# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..xoconstants import Constants, constant, group

class XcollParticleStates(Constants):
    __category__ = "particle_state" # auto-plural -> "particle_states"
    __reverse__  = "unique"         # builds particle_state_names
    __c_prefix__ = "XC"

    LOST_WITHOUT_SPEC       = constant(-300, "Lost in Xcoll but no specific cause recorded.")

    LOST_ON_EVEREST_BLOCK   = constant(-330, "Everest: absorbed by bulk material.")
    LOST_ON_EVEREST_COLL    = constant(-331, "Everest: collimator jaw absorption.")
    LOST_ON_EVEREST_CRYSTAL = constant(-332, "Everest: crystal absorption.")

    LOST_ON_FLUKA_BLOCK     = constant(-333, "FLUKA: absorbed by bulk material.")
    LOST_ON_FLUKA_COLL      = constant(-334, "FLUKA: collimator jaw absorption.")
    LOST_ON_FLUKA_CRYSTAL   = constant(-335, "FLUKA: crystal absorption.")

    LOST_ON_GEANT4_BLOCK    = constant(-336, "Geant4: absorbed by bulk material.")
    LOST_ON_GEANT4_COLL     = constant(-337, "Geant4: collimator jaw absorption.")
    LOST_ON_GEANT4_CRYSTAL  = constant(-338, "Geant4: crystal absorption.")

    LOST_ON_BLACK_ABSORBER  = constant(-340, "Lost on black absorber.")
    LOST_ON_BLACK_CRYSTAL   = constant(-341, "Lost on black crystal.")

    MASSLESS_OR_NEUTRAL     = constant(-350, "Massless or neutral particle.")
    ACC_IONISATION_LOSS     = constant(-351, "Not a real particle: Accumulated ionisation loss.")
    VIRTUAL_ENERGY          = constant(-352, "Not a real particle: Virtual energy deposition.")
    EXCITED_ION_STATE       = constant(-353, "An excited state of an ion (not supported by BDSIM or FLUKA).")

    ERR_INVALID_TRACK       = constant(-390, "Invalid track through Xcoll element.")
    ERR_NOT_IMPLEMENTED     = constant(-391, "Not implemented in Xcoll.")
    ERR_INVALID_XOFIELD     = constant(-392, "Invalid xofield in Xcoll element.")
    ERR                     = constant(-399, "Unknown Xcoll error.")

    HIT_ON_FLUKA_BLOCK      = constant(333, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_FLUKA_COLL       = constant(334, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_FLUKA_CRYSTAL    = constant(335, "Temporary variable to register hits. Should not be present in final states.")

    HIT_ON_GEANT4_BLOCK     = constant(336, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4_COLL      = constant(337, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4_CRYSTAL   = constant(338, "Temporary variable to register hits. Should not be present in final states.")

    # groups
    LOST_ON_BLOCK           = group(LOST_ON_EVEREST_BLOCK, LOST_ON_FLUKA_BLOCK, LOST_ON_GEANT4_BLOCK,
                                    info="Lost on any block")
    LOST_ON_COLL            = group(LOST_ON_EVEREST_COLL, LOST_ON_FLUKA_COLL, LOST_ON_GEANT4_COLL,
                                    info="Lost on any collimator")
    LOST_ON_CRYSTAL         = group(LOST_ON_EVEREST_CRYSTAL, LOST_ON_FLUKA_CRYSTAL, LOST_ON_GEANT4_CRYSTAL,
                                    info="Lost on any crystal")
    LOST_ON_ABSORBER        = group(LOST_ON_BLACK_ABSORBER, LOST_ON_BLACK_CRYSTAL,
                                    info="Lost on any absorber")
    LOST_AS_SPECIAL_STATE   = group(MASSLESS_OR_NEUTRAL, ACC_IONISATION_LOSS, VIRTUAL_ENERGY, EXCITED_ION_STATE,
                                    info="Special states or unsupported particles.")
    USE_IN_LOSSMAP          = group(LOST_ON_BLOCK, LOST_ON_COLL, LOST_ON_CRYSTAL, LOST_ON_ABSORBER, LOST_AS_SPECIAL_STATE,
                                    info="All states that should be used in loss maps.")
