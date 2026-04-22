# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..xoconstants import Constants, constant, group

class XcollParticleStates(Constants):
    _category_ = "particle_state" # auto-plural -> "particle_states"
    _reverse_  = "unique"         # builds particle_state_names
    _c_prefix_ = "XC"

    LOST_WITHOUT_SPEC       = constant(-300, "Lost in Xcoll but no specific cause recorded.")

    LOST_ON_EVEREST         = constant(-330, "Primary loss in an Everest element (block, collimator, crystal).")
    LOST_ON_FLUKA           = constant(-331, "Primary loss in a FLUKA element (block, collimator, crystal).")
    LOST_ON_GEANT4          = constant(-332, "Primary loss in a Geant4 element (block, collimator, crystal).")
    LOST_ON_ABSORBER        = constant(-333, "Primary loss in a black absorber or black crystal.")

    LOST_ON_EVEREST_SEC     = constant(-334, "Secondary loss in an Everest element (block, collimator, crystal).")
    LOST_ON_FLUKA_SEC       = constant(-335, "Secondary loss in a FLUKA element (block, collimator, crystal).")
    LOST_ON_GEANT4_SEC      = constant(-336, "Secondary loss in a Geant4 element (block, collimator, crystal).")
    LOST_ON_ABSORBER_SEC    = constant(-337, "Secondary loss in a black absorber or black crystal.")

    ACC_IONISATION_LOSS     = constant(-350, "Primary loss: Not a real particle: Accumulated ionisation loss.")
    ACC_IONISATION_LOSS_SEC = constant(-351, "Secondary loss: Not a real particle: Accumulated ionisation loss.")
    VIRTUAL_ENERGY          = constant(-352, "Primary loss: Not a real particle: Virtual energy deposition.")
    VIRTUAL_ENERGY_SEC      = constant(-353, "Secondary loss: Not a real particle: Virtual energy deposition.")
    MASSLESS_OR_NEUTRAL     = constant(-354, "Secondary loss: Massless or neutral particle.")
    EXCITED_ION_STATE       = constant(-355, "Secondary loss: An excited state of an ion (not supported by BDSIM or FLUKA).")

    ERR_INVALID_TRACK       = constant(-390, "Invalid track through Xcoll element.")
    ERR_NOT_IMPLEMENTED     = constant(-391, "Not implemented in Xcoll.")
    ERR_INVALID_XOFIELD     = constant(-392, "Invalid xofield in Xcoll element.")
    ERR                     = constant(-399, "Unknown Xcoll error.")

    SECONDARY_PARTICLE      = constant(2,   "The particle has scattered off an Everest/FLUKA/Geant4 element before.")
    HIT_ON_FLUKA            = constant(310, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_FLUKA_SEC        = constant(311, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4           = constant(312, "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4_SEC       = constant(313, "Temporary variable to register hits. Should not be present in final states.")

    # groups
    LOST_AS_SPECIAL_STATE_PRIM = group(ACC_IONISATION_LOSS, VIRTUAL_ENERGY,
                                    info="Special states or unsupported particles (as primary loss).")
    LOST_AS_SPECIAL_STATE_SEC  = group(ACC_IONISATION_LOSS_SEC, VIRTUAL_ENERGY_SEC, MASSLESS_OR_NEUTRAL, EXCITED_ION_STATE,
                                    info="Special states or unsupported particles (as secondary loss).")
    LOST_AS_SPECIAL_STATE      = group(LOST_AS_SPECIAL_STATE_PRIM, LOST_AS_SPECIAL_STATE_SEC,
                                    info="Special states or unsupported particles.")
    USE_IN_LOSSMAP_PRIM        = group(LOST_ON_EVEREST, LOST_ON_FLUKA, LOST_ON_GEANT4, LOST_ON_ABSORBER,
                                       LOST_AS_SPECIAL_STATE_PRIM,
                                       info="All states that should be used in loss maps (as primary loss).")
    USE_IN_LOSSMAP_SEC         = group(LOST_ON_EVEREST_SEC, LOST_ON_FLUKA_SEC, LOST_ON_GEANT4_SEC, LOST_ON_ABSORBER_SEC,
                                       LOST_AS_SPECIAL_STATE_SEC,
                                       info="All states that should be used in loss maps (as secondary loss).")
    USE_IN_LOSSMAP             = group(USE_IN_LOSSMAP_PRIM, USE_IN_LOSSMAP_SEC, info="All states that should be used in loss maps.")
