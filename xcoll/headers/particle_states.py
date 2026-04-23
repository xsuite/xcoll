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

    LOST_ON_MATERIAL        = constant(-330, "Primary loss in an Xcoll material element (block, collimator, crystal).")
    LOST_ON_MATERIAL_SEC    = constant(-331, "Secondary loss in an Xcoll material element (block, collimator, crystal).")

    VIRTUAL_ENERGY          = constant(-350, "Primary loss: Not a real particle: Virtual energy deposition.")
    VIRTUAL_ENERGY_SEC      = constant(-351, "Secondary loss: Not a real particle: Virtual energy deposition.")
    MASSLESS_OR_NEUTRAL     = constant(-352, "Secondary loss: Massless or neutral particle.")
    EXCITED_ION_STATE       = constant(-353, "Secondary loss: An excited state of an ion (not supported by BDSIM or FLUKA).")

    ERR_INVALID_TRACK       = constant(-390, "Invalid track through Xcoll element.")
    ERR_NOT_IMPLEMENTED     = constant(-391, "Not implemented in Xcoll.")
    ERR_INVALID_XOFIELD     = constant(-392, "Invalid xofield in Xcoll element.")
    ERR                     = constant(-399, "Unknown Xcoll error.")

    SECONDARY_PARTICLE      = constant(301,  "The particle has scattered off an Xcoll material element before.")
    HIT_ON_FLUKA            = constant(302,  "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_FLUKA_SEC        = constant(303,  "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4           = constant(304,  "Temporary variable to register hits. Should not be present in final states.")
    HIT_ON_GEANT4_SEC       = constant(305,  "Temporary variable to register hits. Should not be present in final states.")

    # groups
    LOST_AS_SPECIAL_STATE_PRIM = group(VIRTUAL_ENERGY, info="Special states or unsupported particles (as primary loss).")
    LOST_AS_SPECIAL_STATE_SEC  = group(VIRTUAL_ENERGY_SEC, MASSLESS_OR_NEUTRAL, EXCITED_ION_STATE,
                                    info="Special states or unsupported particles (as secondary loss).")
    LOST_AS_SPECIAL_STATE      = group(LOST_AS_SPECIAL_STATE_PRIM, LOST_AS_SPECIAL_STATE_SEC,
                                    info="Special states or unsupported particles.")
    USE_IN_LOSSMAP_PRIM        = group(LOST_ON_MATERIAL, LOST_AS_SPECIAL_STATE_PRIM,
                                       info="All states that should be used in loss maps (as primary loss).")
    USE_IN_LOSSMAP_SEC         = group(LOST_ON_MATERIAL_SEC, LOST_AS_SPECIAL_STATE_SEC,
                                       info="All states that should be used in loss maps (as secondary loss).")
    USE_IN_LOSSMAP             = group(USE_IN_LOSSMAP_PRIM, USE_IN_LOSSMAP_SEC, info="All states that should be used in loss maps.")
