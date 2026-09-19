#ifndef XC_PARTICLE_STATE_H_I0
#define XC_PARTICLE_STATE_H_I0
  #define XC_LOST_WITHOUT_SPEC                -300                          // Lost in Xcoll but no specific cause recorded.
  #define XC_LOST_ON_EVEREST_BLOCK            -330                          // Everest: absorbed by bulk material.
  #define XC_LOST_ON_EVEREST_COLL             -331                          // Everest: collimator jaw absorption.
  #define XC_LOST_ON_EVEREST_CRYSTAL          -332                          // Everest: crystal absorption.
  #define XC_LOST_ON_FLUKA_BLOCK              -333                          // FLUKA: absorbed by bulk material.
  #define XC_LOST_ON_FLUKA_COLL               -334                          // FLUKA: collimator jaw absorption.
  #define XC_LOST_ON_FLUKA_CRYSTAL            -335                          // FLUKA: crystal absorption.
  #define XC_LOST_ON_GEANT4_BLOCK             -336                          // Geant4: absorbed by bulk material.
  #define XC_LOST_ON_GEANT4_COLL              -337                          // Geant4: collimator jaw absorption.
  #define XC_LOST_ON_GEANT4_CRYSTAL           -338                          // Geant4: crystal absorption.
  #define XC_LOST_ON_BLACK_ABSORBER           -340                          // Lost on black absorber.
  #define XC_LOST_ON_BLACK_CRYSTAL            -341                          // Lost on black crystal.
  #define XC_MASSLESS_OR_NEUTRAL              -350                          // Massless or neutral particle.
  #define XC_ACC_IONISATION_LOSS              -351                          // Not a real particle: Accumulated ionisation loss.
  #define XC_VIRTUAL_ENERGY                   -352                          // Not a real particle: Virtual energy deposition.
  #define XC_EXCITED_ION_STATE                -353                          // An excited state of an ion (not supported by BDSIM or FLUKA).
  #define XC_ERR_INVALID_TRACK                -390                          // Invalid track through Xcoll element.
  #define XC_ERR_NOT_IMPLEMENTED              -391                          // Not implemented in Xcoll.
  #define XC_ERR_INVALID_XOFIELD              -392                          // Invalid xofield in Xcoll element.
  #define XC_ERR                              -399                          // Unknown Xcoll error.
  #define XC_HIT_ON_FLUKA_BLOCK               333                           // Temporary variable to register hits. Should not be present in final states.
  #define XC_HIT_ON_FLUKA_COLL                334                           // Temporary variable to register hits. Should not be present in final states.
  #define XC_HIT_ON_FLUKA_CRYSTAL             335                           // Temporary variable to register hits. Should not be present in final states.
  #define XC_HIT_ON_GEANT4_BLOCK              336                           // Temporary variable to register hits. Should not be present in final states.
  #define XC_HIT_ON_GEANT4_COLL               337                           // Temporary variable to register hits. Should not be present in final states.
  #define XC_HIT_ON_GEANT4_CRYSTAL            338                           // Temporary variable to register hits. Should not be present in final states.
#endif /* XC_PARTICLE_STATE_H_I0 */
