# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp
import time


def track(coll, particles):
    import xcoll as xc
    xc.geant4.engine.assert_ready_to_track(self, coll, particles, _necessary_attributes=['geant4_id'])
    track_core(coll, particles)


def track_core(coll, part):
    import xcoll as xc
    xc.geant4.engine.g4link.clearData() # Clear the old data - bunch particles and hits

    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = part._delta.copy()
    npart = part._num_active_particles
    ndead = part._num_lost_particles
    # Using a list allows to package the required coordinates without copying
    coordinates = np.array([part.x, part.y, part.px, part.py,
                   part.zeta, delta_temp, part.chi,
                   part.charge_ratio, part.s,
                   part.pdg_id,part.particle_id, part.state,
                   part.at_element, part.at_turn])
    xc.geant4.engine.g4link.addParticles2(coordinates)
    #xc.geant4.engine.g4link.addParticles(part.x, part.y, part.px, part.py,
    #                    part.zeta, delta_temp, part.chi,
    #                    part.charge_ratio, part.s,
    #                    part.pdg_id,part.particle_id, part.state,
    #                    part.at_element, part.at_turn)

    # The collimators must be defined already in the g4manager
    xc.geant4.engine.g4link.selectCollimator(coll.geant4_id)

    xc.geant4.engine.g4link.collimate() # Performs the physical interaction simulation

    # Modifies the primary coordinates in place and returns a list of arrays for the
    # coordinates of the secondary part.
    #products = g4link.collimateReturn(coordinates)
    #secondaries_x = np.zeros(len(part.x))-9999
    #secondaries_y = np.zeros(len(part.x))-9999
    #secondaries_px = np.zeros(len(part.x))-9999
    #secondaries_py = np.zeros(len(part.x))-9999
    #secondaries_zeta = np.zeros(len(part.x))-9999
    #secondaries_delta = np.zeros(len(part.x))-9999
    #secondaries_charge_ratio = np.zeros(len(part.x))-9999
    #secondaries_s = np.zeros(len(part.x))-9999
    #secondaries_pdg_id = np.zeros(len(part.x))-9999
    #secondaries_parent_particle_id = np.zeros(len(part.x))-9999
    #secondaries_at_element = np.zeros(len(part.x))-9999
    #secondaries_at_turn = np.zeros(len(part.x))-9999
    #secondaries_mass_ratio = np.zeros(len(part.x))-9999
    #secondaries_state = np.zeros(len(part.x))-9999
    t0 = time.time()
    #products = xc.geant4.engine.g4link.collimateReturn(part.x, part.y, part.px, part.py,
    #                    part.zeta, delta_temp, part.chi,
    #                    part.charge_ratio, part.s,
    #                    part.pdg_id,part.particle_id, part.state,
    #                    part.at_element, part.at_turn,
    #                    secondaries_x,secondaries_y,secondaries_px,secondaries_py,secondaries_zeta,
    #                    secondaries_delta,secondaries_charge_ratio,secondaries_s,
    #                    secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element,
    #                    secondaries_at_turn,secondaries_mass_ratio,secondaries_state)
    secondaries_x,secondaries_y,secondaries_px,secondaries_py,secondaries_zeta, \
    secondaries_delta,secondaries_charge_ratio,secondaries_s, \
    secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element, \
    secondaries_at_turn,secondaries_mass_ratio,secondaries_state = xc.geant4.engine.g4link.collimateReturn(part.x,
                                                                                                       part.y, part.px, part.py,
                                                                                                       part.zeta, delta_temp, part.chi,
                                                                                                       part.charge_ratio, part.s,
                                                                                                       part.pdg_id,part.particle_id, part.state,
                                                                                                       part.at_element, part.at_turn)
    t1 = time.time()
    print(f'time for collimatereturn: {t1-t0}')

    # Force the update using the private member _delta
    # as the update_delta method only updates the delta for active particles
    #particles_prev = part.copy()
    temp_secondaries_x = np.asarray(secondaries_x)
    temp_secondaries_y = np.asarray(secondaries_y)
    temp_secondaries_px = np.asarray(secondaries_px)
    temp_secondaries_py = np.asarray(secondaries_py)
    temp_secondaries_zeta = np.asarray(secondaries_zeta)
    temp_secondaries_delta = np.asarray(secondaries_delta)
    temp_secondaries_charge_ratio = np.asarray(secondaries_charge_ratio)
    temp_secondaries_s = np.asarray(secondaries_s)
    temp_secondaries_pdg_id = np.asarray(secondaries_pdg_id)
    temp_secondaries_parent_particle_id = np.asarray(secondaries_parent_particle_id)
    temp_secondaries_at_element = np.asarray(secondaries_at_element)
    temp_secondaries_at_turn = np.asarray(secondaries_at_turn)
    temp_secondaries_mass_ratio = np.asarray(secondaries_mass_ratio)
    temp_secondaries_state = np.asarray(secondaries_state)

    npartsAliveAndDead = npart+ndead
    part.x[:npartsAliveAndDead] = temp_secondaries_x[:npartsAliveAndDead]
    part.y[:npartsAliveAndDead] = temp_secondaries_y[:npartsAliveAndDead]
    part.px[:npartsAliveAndDead] = temp_secondaries_px[:npartsAliveAndDead]
    part.py[:npartsAliveAndDead] = temp_secondaries_py[:npartsAliveAndDead]
    part.zeta[:npartsAliveAndDead] = temp_secondaries_zeta[:npartsAliveAndDead]
    part.charge_ratio[:npartsAliveAndDead] = temp_secondaries_charge_ratio[:npartsAliveAndDead]
    part.s[:npartsAliveAndDead] = temp_secondaries_s[:npartsAliveAndDead]
    part.pdg_id[:npartsAliveAndDead] = temp_secondaries_pdg_id[:npartsAliveAndDead]
    part.parent_particle_id[:npartsAliveAndDead] = temp_secondaries_parent_particle_id[:npartsAliveAndDead]
    part.at_element[:npartsAliveAndDead] = temp_secondaries_at_element[:npartsAliveAndDead]
    part.at_turn[:npartsAliveAndDead] = temp_secondaries_at_turn[:npartsAliveAndDead]
    #part.mass_ratio[:npart] = secondaries_mass_ratio[:npart]
    part.state[:npartsAliveAndDead] = temp_secondaries_state[:npartsAliveAndDead]
    new_delta = part.delta.copy()
    new_delta[:npartsAliveAndDead] = temp_secondaries_delta[:npartsAliveAndDead]
    part.update_delta(new_delta)

    if temp_secondaries_x is not None and temp_secondaries_x[npartsAliveAndDead] != -9999:
        mask = temp_secondaries_state[npartsAliveAndDead:] > -999999
        new_particles = xp.Particles(_context=part._buffer.context,
                p0c = part.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = part.mass0,
                q0 = part.q0,
                s = temp_secondaries_s[npartsAliveAndDead:][mask],
                x = temp_secondaries_x[npartsAliveAndDead:][mask],
                px = temp_secondaries_px[npartsAliveAndDead:][mask],
                y = temp_secondaries_y[npartsAliveAndDead:][mask],
                py = temp_secondaries_py[npartsAliveAndDead:][mask],
                zeta = temp_secondaries_zeta[npartsAliveAndDead:][mask],
                delta = temp_secondaries_delta[npartsAliveAndDead:][mask],
                mass_ratio = temp_secondaries_mass_ratio[npartsAliveAndDead:][mask],
                charge_ratio = temp_secondaries_charge_ratio[npartsAliveAndDead:][mask],
                at_element = temp_secondaries_at_element[npartsAliveAndDead:][mask],
                at_turn = temp_secondaries_at_turn[npartsAliveAndDead:][mask],
                parent_particle_id = temp_secondaries_parent_particle_id[npartsAliveAndDead:][mask],
                pdg_id = temp_secondaries_pdg_id[npartsAliveAndDead:][mask])
        new_particles._init_random_number_generator()
        part.add_particles(new_particles)

    part.reorganize()
