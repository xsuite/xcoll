# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp
import time


def track(coll, particles):
    from ...beam_elements import Geant4Collimator
    if not isinstance(coll, Geant4Collimator):
        raise ValueError("Collimator is not a Geant4Collimator!\nCannot use Geant4 to track.")

    if not coll.active or not coll._tracking:
        coll._equivalent_drift.track(particles)
        return

    npart = particles._num_active_particles
    if npart == 0:
        return

    assert isinstance(particles._buffer.context, xo.ContextCpu)

    # Check the server and whether it's initialised correctly
    from .engine import Geant4Engine
    #g4engine = Geant4Engine().instance
    #if not Geant4Engine.g4link.is_running():
    #    raise ValueError(f"Geant4Engine not yet running!\nPlease do this first, by calling "
    #                   + f"xcoll.Geant4Engine.start().\n(id: {id(g4engine)})")

    #g4link = g4engine.g4link
    Geant4Engine.g4link.clearData() # Clear the old data - bunch particles and hits
    print(f"Processing collimator: {coll.geant4_id}")
    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = particles._delta.copy()
    npart = particles._num_active_particles
    ndead = particles._num_lost_particles
    # Using a list allows to package the required coordinates without copying
    coordinates = np.array([particles.x, particles.y, particles.px, particles.py,
                   particles.zeta, delta_temp, particles.chi,
                   particles.charge_ratio, particles.s,
                   particles.pdg_id,particles.particle_id, particles.state,
                   particles.at_element, particles.at_turn])
    Geant4Engine.g4link.addParticles2(coordinates)
    #Geant4Engine.g4link.addParticles(particles.x, particles.y, particles.px, particles.py,
    #                    particles.zeta, delta_temp, particles.chi,
    #                    particles.charge_ratio, particles.s,
    #                    particles.pdg_id,particles.particle_id, particles.state,
    #                    particles.at_element, particles.at_turn)

    # The collimators must be defined already in the g4manager
    Geant4Engine.g4link.selectCollimator(coll.geant4_id)

    Geant4Engine.g4link.collimate() # Performs the physical interaction simulation

    # Modifies the primary coordinates in place and returns a list of arrays for the
    # coordinates of the secondary particles.
    #products = g4link.collimateReturn(coordinates)
    #secondaries_x = np.zeros(len(particles.x))-9999
    #secondaries_y = np.zeros(len(particles.x))-9999
    #secondaries_px = np.zeros(len(particles.x))-9999
    #secondaries_py = np.zeros(len(particles.x))-9999
    #secondaries_zeta = np.zeros(len(particles.x))-9999
    #secondaries_delta = np.zeros(len(particles.x))-9999
    #secondaries_charge_ratio = np.zeros(len(particles.x))-9999
    #secondaries_s = np.zeros(len(particles.x))-9999
    #secondaries_pdg_id = np.zeros(len(particles.x))-9999
    #secondaries_parent_particle_id = np.zeros(len(particles.x))-9999
    #secondaries_at_element = np.zeros(len(particles.x))-9999
    #secondaries_at_turn = np.zeros(len(particles.x))-9999
    #secondaries_mass_ratio = np.zeros(len(particles.x))-9999
    #secondaries_state = np.zeros(len(particles.x))-9999
    t0 = time.time()
    #products = Geant4Engine.g4link.collimateReturn(particles.x, particles.y, particles.px, particles.py,
    #                    particles.zeta, delta_temp, particles.chi,
    #                    particles.charge_ratio, particles.s,
    #                    particles.pdg_id,particles.particle_id, particles.state,
    #                    particles.at_element, particles.at_turn,
    #                    secondaries_x,secondaries_y,secondaries_px,secondaries_py,secondaries_zeta,
    #                    secondaries_delta,secondaries_charge_ratio,secondaries_s,
    #                    secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element,
    #                    secondaries_at_turn,secondaries_mass_ratio,secondaries_state)
    secondaries_x,secondaries_y,secondaries_px,secondaries_py,secondaries_zeta, \
    secondaries_delta,secondaries_charge_ratio,secondaries_s, \
    secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element, \
    secondaries_at_turn,secondaries_mass_ratio,secondaries_state = Geant4Engine.g4link.collimateReturn(particles.x,
                                                                                                       particles.y, particles.px, particles.py,
                                                                                                       particles.zeta, delta_temp, particles.chi,
                                                                                                       particles.charge_ratio, particles.s,
                                                                                                       particles.pdg_id,particles.particle_id, particles.state,
                                                                                                       particles.at_element, particles.at_turn)
    t1 = time.time()
    print(f'time for collimatereturn: {t1-t0}')

    # Force the update using the private member _delta
    # as the update_delta method only updates the delta for active particles
    #particles_prev = particles.copy()
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
    particles.x[:npartsAliveAndDead] = temp_secondaries_x[:npartsAliveAndDead]
    particles.y[:npartsAliveAndDead] = temp_secondaries_y[:npartsAliveAndDead]
    particles.px[:npartsAliveAndDead] = temp_secondaries_px[:npartsAliveAndDead]
    particles.py[:npartsAliveAndDead] = temp_secondaries_py[:npartsAliveAndDead]
    particles.zeta[:npartsAliveAndDead] = temp_secondaries_zeta[:npartsAliveAndDead]
    particles.charge_ratio[:npartsAliveAndDead] = temp_secondaries_charge_ratio[:npartsAliveAndDead]
    particles.s[:npartsAliveAndDead] = temp_secondaries_s[:npartsAliveAndDead]
    particles.pdg_id[:npartsAliveAndDead] = temp_secondaries_pdg_id[:npartsAliveAndDead]
    particles.parent_particle_id[:npartsAliveAndDead] = temp_secondaries_parent_particle_id[:npartsAliveAndDead]
    particles.at_element[:npartsAliveAndDead] = temp_secondaries_at_element[:npartsAliveAndDead]
    particles.at_turn[:npartsAliveAndDead] = temp_secondaries_at_turn[:npartsAliveAndDead]
    #particles.mass_ratio[:npart] = secondaries_mass_ratio[:npart]
    particles.state[:npartsAliveAndDead] = temp_secondaries_state[:npartsAliveAndDead]
    new_delta = particles.delta.copy()
    new_delta[:npartsAliveAndDead] = temp_secondaries_delta[:npartsAliveAndDead]
    particles.update_delta(new_delta)

    if temp_secondaries_x is None or temp_secondaries_x[npartsAliveAndDead] == -9999:
        particles.reorganize()
    else:
        mask = temp_secondaries_state[npartsAliveAndDead:] > -999999
        new_particles = xp.Particles(_context=particles._buffer.context,
                p0c = particles.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = particles.mass0,
                q0 = particles.q0,
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

        particles.add_particles(new_particles)