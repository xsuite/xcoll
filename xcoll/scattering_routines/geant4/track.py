# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp


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
    secondaries_x = np.zeros(len(particles.x))-9999
    secondaries_y = np.zeros(len(particles.x))-9999
    secondaries_px = np.zeros(len(particles.x))-9999
    secondaries_py = np.zeros(len(particles.x))-9999
    secondaries_zeta = np.zeros(len(particles.x))-9999
    secondaries_delta = np.zeros(len(particles.x))-9999
    secondaries_charge_ratio = np.zeros(len(particles.x))-9999
    secondaries_s = np.zeros(len(particles.x))-9999
    secondaries_pdg_id = np.zeros(len(particles.x))-9999
    secondaries_parent_particle_id = np.zeros(len(particles.x))-9999
    secondaries_at_element = np.zeros(len(particles.x))-9999
    secondaries_at_turn = np.zeros(len(particles.x))-9999
    secondaries_mass_ratio = np.zeros(len(particles.x))-9999
    secondaries_state = np.zeros(len(particles.x))-9999

    products = Geant4Engine.g4link.collimateReturn(particles.x, particles.y, particles.px, particles.py,
                        particles.zeta, delta_temp, particles.chi,
                        particles.charge_ratio, particles.s,
                        particles.pdg_id,particles.particle_id, particles.state,
                        particles.at_element, particles.at_turn,
                        secondaries_x,secondaries_y,secondaries_px,secondaries_py,secondaries_zeta,
                        secondaries_delta,secondaries_charge_ratio,secondaries_s,
                        secondaries_pdg_id,secondaries_parent_particle_id,secondaries_at_element,
                        secondaries_at_turn,secondaries_mass_ratio,secondaries_state)
    # Force the update using the private member _delta
    # as the update_delta method only updates the delta for active particles
    particles.x[:npart] = secondaries_x[:npart]
    particles.y[:npart] = secondaries_y[:npart]
    particles.px[:npart] = secondaries_px[:npart]
    particles.py[:npart] = secondaries_py[:npart]
    particles.zeta[:npart] = secondaries_zeta[:npart]
    particles.charge_ratio[:npart] = secondaries_charge_ratio[:npart]
    particles.s[:npart] = secondaries_s[:npart]
    particles.pdg_id[:npart] = secondaries_pdg_id[:npart]
    particles.parent_particle_id[:npart] = secondaries_parent_particle_id[:npart]
    particles.at_element[:npart] = secondaries_at_element[:npart]
    particles.at_turn[:npart] = secondaries_at_turn[:npart]
    #particles.mass_ratio[:npart] = secondaries_mass_ratio[:npart]
    particles.state[:npart] = secondaries_state[:npart]
    new_delta = particles.delta.copy()
    new_delta[:npart] = secondaries_delta[:npart]
    particles.update_delta(new_delta)
    mask = np.abs(secondaries_parent_particle_id)<100
    if secondaries_x is None or secondaries_x[0] == 0:
        particles.reorganize()
    else:
        mask = secondaries_state[npart:] > -999999
        new_particles = xp.Particles(_context=particles._buffer.context,
                p0c = particles.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = particles.mass0,
                q0 = particles.q0,
                s = secondaries_s[npart:][mask],
                x = secondaries_x[npart:][mask],
                px = secondaries_px[npart:][mask],
                y = secondaries_y[npart:][mask],
                py = secondaries_py[npart:][mask],
                zeta = secondaries_zeta[npart:][mask],
                delta = secondaries_delta[npart:][mask],
                mass_ratio = secondaries_mass_ratio[npart:][mask],
                charge_ratio = secondaries_charge_ratio[npart:][mask],
                at_element = secondaries_at_element[npart:][mask],
                at_turn = secondaries_at_turn[npart:][mask],
                parent_particle_id = secondaries_parent_particle_id[npart:][mask])

        particles.add_particles(new_particles)