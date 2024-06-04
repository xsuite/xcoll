# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp


def track(collimator_id, particles):

    assert isinstance(particles._buffer.context, xo.ContextCpu)
    assert particles._num_active_particles >= 0

    from .engine import Geant4Engine

    g4engine = Geant4Engine() # Get the singleton engine instance
    if not g4engine.connected:
        raise ValueError("Geant4Engine not linked to BDSIM! Cannot track.")

    g4engine.assert_geometry()

    g4link = g4engine.g4link
    g4link.clearData() # Clear the old data - bunch particles and hits

    print(f"Processing collimator: {collimator_id}")
    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = particles._delta.copy()

    # Using a list allows to package the required coordinates without copying
    coordinates = [particles.x, particles.y, particles.px, particles.py,
                   particles.zeta, delta_temp, particles.chi,
                   particles.charge_ratio, particles.s,
                   particles.pdg_id,particles.particle_id, particles.state,
                   particles.at_element, particles.at_turn]



    g4link.addParticles(coordinates)
    # The collimators must be defined already in the g4manager
    g4link.selectCollimator(collimator_id)

    g4link.collimate() # Performs the physical interaction simulation

    # Modifies the primary coordinates in place and returns a list of arrays for the
    # coordinates of the secondary particles.
    products = g4link.collimateReturn(coordinates)

    # Force the update using the private member _delta
    # as the update_delta method only updates the delta for active particles
    particles._delta[:len(delta_temp)] = delta_temp
    particles.update_delta(delta_temp)

    # TODO: This should work also when no products are there
    #       Particles reorganization should still happen

    if products is None or products['x'].size == 0:
        particles.reorganize()
    else:
        new_particles = xp.Particles(_context=particles._buffer.context,
                p0c = particles.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = particles.mass0,
                q0 = particles.q0,
                s = products['s'],
                x = products['x'],
                px = products['px'],
                y = products['y'],
                py = products['py'],
                zeta = products['zeta'],
                delta = products['delta'],
                mass_ratio = products['mass_ratio'],
                charge_ratio = products['charge_ratio'],
                at_element = products['at_element'],
                at_turn = products['at_turn'],
                parent_particle_id = products['parent_particle_id'])

        particles.add_particles(new_particles)
