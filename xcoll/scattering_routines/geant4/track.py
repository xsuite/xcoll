# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp
import time
import io
from rpyc.utils.classic import obtain


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
    Geant4Engine.g4link.clearData() # Clear the old data - bunch particles and hits
    print(f"Processing collimator: {coll.geant4_id}")
    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = particles._delta.copy()
    npart = particles._num_active_particles
    ndead = particles._num_lost_particles

    coordinates = {
        'x': particles.x,
        'y': particles.y,
        'px': particles.px,
        'py': particles.py,
        'zeta': particles.zeta,
        'delta': delta_temp,
        'chi': particles.chi,
        'charge_ratio': particles.charge_ratio,
        's': particles.s,
        'pdg_id': particles.pdg_id,
        'particle_id': particles.particle_id,
        'state': particles.state,
        'at_element': particles.at_element,
        'at_turn': particles.at_turn
    }
    # Use numpy.savez to serialize
    buf = io.BytesIO()
    np.savez(buf, **coordinates)
    buf.seek(0)

    result_blob = Geant4Engine.g4link.add_particles_and_collimate_return(buf.getvalue(), coll.geant4_id)

    # Deserialize
    result_buf = io.BytesIO(result_blob)
    products = np.load(result_buf)

    if False:
        temp_secondaries_x = products['x']
        temp_secondaries_y = products['y']
        temp_secondaries_px = products['px']
        temp_secondaries_py = products['py']
        temp_secondaries_zeta = products['zeta']
        temp_secondaries_delta = products['delta']
        temp_secondaries_charge_ratio = products['charge_ratio']
        temp_secondaries_s = products['s']
        temp_secondaries_pdg_id = products['pdg_id']
        temp_secondaries_parent_particle_id = products['parent_particle_id']
        temp_secondaries_at_element = products['at_element']
        temp_secondaries_at_turn = products['at_turn']
        temp_secondaries_mass_ratio = products['mass_ratio']
        temp_secondaries_state = products['state']

    t5 = time.time()
    npartsAliveAndDead = npart+ndead
    particles.x[:npartsAliveAndDead] = products['x'][:npartsAliveAndDead]
    particles.y[:npartsAliveAndDead] = products['y'][:npartsAliveAndDead]
    particles.px[:npartsAliveAndDead] = products['px'][:npartsAliveAndDead]
    particles.py[:npartsAliveAndDead] = products['py'][:npartsAliveAndDead]
    particles.zeta[:npartsAliveAndDead] = products['zeta'][:npartsAliveAndDead]
    particles.charge_ratio[:npartsAliveAndDead] = products['charge_ratio'][:npartsAliveAndDead]
    particles.s[:npartsAliveAndDead] = products['s'][:npartsAliveAndDead]
    particles.pdg_id[:npartsAliveAndDead] = products['pdg_id'][:npartsAliveAndDead]
    particles.parent_particle_id[:npartsAliveAndDead] = products['parent_particle_id'][:npartsAliveAndDead]
    particles.at_element[:npartsAliveAndDead] = products['at_element'][:npartsAliveAndDead]
    particles.at_turn[:npartsAliveAndDead] = products['at_turn'][:npartsAliveAndDead]
    particles.state[:npartsAliveAndDead] = products['state'][:npartsAliveAndDead]
    new_delta = particles.delta.copy()
    new_delta[:npartsAliveAndDead] = products['delta'][:npartsAliveAndDead]
    particles.update_delta(new_delta)
    print(set(particles.pdg_id))
    if products['x'] is None or products['x'][npartsAliveAndDead] == -9999:
        particles.reorganize()
    else:
        mask = products['state'][npartsAliveAndDead:] > -999999
        new_particles = xp.Particles(_context=particles._buffer.context,
                p0c = particles.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = particles.mass0,
                q0 = particles.q0,
                s = products['s'][npartsAliveAndDead:][mask],
                x = products['x'][npartsAliveAndDead:][mask],
                px = products['px'][npartsAliveAndDead:][mask],
                y = products['y'][npartsAliveAndDead:][mask],
                py = products['py'][npartsAliveAndDead:][mask],
                zeta = products['zeta'][npartsAliveAndDead:][mask],
                delta = products['delta'][npartsAliveAndDead:][mask],
                mass_ratio = products['mass_ratio'][npartsAliveAndDead:][mask],
                charge_ratio = products['charge_ratio'][npartsAliveAndDead:][mask],
                at_element = products['at_element'][npartsAliveAndDead:][mask],
                at_turn = products['at_turn'][npartsAliveAndDead:][mask],
                parent_particle_id = products['parent_particle_id'][npartsAliveAndDead:][mask],
                pdg_id = products['pdg_id'][npartsAliveAndDead:][mask])

        particles.add_particles(new_particles)
    mask = ((particles.pdg_id > 999999999) & ((particles.pdg_id % 10) != 0) & (particles.pdg_id != -999999999))
    particles.state[mask] = -4000