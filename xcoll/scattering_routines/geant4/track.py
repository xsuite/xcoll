# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import xobjects as xo
import xpart as xp
import time
import io
from rpyc.utils.classic import obtain


def track(coll, particles):
    import xcoll as xc
    xc.geant4.engine.assert_ready_to_track(coll, particles, _necessary_attributes=['geant4_id'])
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

    coordinates = {
        'x': part.x,
        'y': part.y,
        'px': part.px,
        'py': part.py,
        'zeta': part.zeta,
        'delta': delta_temp,
        'chi': part.chi,
        'charge_ratio': part.charge_ratio,
        's': part.s,
        'pdg_id': part.pdg_id,
        'particle_id': part.particle_id,
        'state': part.state,
        'at_element': part.at_element,
        'at_turn': part.at_turn
    }
    # Use numpy.savez to serialize
    buf = io.BytesIO()
    np.savez(buf, **coordinates)
    buf.seek(0)

    result_blob = xc.geant4.engine.g4link.add_particles_and_collimate_return(buf.getvalue(), coll.geant4_id)

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

    npartsAliveAndDead = npart+ndead
    part.x[:npartsAliveAndDead] = products['x'][:npartsAliveAndDead]
    part.y[:npartsAliveAndDead] = products['y'][:npartsAliveAndDead]
    part.px[:npartsAliveAndDead] = products['px'][:npartsAliveAndDead]
    part.py[:npartsAliveAndDead] = products['py'][:npartsAliveAndDead]
    part.zeta[:npartsAliveAndDead] = products['zeta'][:npartsAliveAndDead]
    part.charge_ratio[:npartsAliveAndDead] = products['charge_ratio'][:npartsAliveAndDead]
    part.s[:npartsAliveAndDead] = products['s'][:npartsAliveAndDead]
    part.pdg_id[:npartsAliveAndDead] = products['pdg_id'][:npartsAliveAndDead]
    part.parent_particle_id[:npartsAliveAndDead] = products['parent_particle_id'][:npartsAliveAndDead]
    part.at_element[:npartsAliveAndDead] = products['at_element'][:npartsAliveAndDead]
    part.at_turn[:npartsAliveAndDead] = products['at_turn'][:npartsAliveAndDead]
    part.state[:npartsAliveAndDead] = products['state'][:npartsAliveAndDead]
    new_delta = part.delta.copy()
    new_delta[:npartsAliveAndDead] = products['delta'][:npartsAliveAndDead]
    part.update_delta(new_delta)
    print(set(part.pdg_id))
    if products['x'] is None or products['x'][npartsAliveAndDead] == -9999:
        part.reorganize()
    else:
        mask = products['state'][npartsAliveAndDead:] > -999999
        new_particles = xp.Particles(_context=particles._buffer.context,
                p0c = part.p0c[0], # TODO: Should we check that 
                                        #       they are all the same?
                mass0 = part.mass0,
                q0 = part.q0,
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

        part.add_particles(new_particles)
    mask = (part.pdg_id > 999999999) & ((part.pdg_id % 10) != 0)
    part.state[mask] = -4000
