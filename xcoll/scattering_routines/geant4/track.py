# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import xpart as xp
import time
import io

from ...headers.particle_states import LOST_WITHOUT_SPEC, LOST_ON_GEANT4_COLL, EXCITED_ION_STATE


def track(coll, particles):
    import xcoll as xc
    xc.geant4.engine.assert_ready_to_track_or_skip(coll, particles, _necessary_attributes=['geant4_id'])
    track_core(coll, particles)


def track_core(coll, part):
    import xcoll as xc
    xc.geant4.engine._g4link.clearData() # Clear the old data - bunch particles and hits

    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = part._delta.copy()
    npart = part._num_active_particles
    ndead = part._num_lost_particles

    if xc.geant4.engine.reentry_protection_enabled:
        ### remove this part after geant4 bug fixed
        coords = {
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
        buf = io.BytesIO() # Use numpy.savez to serialize
        np.savez(buf, **coords)
        buf.seek(0)
        result_blob = xc.geant4.engine._g4link.add_particles_and_collimate_return(
                                                buf.getvalue(), f'{coll.geant4_id}') # TODO: should geant4_id be a string or an int?
        result_buf = io.BytesIO(result_blob) # Deserialize
        products = np.load(result_buf)
    else:
        coords = [part.x, part.y, part.px, part.py,
                  part.zeta, delta_temp, part.chi,
                  part.charge_ratio, part.s,
                  part.pdg_id,part.particle_id, part.state,
                  part.at_element, part.at_turn]
        xc.geant4.engine._g4link.addParticles(coords)
        xc.geant4.engine._g4link.selectCollimator(f'{coll.geant4_id}')  # TODO: should geant4_id be a string or an int?
        xc.geant4.engine._g4link.collimate()
        products = xc.geant4.engine._g4link.collimateReturn(coords)

    npartsAliveAndDead = npart+ndead
    part.x[:npartsAliveAndDead] = products['x'][:npartsAliveAndDead]
    part.y[:npartsAliveAndDead] = products['y'][:npartsAliveAndDead]
    part.px[:npartsAliveAndDead] = products['px'][:npartsAliveAndDead]
    part.py[:npartsAliveAndDead] = products['py'][:npartsAliveAndDead]
    part.zeta[:npartsAliveAndDead] = products['zeta'][:npartsAliveAndDead]
    part.charge_ratio[:npartsAliveAndDead] = products['charge_ratio'][:npartsAliveAndDead]
    part.s[:npartsAliveAndDead] = products['s'][:npartsAliveAndDead]
    part.pdg_id[:npartsAliveAndDead] = products['pdg_id'][:npartsAliveAndDead]
    # parent particle id should NOT be updated for particles that were sent in
    part.at_element[:npartsAliveAndDead] = products['at_element'][:npartsAliveAndDead]
    part.at_turn[:npartsAliveAndDead] = products['at_turn'][:npartsAliveAndDead]
    part.state[:npartsAliveAndDead] = products['state'][:npartsAliveAndDead]
    new_delta = part.delta.copy()
    new_delta[:npartsAliveAndDead] = products['delta'][:npartsAliveAndDead]
    part.update_delta(new_delta)
    if products['x'] is None or products['x'][npartsAliveAndDead] == -9999:
        part.reorganize()
    else:
        mask = products['state'][npartsAliveAndDead:] > -999999
        new_particles = xp.Particles(_context=part._buffer.context,
                p0c = part.p0c[0],
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

    # Set the state of excited ions - not supported (will fail in BDSIM when resent)
    mask = (part.pdg_id > 999999999) & ((part.pdg_id % 10) != 0)
    part.state[mask] = EXCITED_ION_STATE

    # Set the dead state
    mask = part.state == LOST_WITHOUT_SPEC
    part.state[mask] = LOST_ON_GEANT4_COLL
