"""Particle-based benchmark diagnostics that exclude unused capacity slots.

An Xpart particle array may have a much larger capacity than the number of
primaries initially populated (for example, to accommodate secondaries from
Geant4 or FLUKA). Unallocated slots have negative sentinel states and must
not be mistaken for physical particle losses.
"""

import numpy as np


def count_lost_primaries(particles, n_primaries: int) -> int:
    """Count lost primaries, identified by their original particle IDs.

    Secondary particles are deliberately not counted in the denominator.
    Fail explicitly when the primary ID range is not completely represented,
    rather than publish a misleading loss fraction.
    """
    if n_primaries < 1:
        raise ValueError('n_primaries must be positive')
    state = np.asarray(particles.state)
    ids = np.asarray(particles.particle_id)
    if state.shape != ids.shape:
        raise ValueError('Particle state and ID arrays have different shapes')
    mask = (ids >= 0) & (ids < n_primaries)
    if np.count_nonzero(mask) != n_primaries or np.unique(ids[mask]).size != n_primaries:
        raise ValueError(f'Expected {n_primaries} distinct original primaries; '
                         f'found {np.count_nonzero(mask)} occupied slots with '
                         f'{np.unique(ids[mask]).size} distinct primary IDs')
    return int(np.count_nonzero(mask & (state <= 0)))
