# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from .engine import K2Engine


def _drift(coll, particles, length):
    coll._equivalent_drift.length = length
    coll._equivalent_drift.track(particles)


def track(coll, part):
    if coll.gap is None or not coll._active or not coll._tracking:
        _drift(coll, part, coll.length)
        return

    try:
        from .pyk2f import pyk2_track
    except ImportError:
       raise ValueError("Failed to import pyK2. Cannot track.")

    if not K2Engine.is_running():
        raise ValueError("K2Engine is not running. Start it first with "
                         "`K2Engine.start(line=line)`.")

    k2_id = coll._k2_id
    if k2_id == -1:
        raise ValueError("Collimator not registered in K2Engine!")

    npart = part._num_active_particles
    n_alloc = K2Engine().n_alloc
    if npart > n_alloc:
        raise ValueError(f"Tracking {npart} particles but only {n_alloc} allocated!")

    if npart < 1:
        return

    # Drift to the middle of the collimator (as this is what K2 expects)
    _drift(coll, part, coll.length/2)

    # Initialise arrays for FORTRAN call
    x_particles     = np.zeros(npart, dtype=np.float64)
    xp_particles    = np.zeros(npart, dtype=np.float64)
    y_particles     = np.zeros(npart, dtype=np.float64)
    yp_particles    = np.zeros(npart, dtype=np.float64)
    e_particles     = np.zeros(npart, dtype=np.float64)
    p_particles     = np.zeros(npart, dtype=np.float64)
    delta_particles = np.zeros(npart, dtype=np.float64)
    rvv_particles   = np.zeros(npart, dtype=np.float64)
    rpp_particles   = np.zeros(npart, dtype=np.float64)
    hit             = np.zeros(npart, dtype=np.int32)
    absorbed        = np.zeros(npart, dtype=np.int32)

    x_particles[:npart]     = part.x.copy()
    xp_particles[:npart]    = part.kin_xprime.copy()
    y_particles[:npart]     = part.y.copy()
    yp_particles[:npart]    = part.kin_yprime.copy()
    e_particles[:npart]     = part.energy.copy()
    p_particles[:npart]     = (1 + part.delta)*part.p0c
    delta_particles[:npart] = part.delta.copy()
    rvv_particles[:npart]   = part.rvv.copy()
    rpp_particles[:npart]   = part.rpp.copy()

    # `linside` is an array of logicals in fortran. Beware of the fortran converion:
    # True <=> -1 (https://stackoverflow.com/questions/39454349/numerical-equivalent-of-true-is-1)
    pyk2_track(num_particles=npart,
               x_particles=x_particles,
               xp_particles=xp_particles,
               y_particles=y_particles,
               yp_particles=yp_particles,
               e_particles=e_particles,
               p_particles=p_particles,
               delta_particles=delta_particles,
               rvv_particles=rvv_particles,
               rpp_particles=rpp_particles,
               idx=k2_id,
               hit=hit,
               absorbed=absorbed)

    # Update particles object
    part.add_to_energy(e_particles - part.e[:npart])
    rpp = part.rpp[:npart]
    part.x[:npart]  = x_particles.copy()
    part.px[:npart] = xp_particles / rpp
    part.y[:npart]  = y_particles.copy()
    part.py[:npart] = yp_particles / rpp

    # Masks of hit and survived particles
    mask_lost = absorbed > 0
    mask_hit  = hit > 0
    mask_not_hit = ~mask_hit
    mask_survived_hit = mask_hit & (~mask_lost)

    state_out = part.state[:npart].copy()
    state_out[mask_lost] = -339
    part.state[:npart]   = state_out

    particles.reorganize()

    # Drift to the end of the collimator (as K2 returns particles in the middle)
    _drift(coll, part, coll.length/2)
