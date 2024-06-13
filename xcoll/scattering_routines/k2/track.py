# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #


def track(self, particles):
    from .pyk2 import pyk2_run
    npart = particles._num_active_particles
    if npart > self.k2engine.n_alloc:
        raise ValueError(f"Tracking {npart} particles but only {self.k2engine.n_alloc} allocated!")

    if npart < 1:
        return

    # Drift inactive front
    L = self.inactive_front
    if L > 0:
        rpp = particles.rpp[:npart]
        xp = particles.px[:npart] * rpp
        yp = particles.py[:npart] * rpp
        dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
        particles.x[:npart] += xp * L
        particles.y[:npart] += yp * L
        particles.s[:npart] += L
        particles.zeta[:npart] += dzeta*L

    # Go to collimator reference system (subtract closed orbit)
    x_part = particles.x[:npart] - self.dx
    xp_part = (particles.px[:npart] - self.dpx) * particles.rpp[:npart]
    y_part = particles.y[:npart] - self.dy
    yp_part = (particles.py[:npart] - self.dpy) * particles.rpp[:npart]
    s_part = 0 * x_part
    p_part = particles.energy[:npart] / 1e9 # Energy (not momentum) in GeV

    # Initialise arrays for FORTRAN call
    part_hit_pos = np.zeros(len(x_part), dtype=np.int32)
    part_hit_turn = np.zeros(len(x_part), dtype=np.int32)
    part_abs_pos = np.zeros(len(x_part), dtype=np.int32)
    part_abs_turn = np.zeros(len(x_part), dtype=np.int32)
    part_impact = np.zeros(len(x_part), dtype=np.float64)
    part_indiv = np.zeros(len(x_part), dtype=np.float64)
    part_linteract = np.zeros(len(x_part), dtype=np.float64)
    nhit_stage = np.zeros(len(x_part), dtype=np.int32)
    nabs_type = np.zeros(len(x_part), dtype=np.int32)
    linside = np.zeros(len(x_part), dtype=np.int32)

    # `linside` is an array of logicals in fortran. Beware of the fortran converion:
    # True <=> -1 (https://stackoverflow.com/questions/39454349/numerical-equivalent-of-true-is-1)

    pyk2_run(x_particles=x_part,
                    xp_particles=xp_part,
                    y_particles=y_part,
                    yp_particles=yp_part,
                    s_particles=s_part,
                    p_particles=p_part,              # confusing: this is ENERGY not momentum
                    part_hit_pos=part_hit_pos,       # ignore: sixtrack element of impact
                    part_hit_turn=part_hit_turn,     # ignore: turn of impact
                    part_abs_pos=part_abs_pos,       # ignore: sixtrack element of absorption
                    part_abs_turn=part_abs_turn,     # ignore: turn of absorption
                    part_impact=part_impact,         # impact parameter
                    part_indiv=part_indiv,           # particle divergence
                    part_linteract=part_linteract,   # interaction length
                    nhit_stage=nhit_stage,
                    nabs_type=nabs_type,
                    linside=linside,
                    icoll=self.icoll,
                    ie=1,                            # ignore: structure element index
                    c_length=self.active_length,
                    c_rotation=self.angle/180.*np.pi,
                    c_aperture=2*self.jaw,
                    c_offset=self.offset,
                    c_tilt=np.array([0,0], dtype=np.float64),
                    c_enom=particles.energy0[0]/1e9, # Reference energy
                    onesided=self.onesided,
                    random_generator_seed=-1, # skips rng re-initlization
                    )

    # Masks of hit and survived particles
    mask_lost = part_abs_turn > 0
    mask_hit = part_hit_turn > 0
    mask_not_hit = ~mask_hit
    mask_survived_hit = mask_hit & (~mask_lost)

    state_out = particles.state[:npart].copy()
    state_out[mask_lost] = -333
    particles.state[:npart] = state_out

    # Update particle energy
    ptau_out = particles.ptau[:npart].copy()
    e0 = particles.energy0[:npart]
    beta0 = particles.beta0[:npart]
    ptau_out[mask_survived_hit] = (
            p_part[mask_survived_hit] * 1e9 - e0[mask_survived_hit]
        ) / (
            e0[mask_survived_hit] * beta0[mask_survived_hit]
        )
    particles.ptau[:npart] = ptau_out

    # Update transversal coordinates (moving back from collimator frame)
    x_part_out = particles.x[:npart].copy()
    x_part_out[~mask_lost] = x_part[~mask_lost] + self.dx
    particles.x[:npart] = x_part_out

    y_part_out = particles.y[:npart].copy()
    y_part_out[~mask_lost] = y_part[~mask_lost] + self.dy
    particles.y[:npart] = y_part_out

    rpp_out = particles.rpp[:npart]
    px_out = particles.px[:npart]
    px_out[mask_survived_hit] = xp_part[mask_survived_hit]/rpp_out[mask_survived_hit]
    particles.px[:npart] = px_out

    py_out = particles.py[:npart]
    py_out[mask_survived_hit] = yp_part[mask_survived_hit]/rpp_out[mask_survived_hit]
    particles.py[:npart] = py_out
    
    # Update longitudinal coordinate zeta
    # Note: they have NOT been drifted yet in K2, so we have to do that manually
    rvv_out = particles.rvv[:npart]
    zeta_out = particles.zeta[:npart]
    zeta_out[mask_not_hit] += self.length*(
                        rvv_out[mask_not_hit] - (1 + ( xp_part[mask_not_hit]**2 + yp_part[mask_not_hit]**2)/2 ) 
                    )
    particles.zeta[:npart] = zeta_out
    
    s_out = particles.s[:npart]
    s_out[mask_not_hit] += self.length
    particles.s[:npart] = s_out

    # Drift inactive back
    L = self.inactive_back
    if L > 0:
        rpp = particles.rpp[:npart]
        xp = particles.px[:npart] * rpp
        yp = particles.py[:npart] * rpp
        dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
        particles.x[:npart] += xp * L
        particles.y[:npart] += yp * L
        particles.s[:npart] += L
        particles.zeta[:npart] += dzeta*L

    particles.reorganize()