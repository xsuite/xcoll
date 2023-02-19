import numpy as np


# =================================================================== #
# ===============================  K2  ============================== #
# =================================================================== #
#  Even though K2 has info on the length of the collimator, and even  #
#  uses that info, it is still a thin element in the original code.   #
#  Hence, our thick tracking will always be equivalent to:            #
#    - drift first half length                                        #
#    - scatter                                                        #
#    - drift second half length                                       #
#  To make things worse, the shift to the Closed Orbit happens in the #
#  centre of the collimator...                                        #
# =================================================================== #

def drift_4d(x, y, xp, yp, length):
    x += xp * length
    y += yp * length
    return

def drift_zeta(zeta, rvv, xp, yp, length):
    rv0v = 1./rvv
    dzeta = 1 - rv0v * (1 + (xp**2 + yp**2)/2 )
    zeta += length * dzeta
    return zeta

def drift_6d(particles, length):
    npart = particles._num_active_particles
    rpp = particles.rpp[:npart]
    xp = particles.px[:npart] * rpp
    yp = particles.py[:npart] * rpp
    dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
    particles.x[:npart] += xp * length
    particles.y[:npart] += yp * length
    particles.s[:npart] += length
    particles.zeta[:npart] += dzeta*length
    return

def track(k2collimator, particles):  # TODO: write impacts
        from ...beam_elements import K2Collimator, K2Crystal
        if not isinstance(k2collimator, K2Collimator) and not isinstance(k2collimator, K2Crystal):
            raise ValueError("Collimator is neither a K2Collimator nor a K2Crystal!\nCannot use K2 to track.")

        from .engine import K2Engine
        npart = particles._num_active_particles
        if npart > K2Engine.instance._capacity:
            raise ValueError(f"Tracking {npart} particles but only {K2Engine.instance._capacity} allocated!")
        if npart == 0:
            return
        
        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        if not k2collimator.is_active:
            # Drift full length
            drift_6d(particles, k2collimator.length)
        else:
            drift_6d(particles, k2collimator.inactive_front)
            track_core(k2collimator, particles)
            drift_6d(particles, k2collimator.inactive_back)
            particles.reorganize()
        return



def track_core(k2collimator, particles):
    npart = particles._num_active_particles
    try:
        from .pyk2f import pyk2_startcry, pyk2_run
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

    length = k2collimator.active_length

    # Get coordinates
    x_part  = particles.x[:npart].copy()
    y_part  = particles.y[:npart].copy()
    xp_part = particles.px[:npart] * particles.rpp[:npart]
    yp_part = particles.py[:npart] * particles.rpp[:npart]
    s_part  = 0 * x_part
    e_part  = particles.energy[:npart].copy() / 1e9 # Energy in GeV
    rpp_in  = particles.rpp[:npart].copy()
    rvv_in  = particles.rvv[:npart].copy()
    xp_in   = xp_part.copy()
    yp_in   = yp_part.copy()

    # Move to closed orbit
    x_part  -= k2collimator.dx
    y_part  -= k2collimator.dy

    # Initialise arrays for FORTRAN call
    part_hit       = np.zeros(len(x_part), dtype=np.int32)
    part_abs       = np.zeros(len(x_part), dtype=np.int32)
    part_impact    = np.zeros(len(x_part), dtype=float)
    part_indiv     = np.zeros(len(x_part), dtype=float)
    part_linteract = np.zeros(len(x_part), dtype=float)
    nhit_stage     = np.zeros(len(x_part), dtype=np.int32)
    nabs_type      = np.zeros(len(x_part), dtype=np.int32)
    linside        = np.zeros(len(x_part), dtype=np.int32)

    if k2collimator.jaw_F_L != k2collimator.jaw_B_L or k2collimator.jaw_F_R != k2collimator.jaw_B_R:
        raise NotImplementedError
    opening = k2collimator.jaw_F_L - k2collimator.jaw_F_R
    offset = k2collimator.offset + ( k2collimator.jaw_F_L + k2collimator.jaw_F_R )/2

    # Get material properties
    zatom    = k2collimator.material.Z
    anuc     = k2collimator.material.A
    rho      = k2collimator.material.density
    exenergy = k2collimator.material.excitation_energy
    emr      = k2collimator.material.nuclear_radius
    bnref    = k2collimator.material.nuclear_elastic_slope
    csref0   = k2collimator.material.cross_section[0]
    csref1   = k2collimator.material.cross_section[1]
    csref5   = k2collimator.material.cross_section[5]
    hcut     = k2collimator.material.hcut

    # Get crystal parameters
    from ...beam_elements import K2Crystal
    from xcoll import CrystalMaterial
    if isinstance(k2collimator, K2Crystal):
        if not isinstance(k2collimator.material, CrystalMaterial):
            raise ValueError(f"The collimator material {k2collimator.material.name} cannot be used as a crystal!")
        dlri     = k2collimator.material.crystal_radiation_length
        dlyi     = k2collimator.material.crystal_nuclear_length
        ai       = k2collimator.material.crystal_plane_distance
        eUm      = k2collimator.material.crystal_potential
        collnt   = k2collimator.material.nuclear_collision_length
        is_crystal = True
        radl     = 0
    else:
        radl     = k2collimator.material.radiation_length
        dlri     = 0
        dlyi     = 0
        ai       = 0
        eUm      = 0
        collnt   = 0
        is_crystal = False

    if is_crystal:

        crytilt = k2collimator.align_angle + k2collimator.crytilt

        new_length = np.array(length)
        pyk2_startcry(
            c_length=length,
            new_length=new_length,
            c_rotation=k2collimator.angle/180.*np.pi,
            crytilt=crytilt,
            crybend=k2collimator.bend,
            crythick=k2collimator.thick,
            cryxdim=k2collimator.xdim,
            cryydim=k2collimator.ydim,
            cryorient=k2collimator.orient,
            crymiscut=k2collimator.miscut
        )
        length = new_length

    pyk2_run(x_particles=x_part,
              xp_particles=xp_part,
              y_particles=y_part,
              yp_particles=yp_part,
              s_particles=s_part,
              p_particles=e_part,              # confusing: this is ENERGY not momentum
              part_hit=part_hit,
              part_abs=part_abs,
              part_impact=part_impact,         # impact parameter
              part_indiv=part_indiv,           # particle divergence
              part_linteract=part_linteract,   # interaction length
              nhit_stage=nhit_stage,
              nabs_type=nabs_type,
              linside=linside,
              run_exenergy=exenergy,
              run_anuc=anuc,
              run_zatom=zatom,
              run_emr=emr,
              run_rho=rho,
              run_hcut=hcut,
              run_bnref=bnref,
              run_csref0=csref0,
              run_csref1=csref1,
              run_csref5=csref5,
              run_radl=radl,
              run_dlri=dlri,
              run_dlyi=dlyi,
              run_eum=eUm,
              run_ai=ai,
              run_collnt=collnt,
              is_crystal=is_crystal,
              c_length=length,
              c_rotation=k2collimator.angle/180.*np.pi,
              c_aperture=opening,
              c_offset=offset,
              c_tilt=k2collimator.tilt,
              c_enom=particles.energy0[0]/1e6, # Reference energy
              onesided=k2collimator.onesided
              )

    # Masks of hit and survived particles
    lost = part_abs > 0
    hit = part_hit > 0
    not_hit = ~hit
    not_lost = ~lost
    survived_hit = hit & (~lost)

    # Update energy    ---------------------------------------------------
    # Only particles that hit the jaw and survived need to be updated
    ptau_out = particles.ptau[:npart].copy()
    e0 = particles.energy0[:npart].copy()
    beta0 = particles.beta0[:npart].copy()
    ptau_out[survived_hit] = (
            e_part[survived_hit] * 1e9 - e0[survived_hit]
        ) / (
            e0[survived_hit] * beta0[survived_hit]
        )
    particles.ptau[:npart] = ptau_out
    rpp_out = particles.rpp[:npart].copy()

    # Rescale angles, because K2 did not update angles with new energy!
    # So need to do xp' = xp * p_in / p_out = xp * rpp_out / rpp_in
    # (see collimation.f90 line 1709 and mod_particles.f90 line 210)
    xp_part *= rpp_out/rpp_in
    yp_part *= rpp_out/rpp_in

    # Return from closed orbit
    x_part  += k2collimator.dx
    y_part  += k2collimator.dy

    # Update 4D coordinates    -------------------------------------------
    # Absorbed particles get their coordinates set to the entrance of collimator
    x_out  = particles.x[:npart].copy()
    y_out  = particles.y[:npart].copy()
    px_out = particles.px[:npart].copy()
    py_out = particles.py[:npart].copy()
    # Survived particles get updated coordinates
    x_out[not_lost]  = x_part[not_lost]
    y_out[not_lost]  = y_part[not_lost]
    px_out[not_lost] = xp_part[not_lost]/rpp_out[not_lost]
    py_out[not_lost] = yp_part[not_lost]/rpp_out[not_lost]
    # Write updated coordinates
    particles.x[:npart] = x_out
    particles.y[:npart] = y_out
    particles.px[:npart] = px_out
    particles.py[:npart] = py_out

    # Update longitudinal coordinate zeta    -----------------------------
    # Absorbed particles get coordinates set to the entrance of collimator
    rvv_out = particles.rvv[:npart].copy()
    zeta_out = particles.zeta[:npart].copy()
    # Non-hit particles are just drifting (zeta not yet drifted in K2, so do here)
    zeta_out[not_hit] = drift_zeta(zeta_out[not_hit], rvv_in[not_hit], xp_in[not_hit], yp_in[not_hit], length)
    # Hit and survived particles need correcting:
    # First we drift half the length with the old angles, then half the length with the new angles
    zeta_out[survived_hit] = drift_zeta(zeta_out[survived_hit], rvv_in[survived_hit],
                                        xp_in[survived_hit], yp_in[survived_hit], length/2)
    zeta_out[survived_hit] = drift_zeta(zeta_out[survived_hit], rvv_out[survived_hit],
                                        xp_part[survived_hit], yp_part[survived_hit], length/2)
    # Write updated coordinates
    particles.zeta[:npart] = zeta_out

    # Update s    --------------------------------------------------------
    s_out = particles.s[:npart]
    s_out[not_lost] += length
    particles.s[:npart] = s_out

    # Update state    ----------------------------------------------------
    state_out = particles.state[:npart].copy()
    state_out[lost] = -333
    particles.state[:npart] = state_out

    # TODO update records

    # =================================================================== #


#             print(lost[:40])
#             print(hit[:40])
#             print(survived_hit[:40])
#             print(part_impact[:40])
#             print(part_indiv[:40])
#             print(part_linteract[:40])  #  This is how much of the collimator it traversed -- correct? Or only what was left?
#             print(nabs_type[:40])
#             1:Nuclear-Inelastic,2:Nuclear-Elastic,3:pp-Elastic, 4:Single-Diffractive,5:Coulomb
