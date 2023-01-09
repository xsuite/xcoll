import numpy as np
from ._everest import lib

from .scatter import scatter
from .random import set_rutherford_parameters

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

def track(collimator, particles):  # TODO: write impacts
        from ...beam_elements import Collimator, Crystal
        if not isinstance(collimator, Collimator) and not isinstance(collimator, Crystal):
            raise ValueError("Collimator is neither a Collimator nor a Crystal!\nCannot use Everest to track.")
        if particles._num_active_particles == 0:
            return
        
        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        if not collimator.is_active:
            # Drift full length
            drift_6d(particles, collimator.length)
        else:
            drift_6d(particles, collimator.inactive_front)
            track_core(collimator, particles)
            drift_6d(particles, collimator.inactive_back)
            particles.reorganize()
        return



def track_core(collimator, particles):
    from .scatter_init import calculate_scattering
    from .scatter import scatter
    from ...beam_elements.everest_collimator import Crystal
    from .materials import CrystalMaterial

    npart = particles._num_active_particles
    length = collimator.active_length

    # Get coordinates
    x_part  = particles.x[:npart].copy()
    y_part  = particles.y[:npart].copy()
    xp_part = particles.px[:npart] * particles.rpp[:npart]
    yp_part = particles.py[:npart] * particles.rpp[:npart]
    s_part  = 0 * x_part
    e0_ref  = particles.energy0[0] / 1e9    # Reference energy in GeV
    e_part  = particles.energy[:npart].copy() / 1e9 # Energy in GeV
    rpp_in  = particles.rpp[:npart].copy()
    rvv_in  = particles.rvv[:npart].copy()
    xp_in   = xp_part.copy()
    yp_in   = yp_part.copy()

    # Drift to centre of collimator
    drift_4d(x_part, y_part, xp_part, yp_part, length/2)

    # Move to closed orbit  (dpx = dxp because ref. particle has delta = 0)
    x_part  -= collimator.dx
    y_part  -= collimator.dy
    xp_part -= collimator.dpx
    yp_part -= collimator.dpy

    # Backtrack to start of collimator
    drift_4d(x_part, y_part, xp_part, yp_part, -length/2)

    # Initialise arrays for FORTRAN call
    part_hit       = np.zeros(len(x_part), dtype=np.int32)
    part_abs       = np.zeros(len(x_part), dtype=np.int32)
    part_impact    = np.zeros(len(x_part), dtype=np.float64)
    part_indiv     = np.zeros(len(x_part), dtype=np.float64)
    part_linteract = np.zeros(len(x_part), dtype=np.float64)
    nabs_type      = np.zeros(len(x_part), dtype=np.int32)
    linside        = np.zeros(len(x_part), dtype=np.int32)

    if collimator.jaw_F_L != collimator.jaw_B_L or collimator.jaw_F_R != collimator.jaw_B_R:
        raise NotImplementedError
    opening = collimator.jaw_F_L - collimator.jaw_F_R
    offset = collimator.offset + ( collimator.jaw_F_L + collimator.jaw_F_R )/2

    # Get material properties
    zatom    = collimator.material.Z
    anuc     = collimator.material.A
    rho      = collimator.material.density
    exenergy = collimator.material.excitation_energy
    emr      = collimator.material.nuclear_radius
    bnref    = collimator.material.nuclear_elastic_slope
    csref0   = collimator.material.cross_section[0]
    csref1   = collimator.material.cross_section[1]
    csref5   = collimator.material.cross_section[5]
    hcut     = collimator.material.hcut

    # Get crystal parameters
    if isinstance(collimator, Crystal):
        if not isinstance(collimator.material, CrystalMaterial):
            raise ValueError(f"The collimator material {collimator.material.name} cannot be used as a crystal!")
        dlri     = collimator.material.crystal_radiation_length
        dlyi     = collimator.material.crystal_nuclear_length
        ai       = collimator.material.crystal_plane_distance
        eUm      = collimator.material.crystal_potential
        collnt   = collimator.material.nuclear_collision_length
        is_crystal = True
        radl     = 0
    else:
        radl     = collimator.material.radiation_length
        dlri     = 0
        dlyi     = 0
        ai       = 0
        eUm      = 0
        collnt   = 0
        is_crystal = False

    # Get crystal parameters
    if is_crystal:
        bend = collimator.bend    # bend is bending radius
        cry_tilt = collimator.align_angle + collimator.crytilt
        bend_ang  = length/bend    # temporary value
        if cry_tilt >= -bend_ang:
            length = bend*(np.sin(bend_ang + cry_tilt) - np.sin(cry_tilt))
        else:
            length = bend*(np.sin(bend_ang - cry_tilt) + np.sin(cry_tilt))
    
        cry_rcurv  = bend
        cry_bend = length/cry_rcurv  # final value (with corrected length)
        cry_alayer = collimator.thick
        cry_xmax   = collimator.xdim
        cry_ymax   = collimator.ydim
        cry_orient = collimator.orient
        cry_miscut = collimator.miscut
        cry_cBend  = np.cos(cry_bend)
        cry_sBend  = np.sin(cry_bend)
        cry_cpTilt = np.cos(cry_tilt)
        cry_spTilt = np.sin(cry_tilt)
        cry_cnTilt = np.cos(-cry_tilt)
        cry_snTilt = np.sin(-cry_tilt)

    else:
        cry_tilt = 0
        cry_rcurv  = 0
        cry_bend = 0
        cry_alayer = 0
        cry_xmax   = 0
        cry_ymax   = 0
        cry_orient = 0
        cry_miscut = 0
        cry_cBend  = 0
        cry_sBend  = 0
        cry_cpTilt = 0
        cry_spTilt = 0
        cry_cnTilt = 0
        cry_snTilt = 0

    cprob0, cprob1, cprob2, cprob3, cprob4, cprob5, xintl, bn, ecmsq, xln15s, bpp, csect = calculate_scattering(e0_ref,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref)

    set_rutherford_parameters(zatom=zatom, emr=emr, hcut=hcut)

    # Initilaisation
    p0 = c_enom
    x0 = 0
    xp0 = 0
    nhit   = 0
    nabs   = 0
    fracab = 0
    # Set energy and nucleon change variables as with the coupling
    # ien0, ien1: particle energy entering/leaving the collimator
    # energy in MeV
    nnuc0 = 0
    ien0  = 0
    nnuc1 = 0
    ien1  = 0
    # Crystal tracking parameters
    iProc       = 0
    n_chan      = 0
    n_VR        = 0
    n_amorphous = 0
    s_imp        = 0

    for i in range(npart):

        if (part_abs[i] != 0):
            continue

        x_in[i], xp_in[i], y_in[i], yp_in[i], s_in[i], e_part[i], part_hit[i], part_abs[i], part_impact[i], part_indiv[i],part_linteract[i], nabs_type[i], linside[i], p0, x0, xp0, nhit, nabs, fracab, nnuc0, ien0, nnuc1, ien1, iProc, n_chan, n_VR, n_amorphous, s_imp = scatter(
                x_in[i],
                xp_in[i],
                y_in[i],
                yp_in[i],
                s_in[i],
                e_part[i],                   # confusing: this is ENERGY not momentum
                part_hit[i],
                part_abs[i],
                part_impact[i],         # impact parameter
                part_indiv[i],           # particle divergence
                part_linteract[i],   # interaction length
                nabs_type[i],
                linside[i],
                exenergy,
                anuc,
                zatom,
                emr,
                rho,
                hcut,
                bnref,
                csref0,
                csref1,
                csref5,
                radl,
                dlri,
                dlyi,
                eUm,
                ai,
                collnt,
                cprob0,
                cprob1,
                cprob2,
                cprob3,
                cprob4,
                cprob5,
                xintl,
                bn,
                ecmsq,
                xln15s,
                bpp,
                is_crystal,
                length,
                collimator.angle/180.*np.pi,
                opening,
                offset,
                collimator.tilt[0],
                collimator.tilt[1],
                e0_ref, # Reference energy in MeV
                collimator.onesided,
                csect, 
                cry_tilt,
                cry_rcurv,
                cry_bend,
                cry_alayer,
                cry_xmax,
                cry_ymax,
                cry_orient,
                cry_miscut,
                cry_cBend,
                cry_sBend,
                cry_cpTilt,
                cry_spTilt,
                cry_cnTilt,
                cry_snTilt,
                p0,
                x0,
                xp0,
                nhit,
                nabs,
                fracab,
                nnuc0,
                ien0,
                nnuc1,
                ien1,
                iProc,
                n_chan,
                n_VR,
                n_amorphous,
                s_imp
                )

    # Masks of hit and survived particles
    lost = part_abs > 0
    hit = part_hit > 0
    not_hit = ~hit
    not_lost = ~lost
    survived_hit = hit & (~lost)

    # Backtrack to centre of collimator
    drift_4d(x_part, y_part, xp_part, yp_part, -length/2)

    # Return from closed orbit  (dpx = dxp because ref. particle has delta = 0)
    x_part  += collimator.dx
    y_part  += collimator.dy
    xp_part += collimator.dpx
    yp_part += collimator.dpy

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

    # Drift to end of collimator
    drift_4d(x_part, y_part, xp_part, yp_part, length/2)

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


