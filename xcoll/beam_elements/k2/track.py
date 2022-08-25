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
    return x, y

def drift_zeta(zeta, rvv, xp, yp, length):
    rv0v = 1./rvv
    dzeta = 1 - rv0v * (1 + (xp**2 + yp**2)/2 )
    zeta += length * dzeta
    return zeta


    
def track(k2collimator, particles, npart, reset_seed, is_crystal=False):
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")
    from .scatter_init import calculate_scattering
    from .materials import materials
    from .scatter import scatter

    length = k2collimator.active_length

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
    x_part  -= k2collimator.dx
    y_part  -= k2collimator.dy
    xp_part -= k2collimator.dpx
    yp_part -= k2collimator.dpy

    # Backtrack to start of collimator
    drift_4d(x_part, y_part, xp_part, yp_part, -length/2)

    # Initialise arrays for FORTRAN call
    part_hit       = np.zeros(len(x_part), dtype=np.int32)
    part_abs       = np.zeros(len(x_part), dtype=np.int32)
    part_impact    = np.zeros(len(x_part), dtype=np.float64)
    part_indiv     = np.zeros(len(x_part), dtype=np.float64)
    part_linteract = np.zeros(len(x_part), dtype=np.float64)
    nhit_stage     = np.zeros(len(x_part), dtype=np.int32)
    nabs_type      = np.zeros(len(x_part), dtype=np.int32)
    linside        = np.zeros(len(x_part), dtype=np.int32)

    if k2collimator.jaw_F_L != k2collimator.jaw_B_L or k2collimator.jaw_F_R != k2collimator.jaw_B_R:
        raise NotImplementedError
    opening = k2collimator.jaw_F_L - k2collimator.jaw_F_R
    offset = k2collimator.offset + ( k2collimator.jaw_F_L + k2collimator.jaw_F_R )/2

    # Get material properties
    exenergy = materials[k2collimator.material]['exenergy']
    anuc     = materials[k2collimator.material]['anuc']
    zatom    = materials[k2collimator.material]['zatom']
    emr      = materials[k2collimator.material]['emr']
    rho      = materials[k2collimator.material]['rho']
    hcut     = materials[k2collimator.material]['hcut']
    bnref    = materials[k2collimator.material]['bnref']
    csref0   = materials[k2collimator.material]['csref'][0]
    csref1   = materials[k2collimator.material]['csref'][1]
    csref4   = materials[k2collimator.material]['csref'][4]
    csref5   = materials[k2collimator.material]['csref'][5]
    radl     = materials[k2collimator.material]['radl']
    dlri     = materials[k2collimator.material]['dlri']
    dlyi     = materials[k2collimator.material]['dlyi']
    eUm      = materials[k2collimator.material]['eUm']
    ai       = materials[k2collimator.material]['ai']
    collnt   = materials[k2collimator.material]['collnt']
    nmat     = materials[k2collimator.material]['ID']

    csref=[csref0,csref1,0,0,csref4,csref5]

    # if self.is_crystal and not pyk2.materials[self.material]['can_be_crystal']:
    #  raise ValueError()
    cprob, xintl, bn, ecmsq, xln15s, bpp, csect = calculate_scattering(e0_ref,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref)

    if is_crystal:
        if not materials[k2collimator.material]['can_be_crystal']:
            raise ValueError(f"The collimator material {k2collimator.material} cannot be used as a crystal!")

        crytilt = k2collimator.align_angle + k2collimator.crytilt

        new_length = np.array(length)
        pyk2.pyk2_startcry(
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

    scatter(npart=npart,
            x_part=x_part,
            xp_part=xp_part,
            y_part=y_part,
            yp_part=yp_part,
            s_part=s_part,
            p_part=e_part,              # confusing: this is ENERGY not momentum
            part_hit=part_hit,
            part_abs=part_abs,
            part_impact=part_impact,         # impact parameter
            part_indiv=part_indiv,           # particle divergence
            part_linteract=part_linteract,   # interaction length
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
            run_csref4=csref4,
            run_csref5=csref5,
            run_radl=radl,
            run_dlri=dlri,
            run_dlyi=dlyi,
            run_eum=eUm,
            run_ai=ai,
            run_collnt=collnt,
            run_cprob=cprob,
            run_xintl=xintl,
            run_bn=bn,
            run_ecmsq=ecmsq,
            run_xln15s=xln15s,
            run_bpp=bpp,
            is_crystal=is_crystal,
            c_length=length,
            c_rotation=k2collimator.angle/180.*np.pi,
            c_aperture=opening,
            c_offset=offset,
            c_tilt=k2collimator.tilt,
            c_enom=e0_ref, # Reference energy in MeV
            onesided=k2collimator.onesided,
            length=length,
            material=k2collimator.material, 
            matid=materials[k2collimator.material]['ID'],
            nmat=nmat,
            run_csect=csect
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
    x_part  += k2collimator.dx
    y_part  += k2collimator.dy
    xp_part += k2collimator.dpx
    yp_part += k2collimator.dpy

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


