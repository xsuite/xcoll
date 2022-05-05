import numpy as np

def scatter(*, npart, x_part, xp_part, y_part, yp_part, s_part, p_part, part_hit,
                part_abs, part_impact, part_indiv, part_linteract, nabs_type,
                linside, run_exenergy, run_anuc, run_zatom, run_emr, run_rho,  run_hcut, run_bnref, run_csref0, run_csref1, run_csref4, run_csref5,
                run_radl, run_dlri, run_dlyi, run_eum, run_ai, run_collnt, run_cprob, run_xintl, run_bn, run_ecmsq, run_xln15s, run_bpp, is_crystal, 
                c_length, c_rotation, c_aperture, c_offset, c_tilt, c_enom, onesided, random_generator_seed, length):
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

    cgen = np.zeros(200, dtype=np.float64)
    pyk2.initialise_random(random_generator_seed=random_generator_seed, cgen=cgen, zatom=run_zatom, emr=run_emr, hcut=run_hcut)

    # Initilaisation
    length  = c_length
    p0 = c_enom


    # Initialise scattering processes
    # call k2coll_scatin(p0,coll_anuc,coll_rho,coll_zatom,coll_emr,&
    #                   coll_csref0,coll_csref1,coll_csref5,coll_bnref,cprob,xintl,bn)

    nhit   = 0
    nabs   = 0
    fracab = 0
    mirror = 1

    # Compute rotation factors for collimator rotation
    cRot   = np.cos(c_rotation)
    sRot   = np.sin(c_rotation)
    cRRot  = np.cos(-c_rotation)
    sRRot  = np.sin(-c_rotation)

    # Set energy and nucleon change variables as with the coupling
    # ien0,ien1: particle energy entering/leaving the collimator
    # energy in MeV
    nnuc0 = 0
    ien0  = 0
    nnuc1 = 0
    ien1  = 0

    for i in range(npart):

        x_in = np.array(x_part[i])
        xp_in = np.array(xp_part[i])
        y_in = np.array(y_part[i])
        yp_in = np.array(yp_part[i])
        s_in = np.array(s_part[i])
        p_in = np.array(p_part[i])
        val_part_hit = np.array(part_hit[i])
        val_part_abs = np.array(part_abs[i])
        val_part_impact = np.array(part_impact[i])
        val_part_indiv = np.array(part_indiv[i])
        val_part_linteract = np.array(part_linteract[i])
        val_nabs_type = np.array(nabs_type[i])
        val_linside = np.array(linside[i])

        pyk2.pyk2_run(
                    x_in=x_in,
                    xp_in=xp_in,
                    y_in=y_in,
                    yp_in=yp_in,
                    s_in=s_in,
                    p_in=p_in,
                    val_part_hit=val_part_hit,
                    val_part_abs=val_part_abs,
                    val_part_impact=val_part_impact,
                    val_part_indiv=val_part_indiv,
                    val_part_linteract=val_part_linteract,
                    val_nabs_type=val_nabs_type,
                    val_linside=val_linside,
                    run_exenergy=run_exenergy,
                    run_anuc=run_anuc,
                    run_zatom=run_zatom,
                    run_emr=run_emr,
                    run_rho=run_rho,
                    run_hcut=run_hcut,
                    run_bnref=run_bnref,
                    run_csref0=run_csref0,
                    run_csref1=run_csref1,
                    run_csref4=run_csref4,
                    run_csref5=run_csref5,
                    run_radl=run_radl,
                    run_dlri=run_dlri,
                    run_dlyi=run_dlyi,
                    run_eum=run_eum,
                    run_ai=run_ai,
                    run_collnt=run_collnt,
                    run_cprob=run_cprob,
                    run_xintl=run_xintl,
                    run_bn=run_bn,
                    run_ecmsq=run_ecmsq,
                    run_xln15s=run_xln15s,
                    run_bpp=run_bpp,
                    run_cgen=cgen,
                    is_crystal=is_crystal,
                    c_length=c_length,
                    c_aperture=c_aperture,
                    c_offset=c_offset,
                    c_tilt=c_tilt,
                    onesided=onesided,
                    length=length,
                    p0=p0,
                    nhit=nhit,
                    nabs=nabs,
                    fracab=fracab,
                    mirror=mirror,
                    crot=cRot,
                    srot=sRot,
                    crrot=cRRot,
                    srrot=sRRot,
                    nnuc0=nnuc0,
                    nnuc1=nnuc1,
                    ien0=ien0,
                    ien1=ien1
                    )

        x_part[i] = x_in
        xp_part[i] = xp_in
        y_part[i] = y_in
        yp_part[i] = yp_in
        s_part[i] = s_in
        p_part[i] = p_in
        part_hit[i] = val_part_hit
        part_abs[i] = val_part_abs
        part_impact[i] = val_part_impact
        part_indiv[i] = val_part_indiv
        part_linteract[i] = val_part_linteract
        nabs_type[i] = val_nabs_type
        linside[i] = val_linside
