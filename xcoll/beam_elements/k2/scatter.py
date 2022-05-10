import numpy as np

def scatter(*, npart, x_part, xp_part, y_part, yp_part, s_part, p_part, part_hit,
                part_abs, part_impact, part_indiv, part_linteract, nabs_type,linside, run_exenergy, run_anuc, run_zatom,
                run_emr, run_rho,  run_hcut, run_bnref, run_csref0, run_csref1, run_csref4, run_csref5,run_radl, run_dlri, 
                run_dlyi, run_eum, run_ai, run_collnt, run_cprob, run_xintl, run_bn, run_ecmsq, run_xln15s, run_bpp, is_crystal, 
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

        val_part_hit = np.array(part_hit[i])
        val_part_abs = np.array(part_abs[i])
        val_part_impact = np.array(part_impact[i])
        val_part_indiv = np.array(part_indiv[i])
        val_part_linteract = np.array(part_linteract[i])
        val_nabs_type = np.array(nabs_type[i])
        val_linside = np.array(linside[i])

        if (val_part_abs != 0):
            continue

        x_in = x_part[i]
        xp_in = xp_part[i]
        y_in = y_part[i]
        yp_in = yp_part[i]
        s_in = s_part[i]
        p_in = p_part[i]
        
        val_part_impact = -1
        val_part_linteract = -1
        val_part_indiv = -1

        x = x_in
        xp = xp_in
        xp_in0 = xp_in
        z = y_in
        zp = yp_in
        p = p_in
        sp = 0
        dpop = (p - p0)/p0

        # Transform particle coordinates to get into collimator coordinate  system
        # First do rotation into collimator frame
        x  =  x_in*cRot + sRot*y_in
        z  =  y_in*cRot - sRot*x_in
        xp = xp_in*cRot + sRot*yp_in
        zp = yp_in*cRot - sRot*xp_in
        
        # For one-sided collimators consider only positive X. For negative X jump to the next particle
        if (onesided & (x < 0)):
            continue

        # Log input energy + nucleons as per the FLUKA coupling
        nnuc0   = nnuc0 + 1
        ien0    = ien0 + p_in * 1.0e3

        # Now mirror at the horizontal axis for negative X offset
        if (x < 0):
            mirror    = -1
            tiltangle = -1*c_tilt[1]
        else:
            mirror    = 1
            tiltangle = c_tilt[0]
    
        x  = mirror*x
        xp = mirror*xp

        # Shift with opening and offset
        x = (x - c_aperture/2) - mirror*c_offset

        # Include collimator tilt
        if(tiltangle > 0):
            xp = xp - tiltangle
        
        if(tiltangle < 0):
            x  = x + np.sin(tiltangle) * c_length
            xp = xp - tiltangle


        # CRY Only: x_in0 has to be assigned after the change of reference frame
        x_in0 = x

        # particle passing above the jaw are discarded => take new event
        # entering by the face, shorten the length (zlm) and keep track of
        # entrance longitudinal coordinate (keeps) for histograms

        # The definition is that the collimator jaw is at x>=0.

        # 1) Check whether particle hits the collimator
        isimp = False
        s     = 0
        zlm = -1*length


        if (is_crystal):

            ##########################################

            val_part_hit = np.array(val_part_hit, dtype=np.int64)
            val_part_abs = np.array(val_part_abs, dtype=np.int64)
            val_part_impact = np.array(val_part_impact, dtype=np.float64)
            val_part_indiv = np.array(val_part_indiv, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            run_bn = np.array(run_bn, dtype=np.float64)
            nhit = np.array(nhit, dtype=np.int64)
            nabs = np.array(nabs, dtype=np.int64)
            isimp = np.array(isimp)
            s = np.array(s, dtype=np.float64)
            zlm = np.array(zlm, dtype=np.float64)
            x = np.array(x, dtype=np.float64)
            xp = np.array(xp, dtype=np.float64)
            xp_in0 = np.array(xp_in0, dtype=np.float64)
            z = np.array(z, dtype=np.float64)
            zp = np.array(zp, dtype=np.float64)
            p = np.array(p, dtype=np.float64)
            x_in = np.array(x_in, dtype=np.float64)

            ##########################################

            pyk2.pyk2_crystal(
                    val_part_hit=val_part_hit,
                    val_part_abs=val_part_abs,
                    val_part_impact=val_part_impact,
                    val_part_indiv=val_part_indiv,
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
                    run_dlri=run_dlri,
                    run_dlyi=run_dlyi,
                    run_eum=run_eum,
                    run_ai=run_ai,
                    run_collnt=run_collnt,
                    run_bn=run_bn,
                    c_length=c_length,
                    nhit=nhit,
                    nabs=nabs,
                    isimp=isimp,
                    s=s,
                    zlm=zlm,
                    x=x,
                    xp=xp,
                    xp_in0=xp_in0,
                    z=z,
                    zp=zp,
                    p=p,
                    x_in0=x_in0
                    )

            if (nabs != 0):
                val_part_abs = 1
                val_part_linteract = zlm

            simp  = (s - c_length) + simp
            sout  = s
            xout  = x
            xpout = xp
            yout  = z
            ypout = zp

##################################################################################################

        else:

            if (x >= 0):
            # Particle hits collimator and we assume interaction length ZLM equal
            # to collimator length (what if it would leave collimator after
            # small length due to angle???)
                zlm  = length
                val_part_impact = x
                val_part_indiv  = xp
            
            elif (xp <= 0):
            # Particle does not hit collimator. Interaction length ZLM is zero.
                zlm = 0
            
            else:
            # Calculate s-coordinate of interaction point
                s = (-1*x)/xp
                # if(s <= 0):
                # write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
                # call prror
            
                if (s < length):
                    zlm = length - s
                    val_part_impact = 0
                    val_part_indiv  = xp
                else:
                    zlm = 0
           

            # First do the drift part
            # DRIFT PART
            drift_length = length - zlm
            if (drift_length > 0):
                x  = x  + xp * drift_length
                z  = z  + zp * drift_length
                sp = sp + drift_length

            # Now do the scattering part
            if (zlm > 0):
                if (not val_linside):
                # first time particle hits collimator: entering jaw
                    val_linside = True

                s_impact = sp
                nhit = nhit + 1

                ##########################################

               
                val_part_hit = np.array(val_part_hit, dtype=np.int64)
                val_part_abs = np.array(val_part_abs, dtype=np.int64)
                val_part_impact = np.array(val_part_impact, dtype=np.float64)
                val_part_indiv = np.array(val_part_indiv, dtype=np.float64)
                val_part_linteract = np.array(val_part_linteract, dtype=np.float64)
                val_nabs_type = np.array(val_nabs_type, dtype=np.int64)
                val_linside = np.array(val_linside)
                run_exenergy = np.array(run_exenergy, dtype=np.float64)
                run_bn = np.array(run_bn, dtype=np.float64)
                length = np.array(length, dtype=np.float64)
                p0 = np.array(p0, dtype=np.float64)
                nhit = np.array(nhit, dtype=np.int64)
                nabs = np.array(nabs, dtype=np.int64)
                fracab = np.array(fracab, dtype=np.float64)
                isimp = np.array(isimp)
                s = np.array(s, dtype=np.float64)
                zlm = np.array(zlm, dtype=np.float64)
                x = np.array(x, dtype=np.float64)
                xp = np.array(xp, dtype=np.float64)
                z = np.array(z, dtype=np.float64)
                zp = np.array(zp, dtype=np.float64)
                sp = np.array(sp, dtype=np.float64)
                dpop = np.array(dpop, dtype=np.float64)

                ##########################################

                pyk2.pyk2_jaw(
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
                    run_rho=run_rho,
                    run_radl=run_radl,
                    run_cprob=run_cprob,
                    run_xintl=run_xintl,
                    run_bn=run_bn,
                    run_ecmsq=run_ecmsq,
                    run_xln15s=run_xln15s,
                    run_bpp=run_bpp,
                    run_cgen=cgen,
                    length=length,
                    p0=p0,
                    nhit=nhit,
                    nabs=nabs,
                    fracab=fracab,
                    isimp=isimp,
                    s=s,
                    zlm=zlm,
                    x=x,
                    xp=xp,
                    z=z,
                    zp=zp,
                    sp=sp,
                    dpop=dpop
                    )

                val_nabs_type = nabs
                val_part_hit  = 1

                isimp = True
                simp  = s_impact
                sout  = (s+sp)
                xout  = x
                xpout = xp
                yout  = z
                ypout = zp

                # Writeout should be done for both inelastic and single diffractive. doing all transformations
                # in x_flk and making the set to 99.99 mm conditional for nabs=1
                if (nabs == 1 or nabs == 4):
                # Transform back to lab system for writeout.
                # keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

                # Finally, the actual coordinate change to 99 mm
                    if (nabs == 1):
                        fracab = fracab + 1
                        x = 99.99e-3
                        z = 99.99e-3
                        val_part_linteract = zlm
                        val_part_abs = 1
                    # Collimator jaw interaction

            if (nabs != 1 and zlm > 0):
            # Do the rest drift, if particle left collimator early
                drift_length = (length-(s+sp))

                if (drift_length > 1.0e-15):
                    val_linside = False
                    x  = x  + xp * drift_length
                    z  = z  + zp * drift_length
                    sp = sp + drift_length
            
                val_part_linteract = zlm - drift_length
            
####################################################################################################################

        # Transform back to particle coordinates with opening and offset
        if(x < 99.0e-3):
            # Include collimator tilt
            if(tiltangle > 0):
                x  = x  + tiltangle*c_length
                xp = xp + tiltangle
            elif(tiltangle < 0):
                x  = x  + tiltangle*c_length
                xp = xp + tiltangle
                x  = x  - np.sin(tiltangle) * c_length

            # Transform back to particle coordinates with opening and offset
            z00 = z
            x00 = x + mirror*c_offset
            x   = (x + c_aperture/2) + mirror*c_offset

            # Now mirror at the horizontal axis for negative X offset
            x  = mirror * x
            xp = mirror * xp

            # Last do rotation into collimator frame
            x_in  =  x*cRRot +  z*sRRot
            y_in  =  z*cRRot -  x*sRRot
            xp_in = xp*cRRot + zp*sRRot
            yp_in = zp*cRRot - xp*sRRot

            # Log output energy + nucleons as per the FLUKA coupling
            # Do not log dead particles
            nnuc1       = nnuc1 + 1                           # outcoming nucleons
            ien1        = ien1  + p_in * 1e3                 # outcoming energy

            if(is_crystal):
                p_in = p
                s_in = s_in + s
            else:
                p_in = (1 + dpop) * p0
                s_in = sp

        else:
            x_in = x
            y_in = z

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
