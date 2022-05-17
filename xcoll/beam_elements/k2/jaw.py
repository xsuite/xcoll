import numpy as np

def jaw(*, run_exenergy, run_anuc, run_zatom, run_rho, run_radl, run_cprob, run_xintl, run_bn, run_ecmsq, run_xln15s, run_bpp, cgen, p0, nabs, s, zlm, x, xp, z, zp, dpop):
    
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

    # Note that the input parameter is dpop. Here the momentum p is constructed out of this input.
    p    = p0*(1+dpop)
    nabs = 0
      
    # Initialize the interaction length to input interaction length
    rlen = zlm
    m_dpodx = 0.
    tx = 0.
    tz = 0.

    # Do a step for a point-like interaction.
    # Get monte-carlo interaction length.
    while (True):

        run_zlm1 = (-1*run_xintl)*np.log(pyk2.pyk2_rand())
                        
        # If the monte-carlo interaction length is longer than the
        # remaining collimator length, then put it to the remaining
        # length, do multiple coulomb scattering and return.
        # LAST STEP IN ITERATION LOOP
        if(run_zlm1 > rlen):
            
            run_zlm1 = rlen
            
            s = np.array(s, dtype=np.float64)
            p0 = np.array(p0, dtype=np.float64)
            x = np.array(x, dtype=np.float64)
            xp = np.array(xp, dtype=np.float64)
            z = np.array(z, dtype=np.float64)
            zp = np.array(zp, dtype=np.float64)
            dpop = np.array(dpop, dtype=np.float64)
            ####################################################
            pyk2.pyk2_mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop)
            ####################################################

            s = (zlm-rlen)+s

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            #############################################################################
            pyk2.pyk2_calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)  # DM routine to include tail
            #############################################################################
            p = p-m_dpodx*s
                    
            dpop = (p-p0)/p0
            break
        # Otherwise do multi-coulomb scattering.
        # REGULAR STEP IN ITERATION LOOP

        s = np.array(s, dtype=np.float64)
        p0 = np.array(p0, dtype=np.float64)
        x = np.array(x, dtype=np.float64)
        xp = np.array(xp, dtype=np.float64)
        z = np.array(z, dtype=np.float64)
        zp = np.array(zp, dtype=np.float64)
        dpop = np.array(dpop, dtype=np.float64)
        ####################################################
        pyk2.pyk2_mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop)
        ####################################################
        # Check if particle is outside of collimator (X.LT.0) after
        # MCS. If yes, calculate output longitudinal position (s),
        # reduce momentum (output as dpop) and return.
        # PARTICLE LEFT COLLIMATOR BEFORE ITS END.

        if(x <= 0):

            s = (zlm-rlen)+s

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            #############################################################################
            pyk2.pyk2_calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)
            #############################################################################

            p = p-m_dpodx*s
            dpop = (p-p0)/p0
            break

        # Check whether particle is absorbed. If yes, calculate output
        # longitudinal position (s), reduce momentum (output as dpop)
        # and return.
        # PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

        ######################################
        inter    = pyk2.pyk2_ichoix(run_cprob)
        ######################################

        nabs     = inter

        if(inter == 1):

            s = (zlm-rlen)+run_zlm1

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            #############################################################################
            pyk2.pyk2_calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)
            #############################################################################

            p = p-m_dpodx*s
            dpop = (p-p0)/p0

            break


        # Now treat the other types of interaction, as determined by ICHOIX:

        # Nuclear-Elastic:          inter = 2
        # pp Elastic:               inter = 3
        # Single-Diffractive:       inter = 4    (changes momentum p)
        # Coulomb:                  inter = 5

        # Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.

        p = np.array(p, dtype=np.float64)
        ######################################################################
        t = pyk2.pyk2_gettran(inter,p,run_bn,cgen,run_ecmsq,run_xln15s,run_bpp)
        ######################################################################

        # Tetat calculates from the rms transverse momentum transfer in
        # monte-carlo fashion the angle changes for x and z planes. The
        # angle change is proportional to SQRT(t) and 1/p, as expected.

        tx = np.array(tx, dtype=np.float64)
        tz = np.array(tz, dtype=np.float64)
        ##########################
        pyk2.pyk2_tetat(t,p,tx,tz)
        ##########################

        # Apply angle changes
        xp = xp+tx
        zp = zp+tz

        # Treat single-diffractive scattering.
        if(inter == 4):
            # added update for s
            s    = (zlm-rlen)+run_zlm1

            # Add this code to get the momentum transfer also in the calling routine
            dpop = (p-p0)/p0

        # Calculate the remaining interaction length and close the iteration loop.
        rlen = rlen-run_zlm1
                  
    return run_exenergy, run_bn, p0, nabs, s, zlm, x, xp, z, zp, dpop