import numpy as np

def jaw(*, run_exenergy, run_anuc, run_zatom, run_rho, run_radl, run_cprob, run_xintl, run_bn, run_ecmsq, run_xln15s, run_bpp, cgen, p0, nabs, s, zlm, x, xp, z, zp, dpop):
    
    from .k2_random import get_random

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

        run_zlm1 = (-1*run_xintl)*np.log(get_random())
                        
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
            ##################################################################
            s, x, xp, z, zp, dpop = mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop)
            ##################################################################

            s = (zlm-rlen)+s

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            ###########################################################################################
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)  # DM routine to include tail
            ###########################################################################################
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
        ##################################################################
        s, x, xp, z, zp, dpop = mcs(s,run_radl,run_zlm1,p0,x,xp,z,zp,dpop)
        ##################################################################
        # Check if particle is outside of collimator (X.LT.0) after
        # MCS. If yes, calculate output longitudinal position (s),
        # reduce momentum (output as dpop) and return.
        # PARTICLE LEFT COLLIMATOR BEFORE ITS END.

        if(x <= 0):

            s = (zlm-rlen)+s

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            ###########################################################################################
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)
            ###########################################################################################

            p = p-m_dpodx*s
            dpop = (p-p0)/p0
            break

        # Check whether particle is absorbed. If yes, calculate output
        # longitudinal position (s), reduce momentum (output as dpop)
        # and return.
        # PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

        ###################################
        inter = ichoix(run_cprob)
        ###################################

        nabs = inter

        if(inter == 1):

            s = (zlm-rlen)+run_zlm1

            m_dpodx = np.array(m_dpodx, dtype=np.float64)
            run_exenergy = np.array(run_exenergy, dtype=np.float64)
            ###########################################################################################
            run_exenergy, m_dpodx = calcionloss(p,rlen,run_exenergy,run_anuc,run_zatom,run_rho,m_dpodx)
            ###########################################################################################

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
        ################################################################
        t, p = gettran(inter,p,run_bn,cgen,run_ecmsq,run_xln15s,run_bpp)
        ################################################################

        # Tetat calculates from the rms transverse momentum transfer in
        # monte-carlo fashion the angle changes for x and z planes. The
        # angle change is proportional to SQRT(t) and 1/p, as expected.

        tx = np.array(tx, dtype=np.float64)
        tz = np.array(tz, dtype=np.float64)
        ##########################
        tx, tz = tetat(t,p,tx,tz)
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


   
def calcionloss(p,rlen,il_exenergy,il_anuc,il_zatom,il_rho,enlo):

    from .k2_random import get_random

    mom    = p*1.0e3                     #[GeV/c] -> [MeV/c]
    enr    = (mom*mom + 938.271998*938.271998)**0.5 #[MeV]
    gammar = enr/938.271998
    betar  = mom/enr
    bgr    = betar*gammar
    kine   = ((2*0.510998902)*bgr)*bgr
    k = 0.307075

    # Mean excitation energy
    exEn = il_exenergy*1.0e3 # [MeV]

    # tmax is max energy loss from kinematics
    tmax = kine/(1 + (2*gammar)*(0.510998902/938.271998) + (0.510998902/938.271998)**2) # [MeV]

    # Plasma energy - see PDG 2010 table 27.1
    plen = (((il_rho*il_zatom)/il_anuc)**0.5)*28.816e-6 # [MeV]

    # Calculate threshold energy
    # Above this threshold, the cross section for high energy loss is calculated and then
    # a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
    # below threshold, only the standard Bethe-Bloch is used (all particles get average energy loss)

    # thl is 2*width of Landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
    thl = ((((4*(k*il_zatom))*rlen)*1.0e2)*il_rho)/(il_anuc*betar**2) # [MeV]

    # Bethe-Bloch mean energy loss
    enlo = ((k*il_zatom)/(il_anuc*betar**2)) * (0.5*np.log((kine*tmax)/(exEn*exEn)) - betar**2 - np.log(plen/exEn) - np.log(bgr) + 0.5)
    enlo = ((enlo*il_rho)*1.0e-1)*rlen # [GeV]

    # Threshold Tt is Bethe-Bloch + 2*width of Landau distribution
    Tt = enlo*1.0e3 + thl # [MeV]

    # Cross section - see Alfredo's presentation for derivation
    cs_tail = ((k*il_zatom)/(il_anuc*betar**2)) * (0.5*((1/Tt)-(1/tmax)) - (np.log(tmax/Tt)*betar**2)/(2*tmax) + (tmax-Tt)/((4*gammar**2)*938.271998**2))

    # Probability of being in tail: cross section * density * path length
    prob_tail = ((cs_tail*il_rho)*rlen)*1.0e2

    # Determine based on random number if tail energy loss occurs.
    if (get_random() < prob_tail):
        enlo = ((k*il_zatom)/(il_anuc*betar**2)) * (0.5*np.log((kine*tmax)/(exEn*exEn)) - betar**2 - np.log(plen/exEn) - np.log(bgr) + 0.5 + tmax**2/((8*gammar**2)*938.271998**2))
        enlo = (enlo*il_rho)*1.0e-1 # [GeV/m]
    else:
        # If tail energy loss does not occur, just use the standard Bethe-Bloch
        enlo = enlo/rlen  # [GeV/m]

    return il_exenergy, enlo




def gettran(inter,p,tt_bn,tt_cgen,tt_ecmsq,tt_xln15s,tt_bpp):

    from .k2_random import get_random, get_random_ruth

    # Neither if-statements below have an else, so defaulting function return to zero.
    result = 0

    if (inter==2): # Nuclear Elastic
        result = (-1*np.log(get_random()))/tt_bn
    
    elif (inter==3): # pp Elastic
        result = (-1*np.log(get_random()))/tt_bpp

    elif (inter==4): # Single Diffractive
        xm2 = np.exp(get_random() * tt_xln15s)
        p   = p * (1 - xm2/tt_ecmsq)
    
        if (xm2 < 2):
            bsd = 2 * tt_bpp

        elif ((xm2 >= 2) & (xm2 <= 5)):
            bsd = ((106.0 - 17.0*xm2)*tt_bpp)/36.0

        else:
            bsd = (7*tt_bpp)/12.0
   
        result = (-1*np.log(get_random()))/bsd

    elif (inter==5): # Coulomb
        result = get_random_ruth(tt_cgen)

    return result, p



def tetat(t,p,tx,tz):

    from .k2_random import get_random

    teta = np.sqrt(t)/p
    while (True):
        va  = 2*get_random() - 1
        vb  = get_random()
        va2 = va**2
        vb2 = vb**2
        r2  = va2 + vb2

        if(r2 < 1):
            break
        
    tx  = (teta*((2*va)*vb))/r2
    tz  = (teta*(va2 - vb2))/r2

    return tx, tz



def mcs(s, mc_radl, mc_zlm1, mc_p0, mc_x, mc_xp, mc_z, mc_zp, mc_dpop):

    theta    = 13.6e-3/(mc_p0*(1+mc_dpop)) # dpop   = (p - p0)/p0
    h   = 0.001
    dh  = 0.0001
    bn0 = 0.4330127019

    mc_x     = (mc_x/theta)/mc_radl
    mc_xp    = mc_xp/theta
    mc_z     = (mc_z/theta)/mc_radl
    mc_zp    = mc_zp/theta
    rlen0 = mc_zlm1/mc_radl
    rlen  = rlen0

    while (True): #10
        ae = bn0*mc_x
        be = bn0*mc_xp

        # #######################################
        # ae = np.array(ae, dtype=np.float64)
        # be = np.array(be, dtype=np.float64)
        # dh = np.array(dh, dtype=np.float64)
        # rlen = np.array(rlen, dtype=np.float64)
        # s = np.array(s, dtype=np.float64)
        # #######################################
        s = soln3(ae,be,dh,rlen,s)

        if (s < h):
            s = h

        mc_x, mc_xp = scamcs(mc_x,mc_xp,s)

        if (mc_x <= 0):
            s = (rlen0-rlen)+s
            break # go to 20

        if ((s+dh) >= rlen):
            s = rlen0
            break # go to 20
        # go to 10
        rlen = rlen-s

    mc_z, mc_zp = scamcs(mc_z,mc_zp,s)

    s  = s*mc_radl
    mc_x  = (mc_x*theta)*mc_radl
    mc_xp = mc_xp*theta
    mc_z  = (mc_z*theta)*mc_radl
    mc_zp = mc_zp*theta

    return s, mc_x, mc_xp, mc_z, mc_zp, mc_dpop



def scamcs(xx, xxp, s):

    from .k2_random import get_random

    x0  = xx
    xp0 = xxp

    while (True):
        v1 = 2*get_random() - 1
        v2 = 2*get_random() - 1
        r2 = v1**2 + v2**2

        if(r2 < 1):
            break

    a   = np.sqrt((-2*np.log(r2))/r2)
    z1  = v1*a
    z2  = v2*a
    ss  = np.sqrt(s)
    sss = 1 + 0.038*np.log(s)
    xx  = x0  + s*(xp0 + ((0.5*ss)*sss)*(z2 + z1*0.577350269))
    xxp = xp0 + (ss*z2)*sss

    return xx, xxp



def soln3(a, b, dh, smax, s):

    if(b == 0):
        s = a**0.6666666666666667
        # s = a**(two/three)

        if (s > smax): 
            s = smax
        return s

    if(a == 0):     
        if(b > 0):
            s = b**2
        else:
            s = 0
        
        if (s > smax):
            s = smax
        return s
        
    if (b > 0):

        if (smax**3 <= (a + b*smax)**2):
            s = smax
            return s
      
        else:
            s = smax*0.5
            s = iterat(a,b,dh,s)
    
    else:
        c = (-1*a)/b
        if (smax < c):
            if (smax**3 <= (a + b*smax)**2):
                s = smax
                return s

            else:
                s = smax*0.5
                s = iterat(a,b,dh,s)
    
        else:
            s = c*0.5
            s = iterat(a,b,dh,s)
   
    return s



def iterat(a, b, dh, s):

    ds = s

    while (True):
        ds = ds*0.5

        if (s**3 < (a+b*s)**2):
            s = s+ds
        else:
            s = s-ds

        if (ds < dh):
            break

        else: 
            continue

    return s 



def ichoix(ich_cprob):

    from .k2_random import get_random

    aran = get_random()
    
    for i in range(5):
        i += 1
        if(aran <= ich_cprob[i]):
            break
    return i

