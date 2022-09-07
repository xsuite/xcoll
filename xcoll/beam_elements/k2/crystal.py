from calendar import c
from multiprocessing.resource_sharer import stop
import numpy as np

from xcoll.beam_elements.k2.k2_random import get_random_gauss

# Rutherford Scatter
tlcut_cry = 0.0009982
# cgen_cry(200,1)
# emr_curr_cry
# zatom_curr_cry

aTF = 0.194e-10 # Screening function [m]
dP  = 1.920e-10 # Distance between planes (110) [m]
u1  = 0.075e-10 # Thermal vibrations amplitude

# pp cross-sections and parameters for energy dependence
pptref_cry = 0.040
freeco_cry = 1.618

# Processes
proc_out         =  -1     # Crystal not hit
proc_AM          =   1     # Amorphous
proc_VR          =   2     # Volume reflection
proc_CH          =   3     # Channeling
proc_VC          =   4     # Volume capture
proc_absorbed    =   5     # Absorption
proc_DC          =   6     # Dechanneling
proc_pne         =   7     # Proton-neutron elastic interaction
proc_ppe         =   8     # Proton-proton elastic interaction
proc_diff        =   9     # Single diffractive
proc_ruth        =  10     # Rutherford scattering
proc_ch_absorbed =  15     # Channeling followed by absorption
proc_ch_pne      =  17     # Channeling followed by proton-neutron elastic interaction
proc_ch_ppe      =  18     # Channeling followed by proton-proton elastic interaction
proc_ch_diff     =  19     # Channeling followed by single diffractive
proc_ch_ruth     =  20     # Channeling followed by Rutherford scattering
proc_TRVR        = 100     # Volume reflection in VR-AM transition region
proc_TRAM        = 101     # Amorphous in VR-AM transition region

temp = 0

def crystal(*,x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isimp,val_part_hit,val_part_abs,val_part_impact,val_part_indiv,c_length,exenergy,rho,anuc,zatom,emr,
            dlri,dlyi,ai,eum,collnt,hcut,bnref,csref0,csref1,csref4,csref5,csect,nhit,nabs,
            cry_tilt, cry_rcurv, cry_bend, cry_alayer, cry_xmax, cry_ymax, cry_orient, cry_miscut, cry_cBend, 
            cry_sBend, cry_cpTilt, cry_spTilt, cry_cnTilt, cry_snTilt, iProc, n_chan, n_VR, n_amorphous):


    from .k2_random import get_random
    
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.") 

    val_part_hit = np.array(val_part_hit)
    val_part_abs = np.array(val_part_abs)
    val_part_impact = np.array(val_part_impact)
    val_part_indiv = np.array(val_part_indiv)

    exenergy = np.array(exenergy)
    anuc = np.array(anuc)
    zatom = np.array(zatom)
    emr = np.array(emr)
    rho = np.array(rho)
    hcut = np.array(hcut)
    bnref = np.array(bnref)
    s_imp = np.array(s_imp)

    csref0 = np.array(csref0)
    csref1 = np.array(csref1)
    csref4 = np.array(csref4)
    csref5 = np.array(csref5)

    csect = np.array(csect)
    dlri = np.array(dlri)
    dlyi = np.array(dlyi)
    eum = np.array(eum)
    ai = np.array(ai)
    collnt = np.array(collnt)
    c_length = np.array(c_length)

    nhit = np.array(nhit)
    nabs = np.array(nabs)
    isimp = np.array(isimp)
    s = np.array(s)
    zlm = np.array(zlm)
    x0 = np.array(x0)
    xp0 = np.array(xp0)

    x = np.array(x)
    xp = np.array(xp)
    z = np.array(z)
    zp = np.array(zp)
    p = np.array(p)
    csect = np.array(csect)

    iProc=np.array(iProc)
    n_chan=np.array(n_chan)
    n_VR=np.array(n_VR)
    n_amorphous=np.array(n_amorphous)

    s_temp     = 0
    s_shift    = 0
    s_rot      = 0
    s_int      = 0
    x_temp     = 0
    x_shift    = 0
    x_rot      = 0
    x_int      = 0
    xp_temp    = 0
    xp_shift   = 0
    xp_rot     = 0
    xp_int     = 0
    xp_tangent = 0
    tilt_int   = 0
    shift      = 0
    delta      = 0
    a_eq       = 0
    b_eq       = 0
    c_eq       = 0
    s_imp      = 0

    iProc       = proc_out
    cry_length  = c_length

    # Transform in the crystal reference system
    # 1st transformation: shift of the center of the reference frame
    if (cry_tilt < 0):
        s_shift = s
        shift   = cry_rcurv*(1 - cry_cpTilt)

        if (cry_tilt < -cry_bend):
            shift = cry_rcurv*(cry_cnTilt - np.cos(cry_bend - cry_tilt))

        x_shift = x - shift

    else:
        s_shift = s
        x_shift = x

  
    # 2nd transformation: rotation
    s_rot  = x_shift*cry_spTilt + s_shift*cry_cpTilt
    x_rot  = x_shift*cry_cpTilt - s_shift*cry_spTilt
    xp_rot = xp - cry_tilt

    # 3rd transformation: drift to the new coordinate s=0
    xp = xp_rot
    x  = x_rot - xp_rot*s_rot
    z  = z - zp*s_rot
    s  = 0
  
    # Check that particle hit the crystal
    if (x >= 0 and x < cry_xmax):
        # MISCUT first step: P coordinates (center of curvature of crystalline planes)
        s_P = (cry_rcurv-cry_xmax)*np.sin(-cry_miscut)
        x_P = cry_xmax + (cry_rcurv-cry_xmax)*np.cos(-cry_miscut)

        x=np.array(x)
        xp=np.array(xp)
        z=np.array(z)
        zp=np.array(zp)
        p=np.array(p)
        iProc=np.array(iProc)

        x,xp,z,zp,p,iProc = interact(x,xp,z,zp,p,cry_length,s_P,x_P,exenergy,rho,anuc,zatom,emr,dlri,dlyi,
                        ai,eum,collnt,hcut,csref0,csref1,csref4,csref5,bnref,csect,cry_tilt,
                        cry_rcurv,cry_alayer,cry_xmax,cry_ymax,cry_orient,cry_miscut,cry_bend,cry_cBend,
                        cry_sBend,cry_cpTilt,cry_spTilt,cry_cnTilt,cry_snTilt,iProc)

                
        s   = cry_rcurv*cry_sBend
        zlm = cry_rcurv*cry_sBend

        if (iProc != proc_out):
            isimp    = True
            nhit     = nhit + 1
            val_part_hit     = 1
            val_part_impact   = x0
            val_part_indiv    = xp0

    else:

        if (x < 0): # Crystal can be hit from below
            xp_tangent = np.sqrt((-(2*x)*cry_rcurv + x**2)/(cry_rcurv**2))

        else:              # Crystal can be hit from above
            xp_tangent = np.arcsin((cry_rcurv*(1 - cry_cBend) - x)/np.sqrt(((2*cry_rcurv)*(cry_rcurv - x))*(1 - cry_cBend) + x**2))

        # If the hit is below, the angle must be greater or equal than the tangent,
        # or if the hit is above, the angle must be smaller or equal than the tangent
        if ((x < 0 and xp >= xp_tangent) or (x >= 0 and xp <= xp_tangent)):

            # If it hits the crystal, calculate in which point and apply the transformation and drift to that point
            a_eq  = 1 + xp**2
            b_eq  = (2*xp)*(x - cry_rcurv)
            c_eq  = -(2*x)*cry_rcurv + x**2
            delta = b_eq**2 - 4*(a_eq*c_eq)
            s_int = (-b_eq - np.sqrt(delta))/(2*a_eq)
            s_imp = s_int

            # MISCUT first step: P coordinates (center of curvature of crystalline planes)
            s_P_tmp = (cry_rcurv-cry_xmax)*np.sin(-cry_miscut)
            x_P_tmp = cry_xmax + (cry_rcurv-cry_xmax)*np.cos(-cry_miscut)

            if (s_int < cry_rcurv*cry_sBend):
                # Transform to a new reference system: shift and rotate
                x_int  = xp*s_int + x
                xp_int = xp
                z      = z + zp*s_int
                x      = 0
                s      = 0

                tilt_int = s_int/cry_rcurv
                xp    = xp-tilt_int

                # MISCUT first step (bis): transform P in new reference system
                # Translation
                s_P_tmp = s_P_tmp - s_int
                x_P_tmp = x_P_tmp - x_int
                # Rotation
                s_P = s_P_tmp*np.cos(tilt_int) + x_P_tmp*np.sin(tilt_int)
                x_P = -s_P_tmp*np.sin(tilt_int) + x_P_tmp*np.cos(tilt_int)

                x=np.array(x)
                xp=np.array(xp)
                z=np.array(z)
                zp=np.array(zp)
                p=np.array(p)
                iProc=np.array(iProc)

                x,xp,z,zp,p,iProc = interact(x,xp,z,zp,p,cry_length-(tilt_int*cry_rcurv),s_P,x_P,exenergy,rho,anuc,
                                zatom,emr,dlri,dlyi,ai,eum,collnt,hcut,csref0,csref1,csref4,csref5,bnref,
                                csect,cry_tilt,cry_rcurv,cry_alayer,cry_xmax,cry_ymax,cry_orient,cry_miscut,
                                cry_bend,cry_cBend,cry_sBend,cry_cpTilt,cry_spTilt,cry_cnTilt,cry_snTilt,iProc)

                s   = cry_rcurv*np.sin(cry_bend - tilt_int)
                zlm = cry_rcurv*np.sin(cry_bend - tilt_int)
                
                if (iProc != proc_out):
                    x_rot    = x_int
                    s_rot    = s_int
                    xp_rot   = xp_int
                    s_shift  =  s_rot*cry_cnTilt + x_rot*cry_snTilt
                    x_shift  = -s_rot*cry_snTilt + x_rot*cry_cnTilt
                    xp_shift = xp_rot + cry_tilt

                    if (cry_tilt < 0):
                        x0  = x_shift + shift
                        xp0 = xp_shift

                    else:
                        x0  = x_shift
                        xp0 = xp_shift

                    isimp     = True
                    nhit      = nhit + 1
                    val_part_hit    = 1
                    val_part_impact = x0
                    val_part_indiv  = xp0
                    #

                # un-rotate
                x_temp  = x
                s_temp  = s
                xp_temp = xp
                s       =  s_temp*np.cos(-tilt_int) + x_temp*np.sin(-tilt_int)
                x    = -s_temp*np.sin(-tilt_int) + x_temp*np.cos(-tilt_int)
                xp   = xp_temp + tilt_int

                # 2nd: shift back the 2 axis
                x = x + x_int
                s = s + s_int

            else:
                s = cry_rcurv*np.sin(cry_length/cry_rcurv)
                x = x + s*xp
                z = z + s*zp

        else:
            s = cry_rcurv*np.sin(cry_length/cry_rcurv)
            x = x + s*xp
            z = z + s*zp

    # transform back from the crystal to the collimator reference system
    # 1st: un-rotate the coordinates
    x_rot  = x
    s_rot  = s
    xp_rot = xp

    s_shift  =  s_rot*cry_cnTilt + x_rot*cry_snTilt
    x_shift  = -s_rot*cry_snTilt + x_rot*cry_cnTilt
    xp_shift = xp_rot + cry_tilt

    # 2nd: shift back the reference frame
    if (cry_tilt < 0):
        s  = s_shift
        x  = x_shift + shift
        xp = xp_shift

    else:
        x  = x_shift
        s  = s_shift
        xp = xp_shift
    #

    # 3rd: shift to new S=Length position
    x = xp*(c_length - s) + x
    z = zp*(c_length - s) + z
    s = c_length

    nabs = 0
  
    if (iProc == proc_AM):
        n_amorphous = n_amorphous + 1

    elif (iProc == proc_VR):
        n_VR = n_VR + 1

    elif (iProc == proc_CH):
        n_chan = n_chan + 1

    elif (iProc == proc_absorbed):
        nabs = 1   # TODO: do we need to set part_abs_pos etc?

    elif (iProc == proc_ch_absorbed):
        nabs = 1
    #
    
    return val_part_hit, val_part_abs, val_part_impact, val_part_indiv, nhit, nabs, s_imp, isimp, s, zlm, x, xp, x0, xp0, z, zp, p, iProc, n_chan, n_VR, n_amorphous



def interact(x,xp,y,yp,pc,length,s_P,x_P,exenergy,rho,anuc,zatom,emr,dlri,dlyi,ai,eUm,collnt,hcut,csref0,csref1,csref4,
            csref5,bnref,csect,cry_tilt,cry_rcurv,cry_alayer,cry_xmax,cry_ymax,cry_orient,cry_miscut,cry_bend,cry_cBend,
            cry_sBend,cry_cpTilt,cry_spTilt,cry_cnTilt,cry_snTilt,iProc):

    from .k2_random import get_random, get_random_gauss
    
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

    dest = 0.
    pmap = 938.271998
    pmae = 0.51099890
    crade = 2.817940285e-15

    c_v1 =  1.7   # Fitting coefficient
    c_v2 = -1.5   # Fitting coefficient

    nam = 1 # Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
    zn  = 1

    # dE/dX and dechanneling length calculation
    mom    = pc*1.0e3                # [GeV]
    enr    = np.sqrt(mom**2 + pmap**2) # [MeV]
    gammar = enr/pmap
    betar  = mom/enr
    bgr    = betar*gammar
    mep    = pmae/pmap  # Electron/proton

    tmax = (2*pmae*bgr**2)/(1 + 2*gammar*mep + mep**2)  # [MeV]
    plen = np.sqrt((rho*zatom)/anuc)*28.816e-6 # [MeV]

    const_dech = ((256.0/(9*np.pi**2)) * (1/(np.log(((2*pmae)*gammar)/(exenergy*1.0e3)) - 1))) * ((aTF*dP)/(crade*pmae)) # [m/MeV]
    const_dech = const_dech*1.0e3 # [m/GeV]

    s        = 0
    s_length = cry_rcurv*np.sin(length/cry_rcurv)
    L_chan   = length

    # MISCUT second step: fundamental coordinates (crystal edges and plane curvature radius)
    s_K = cry_rcurv*np.sin(length/cry_rcurv)
    x_K = cry_rcurv*(1-np.cos(length/cry_rcurv))
    s_M = (cry_rcurv-cry_xmax)*np.sin(length/cry_rcurv)
    x_M = cry_xmax + (cry_rcurv-cry_xmax)*(1-np.cos(length/cry_rcurv))
    r   = np.sqrt(s_P**2 + (x-x_P)**2)

    # MISCUT third step: F coordinates (exit point) on crystal exit face
    A_F = (np.tan(length/cry_rcurv))**2 + 1
    B_F = ((-2)*(np.tan(length/cry_rcurv))**2)*cry_rcurv + (2*np.tan(length/cry_rcurv))*s_P - 2*x_P
    C_F = ((np.tan(length/cry_rcurv))**2)*(cry_rcurv**2) - ((2*np.tan(length/cry_rcurv))*s_P)*cry_rcurv + s_P**2 + x_P**2 - r**2

    x_F = (-B_F-np.sqrt(B_F**2-4*(A_F*C_F)))/(2*A_F)
    s_F = (-np.tan(length/cry_rcurv))*(x_F-cry_rcurv)

    if (x_F < x_K or x_F > x_M or s_F < s_M or s_F > s_K):
        
        if (cry_miscut == 0 and abs(x_F-x_K) <= 1.0e-13 and abs(s_F-s_K) <= 1.0e3):
        # no miscut, entrance from below: exit point is K (lower edge)
            x_F = x_K
            s_F = s_K

        elif (cry_miscut == 0 and abs(x_F-x_M) <= 1.0e3 and abs(s_F-s_M) <= 1.0e3):
        # no miscut, entrance from above: exit point is M (upper edge)
            x_F = x_M
            s_F = s_M

        else:
        # MISCUT Third step (bis): F coordinates (exit point)  on bent side
            if (cry_miscut < 0):
            # Intersect with bottom side
                alpha_F = (cry_rcurv-x_P)/x_P
                beta_F = -(r**2-s_P**2-x_P**2)/(2*s_P)
                A_F = alpha_F**2 + 1
                B_F = 2*(alpha_F*beta_F) - 2*cry_rcurv
                C_F = beta_F**2

            else:
            # Intersect with top side
                alpha_F = (cry_rcurv-x_P)/s_P
                beta_F = -(r**2+cry_xmax*(cry_xmax-(2*cry_rcurv))-s_P**2-x_P**2)/(2*s_P)
                A_F = alpha_F**2 + 1
                B_F = 2*(alpha_F*beta_F) - 2*cry_rcurv
                C_F = beta_F**2 - cry_xmax*(cry_xmax-2*cry_rcurv)
            
            x_F = (-B_F-np.sqrt(B_F**2-4*(A_F*C_F)))/(2*A_F)
            s_F = alpha_F*x_F + beta_F

    # MISCUT 4th step: deflection and length calculation
    a = np.sqrt(s_F**2+(x-x_F)**2)
    tP = np.arccos((2*(r**2)-a**2)/(2*(r**2)))
    tdefl = np.arcsin((s_F-s_P)/r)
    L_chan = r*tP

    xp_rel = xp - cry_miscut

    ymin = -cry_ymax/2
    ymax =  cry_ymax/2

    # FIRST CASE: p don't interact with crystal
    if (y < ymin or y > ymax or x > cry_xmax):
        x  = x + xp*s_length
        y     = y + yp*s_length
        iProc = proc_out
        return x, xp, y, yp, pc, iProc

    # SECOND CASE: p hits the amorphous layer
    elif (x < cry_alayer or y-ymin < cry_alayer or ymax-y < cry_alayer):
        x0    = x
        y0    = y
        a_eq  = 1 + xp**2
        b_eq  = (2*x)*xp - (2*xp)*cry_rcurv
        c_eq  = x**2 - (2*x)*cry_rcurv
        delta = b_eq**2 - (4*a_eq)*c_eq
        s     = (-b_eq+np.sqrt(delta))/(2*a_eq)
        if (s >= s_length):
            s = s_length
        
        x   =  xp*s + x0
        len_xs = np.sqrt((x-x0)**2 + s**2)

        if (yp >= 0 and y + yp*s <= ymax):
            len_ys = yp*len_xs

        elif(yp < 0 and y + yp*s >= ymin):
            len_ys = yp*len_xs

        else:
            s      = (ymax-y)/yp
            len_ys = np.sqrt((ymax-y)**2 + s**2)
            x   = x0 + xp*s
            len_xs = np.sqrt((x-x0)**2 + s**2)
        
        am_len = np.sqrt(len_xs**2 + len_ys**2)
        s     = s/2
        x  = x0 + xp*s
        y     = y0 + yp*s
        iProc = proc_AM

        dest=np.array(dest)
        dest = calcionloss(am_len,dest,betar,bgr,gammar,tmax,plen,
                            exenergy,zatom,rho,anuc)
        
        xp=np.array(xp)
        yp=np.array(yp)
        pc=np.array(pc)
        iProc=np.array(iProc) 
        pyk2.pyk2_crymoveam(nam,am_len,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,csref0,
                        csref1,csref4,csref5,collnt,iProc)

        x = x + xp*(s_length-s)
        y = y + yp*(s_length-s)
        return x, xp, y, yp, pc, iProc

    elif (x > cry_xmax-cry_alayer and x < cry_xmax):
        iProc = proc_AM
        
        dest=np.array(dest)  
        dest = dest = calcionloss(s_length,dest,betar,bgr,gammar,tmax,plen,
                            exenergy,zatom,rho,anuc)
    
        xp=np.array(xp)
        yp=np.array(yp)
        pc=np.array(pc)
        iProc=np.array(iProc)  
        pyk2.pyk2_crymoveam(nam,s_length,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,csref0,
                    csref1,csref4,csref5,collnt,iProc)

        return x, xp, y, yp, pc, iProc
    

    #THIRD CASE: the p interacts with the crystal.
    #Define typical angles/probabilities for orientation 110
    xpcrit0 = np.sqrt((2.0e-9*eUm)/pc)    #Critical angle (rad) for straight crystals
    Rcrit   = (pc/(2.0e-6*eUm))*ai #Critical curvature radius [m]

    #If R>Rcritical=>no channeling is possible (ratio<1)
    ratio  = cry_rcurv/Rcrit
    xpcrit = (xpcrit0*(cry_rcurv-Rcrit))/cry_rcurv #Critical angle for curved crystal

    if (ratio <= 1): #no possibile channeling
        Ang_rms = ((c_v1*0.42)*xpcrit0)*np.sin(1.4*ratio) #RMS scattering
        Ang_avr = ((c_v2*xpcrit0)*5.0e-2)*ratio                         #Average angle reflection
        Vcapt   = 0                                                #Probability of VC

    elif (ratio <= 3): #Strongly bent crystal
        Ang_rms = ((c_v1*0.42)*xpcrit0)*np.sin(0.4713*ratio + 0.85) #RMS scattering
        Ang_avr = (c_v2*xpcrit0)*(0.1972*ratio - 0.1472)                  #Average angle reflection
        Vcapt   = 7.0e-4*(ratio - 0.7)/pc**2.0e-1                           #Correction by sasha drozdin/armen
        #K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

    else: #Rcry >> Rcrit
        Ang_rms = (c_v1*xpcrit0)*(1/ratio)                #RMS scattering
        Ang_avr = (c_v2*xpcrit0)*(1 - 1.6667/ratio) #Average angle for VR
        Vcapt   = 7.0e-4*(ratio - 0.7)/pc**2.0e-1 #Probability for VC correction by sasha drozdin/armen
        #K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

    if (cry_orient == 2):
        Ang_avr = Ang_avr*0.93
        Ang_rms = Ang_rms*1.05
        xpcrit  = xpcrit*0.98

    if (np.abs(xp_rel) < xpcrit):
        alpha  = xp_rel/xpcrit
        Chann  = np.sqrt(0.9*(1 - alpha**2))*np.sqrt(1-(1/ratio)) #Saturation at 95%
        N_atom = 1.0e-1

        #if they can channel: 2 options
        if (get_random() <= Chann): #option 1:channeling
            TLdech1 = (const_dech*pc)*(1-1/ratio)**2 #Updated calculate typical dech. length(m)

            if(get_random() <= N_atom):
                TLdech1 = ((const_dech/2.0e2)*pc)*(1-1/ratio)**2  #Updated dechanneling length (m)      

            Dechan = -np.log(get_random()) #Probability of dechanneling
            Ldech  = TLdech1*Dechan   #Actual dechan. length

            #careful: the dechanneling lentgh is along the trajectory
            #of the particle -not along the longitudinal coordinate...
            if (Ldech < L_chan):
                iProc = proc_DC
                Dxp   = Ldech/r #Change angle from channeling [mrad]
                Sdech = Ldech*np.cos(cry_miscut + 0.5*Dxp)
                x  = x  + Ldech*(np.sin(0.5*Dxp+cry_miscut)) #Trajectory at channeling exit
                xp    = xp + Dxp + (2*(get_random()-0.5))*xpcrit
                y     = y  + yp * Sdech

                dest=np.array(dest)
                dest = calcionloss(Ldech,dest,betar,bgr,gammar,tmax,plen,
                                    exenergy,zatom,rho,anuc)
                pc = pc - 0.5*dest*Ldech #Energy loss to ionization while in CH [GeV]
                x  = x  + (0.5*(s_length-Sdech))*xp
                y  = y  + (0.5*(s_length-Sdech))*yp

                dest=np.array(dest)
                dest = calcionloss(s_length-Sdech,dest,betar,bgr,gammar,tmax,plen,
                                    exenergy,zatom,rho,anuc)
                xp=np.array(xp)
                yp=np.array(yp)
                pc=np.array(pc)
                iProc=np.array(iProc) 
                pyk2.pyk2_crymoveam(nam,s_length-Sdech,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,
                                csref0,csref1,csref4,csref5,collnt,iProc)
                x = x + (0.5*(s_length-Sdech))*xp
                y = y + (0.5*(s_length-Sdech))*yp

            else:
                iProc = proc_CH
                xpin  = xp
                ypin  = yp

                #check if a nuclear interaction happen while in CH
                x=np.array(x)
                xp=np.array(xp)
                yp=np.array(yp)
                pc=np.array(pc)
                iProc=np.array(iProc) 
                pyk2.pyk2_crymovech(nam,L_chan,x,xp,yp,pc,cry_rcurv,Rcrit,rho,anuc,zatom,emr,hcut,bnref,csect,
                                csref0,csref1,csref4,csref5,eUm,collnt,iProc)

                if (iProc != proc_CH):
                    #if an nuclear interaction happened, move until the middle with initial xp,yp:
                    #propagate until the "crystal exit" with the new xp,yp accordingly with the rest
                    #of the code in "thin lens approx"
                    x = x + (0.5*L_chan)*xpin
                    y = y + (0.5*L_chan)*ypin
                    x = x + (0.5*L_chan)*xp
                    y = y + (0.5*L_chan)*yp

                    dest=np.array(dest)
                    dest = calcionloss(length,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    pc = pc - dest*length #energy loss to ionization [GeV]

                else:
                    Dxp = tdefl + (0.5*get_random_gauss(0))*xpcrit #Change angle[rad]
                    
                    xp  = Dxp
                    x = x + L_chan*(np.sin(0.5*Dxp)) #Trajectory at channeling exit
                    y   = y + s_length * yp

                    dest=np.array(dest)
                    dest = calcionloss(length,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    pc = pc - (0.5*dest)*length #energy loss to ionization [GeV]      

        else: #Option 2: VR
            #good for channeling but don't channel (1-2)
            iProc = proc_VR

            xp = xp + (0.45*(xp_rel/xpcrit + 1))*Ang_avr
            x  = x  + (0.5*s_length)*xp
            y  = y  + (0.5*s_length)*yp

            dest=np.array(dest)
            dest = calcionloss(s_length,dest,betar,bgr,gammar,tmax,plen,
                                exenergy,zatom,rho,anuc)
            xp=np.array(xp)
            yp=np.array(yp)
            pc=np.array(pc)
            iProc=np.array(iProc) 
            pyk2.pyk2_crymoveam(nam,s_length,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,csref0,
                            csref1,csref4,csref5,collnt,iProc)

            x = x + (0.5*s_length)*xp
            y = y + (0.5*s_length)*yp  

    else: #case 3-2: no good for channeling. check if the  can VR
        Lrefl = xp_rel*r #Distance of refl. point [m]
        Srefl = np.sin(xp_rel/2 + cry_miscut)*Lrefl

        if(Lrefl > 0 and Lrefl < L_chan): #VR point inside

        #2 options: volume capture and volume reflection

            if(get_random() > Vcapt or zn == 0): #Option 1: VR
                iProc = proc_VR
                x  = x + xp*Srefl
                y     = y + yp*Srefl
                Dxp   = Ang_avr
                xp    = xp + Dxp + Ang_rms*get_random_gauss(0)
                x  = x  + (0.5*xp)*(s_length - Srefl)
                y     = y  + (0.5*yp)*(s_length - Srefl)

                dest=np.array(dest)
                dest = calcionloss(s_length-Srefl,dest,betar,bgr,gammar,tmax,plen,
                                    exenergy,zatom,rho,anuc)
                xp=np.array(xp)
                yp=np.array(yp)
                pc=np.array(pc)
                iProc=np.array(iProc) 
                pyk2.pyk2_crymoveam(nam,s_length-Srefl,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,
                                csref0,csref1,csref4,csref5,collnt,iProc)
                x = x + (0.5*xp)*(s_length - Srefl)
                y = y + (0.5*yp)*(s_length - Srefl)

            else: #Option 2: VC
                x = x + xp*Srefl
                y = y + yp*Srefl

                TLdech2 = (const_dech/1.0e1)*pc*(1-1/ratio)**2          #Updated typical dechanneling length(m)
                Ldech   = TLdech2*(np.sqrt(1.0e-2 - np.log(get_random())) - 1.0e-1)**2 #Updated DC length
                tdech   = Ldech/cry_rcurv
                Sdech   = Ldech*np.cos(xp + 0.5*tdech)

                if (Ldech < length-Lrefl):
                    iProc = proc_DC
                    Dxp   = Ldech/cry_rcurv + (0.5*get_random_gauss(0))*xpcrit
                    x  = x + Ldech*(np.sin(0.5*Dxp+xp)) #Trajectory at channeling exit
                    y     = y + Sdech*yp
                    xp    =  Dxp
                    Red_S = (s_length - Srefl) - Sdech
                    x  = x + (0.5*xp)*Red_S
                    y     = y + (0.5*yp)*Red_S

                    dest=np.array(dest)
                    dest = calcionloss(Srefl,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)

                    pc = pc - dest*Srefl #"added" energy loss before capture

                    dest=np.array(dest)
                    dest = calcionloss(Sdech,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    pc = pc - (0.5*dest)*Sdech #"added" energy loss while captured

                    dest=np.array(dest)
                    dest = calcionloss(Red_S,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    xp=np.array(xp)
                    yp=np.array(yp)
                    pc=np.array(pc)
                    iProc=np.array(iProc) 
                    pyk2.pyk2_crymoveam(nam,Red_S,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,csref0,
                                    csref1,csref4,csref5,collnt,iProc)
                    x = x + (0.5*xp)*Red_S
                    y = y + (0.5*yp)*Red_S

                else:
                    iProc   = proc_VC
                    Rlength = length - Lrefl
                    tchan   = Rlength/cry_rcurv
                    Red_S   = Rlength*np.cos(xp + 0.5*tchan)

                    dest=np.array(dest)
                    dest = calcionloss(Lrefl,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    pc   = pc - dest*Lrefl #"added" energy loss before capture
                    xpin = xp
                    ypin = yp

                    #Check if a nuclear interaction happen while in ch
                    x=np.array(x)
                    xp=np.array(xp)
                    yp=np.array(yp)
                    pc=np.array(pc)
                    iProc=np.array(iProc) 
                    pyk2.pyk2_crymovech(nam,Rlength,x,xp,yp,pc,cry_rcurv,Rcrit,rho,anuc,zatom,emr,hcut,bnref,csect,
                                    csref0,csref1,csref4,csref5,eUm,collnt,iProc)
                                    
                    if (iProc != proc_VC):
                        #if an nuclear interaction happened, move until the middle with initial xp,yp: propagate until
                        #the "crystal exit" with the new xp,yp aciordingly with the rest of the code in "thin lens approx"
                        x = x + (0.5*Rlength)*xpin
                        y = y + (0.5*Rlength)*ypin
                        x = x + (0.5*Rlength)*xp
                        y = y + (0.5*Rlength)*yp

                        dest=np.array(dest)
                        dest = calcionloss(Rlength,dest,betar,bgr,gammar,tmax,plen,
                                            exenergy,zatom,rho,anuc)
                        pc = pc - dest*Rlength

                    else:
                        Dxp = (length-Lrefl)/cry_rcurv
                        x = x + np.sin(0.5*Dxp+xp)*Rlength #Trajectory at channeling exit
                        y   = y + Red_S*yp
                        xp  = tdefl + (0.5*get_random_gauss(0))*xpcrit #[mrad]

                        dest=np.array(dest)
                        dest = calcionloss(Rlength,dest,betar,bgr,gammar,tmax,plen,
                                            exenergy,zatom,rho,anuc)
                        pc = pc - (0.5*dest)*Rlength  #"added" energy loss once captured

        else:
            #Case 3-3: move in amorphous substance (big input angles)
            #Modified for transition vram daniele
            if (xp_rel > tdefl-cry_miscut + 2*xpcrit or xp_rel < -xpcrit):
                iProc = proc_AM
                x  = x + (0.5*s_length)*xp
                y     = y + (0.5*s_length)*yp
                if(zn > 0):
                    dest=np.array(dest)
                    dest = calcionloss(s_length,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    xp=np.array(xp)
                    yp=np.array(yp)
                    pc=np.array(pc)
                    iProc=np.array(iProc) 
                    pyk2.pyk2_crymoveam(nam,s_length,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,csref0,
                                    csref1,csref4,csref5,collnt,iProc)
            
                x = x + (0.5*s_length)*xp
                y = y + (0.5*s_length)*yp

            else:
                Pvr = (xp_rel-(tdefl-cry_miscut))/(2*xpcrit)
                if(get_random() > Pvr):

                    iProc = proc_TRVR
                    x  = x + xp*Srefl
                    y     = y + yp*Srefl

                    Dxp = (((-3*Ang_rms)*xp_rel)/(2*xpcrit) + Ang_avr) + ((3*Ang_rms)*(tdefl-cry_miscut))/(2*xpcrit)
                    xp  = xp + Dxp
                    x = x + (0.5*xp)*(s_length-Srefl)
                    y   = y + (0.5*yp)*(s_length-Srefl)

                    dest=np.array(dest)
                    dest = calcionloss(s_length-Srefl,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    xp=np.array(xp)
                    yp=np.array(yp)
                    pc=np.array(pc)
                    iProc=np.array(iProc) 
                    pyk2.pyk2_crymoveam(nam,s_length-Srefl,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,
                                    csref0,csref1,csref4,csref5,collnt,iProc)
                    x = x + (0.5*xp)*(s_length - Srefl)
                    y = y + (0.5*yp)*(s_length - Srefl)

                else:
                    iProc = proc_TRAM
                    x = x + xp*Srefl
                    y = y + yp*Srefl
                    Dxp = ((((-1*(13.6/pc))*np.sqrt(s_length/dlri))*1.0e-3)*xp_rel)/(2*xpcrit) + (((13.6/pc)*np.sqrt(s_length/dlri))*1.0e-3)*(1+(tdefl-cry_miscut)/(2*xpcrit))
                    xp = xp+Dxp
                    x  = x + (0.5*xp)*(s_length-Srefl)
                    y  = y + (0.5*yp)*(s_length-Srefl)

                    dest=np.array(dest)
                    dest = calcionloss(s_length-Srefl,dest,betar,bgr,gammar,tmax,plen,
                                        exenergy,zatom,rho,anuc)
                    xp=np.array(xp)
                    yp=np.array(yp)
                    pc=np.array(pc)
                    iProc=np.array(iProc) 
                    pyk2.pyk2_crymoveam(nam,s_length-Srefl,dest,dlyi,dlri,xp,yp,pc,anuc,zatom,emr,hcut,bnref,
                                    csref0,csref1,csref4,csref5,collnt,iProc)
                    x = x + (0.5*xp)*(s_length - Srefl)
                    y = y + (0.5*yp)*(s_length - Srefl)
                
    return x, xp, y, yp, pc, iProc



def calcionloss(dz,EnLo,betar,bgr,gammar,tmax,plen,exenergy,zatom,rho,anuc):

    from .k2_random import get_random

    k = 0.307075 # Constant in front bethe-bloch [mev g^-1 cm^2]
    pmae = 0.51099890
    pmap = 938.271998

    thl  = (((((4*k)*zatom)*dz)*1.0e2)*rho)/(anuc*betar**2) # [MeV]
    EnLo = ((k*zatom)/(anuc*betar**2)) * (
        0.5*np.log(((((2*pmae)*bgr)*bgr)*tmax)/(1.0e6*exenergy**2)) -
        betar**2 - np.log(plen/(exenergy*1.0e3)) - np.log(bgr) + 0.5   
    )
    EnLo = ((EnLo*rho)*1.0e-1)*dz # [GeV]
    Tt   = (EnLo*1.0e3)+thl          # [MeV]

    cs_tail   = ((k*zatom)/(anuc*betar**2)) * ((0.5*((1/Tt)-(1/tmax))) -
        (np.log(tmax/Tt)*(betar**2)/(2*tmax)) + ((tmax-Tt)/((4*(gammar**2))*(pmap**2))))
    prob_tail = ((cs_tail*rho)*dz)*1.0e2

    if (get_random() < prob_tail):
        EnLo = ((k*zatom)/(anuc*betar**2)) * (
        0.5*np.log((2*pmae*bgr*bgr*tmax)/(1.0e6*exenergy**2)) -     
        betar**2 - np.log(plen/(exenergy*1.0e3)) - np.log(bgr) + 0.5 +
        tmax**2/(8*(gammar**2)*(pmap**2))
        )
        EnLo = (EnLo*rho)*1.0e-1 # [GeV/m]

    else:
        EnLo = EnLo/dz # [GeV/m]

    return EnLo