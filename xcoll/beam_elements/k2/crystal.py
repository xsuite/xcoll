import numpy as np

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


def crystal(*,x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isimp,val_part_hit,val_part_abs,val_part_impact,val_part_indiv,c_length,exenergy,rho,anuc,zatom,emr,
            dlri,dlyi,ai,eum,collnt,hcut,bnref,csref0,csref1,csref4,csref5,csect,nhit,nabs,
            cry_tilt, cry_rcurv, cry_bend, cry_alayer, cry_xmax, cry_ymax, cry_orient, cry_miscut, cry_cBend, 
            cry_sBend, cry_cpTilt, cry_spTilt, cry_cnTilt, cry_snTilt, iProc, n_chan, n_VR, n_amorphous
            ):


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

    # c_rcurv  = 0
    # c_alayer = 0
    # c_xmax   = 0
    # c_ymax   = 0
    # c_orient = 0
    # c_miscut = 0
    # cry_bend = 0
    # c_cBend  = 0
    # c_sBend  = 0
    # c_cpTilt = 0
    # c_spTilt = 0
    # c_cnTilt = 0
    # c_snTilt = 0


    # s_temp     = 0
    # s_shift    = 0
    # s_rot      = 0
    # s_int      = 0
    # x_temp     = 0
    # x_shift    = 0
    # x_rot      = 0
    # x_int      = 0
    # xp_temp    = 0
    # xp_shift   = 0
    # xp_rot     = 0
    # xp_int     = 0
    # xp_tangent = 0
    # tilt_int   = 0
    # shift      = 0
    # delta      = 0
    # a_eq       = 0
    # b_eq       = 0
    # c_eq       = 0
    # # simp      = 0

    # iProc       = proc_out

    # # Transform in the crystal reference system
    # # 1st transformation: shift of the center of the reference frame
    # if (cry_tilt < 0):
    #     s_shift = s
    #     shift   = c_rcurv*(1 - c_cpTilt)

    #     if (cry_tilt < -cry_bend):
    #         shift = c_rcurv*(c_cnTilt - np.cos(cry_bend - cry_tilt))
    #     # end if
    #     x_shift = x - shift
    # else:
    #     s_shift = s
    #     x_shift = x
    # # end if
  
    # # 2nd transformation: rotation
    # s_rot  = x_shift*c_spTilt + s_shift*c_cpTilt
    # x_rot  = x_shift*c_cpTilt - s_shift*c_spTilt
    # xp_rot = xp - cry_tilt

    # # 3rd transformation: drift to the new coordinate s=0
    # xp = xp_rot
    # x  = x_rot - xp_rot*s_rot
    # z  = z - zp*s_rot
    # s  = 0
  
    # # Check that particle hit the crystal
    # if (x >= 0 and x < c_xmax):
    #     # MISCUT first step: P coordinates (center of curvature of crystalline planes)
    #     s_P = (c_rcurv-c_xmax)*np.sin(-c_miscut)
    #     x_P = c_xmax + (c_rcurv-c_xmax)*np.cos(-c_miscut)

    #     x=np.array(x)
    #     xp=np.array(xp)
    #     z=np.array(z)
    #     zp=np.array(zp)
    #     p=np.array(p)

    #     pyk2.pyk2_cryinteract(x,xp,z,zp,p,cry_length,s_P,x_P,exenergy,rho,anuc,zatom,emr,dlri,dlyi,
    #                     ai,eum,collnt,hcut,csref0,csref1,csref4,csref5,bnref,csect)

    #     s   = c_rcurv*c_sBend
    #     zlm = c_rcurv*c_sBend

    #     if (iProc != proc_out):
    #         isimp    = True
    #         nhit     = nhit + 1
    #         val_part_hit     = 1
    #         val_part_impact   = x0
    #         val_part_indiv    = xp0
    #     # end if

    # else:

    #     if (x < 0): # Crystal can be hit from below
    #         xp_tangent = np.sqrt((-(2*x)*c_rcurv + x**2)/(c_rcurv**2))
    #     else:              # Crystal can be hit from above
    #         xp_tangent = np.sin((c_rcurv*(1 - c_cBend) - x)/np.sqrt(((2*c_rcurv)*(c_rcurv - x))*(1 - c_cBend) + x**2))
    #     # end if

    #     # If the hit is below, the angle must be greater or equal than the tangent,
    #     # or if the hit is above, the angle must be smaller or equal than the tangent
    #     if ((x < 0 and xp >= xp_tangent) or (x >= 0 and xp <= xp_tangent)):

    #         # If it hits the crystal, calculate in which point and apply the transformation and drift to that point
    #         a_eq  = 1 + xp**2
    #         b_eq  = (2*xp)*(x - c_rcurv)
    #         c_eq  = -(2*x)*c_rcurv + x**2
    #         delta = b_eq**2 - 4*(a_eq*c_eq)
    #         s_int = (-b_eq - np.sqrt(delta))/(2*a_eq)
    #         simp = s_int

    #         # MISCUT first step: P coordinates (center of curvature of crystalline planes)
    #         s_P_tmp = (c_rcurv-c_xmax)*np.sin(-c_miscut)
    #         x_P_tmp = c_xmax + (c_rcurv-c_xmax)*np.cos(-c_miscut)

    #         if (s_int < c_rcurv*c_sBend):
    #             # Transform to a new reference system: shift and rotate
    #             x_int  = xp*s_int + x
    #             xp_int = xp
    #             z      = z + zp*s_int
    #             x      = 0
    #             s      = 0

    #             tilt_int = s_int/c_rcurv
    #             xp    = xp-tilt_int

    #             # MISCUT first step (bis): transform P in new reference system
    #             # Translation
    #             s_P_tmp = s_P_tmp - s_int
    #             x_P_tmp = x_P_tmp - x_int
    #             # Rotation
    #             s_P = s_P_tmp*np.cos(tilt_int) + x_P_tmp*np.sin(tilt_int)
    #             x_P = -s_P_tmp*np.sin(tilt_int) + x_P_tmp*np.cos(tilt_int)

    #             x=np.array(x)
    #             xp=np.array(xp)
    #             z=np.array(z)
    #             zp=np.array(zp)
    #             p=np.array(p)

    #             pyk2.pyk2_cryinteract(x,xp,z,zp,p,cry_length-(tilt_int*c_rcurv),s_P,x_P,exenergy,rho,anuc,
    #                             zatom,emr,dlri,dlyi,ai,eum,collnt,hcut,csref0,csref1,
    #                             csref4,csref5,bnref,csect)
    #             s   = c_rcurv*np.sin(cry_bend - tilt_int)
    #             zlm = c_rcurv*np.sin(cry_bend - tilt_int)
                
    #             if (iProc != proc_out):
    #                 x_rot    = x_int
    #                 s_rot    = s_int
    #                 xp_rot   = xp_int
    #                 s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
    #                 x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
    #                 xp_shift = xp_rot + cry_tilt

    #                 if (cry_tilt < 0):
    #                     x0  = x_shift + shift
    #                     xp0 = xp_shift
    #                 else:
    #                     x0  = x_shift
    #                     xp0 = xp_shift
    #                 # end if

    #                 isimp     = True
    #                 nhit      = nhit + 1
    #                 val_part_hit      = 1
    #                 val_part_impact    = x0
    #                 val_part_indiv     = xp0
    #                 # end if

    #             # un-rotate
    #             x_temp  = x
    #             s_temp  = s
    #             xp_temp = xp
    #             s       =  s_temp*np.cos(-tilt_int) + x_temp*np.sin(-tilt_int)
    #             x    = -s_temp*np.sin(-tilt_int) + x_temp*np.cos(-tilt_int)
    #             xp   = xp_temp + tilt_int

    #             # 2nd: shift back the 2 axis
    #             x = x + x_int
    #             s = s + s_int

    #         else:

    #             s = c_rcurv*np.sin(cry_length/c_rcurv)
    #             x = x + s*xp
    #             z = z + s*zp

    #         # end if

    #     else:

    #         s = c_rcurv*np.sin(cry_length/c_rcurv)
    #         x = x + s*xp
    #         z = z + s*zp

    #     # end if

    # # end if

    # # transform back from the crystal to the collimator reference system
    # # 1st: un-rotate the coordinates
    # x_rot  = x
    # s_rot  = s
    # xp_rot = xp

    # s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
    # x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
    # xp_shift = xp_rot + cry_tilt

    # # 2nd: shift back the reference frame
    # if (cry_tilt < 0):
    #     s  = s_shift
    #     x  = x_shift + shift
    #     xp = xp_shift
    # else:
    #     x  = x_shift
    #     s  = s_shift
    #     xp = xp_shift
    # # end if

    # # 3rd: shift to new S=Length position
    # x = xp*(c_length - s) + x
    # z = zp*(c_length - s) + z
    # s = c_length

    # nabs = 0
  
    # if (iProc == proc_AM):
    #     n_amorphous = n_amorphous + 1

    # elif (iProc == proc_VR):
    #     n_VR = n_VR + 1

    # elif (iProc == proc_CH):
    #     n_chan = n_chan + 1

    # elif (iProc == proc_absorbed):
    #     nabs = 1   # TODO: do we need to set part_abs_pos etc?

    # elif (iProc == proc_ch_absorbed):
    #     nabs = 1
    # # end if


    ########################################################

    pyk2.pyk2_docrystal(x=x,
                        xp=xp,
                        z=z,
                        zp=zp,
                        s=s,
                        p=p,
                        x0=x0,
                        xp0=xp0, 
                        zlm=zlm,
                        s_imp=s_imp,
                        isimp=isimp,
                        nhit=nhit,
                        nabs=nabs,
                        lhit=val_part_hit,
                        part_abs=val_part_abs,
                        impact=val_part_impact,
                        indiv=val_part_indiv,
                        c_length=c_length,
                        exenergy=exenergy,
                        rho=rho,
                        anuc=anuc,
                        zatom=zatom,
                        emr=emr,
                        dlri=dlri,
                        dlyi=dlyi,
                        ai=ai,
                        eum=eum,
                        collnt=collnt,
                        hcut=hcut,
                        csref0=csref0,
                        csref1=csref1,
                        csref4=csref4,
                        csref5=csref5,
                        bnref=bnref,
                        csect=csect,
                        cry_tilt=cry_tilt,
                        c_rcurv=cry_rcurv,
                        c_alayer=cry_alayer,
                        c_xmax=cry_xmax,
                        c_ymax=cry_ymax,
                        c_orient=cry_orient,
                        c_miscut=cry_miscut,
                        cry_bend=cry_bend,
                        c_cbend=cry_cBend,
                        c_sbend=cry_sBend,
                        c_cptilt=cry_cpTilt,
                        c_sptilt=cry_spTilt,
                        c_cntilt=cry_cnTilt,
                        c_sntilt=cry_snTilt,
                        iproc=iProc,
                        n_chan=n_chan,
                        n_vr=n_VR,
                        n_amorphous=n_amorphous
                        )

    ########################################################
    
    return val_part_hit, val_part_abs, val_part_impact, val_part_indiv, nhit, nabs, s_imp, isimp, s, zlm, x, xp, x0, xp0, z, zp, p, iProc, n_chan, n_VR, n_amorphous
