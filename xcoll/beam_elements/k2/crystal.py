import numpy as np

def crystal(*, val_part_hit, val_part_abs, val_part_impact, val_part_indiv, run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, run_csref0, run_csref1, run_csref4, run_csref5, run_dlri, run_dlyi, run_eum, run_ai, run_collnt, run_bn, c_length, nhit, nabs, isimp, s, zlm, x, xp, xp_in0, z, zp, p, x_in0):

    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.") 

    val_part_hit = np.array(val_part_hit, dtype=np.int64)
    val_part_abs = np.array(val_part_abs, dtype=np.int64)
    val_part_impact = np.array(val_part_impact, dtype=np.int64)
    val_part_indiv = np.array(val_part_indiv, dtype=np.int64)

    run_exenergy = np.array(run_exenergy, dtype=np.int64)
    run_anuc = np.array(run_anuc, dtype=np.int64)
    zatom = np.array(zatom, dtype=np.int64)
    run_emr = np.array(run_emr, dtype=np.int64)
    run_rho = np.array(run_rho, dtype=np.int64)
    run_hcut = np.array(run_hcut, dtype=np.int64)
    run_bnref = np.array(run_bnref, dtype=np.int64)

    run_csref0 = np.array(run_csref0, dtype=np.int64)
    run_csref1 = np.array(run_csref1, dtype=np.int64)
    run_csref4 = np.array(run_csref4, dtype=np.int64)
    run_csref5 = np.array(run_csref5, dtype=np.int64)

    run_dlri = np.array(run_dlri, dtype=np.int64)
    run_dlyi = np.array(run_dlyi, dtype=np.int64)
    run_eUm = np.array(run_eUm, dtype=np.int64)
    run_ai = np.array(run_ai, dtype=np.int64)
    run_collnt = np.array(run_collnt, dtype=np.int64)
    run_bn = np.array(run_bn, dtype=np.int64)
    c_length = np.array(c_length, dtype=np.int64)

    nhit = np.array(nhit, dtype=np.int64)
    nabs = np.array(nabs, dtype=np.int64)
    lhit = np.array(lhit, dtype=np.int64)
    isImp = np.array(isImp, dtype=np.int64)
    s = np.array(s, dtype=np.int64)
    zlm = np.array(zlm, dtype=np.int64)

    x = np.array(x, dtype=np.int64)
    xp = np.array(xp, dtype=np.int64)
    xp_in0 = np.array(xp_in0, dtype=np.int64)
    z = np.array(z, dtype=np.int64)
    zp = np.array(zp, dtype=np.int64)
    p = np.array(p, dtype=np.int64)
    x_in0 = np.array(x_in0, dtype=np.int64)

    ########################################################

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

    ########################################################

    iProc       = 0
    n_chan      = 0
    n_VR        = 0
    n_amorphous = 0

    # Shared settings for the currently active crystal
    c_orient   = 0    # Crystal orientation [0-2]
    c_rcurv    = 0 # Crystal geometrical parameters [m]
    c_xmax     = 0 # Crystal geometrical parameters [m]
    c_ymax     = 0 # Crystal geometrical parameters [m]
    c_alayer   = 0 # Crystal amorphous layer [mm]
    c_miscut   = 0 # Crystal miscut angle in rad
    c_cpTilt   = 0 # Cosine of positive crystal tilt
    c_spTilt   = 0 # Sine of positive crystal tilt
    c_cnTilt   = 0 # Cosine of negative crystal tilt
    c_snTilt   = 0 # Sine of negative crystal tilt
    c_cBend    = 0 # Cosine of crystal bend
    c_sBend    = 0 # Sine of crystal bend
    cry_tilt   = 0 # Crystal tilt angle in rad
    cry_length = 0 # Crystal length [m]
    cry_bend   = 0 # Crystal bending angle in rad

    # Rutherford Scatter
    tlcut_cry = 0.0009982
    cgen_cry = np.empty(200,dtype=float64)
    # mcurr_cry
    zatom_curr_cry # Current zatom, used for Rutherford scattering integration
    emr_curr_cry # Current emr, used for Rutherford scattering integration
    
    # enr
    # mom
    # betar
    # gammar
    # bgr
    # tmax
    # plen

    aTF = 0.194e-10 # Screening function [m]
    dP  = 1.920e-10 # Distance between planes (110) [m]
    u1  = 0.075e-10 # Thermal vibrations amplitude

    # pp cross-sections and parameters for energy dependence
    pptref_cry = 0.040
    freeco_cry = 1.618

    # Crystal Specific Material Arrays
    # logical,         save :: validMat(nmat) = .false. # True for materials the crystal module supports
    # real(kind=fPrec),save :: dlri     = zero
    # real(kind=fPrec),save :: dlyi     = zero
    # real(kind=fPrec),save :: ai       = zero
    # real(kind=fPrec),save :: eUm      = zero
    # real(kind=fPrec),save :: collnt   = zero    # Nuclear Collision length [m]

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

    # ================================================================================================ #
    #  Initialise the crystal module
    # ================================================================================================ #
    
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

    # Determining if particle previously interacted with a crystal and storing the process ID
    if(cry_proc_tmp != proc_out):
        cry_proc_prev = cry_proc_tmp
        
    iProc       = proc_out
    cry_proc    = proc_out

  # Transform in the crystal reference system
  # 1st transformation: shift of the center of the reference frame
    if(cry_tilt < 0):
        s_shift = s
        shift   = c_rcurv*(1 - c_cpTilt)
        if(cry_tilt < -cry_bend):
            shift = c_rcurv*(c_cnTilt - np.cos(cry_bend - cry_tilt))
            
        x_shift = x - shift

    else:
        s_shift = s
        x_shift = x

    # 2nd transformation: rotation
    s_rot  = x_shift*c_spTilt + s_shift*c_cpTilt
    x_rot  = x_shift*c_cpTilt - s_shift*c_spTilt
    xp_rot = xp - cry_tilt

    # 3rd transformation: drift to the new coordinate s=0
    xp = xp_rot
    x  = x_rot - xp_rot*s_rot
    z  = z - zp*s_rot
    s  = 0

    # Check that particle hit the crystal
    if(x >= 0 & x < c_xmax):

    # MISCUT first step: P coordinates (center of curvature of crystalline planes)
        s_P = (c_rcurv-c_xmax)*np.sin(-c_miscut)
        x_P = c_xmax + (c_rcurv-c_xmax)*np.cos(-c_miscut)

        call cry_interact(x,xp,z,zp,p,cry_length,s_P,x_P,run_exenergy,run_anuc,run_zatom,run_emr,run_rho,run_hcut,
                      run_bnref,run_csref0,run_csref1,run_csref4,run_csref5,run_dlri,run_dlyi,run_eum,run_ai,run_collnt,run_bn)

        s   = c_rcurv*c_sBend
        zlm = c_rcurv*c_sBend

        if(iProc != proc_out):
            
            isImp  = True
            nhit   = nhit + 1
            lhit   = 1
            val_part_impact = x0
            val_part_indiv  = xp0

    else:

        if(x < 0): # Crystal can be hit from below
            xp_tangent = np.sqrt((-(2*x)*c_rcurv + x**2)/(c_rcurv**2))

        else:         # Crystal can be hit from above
            xp_tangent = np.arcsin((c_rcurv*(1 - c_cBend) - x)/np.sqrt(((2*c_rcurv)*(c_rcurv - x))*(1 - c_cBend) + x**2))

        # If the hit is below, the angle must be greater or equal than the tangent,
        # or if the hit is above, the angle must be smaller or equal than the tangent
        if((x < 0 & xp >= xp_tangent) | (x >= 0 & xp <= xp_tangent)):

            # If it hits the crystal, calculate in which point and apply the transformation and drift to that point
            a_eq  = 1 + xp**2
            b_eq  = (2*xp)*(x - c_rcurv)
            c_eq  = -(2*x)*c_rcurv + x**2
            delta = b_eq**2 - four*(a_eq*c_eq)
            s_int = (-b_eq - np.sqrt(delta))/(2*a_eq)
            s_imp = s_int

            # MISCUT first step: P coordinates (center of curvature of crystalline planes)
            s_P_tmp = (c_rcurv-c_xmax)*np.sin(-c_miscut)
            x_P_tmp = c_xmax + (c_rcurv-c_xmax)*np.cos(-c_miscut)


            if(s_int < c_rcurv*c_sBend):
                # Transform to a new reference system: shift and rotate
                x_int  = xp*s_int + x
                xp_int = xp
                z      = z + zp*s_int
                x      = 0
                s      = 0

                tilt_int = s_int/c_rcurv
                xp       = xp-tilt_int

                # MISCUT first step (bis): transform P in new reference system
                # Translation
                s_P_tmp = s_P_tmp - s_int
                x_P_tmp = x_P_tmp - x_int
                # Rotation
                s_P = s_P_tmp*np.cos(tilt_int) + x_P_tmp*np.sin(tilt_int)
                x_P = -s_P_tmp*np.sin(tilt_int) + x_P_tmp*np.cos(tilt_int)

                call cry_interact(x,xp,z,zp,p,cry_length-(tilt_int*c_rcurv),s_P,x_P,run_exenergy,run_anuc,
                        run_zatom,run_emr,run_rho,run_hcut,run_bnref,run_csref0,run_csref1,run_csref4,run_csref5,
                        run_dlri,run_dlyi,run_eum,run_ai,run_collnt,run_bn)

                s   = c_rcurv*np.sin(cry_bend - tilt_int)
                zlm = c_rcurv*np.sin(cry_bend - tilt_int)
                
                if(iProc != proc_out):
                    x_rot    = x_int
                    s_rot    = s_int
                    xp_rot   = xp_int
                    s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
                    x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
                    xp_shift = xp_rot + cry_tilt

                    if(cry_tilt < 0):
                        x0  = x_shift + shift
                        xp0 = xp_shift

                    else:
                        x0  = x_shift
                        xp0 = xp_shift

                    isImp        = .true.
                    nhit         = nhit + 1
                    lhit      = 1
                    val_part_impact    = x0
                    val_part_indiv     = xp0

                # un-rotate
                x_temp  = x
                s_temp  = s
                xp_temp = xp
                s       =  s_temp*np.cos(-tilt_int) + x_temp*np.sin(-tilt_int)
                x       = -s_temp*np.sin(-tilt_int) + x_temp*np.cos(-tilt_int)
                xp      = xp_temp + tilt_int

                # 2nd: shift back the 2 axis
                x = x + x_int
                s = s + s_int

            else:

                s = c_rcurv*np.sin(cry_length/c_rcurv)
                x = x + s*xp
                z = z + s*zp

        else:

            s = c_rcurv*np.sin(cry_length/c_rcurv)
            x = x + s*xp
            z = z + s*zp

    # transform back from the crystal to the collimator reference system
    # 1st: un-rotate the coordinates
    x_rot  = x
    s_rot  = s
    xp_rot = xp

    s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
    x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
    xp_shift = xp_rot + cry_tilt

    # 2nd: shift back the reference frame
    if(cry_tilt < 0):
        s  = s_shift
        x  = x_shift + shift
        xp = xp_shift
    else:
        x  = x_shift
        s  = s_shift
        xp = xp_shift
    
    # 3rd: shift to new S=Length position
    x = xp*(c_length - s) + x
    z = zp*(c_length - s) + z
    s = c_length

    nabs = 0
    cry_proc = iProc
    if(iProc == proc_AM):
        n_amorphous = n_amorphous + 1
    else if(iProc == proc_VR):
        n_VR = n_VR + 1
    else if(iProc == proc_CH):
        n_chan = n_Chan + 1
    else if(iProc == proc_absorbed):
        nabs = 1   # todo: do we need to set part_abs_pos etc?
    else if(iProc == proc_ch_absorbed):
        nabs = 1

    # Storing the process ID for the next interaction
    cry_proc_tmp = cry_proc


    return val_part_hit, val_part_abs, val_part_impact, val_part_indiv, run_exenergy, run_bn, nhit, nabs, isimp, s, zlm, x, xp, xp_in0, z, zp, p, x_in