import numpy as np
from .materials import materials

def k2_track(*, material, particles, closed_orbit, angle, jaws, offset, npart, length, is_crystal, onesided):
    from .pyk2f import pyk2_run, pyk2_start_run, pyk2_per_particle, pyk2_finish

    # Go to collimator reference system (subtract closed orbit)
    dx = closed_orbit[0]
    dy = closed_orbit[1]
    dpx = closed_orbit[2]
    dpy = closed_orbit[3]
    x_part = particles.x[:npart] - dx
    xp_part = (particles.px[:npart] - dpx) * particles.rpp[:npart]
    y_part = particles.y[:npart] - dy
    yp_part = (particles.py[:npart] - dpy) * particles.rpp[:npart]
    s_part = 0 * x_part
    p_part = particles.energy[:npart] / 1e9 # Energy (not momentum) in GeV

    p0 = particles.energy0[0]/1e9 # Reference energy

    # `linside` is an array of logicals in fortran. Beware of the fortran converion:
    # True <=> -1 (https://stackoverflow.com/questions/39454349/numerical-equivalent-of-true-is-1)

    jaw_F_L = jaws[0]
    jaw_F_R = jaws[1]
    jaw_B_L = jaws[2]
    jaw_B_R = jaws[3]
    if jaw_F_L != jaw_B_L or jaw_F_R != jaw_B_R:
        raise NotImplementedError
    opening = jaw_F_L - jaw_F_R
    offset = offset + ( jaw_F_L + jaw_F_R )/2

    # Initialise arrays for FORTRAN call
    part_hit = np.zeros(len(x_part), dtype=np.int32)
    part_abs = np.zeros(len(x_part), dtype=np.int32)
    part_impact = np.zeros(len(x_part), dtype=np.float)
    part_indiv = np.zeros(len(x_part), dtype=np.float)
    part_linteract = np.zeros(len(x_part), dtype=np.float)
    nhit_stage = np.zeros(len(x_part), dtype=np.int32)
    nabs_type = np.zeros(len(x_part), dtype=np.int32)
    linside = np.zeros(len(x_part), dtype=np.int32)

    # Get material parameters
    exenergy = materials[material]['exenergy']
    anuc     = materials[material]['anuc']
    zatom    = materials[material]['zatom']
    emr      = materials[material]['emr']
    rho      = materials[material]['rho']
    hcut     = materials[material]['hcut']
    bnref    = materials[material]['bnref']
    csref0   = materials[material]['csref'][0]
    csref1   = materials[material]['csref'][1]
    csref4   = materials[material]['csref'][4]
    csref5   = materials[material]['csref'][5]
    radl     = materials[material]['radl']
    dlri     = materials[material]['dlri']
    dlyi     = materials[material]['dlyi']
    eUm      = materials[material]['eUm']
    ai       = materials[material]['ai']
    collnt   = materials[material]['collnt']
    # if self.is_crystal and not pyk2.materials[self.material]['can_be_crystal']:
    #   raise ValueError()

    # Prepare random generator in FORTRAN
    pyk2_start_run(num_particles=npart, hcut=hcut, zatom=zatom, emr=emr,
                       x_part=x_part, xp_part=xp_part, y_part=y_part, yp_part=yp_part, p_part=p_part, s_part=s_part)

    # Prepare scattering parameters (was k2coll_scatin)
    cprob, xintl, bn, ecmsq, xln15s, bpp = calculate_scattering(p0,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref)

    # Initialise number of hits and absorptions
    nhit   = 0
    nabs   = 0
    fracab = 0.
    mirror = 1.

    # Compute rotation factors for collimator rotation
    cRot   = np.cos(angle/180.*np.pi)
    sRot   = np.sin(angle/180.*np.pi)
    cRRot  = np.cos(-angle/180.*np.pi)
    sRRot  = np.sin(-angle/180.*np.pi)

    # Set energy and nucleon change variables as with the coupling
    # ien0,ien1: particle energy entering/leaving the collimator
    nnuc0 = 0
    ien0  = 0.
    nnuc1 = 0
    ien1  = 0.

    # Main loop over particles
    for this_part in range(1,npart):
        print(this_part)
        pyk2_per_particle(
            j=this_part,
            num_particles=npart,
            thisp0=p0,
            nhit=nhit,
            nabs=nabs,
            fracab=fracab,
            mirror=mirror,
            crot=cRot,
            srot=sRot,
            crrot=cRRot,
            srrot=sRRot,
            nnuc0=nnuc0,
            ien0=ien0,
            nnuc1=nnuc1,
            ien1=ien1,
            coll_exenergy=exenergy,
            coll_anuc=anuc,
            coll_zatom=zatom,
            coll_emr=emr,
            coll_rho=rho,
            coll_hcut=hcut,
            coll_bnref=bnref,
            coll_csref0=csref0,
            coll_csref1=csref1,
            coll_csref4=csref4,
            coll_csref5=csref5,
            coll_radl=radl,
            coll_dlri=dlri,
            coll_dlyi=dlyi,
            coll_eum=eUm,
            coll_ai=ai,
            coll_collnt=collnt,
            coll_cprob=cprob,
            coll_xintl=xintl,
            coll_bn=bn,
            coll_ecmsq=ecmsq,
            coll_xln15s=xln15s,
            coll_bpp=bpp,
            is_crystal=is_crystal,
            c_length=length,
            c_aperture=opening,
            c_offset=offset,
            c_tilt=np.array([0,0], dtype=np.float64),
            x_in=x_part,
            xp_in=xp_part,
            y_in=y_part,
            yp_in=yp_part,
            p_in=p_part,
            s_in=s_part,
            lhit=part_hit,
            part_abs_local=part_abs,
            impact=part_impact,
            indiv=part_indiv,
            lint=part_linteract,
            onesided=onesided,
            nhit_stage=nhit_stage,
            nabs_type=nabs_type,
            linside=linside
        )

    # Retrieve results
    pyk2_finish(num_particles=npart, x_part=x_part, xp_part=xp_part, y_part=y_part, yp_part=yp_part, p_part=p_part, s_part=s_part)

    return x_part, xp_part, y_part, yp_part, p_part, s_part, part_hit, part_abs




#         #... preparation (coll_k2 line 176 - 258)
#         if (part_abs_local[this_part]) != 0):
#             continue

#         impact[this_part] = -1
#         lint[this_part]   = -1
#         indiv[this_part]  = -1

#         x      = x_in[this_part]
#         xp     = xp_in[this_part]
#         xp_in0 = xp_in[this_part]
#         z      = y_in[this_part]
#         zp     = yp_in[this_part]
#         p      = p_in[this_part]
#         sp     = 0
#         dpop   = (p - p0)/p0
#         x_flk  = 0
#         y_flk  = 0
#         xp_flk = 0
#         yp_flk = 0

#         # ! Transform particle coordinates to get into collimator coordinate  system
#         # ! First do rotation into collimator frame
#         x  =  x_in[this_part]*cRot + sRot*y_in[this_part]
#         z  =  y_in[this_part]*cRot - sRot*x_in[this_part]
#         xp = xp_in[this_part]*cRot + sRot*yp_in[this_part]
#         zp = yp_in[this_part]*cRot - sRot*xp_in[this_part]

#         #! For 1-sided collimators consider only positive X. For negative X this_partump to the next particle
#         if (onesided && x < 0):
#             continue

#         #! Log input energy + nucleons as per the FLUKA coupling
#         nnuc0   = nnuc0 + naa[this_part]
#         ien0    = ien0 + rcp[this_part] * 1.0e3


#         #! Now mirror at the horizontal axis for negative X offset
#         if (x < 0):
#             mirror    = -1
#             tiltangle = -1*c_tilt[2]
#         else:
#             mirror    = 1
#             tiltangle = c_tilt[1]
    
#         x  = mirror*x
#         xp = mirror*xp

#         #! Shift with opening and offset
#         x = (x - c_aperture/2) - mirror*c_offset

#         #! Include collimator tilt
#         if (tiltangle > 0):
#             xp = xp - tiltangle

#         if (tiltangle < 0):
#             x  = x + np.sin(tiltangle) * c_length
#             xp = xp - tiltangle
#         end if

#         #! CRY Only: x_in0 has to be assigned after the change of reference frame
#         x_in0 = x

#     # ! After finishing the coordinate transformation, or the coordinate manipulations in case of pencil beams,
#     # ! save the initial coordinates of the impacting particles
#     xIn  = x
#     xpIn = xp
#     yIn  = z
#     ypIn = zp

#     # ! particle passing above the jaw are discarded => take new event
#     # ! entering by the face, shorten the length (zlm) and keep track of
#     # ! entrance longitudinal coordinate (keeps) for histograms

#     # ! The definition is that the collimator jaw is at x>=0.

#     # ! 1) Check whether particle hits the collimator
#     isImp = False
#     s     = 0
#     keeps = 0
#     zlm   = -1*length

#     if (is_crystal): #! This is a crystal collimator

#         x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,nhit,nabs,lhit[this_part],part_abs[this_part],impact[this_part],indiv[this_part],exenergy,bn,isImp = pyk2.pyk2_doCrystal(j,x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit,
#         part_abs_local,impact,indiv,c_length,exenergy,anuc,zatom,emr,rho,
#         hcut,bnref,csref0,csref1,csref4,csref5,dlri,dlyi,
#         eUm,ai,collnt,bn)

#         if(nabs != 0):
#             part_abs_local[this_part] = 1
#             lint[this_part] = zlm

#         sImp  = (s - c_length) + sImp
#         sOut  = s
#         xOut  = x
#         xpOut = xp
#         yOut  = z
#         ypOut = zp

#     else:
#            # ... stuff (coll_k2 line 281 - 321)
#         if(x >= 0):
#         # ! Particle hits collimator and we assume interaction length ZLM equal
#         # ! to collimator length (what if it would leave collimator after
#         # ! small length due to angle???)
#             zlm = length
#             impact[this_part] = x
#             indiv[this_part] = xp

#         elif (xp <= 0):
#         # ! Particle does not hit collimator. Interaction length ZLM is zero.
#             zlm = 0

#         else:
#         # ! Calculate s-coordinate of interaction point
#             s = (-1*x)/xp
#             if(s <= 0):
#                 write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
#            # !call prror
#             if(s < length):
#                 zlm       = length - s
#                 impact(j) = zero
#                 indiv(j)  = xp
#             else:
#                 zlm = 0

#     #   ! First do the drift part
#     #   ! DRIFT PART
#         drift_length = length - zlm

#         if(drift_length > 0):
#             if(iexact):
#                 zpj = np.sqrt(1-xp**2-zp**2)
#                 x   = x  + drift_length*(xp/zpj)
#                 z   = z  + drift_length*(zp/zpj)
#                 sp  = sp + drift_length
#             else:
#                 x  = x  + xp* drift_length
#                 z  = z  + zp * drift_length
#                 sp = sp + drift_length

#         if (zlm > 0): # Collimator jaw interaction            
#                 #... stuff (coll_k2 line 325 - 342)
#             if (not linside[this_part]):
#             #! first time particle hits collimator: entering jaw
#                 linside[this_part] = True

#                 if (dowrite_impact):
#                     if (tiltangle > 0):
#                          x_Dump = (x + c_aperture/2 + tiltangle*sp)*mirror + c_offset
#                     else:
#                         x_Dump = (x + c_aperture/2 + tiltangle*(sp - c_length))*mirror + c_offset
            
#                 xpDump = (xp + tiltangle)*mirror
#                 y_Dump = z
#                 ypDump = zp
#                 s_Dump = sp+real(j_slices-1,fPrec)*c_length

#             s_impact = sp
#             nhit = nhit + 1

#             pyk2.pyk2_jaw( ... )
#             #... stuff (coll_k2 line 345 - 391)
#             nabs_type[this_part] = nabs
#             lhit[this_part]  = 1

#             isImp = True
#             sImp  = s_impact+(real(j_slices,fPrec)-1)*c_length
#             sOut  = (s+sp)+(real(j_slices,fPrec)-1)*c_length
#             xOut  = x
#             xpOut = xp
#             yOut  = z
#             ypOut = zp

#             #! Writeout should be done for both inelastic and single diffractive. doing all transformations
#             #! in x_flk and making the set to 99.99 mm conditional for nabs=1
#             if (dowrite_impact or nabs == 1 or nabs == 4):
#             #! Transform back to lab system for writeout.
#             #! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

#                 x_flk  = xInt
#                 xp_flk = xpInt

#                 if (tiltangle > 0):
#                     x_flk  = x_flk  + tiltangle*(sInt+sp)
#                     xp_flk = xp_flk + tiltangle
                
#                 elif (tiltangle < 0):
#                     xp_flk = xp_flk + tiltangle
#                     x_flk  = x_flk  - sin_mb(tiltangle) * (length-(sInt+sp))
    

#                 x_flk  = (x_flk + c_aperture/2) + mirror*c_offset
#                 x_flk  = mirror*x_flk
#                 xp_flk = mirror*xp_flk
#                 y_flk  = (  yInt*cRRot -  x_flk*sRRot)*1.0e3
#                 yp_flk = ( ypInt*cRRot - xp_flk*sRRot)*1.0e3
#                 x_flk  = ( x_flk*cRRot +   yInt*sRRot)*1.0e3
#                 xp_flk = (xp_flk*cRRot +  ypInt*sRRot)*1.0e3
#                 s_flk  = (sInt+sp)+(real(j_slices,fPrec)-1)*c_length

#             #! Finally, the actual coordinate change to 99 mm
#             if (nabs == 1):
#                 fracab  = fracab + 1
#                 x       = 99.99e-3
#                 z       = 99.99e-3
#                 lint[this_part] = zlm
#                 part_abs_local[this_part] = 1
          
#         # ! Collimator jaw interaction
#         if (nabs != 1 && zlm > 0):
#         #! Do the rest drift, if particle left collimator early
#             drift_length = (length-(s+sp))

#         if (drift_length > 1.0e-15):
#             linside[this_part] = False

#             if (dowrite_impact):
#                 if(tiltangle > 0):
#                     x_Dump = (x + c_aperture/2 + tiltangle*(s+sp))*mirror + c_offset
#                 else:
#                     x_Dump = (x + c_aperture/2 + tiltangle*(s+sp-c_length))*mirror + c_offset
         
#                 xpDump = (xp+tiltangle)*mirror
#                 y_Dump = z
#                 ypDump = zp
#                 s_Dump = s+sp+real(j_slices-1,fPrec)*c_length

#             if (iexact):
#                 zpj = sqrt(one-xp**2-zp**2)
#                 x   = x  + drift_length*(xp/zpj)
#                 z   = z  + drift_length*(zp/zpj)
#                 sp  = sp + drift_length
#             else:
#                 x  = x  + xp * drift_length
#                 z  = z  + zp * drift_length
#                 sp = sp + drift_length
          
#         lint[this_part] = zlm - drift_length
     
#         #     ... stuff (coll_k2 line 393 - 422)

#         # ... post-process (coll_k2 line 426 - 468)

#     if (x < 99.0e-3):

#       #! Include collimator tilt
#         if(tiltangle > 0):
#             x  = x  + tiltangle*c_length
#             xp = xp + tiltangle

#         elif (tiltangle < 0):
#             x  = x  + tiltangle*c_length
#             xp = xp + tiltangle
#             x  = x  - np.sin(tiltangle) * c_length


#     #! Transform back to particle coordinates with opening and offset
#         z00 = z
#         x00 = x + mirror*c_offset
#         x   = (x + c_aperture/2) + mirror*c_offset

#     #! Now mirror at the horizontal axis for negative X offset
#         x  = mirror * x
#         xp = mirror * xp

#     #! Last do rotation into collimator frame
#         x_in(j)  =  x*cRRot +  z*sRRot
#         y_in(j)  =  z*cRRot -  x*sRRot
#         xp_in(j) = xp*cRRot + zp*sRRot
#         yp_in(j) = zp*cRRot - xp*sRRot

# # ! Log output energy + nucleons as per the FLUKA coupling
# # ! Do not log dead particles
#         nnuc1       = nnuc1 + naa[this_part]                          #! outcoming nucleons
#         ien1        = ien1  + rcp[this_part] * 1.0e3                  #! outcoming energy

#         if(is_crystal):
#             p_in(j) = p
#             s_in(j) = s_in[this_[part] + s
#         else:
#             p_in(j) = (1 + dpop) * p0
#             s_in(j) = sp + (real(j_slices,fPrec)-1) * c_length

#     else:
#         x_in(j) = x
#         y_in(j) = z
    

# #! End of loop over all particles





def calculate_scattering(p0,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref):
    
    # Output parameters
    cprob = np.array([0,0,0,0,0,0], dtype=np.float64)
    csect = np.array([0,0,0,0,0,0], dtype=np.float64) # Cross section
    
    # Constants 
    pptref = 0.04
    freeco = 1.618
    pmap  = 938.271998
    fnavo  = 6.02214076e23

    ecmsq = (2*(pmap*1.0e-3)) * p0
    xln15s = np.log(0.15*ecmsq)
    # Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    pptot = 0.041084 - 0.0023302*np.log(ecmsq) + 0.00031514*np.log(ecmsq)**2
    # Claudia used the fit from TOTEM for ppel (in barn)
    ppel = (11.7-1.59*np.log(ecmsq)+0.134*np.log(ecmsq)**2)/1e3
    # Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    ppsd = (4.3+0.3*np.log(ecmsq))/1e3

    # Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    bpp = 7.156 + 1.439*np.log(np.sqrt(ecmsq))

    # freep: number of nucleons involved in single scattering
    freep = freeco * anuc**(1/3)

    # compute pp and pn el+single diff contributions to cross-section
    # (both added : quasi-elastic or qel later)
    csect[3] = freep * ppel
    csect[4] = freep * ppsd

    # correct TOT-CSec for energy dependence of qel
    # TOT CS is here without a Coulomb contribution
    csect[0] = csref0 + freep * (pptot - pptref)
    bn = (bnref * csect[0]) / csref0

    # also correct inel-CS
    csect[1] = (csref1 * csect[0]) / csref0

    # Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect[2] = ((csect[0] - csect[1]) - csect[3]) - csect[4]
    csect[5] = csref5

    # Now add Coulomb
    csect[0] = csect[0] + csect[5]

    # Interaction length in meter
    xintl = (1.0e-2*anuc)/(((fnavo * rho)*csect[0])*1e-24)

    # Filling CProb with cumulated normalised Cross-sections
    cprob[5] = 1
    for i in range(1,5,1):
        cprob[i] = cprob[i-1] + csect[i]/csect[0]

    return cprob, xintl, bn, ecmsq, xln15s, bpp



# def calculate_ion_loss(PC,DZ,exenergy,anuc,zatom,rho):

#     # Constants
#     mom    = PC*1.0e3                         # [GeV/c] -> [MeV/c]
#     enr    = (mom*mom + 938.271998**2)**(1/2) # [MeV]
#     gammar = enr/938.271998
#     betar  = mom/enr
#     bgr    = betar*gammar
#     kine   = ((2*0.510998902)*bgr)*bgr
#     k = 0.307075

#     # Mean excitation energy
#     exEn = exenergy*1.0e3 # [MeV]

#     # Tmax is max energy loss from kinematics
#     Tmax = kine/(1 + (2*gammar)*(0.510998902/938.271998) + (0.510998902/938.271998)**2) # [MeV]

#     # Plasma energy - see PDG 2010 table 27.1
#     plen = (((rho*zatom)/anuc)**(1/2))*28.816e-6 # [MeV]

#     # Calculate threshold energy
#     # Above this threshold, the cross section for high energy loss is calculated and then
#     # a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
#     # below threshold, only the standard Bethe-Bloch is used (all particles get average energy loss)

#     # thl is 2*width of Landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
#     thl = ((((4*(k*zatom))*DZ)*1.0e2)*rho)/(anuc*betar**2) # [MeV]

#     # Bethe-Bloch mean energy loss
#     EnLo = ((k*zatom)/(anuc*betar**2)) * ((1/2)*np.log((kine*Tmax)/(exEn**2)) - betar**2 - np.log(plen/exEn) - np.log(bgr) + (1/2))
#     EnLo = ((EnLo*rho)*1.0e-1)*DZ # [GeV]

#     # Threshold Tt is Bethe-Bloch + 2*width of Landau distribution
#     Tt = EnLo*1.0e3 + thl # [MeV]

#     # Cross section - see Alfredo's presentation for derivation
#     cs_tail = ((k*zatom)/(anuc*betar**2)) * ((1/2)*((1/Tt)-(1/Tmax)) - (np.log(Tmax/Tt)*betar**2)/(2*Tmax) + (Tmax-Tt)/((4*gammar**2)*938.271998**2))

#     # Probability of being in tail: cross section * density * path length
#     prob_tail = ((cs_tail*rho)*DZ)*1.0e2

#     # Determine based on random number if tail energy loss occurs.
#     if (coll_rand() < prob_tail):
#         EnLo = ((k*zatom)/(anuc*betar**2)) * ((1/2)*np.log((kine*Tmax)/(exEn**2)) - betar**2 - np.log(plen/exEn) - np.log(bgr) + (1/2)+ Tmax**2/((8*gammar**2)*938.271998**2))
#         EnLo = (EnLo*rho)*1.0e-1 # [GeV/m]
#     else:
#     # If tail energy loss does not occur, just use the standard Bethe-Bloch
#         EnLo = EnLo/DZ  # [GeV/m]
    
#     return EnLo

