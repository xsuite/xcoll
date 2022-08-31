

! ================================================================================================ !
!
!  Crystal Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  Written by: Igor Yazynin, Valentina Previtali and Daniele Mirarchi, BE-ABP-HSS
!  Re-written for SixTrack 5 by: Marco D'Andrea and Veronica K. Berglyd Olsen, BE-ABP-HSS (2019)
!
!  Last modified: 2021-12-03
!
! ================================================================================================ !
module coll_crystal

  use floatPrecision
  use numerical_constants

  implicit none

  ! integer,          private, save :: iProc       = 0
  ! integer,          private, save :: n_chan      = 0
  ! integer,          private, save :: n_VR        = 0
  ! integer,          private, save :: n_amorphous = 0

  ! Shared settings for the currently active crystal
  ! integer,          private, save :: c_orient   = 0    ! Crystal orientation [0-2]
  ! real(kind=fPrec), private, save :: c_rcurv    = zero ! Crystal geometrical parameters [m]
  ! real(kind=fPrec), private, save :: c_xmax     = zero ! Crystal geometrical parameters [m]
  ! real(kind=fPrec), private, save :: c_ymax     = zero ! Crystal geometrical parameters [m]
  ! real(kind=fPrec), private, save :: c_alayer   = zero ! Crystal amorphous layer [mm]
  ! real(kind=fPrec), private, save :: c_miscut   = zero ! Crystal miscut angle in rad
  ! real(kind=fPrec), private, save :: c_cpTilt   = zero ! Cosine of positive crystal tilt
  ! real(kind=fPrec), private, save :: c_spTilt   = zero ! Sine of positive crystal tilt
  ! real(kind=fPrec), private, save :: c_cnTilt   = zero ! Cosine of negative crystal tilt
  ! real(kind=fPrec), private, save :: c_snTilt   = zero ! Sine of negative crystal tilt
  ! real(kind=fPrec), private, save :: c_cBend    = zero ! Cosine of crystal bend
  ! real(kind=fPrec), private, save :: c_sBend    = zero ! Sine of crystal bend
  ! real(kind=fPrec), private, save :: cry_tilt   = zero ! Crystal tilt angle in rad
  ! real(kind=fPrec), private, save :: cry_length = zero ! Crystal length [m]
  ! real(kind=fPrec), private, save :: cry_bend   = zero ! Crystal bending angle in rad

  ! Rutherford Scatter
  real(kind=fPrec), parameter     :: tlcut_cry = 0.0009982_fPrec
  real(kind=fPrec), private, save :: cgen_cry(200,1)
  real(kind=fPrec), private, save :: emr_curr_cry
  real(kind=fPrec), private, save :: zatom_curr_cry

  real(kind=fPrec), parameter :: aTF = 0.194e-10_fPrec ! Screening function [m]
  real(kind=fPrec), parameter :: dP  = 1.920e-10_fPrec ! Distance between planes (110) [m]
  real(kind=fPrec), parameter :: u1  = 0.075e-10_fPrec ! Thermal vibrations amplitude

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref_cry = 0.040_fPrec
  real(kind=fPrec), parameter :: freeco_cry = 1.618_fPrec

  ! Processes
  integer, parameter :: proc_out         =  -1     ! Crystal not hit
  integer, parameter :: proc_AM          =   1     ! Amorphous
  integer, parameter :: proc_VR          =   2     ! Volume reflection
  integer, parameter :: proc_CH          =   3     ! Channeling
  integer, parameter :: proc_VC          =   4     ! Volume capture
  integer, parameter :: proc_absorbed    =   5     ! Absorption
  integer, parameter :: proc_DC          =   6     ! Dechanneling
  integer, parameter :: proc_pne         =   7     ! Proton-neutron elastic interaction
  integer, parameter :: proc_ppe         =   8     ! Proton-proton elastic interaction
  integer, parameter :: proc_diff        =   9     ! Single diffractive
  integer, parameter :: proc_ruth        =  10     ! Rutherford scattering
  integer, parameter :: proc_ch_absorbed =  15     ! Channeling followed by absorption
  integer, parameter :: proc_ch_pne      =  17     ! Channeling followed by proton-neutron elastic interaction
  integer, parameter :: proc_ch_ppe      =  18     ! Channeling followed by proton-proton elastic interaction
  integer, parameter :: proc_ch_diff     =  19     ! Channeling followed by single diffractive
  integer, parameter :: proc_ch_ruth     =  20     ! Channeling followed by Rutherford scattering
  integer, parameter :: proc_TRVR        = 100     ! Volume reflection in VR-AM transition region
  integer, parameter :: proc_TRAM        = 101     ! Amorphous in VR-AM transition region

contains


! ! ================================================================================================ !
! !  Initialise a crystal collimator
! ! ================================================================================================ !
! subroutine cry_startElement(c_length, new_length, cryTilt, cryBend, cryThick, cryXDim, &
!                             cryYDim, cryOrient, cryMiscut)
!   use mathlib_bouncer

!   real(kind=fPrec), intent(in)    :: c_length     ! Collimator length in m
!   real(kind=fPrec), intent(inout) :: new_length
!   real(kind=fPrec), intent(in)    :: cryTilt
!   real(kind=fPrec), intent(in)    :: cryBend
!   real(kind=fPrec), intent(in)    :: cryThick
!   real(kind=fPrec), intent(in)    :: cryXDim
!   real(kind=fPrec), intent(in)    :: cryYDim
!   integer,          intent(in)    :: cryOrient
!   real(kind=fPrec), intent(in)    :: cryMiscut

!   real(kind=fPrec) bendAng
  
!   cry_tilt = cryTilt       ! global variable

!   bendAng  = c_length/cryBend     ! cryBend is bending radius
!   if(cry_tilt >= (-one)*bendAng) then
!     cry_length = cryBend*(sin_mb(bendAng + cry_tilt) - sin_mb(cry_tilt))
!   else
!     cry_length = cryBend*(sin_mb(bendAng - cry_tilt) + sin_mb(cry_tilt))
!   end if

!   c_rcurv  = cryBend
!   c_alayer = cryThick
!   c_xmax   = cryXDim
!   c_ymax   = cryYDim
!   c_orient = cryOrient
!   c_miscut = cryMiscut
!   cry_bend = cry_length/c_rcurv
!   c_cBend  = cos_mb(cry_bend)
!   c_sBend  = sin_mb(cry_bend)
!   c_cpTilt = cos_mb(cry_tilt)
!   c_spTilt = sin_mb(cry_tilt)
!   c_cnTilt = cos_mb(-cry_tilt)
!   c_snTilt = sin_mb(-cry_tilt)

!   n_chan      = 0
!   n_VR        = 0
!   n_amorphous = 0

!   new_length = cry_length


! end subroutine cry_startElement

! ================================================================================================ !
!  Subroutine for checking for interactions with crystal
! ================================================================================================ !
subroutine cry_doCrystal(cd_x,cd_xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs,lhit,part_abs,impact,indiv,&
                        c_length,cd_exenergy,cd_rho,cd_anuc,cd_zatom,cd_emr,cd_dlri,cd_dlyi,cd_ai,cd_eUm,cd_collnt,&
                        cd_hcut,cd_bnref,cd_csref0,cd_csref1,cd_csref4,cd_csref5,cd_csect, cry_tilt, c_rcurv, &
                        c_alayer, c_xmax, c_ymax, c_orient, c_miscut, cry_bend, c_cBend, c_sBend, c_cpTilt, &
                        c_spTilt, c_cnTilt, c_snTilt, iProc, n_chan, n_VR, n_amorphous)
                 

  use parpro
  use coll_common, only : cry_proc, cry_proc_prev, cry_proc_tmp
  use mathlib_bouncer

  ! integer,          intent(in)    :: mat

  real(kind=fPrec), intent(inout) :: cd_x,cd_xp
  real(kind=fPrec), intent(inout) :: z,zp
  real(kind=fPrec), intent(inout) :: s,p
  real(kind=fPrec), intent(inout) :: x0,xp0
  real(kind=fPrec), intent(inout) :: zlm,s_imp
  integer,          intent(inout) :: nhit,nabs
  integer,          intent(inout) :: lhit
  integer,          intent(inout) :: part_abs
  real(kind=fPrec), intent(inout) :: impact
  real(kind=fPrec), intent(inout) :: indiv
  real(kind=fPrec), intent(in)    :: c_length

  real(kind=fPrec), intent(in)    :: cd_exenergy
  real(kind=fPrec), intent(in)    :: cd_rho
  real(kind=fPrec), intent(in)    :: cd_anuc
  real(kind=fPrec), intent(in)    :: cd_zatom
  real(kind=fPrec), intent(in)    :: cd_emr

  real(kind=fPrec), intent(in)    :: cd_dlri
  real(kind=fPrec), intent(in)    :: cd_dlyi
  real(kind=fPrec), intent(in)    :: cd_ai
  real(kind=fPrec), intent(in)    :: cd_eUm
  real(kind=fPrec), intent(in)    :: cd_collnt

  real(kind=fPrec), intent(in)    :: cd_hcut
  real(kind=fPrec), intent(in)    :: cd_bnref
  real(kind=fPrec), intent(in)    :: cd_csect

  real(kind=fPrec), intent(in)    :: cd_csref0
  real(kind=fPrec), intent(in)    :: cd_csref1
  real(kind=fPrec), intent(in)    :: cd_csref4
  real(kind=fPrec), intent(in)    :: cd_csref5

  real(kind=fPrec), intent(in)    :: cry_tilt
  real(kind=fPrec), intent(in)    :: c_rcurv
  real(kind=fPrec), intent(in)    :: c_alayer
  real(kind=fPrec), intent(in)    :: c_xmax
  real(kind=fPrec), intent(in)    :: c_ymax
  integer,          intent(in)    :: c_orient
  real(kind=fPrec), intent(in)    :: c_miscut
  real(kind=fPrec), intent(in)    :: cry_bend
  real(kind=fPrec), intent(in)    :: c_cBend
  real(kind=fPrec), intent(in)    :: c_sBend
  real(kind=fPrec), intent(in)    :: c_cpTilt
  real(kind=fPrec), intent(in)    :: c_spTilt
  real(kind=fPrec), intent(in)    :: c_cnTilt
  real(kind=fPrec), intent(in)    :: c_snTilt
  integer,          intent(inout) :: iProc
  integer,          intent(inout) :: n_chan
  integer,          intent(inout) :: n_VR
  integer,          intent(inout) :: n_amorphous

  logical,          intent(inout) :: isImp

  real(kind=fPrec) s_temp,s_shift,s_rot,s_int
  real(kind=fPrec) x_temp,x_shift,x_rot,x_int
  real(kind=fPrec) xp_temp,xp_shift,xp_rot,xp_int,xp_tangent
  real(kind=fPrec) tilt_int,shift,delta,a_eq,b_eq,c_eq
  real(kind=fPrec) s_P, x_P, s_P_tmp, x_P_tmp
  real(kind=fPrec) cry_length
  
  cry_length = c_length ! hack

  ! PYTHON
  ! cry_tilt = cryTilt       ! global variable

  ! bendAng  = c_length/cryBend     ! cryBend is bending radius
  ! if(cry_tilt >= (-one)*bendAng) then
  !   cry_length = cryBend*(sin_mb(bendAng + cry_tilt) - sin_mb(cry_tilt))
  ! else
  !   cry_length = cryBend*(sin_mb(bendAng - cry_tilt) + sin_mb(cry_tilt))
  ! end if

  ! c_rcurv  = cryBend
  ! c_alayer = cryThick
  ! c_xmax   = cryXDim
  ! c_ymax   = cryYDim
  ! c_orient = cryOrient
  ! c_miscut = cryMiscut
  ! cry_bend = cry_length/c_rcurv
  ! c_cBend  = cos_mb(cry_bend)
  ! c_sBend  = sin_mb(cry_bend)
  ! c_cpTilt = cos_mb(cry_tilt)
  ! c_spTilt = sin_mb(cry_tilt)
  ! c_cnTilt = cos_mb(-cry_tilt)
  ! c_snTilt = sin_mb(-cry_tilt)

  ! iProc       = 0
  ! n_chan      = 0
  ! n_VR        = 0
  ! n_amorphous = 0

  ! new_length = cry_length

  s_temp     = zero
  s_shift    = zero
  s_rot      = zero
  s_int      = zero
  x_temp     = zero
  x_shift    = zero
  x_rot      = zero
  x_int      = zero
  xp_temp    = zero
  xp_shift   = zero
  xp_rot     = zero
  xp_int     = zero
  xp_tangent = zero
  tilt_int   = zero
  shift      = zero
  delta      = zero
  a_eq       = zero
  b_eq       = zero
  c_eq       = zero
  s_imp      = zero

  iProc       = proc_out
  
  ! Transform in the crystal reference system
  ! 1st transformation: shift of the center of the reference frame
  if(cry_tilt < zero) then
    s_shift = s
    shift   = c_rcurv*(one - c_cpTilt)
    if(cry_tilt < -cry_bend) then
      shift = c_rcurv*(c_cnTilt - cos_mb(cry_bend - cry_tilt))
    end if
    x_shift = cd_x - shift
  else
    s_shift = s
    x_shift = cd_x
  end if
  
  ! 2nd transformation: rotation
  s_rot  = x_shift*c_spTilt + s_shift*c_cpTilt
  x_rot  = x_shift*c_cpTilt - s_shift*c_spTilt
  xp_rot = cd_xp - cry_tilt

  ! 3rd transformation: drift to the new coordinate s=0
  cd_xp = xp_rot
  cd_x  = x_rot - xp_rot*s_rot
  z  = z - zp*s_rot
  s  = zero
  
  ! Check that particle hit the crystal
  if(cd_x >= zero .and. cd_x < c_xmax) then
    ! MISCUT first step: P coordinates (center of curvature of crystalline planes)
    s_P = (c_rcurv-c_xmax)*sin_mb(-c_miscut)
    x_P = c_xmax + (c_rcurv-c_xmax)*cos_mb(-c_miscut)

    call cry_interact(cd_x,cd_xp,z,zp,p,cry_length,s_P,x_P,cd_exenergy,cd_rho,cd_anuc,cd_zatom,cd_emr,cd_dlri,cd_dlyi,&
                      cd_ai,cd_eUm,cd_collnt,cd_hcut,cd_csref0,cd_csref1,cd_csref4,cd_csref5,cd_bnref,cd_csect,cry_tilt,&
                      c_rcurv,c_alayer,c_xmax,c_ymax,c_orient,c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,c_spTilt,&
                      c_cnTilt,c_snTilt, iProc)

    s   = c_rcurv*c_sBend
    zlm = c_rcurv*c_sBend
    if(iProc /= proc_out) then
      isImp    = .true.
      nhit     = nhit + 1
      lhit     = 1
      impact   = x0
      indiv    = xp0
    end if

  else

    if(cd_x < zero) then ! Crystal can be hit from below
      xp_tangent = sqrt((-(two*cd_x)*c_rcurv + cd_x**2)/(c_rcurv**2))
    else              ! Crystal can be hit from above
      xp_tangent = asin_mb((c_rcurv*(one - c_cBend) - cd_x)/sqrt(((two*c_rcurv)*(c_rcurv - cd_x))*(one - c_cBend) + cd_x**2))
    end if

    ! If the hit is below,the angle must be greater or equal than the tangent,
    ! or if the hit is above, the angle must be smaller or equal than the tangent
    if((cd_x < zero .and. cd_xp >= xp_tangent) .or. (cd_x >= zero .and. cd_xp <= xp_tangent)) then

      ! If it hits the crystal, calculate in which point and apply the transformation and drift to that point
      a_eq  = one + cd_xp**2
      b_eq  = (two*cd_xp)*(cd_x - c_rcurv)
      c_eq  = -(two*cd_x)*c_rcurv + cd_x**2
      delta = b_eq**2 - four*(a_eq*c_eq)
      s_int = (-b_eq - sqrt(delta))/(two*a_eq)
      s_imp = s_int

      ! MISCUT first step: P coordinates (center of curvature of crystalline planes)
      s_P_tmp = (c_rcurv-c_xmax)*sin_mb(-c_miscut)
      x_P_tmp = c_xmax + (c_rcurv-c_xmax)*cos_mb(-c_miscut)

      if(s_int < c_rcurv*c_sBend) then
        ! Transform to a new reference system: shift and rotate
        x_int  = cd_xp*s_int + cd_x
        xp_int = cd_xp
        z      = z + zp*s_int
        cd_x   = zero
        s      = zero

        tilt_int = s_int/c_rcurv
        cd_xp    = cd_xp-tilt_int

        ! MISCUT first step (bis): transform P in new reference system
        ! Translation
        s_P_tmp = s_P_tmp - s_int
        x_P_tmp = x_P_tmp - x_int
        ! Rotation
        s_P = s_P_tmp*cos_mb(tilt_int) + x_P_tmp*sin_mb(tilt_int)
        x_P = -s_P_tmp*sin_mb(tilt_int) + x_P_tmp*cos_mb(tilt_int)

        call cry_interact(cd_x,cd_xp,z,zp,p,cry_length-(tilt_int*c_rcurv),s_P,x_P,cd_exenergy,cd_rho,cd_anuc,&
                          cd_zatom,cd_emr,cd_dlri,cd_dlyi,cd_ai,cd_eUm,cd_collnt,cd_hcut,cd_csref0,cd_csref1,&
                          cd_csref4,cd_csref5,cd_bnref,cd_csect,cry_tilt,c_rcurv,c_alayer,c_xmax,c_ymax,c_orient,&
                          c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,c_spTilt,c_cnTilt,c_snTilt,iProc)
        s   = c_rcurv*sin_mb(cry_bend - tilt_int)
        zlm = c_rcurv*sin_mb(cry_bend - tilt_int)
        if(iProc /= proc_out) then
          x_rot    = x_int
          s_rot    = s_int
          xp_rot   = xp_int
          s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
          x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
          xp_shift = xp_rot + cry_tilt

          if(cry_tilt < zero) then
            x0  = x_shift + shift
            xp0 = xp_shift
          else
            x0  = x_shift
            xp0 = xp_shift
          end if

          isImp     = .true.
          nhit      = nhit + 1
          lhit      = 1
          impact    = x0
          indiv     = xp0
        end if

        ! un-rotate
        x_temp  = cd_x
        s_temp  = s
        xp_temp = cd_xp
        s       =  s_temp*cos_mb(-tilt_int) + x_temp*sin_mb(-tilt_int)
        cd_x    = -s_temp*sin_mb(-tilt_int) + x_temp*cos_mb(-tilt_int)
        cd_xp   = xp_temp + tilt_int

        ! 2nd: shift back the 2 axis
        cd_x = cd_x + x_int
        s = s + s_int

      else

        s = c_rcurv*sin_mb(cry_length/c_rcurv)
        cd_x = cd_x + s*cd_xp
        z = z + s*zp

      end if

    else

      s = c_rcurv*sin_mb(cry_length/c_rcurv)
      cd_x = cd_x + s*cd_xp
      z = z + s*zp

    end if

  end if

  ! transform back from the crystal to the collimator reference system
  ! 1st: un-rotate the coordinates
  x_rot  = cd_x
  s_rot  = s
  xp_rot = cd_xp

  s_shift  =  s_rot*c_cnTilt + x_rot*c_snTilt
  x_shift  = -s_rot*c_snTilt + x_rot*c_cnTilt
  xp_shift = xp_rot + cry_tilt

  ! 2nd: shift back the reference frame
  if(cry_tilt < zero) then
    s  = s_shift
    cd_x  = x_shift + shift
    cd_xp = xp_shift
  else
    cd_x  = x_shift
    s  = s_shift
    cd_xp = xp_shift
  end if

  ! 3rd: shift to new S=Length position
  cd_x = cd_xp*(c_length - s) + cd_x
  z = zp*(c_length - s) + z
  s = c_length

  nabs = 0
  
  if(iProc == proc_AM) then
    n_amorphous = n_amorphous + 1
  else if(iProc == proc_VR) then
    n_VR = n_VR + 1
  else if(iProc == proc_CH) then
    n_chan = n_Chan + 1
  else if(iProc == proc_absorbed) then
    nabs = 1   ! TODO: do we need to set part_abs_pos etc?
  else if(iProc == proc_ch_absorbed) then
    nabs = 1
  end if

  ! Storing the process ID for the next interaction

end subroutine cry_doCrystal

! ================================================================================================ !
!  Subroutine for the movements of the particles in the crystal
!  Simple tranport protons in crystal 2
! ================================================================================================ !
subroutine cry_interact(ci_x,xp,y,yp,pc,length,s_P,x_P,ci_exenergy,ci_rho,ci_anuc,ci_zatom,ci_emr,&
                      ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt,ci_hcut,ci_csref0,ci_csref1,ci_csref4,&
                      ci_csref5,ci_bnref,ci_csect,ci_tilt,ci_rcurv,ci_alayer,ci_xmax,ci_ymax,ci_orient,&
                      ci_miscut,ci_bend,ci_cBend,ci_sBend,ci_cpTilt,ci_spTilt,ci_cnTilt,ci_snTilt,ci_iProc)

  use mod_ranlux
  use mod_funlux
  use mod_common_main
  use floatPrecision
  ! use coll_materials, only : zatom, exenergy, rho, anuc
  use mathlib_bouncer
  use physical_constants

  ! integer,          intent(in)    :: is  ! Material number
  real(kind=fPrec), intent(inout) :: ci_x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: y
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: length
  real(kind=fPrec), intent(in)    :: s_P
  real(kind=fPrec), intent(in)    :: x_P

  real(kind=fPrec), intent(in)    :: ci_exenergy
  real(kind=fPrec), intent(in)    :: ci_rho
  real(kind=fPrec), intent(in)    :: ci_anuc
  real(kind=fPrec), intent(in)    :: ci_zatom
  real(kind=fPrec), intent(in)    :: ci_emr
  
  real(kind=fPrec), intent(in)    :: ci_dlri
  real(kind=fPrec), intent(in)    :: ci_dlyi
  real(kind=fPrec), intent(in)    :: ci_ai
  real(kind=fPrec), intent(in)    :: ci_eUm
  real(kind=fPrec), intent(in)    :: ci_collnt
  real(kind=fPrec), intent(in)    :: ci_hcut
  real(kind=fPrec), intent(in)    :: ci_csref0
  real(kind=fPrec), intent(in)    :: ci_csref1
  real(kind=fPrec), intent(in)    :: ci_csref4
  real(kind=fPrec), intent(in)    :: ci_csref5

  real(kind=fPrec), intent(in)    :: ci_bnref
  real(kind=fPrec), intent(in)    :: ci_csect

  real(kind=fPrec), intent(in)    :: ci_tilt
  real(kind=fPrec), intent(in)    :: ci_rcurv
  real(kind=fPrec), intent(in)    :: ci_alayer
  real(kind=fPrec), intent(in)    :: ci_xmax
  real(kind=fPrec), intent(in)    :: ci_ymax
  integer,          intent(in)    :: ci_orient
  real(kind=fPrec), intent(in)    :: ci_miscut
  real(kind=fPrec), intent(in)    :: ci_bend
  real(kind=fPrec), intent(in)    :: ci_cBend
  real(kind=fPrec), intent(in)    :: ci_sBend
  real(kind=fPrec), intent(in)    :: ci_cpTilt
  real(kind=fPrec), intent(in)    :: ci_spTilt
  real(kind=fPrec), intent(in)    :: ci_cnTilt
  real(kind=fPrec), intent(in)    :: ci_snTilt

  integer,          intent(inout) :: ci_iProc

  integer nam,zn                        ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  real(kind=fPrec) ymax,ymin            ! Crystal geometrical parameters
  real(kind=fPrec) s_length             ! Element length along s
  real(kind=fPrec) DESt                 ! Changed energy loss by ionization now calculated and not tabulated
  real(kind=fPrec) x0,y0                ! Coordinates of the particle [m]
  real(kind=fPrec) s                    ! Long coordinates of the particle [m]
  real(kind=fPrec) a_eq,b_eq,c_eq,delta ! Second order equation param.
  real(kind=fPrec) Ang_rms, Ang_avr     ! Volume reflection mean angle [rad]
  real(kind=fPrec) Dechan               ! Probability for dechanneling
  real(kind=fPrec) Lrefl,Srefl          ! Distance of the reflection point [m]
  real(kind=fPrec) Vcapt                ! Volume capture probability
  real(kind=fPrec) Chann                ! Channeling probability
  real(kind=fPrec) N_atom               ! Probability for entering channeling near atomic planes
  real(kind=fPrec) Dxp                  ! Variation in angle
  real(kind=fPrec) xpcrit               ! Critical angle for curved crystal[rad]
  real(kind=fPrec) xpcrit0              ! Critical angle for str. crystal [rad]
  real(kind=fPrec) Rcrit                ! Critical curvature radius [m]
  real(kind=fPrec) ratio                ! X=c_rcurv/Rcrit
  real(kind=fPrec) TLdech1,TLdech2      ! Typical dechanneling length(2) [m]
  real(kind=fPrec) tdech,Ldech,Sdech    ! Angle, lenght, and S coordinate of dechanneling point
  real(kind=fPrec) Rlength,Red_S        ! Reduced length/s coordinate (in case of dechanneling)
  real(kind=fPrec) am_len               ! Amorphous length
  real(kind=fPrec) len_xs,len_ys        ! Amorphous length
  real(kind=fPrec) xp_rel               ! Xp-c_miscut angle in mrad
  real(kind=fPrec) alpha                ! Par for new chann prob
  real(kind=fPrec) Pvr                  ! Prob for VR->AM transition

  ! Quantities for length and deflection calculation
  real(kind=fPrec) const_dech,xpin,ypin,tchan,tdefl,tP,L_chan,mep
  real(kind=fPrec) s_K,x_K,s_M,x_M,s_F,x_F,r,a
  real(kind=fPrec) A_F,B_F,C_F,alpha_F,beta_F

  real(kind=fPrec), parameter :: c_v1 =  1.7_fPrec ! Fitting coefficient
  real(kind=fPrec), parameter :: c_v2 = -1.5_fPrec ! Fitting coefficient

  real(kind=fPrec) enr
  real(kind=fPrec) mom
  real(kind=fPrec) betar
  real(kind=fPrec) bgr
  real(kind=fPrec) gammar
  real(kind=fPrec) tmax
  real(kind=fPrec) plen

  nam = 1 ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
  zn  = 1

  ! dE/dX and dechanneling length calculation
  mom    = pc*c1e3                ! [GeV]
  enr    = sqrt(mom**2 + pmap**2) ! [MeV]
  gammar = enr/pmap
  betar  = mom/enr
  bgr = betar*gammar
  mep    = pmae/pmap  ! Electron/proton

  tmax = (two*pmae*bgr**2)/(one + two*gammar*mep + mep**2)  ! [MeV]
  plen = sqrt((ci_rho*ci_zatom)/ci_anuc)*28.816e-6_fPrec ! [MeV]

  const_dech = ((256.0_fPrec/(nine*pi**2)) * &
    (one/(log_mb(((two*pmae)*gammar)/(ci_exenergy*c1e3)) - one))) * ((aTF*dP)/(crade*pmae)) ! [m/MeV]
  const_dech = const_dech*c1e3 ! [m/GeV]

  s        = zero
  s_length = ci_rcurv*sin_mb(length/ci_rcurv)
  L_chan   = length

  ! MISCUT second step: fundamental coordinates (crystal edges and plane curvature radius)
  s_K = ci_rcurv*sin_mb(length/ci_rcurv)
  x_K = ci_rcurv*(1-cos_mb(length/ci_rcurv))
  s_M = (ci_rcurv-ci_xmax)*sin_mb(length/ci_rcurv)
  x_M = ci_xmax + (ci_rcurv-ci_xmax)*(1-cos_mb(length/ci_rcurv))
  r   = sqrt(s_P**2 + (ci_x-x_P)**2)

  ! MISCUT third step: F coordinates (exit point) on crystal exit face
  A_F = (tan_mb(length/ci_rcurv))**2 + one
  B_F = ((-two)*(tan_mb(length/ci_rcurv))**2)*ci_rcurv + (two*tan_mb(length/ci_rcurv))*s_P - two*x_P
  C_F = ((tan_mb(length/ci_rcurv))**2)*(ci_rcurv**2) - ((two*tan_mb(length/ci_rcurv))*s_P)*ci_rcurv + s_P**2 + x_P**2 - r**2

  x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
  s_F = (-tan_mb(length/ci_rcurv))*(x_F-ci_rcurv)

  if(x_F >= x_K .and. x_F <= x_M .and. s_F >= s_M .and. s_F <= s_K) then
    ! No additional steps required for miscut
  else if (ci_miscut == 0 .and. abs(x_F-x_K) <= c1m13 .and. abs(s_F-s_K) <= c1m13) then
    ! no miscut, entrance from below: exit point is K (lower edge)
    x_F = x_K
    s_F = s_K
  else if (ci_miscut == 0 .and. abs(x_F-x_M) <= c1m13 .and. abs(s_F-s_M) <= c1m13) then
    ! no miscut, entrance from above: exit point is M (upper edge)
    x_F = x_M
    s_F = s_M
  else

    ! MISCUT Third step (bis): F coordinates (exit point)  on bent side
    if(ci_miscut < 0) then
      ! Intersect with bottom side
      alpha_F = (ci_rcurv-x_P)/x_P
      beta_F = -(r**2-s_P**2-x_P**2)/(two*s_P)
      A_F = alpha_F**2 + one
      B_F = two*(alpha_F*beta_F) - two*ci_rcurv
      C_F = beta_F**2
    else
      ! Intersect with top side
      alpha_F = (ci_rcurv-x_P)/s_P
      beta_F = -(r**2+ci_xmax*(ci_xmax-(two*ci_rcurv))-s_P**2-x_P**2)/(two*s_P)
      A_F = alpha_F**2 + one
      B_F = two*(alpha_F*beta_F) - two*ci_rcurv
      C_F = beta_F**2 - ci_xmax*(ci_xmax-two*ci_rcurv)
    endif
    
    x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
    s_F = alpha_F*x_F + beta_F
    
  endif

  ! MISCUT fourth step: deflection and length calculation
  a = sqrt(s_F**2+(ci_x-x_F)**2)
  tP = acos_mb((2*(r**2)-a**2)/(2*(r**2)))
  tdefl = asin_mb((s_f-s_P)/r)
  L_chan = r*tP

  xp_rel = xp - ci_miscut

  ymin = -ci_ymax/two
  ymax =  ci_ymax/two

  ! FIRST CASE: p don't interact with crystal
  if(y < ymin .or. y > ymax .or. ci_x > ci_xmax) then
    ci_x  = ci_x + xp*s_length
    y     = y + yp*s_length
    ci_iProc = proc_out
    return

  ! SECOND CASE: p hits the amorphous layer
  else if(ci_x < ci_alayer .or. y-ymin < ci_alayer .or. ymax-y < ci_alayer) then
    x0    = ci_x
    y0    = y
    a_eq  = one + xp**2
    b_eq  = (two*ci_x)*xp - (two*xp)*ci_rcurv
    c_eq  = ci_x**2 - (two*ci_x)*ci_rcurv
    delta = b_eq**2 - (four*a_eq)*c_eq
    s     = (-b_eq+sqrt(delta))/(two*a_eq)
    if(s >= s_length) then
      s = s_length
    end if
    ci_x   =  xp*s + x0
    len_xs = sqrt((ci_x-x0)**2 + s**2)
    if(yp >= zero .and. y + yp*s <= ymax) then
      len_ys = yp*len_xs
    else if(yp < zero .and. y + yp*s >= ymin) then
      len_ys = yp*len_xs
    else
      s      = (ymax-y)/yp
      len_ys = sqrt((ymax-y)**2 + s**2)
      ci_x   = x0 + xp*s
      len_xs = sqrt((ci_x-x0)**2 + s**2)
    end if
    am_len = sqrt(len_xs**2 + len_ys**2)
    s     = s/two
    ci_x  = x0 + xp*s
    y     = y0 + yp*s
    ci_iProc = proc_AM
    call cry_calcIonLoss(pc,am_len,dest,betar,bgr,gammar,tmax,plen,&
                         ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    call cry_moveAM(nam,am_len,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
                    ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    ci_x = ci_x + xp*(s_length-s)
    y = y + yp*(s_length-s)
    return

  else if(ci_x > ci_xmax-ci_alayer .and. ci_x < ci_xmax) then
    ci_iProc = proc_AM
    call cry_calcIonLoss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
                        ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    call cry_moveAM(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
                   ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    return

  end if

  ! THIRD CASE: the p interacts with the crystal.
  ! Define typical angles/probabilities for orientation 110
  xpcrit0 = sqrt((c2m9*ci_eUm)/pc)    ! Critical angle (rad) for straight crystals
  Rcrit   = (pc/(c2m6*ci_eUm))*ci_ai ! Critical curvature radius [m]

  ! If R>Rcritical=>no channeling is possible (ratio<1)
  ratio  = ci_rcurv/Rcrit
  xpcrit = (xpcrit0*(ci_rcurv-Rcrit))/ci_rcurv ! Critical angle for curved crystal

  if(ratio <= one) then ! no possibile channeling
    Ang_rms = ((c_v1*0.42_fPrec)*xpcrit0)*sin_mb(1.4_fPrec*ratio) ! RMS scattering
    Ang_avr = ((c_v2*xpcrit0)*c5m2)*ratio                         ! Average angle reflection
    Vcapt   = zero                                                ! Probability of VC

  else if(ratio <= three) then ! Strongly bent crystal
    Ang_rms = ((c_v1*0.42_fPrec)*xpcrit0)*sin_mb(0.4713_fPrec*ratio + 0.85_fPrec) ! RMS scattering
    Ang_avr = (c_v2*xpcrit0)*(0.1972_fPrec*ratio - 0.1472_fPrec)                  ! Average angle reflection
    Vcapt   = 7.0e-4_fPrec*(ratio - 0.7_fPrec)/pc**c2m1                           ! Correction by sasha drozdin/armen
    ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

  else ! Rcry >> Rcrit
    Ang_rms = (c_v1*xpcrit0)*(one/ratio)                ! RMS scattering
    Ang_avr = (c_v2*xpcrit0)*(one - 1.6667_fPrec/ratio) ! Average angle for VR
    Vcapt   = 7.0e-4_fPrec*(ratio - 0.7_fPrec)/pc**c2m1 ! Probability for VC correction by sasha drozdin/armen
    ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)

  end if

  if(ci_orient == 2) then
    Ang_avr = Ang_avr*0.93_fPrec
    Ang_rms = Ang_rms*1.05_fPrec
    xpcrit  = xpcrit*0.98_fPrec
  end if

  if(abs(xp_rel) < xpcrit) then
    alpha  = xp_rel/xpcrit
    Chann  = sqrt(0.9_fPrec*(one - alpha**2))*sqrt(one-(one/ratio)) ! Saturation at 95%
    N_atom = c1m1

    ! if they can channel: 2 options
    if(coll_rand() <= chann) then ! option 1:channeling

      TLdech1 = (const_dech*pc)*(one-one/ratio)**2 ! Updated calculate typical dech. length(m)
      if(coll_rand() <= n_atom) then
        TLdech1 = ((const_dech/c2e2)*pc)*(one-one/ratio)**2  ! Updated dechanneling length (m)
      end if

      Dechan = -log_mb(coll_rand()) ! Probability of dechanneling
      Ldech  = TLdech1*Dechan   ! Actual dechan. length

      ! careful: the dechanneling lentgh is along the trajectory
      ! of the particle -not along the longitudinal coordinate...
      if(ldech < l_chan) then
        ci_iProc = proc_DC
        Dxp   = Ldech/r ! Change angle from channeling [mrad]
        Sdech = Ldech*cos_mb(ci_miscut + half*Dxp)
        ci_x  = ci_x  + Ldech*(sin_mb(half*Dxp+ci_miscut)) ! Trajectory at channeling exit
        xp    = xp + Dxp + (two*(coll_rand()-half))*xpcrit
        y     = y  + yp * Sdech

        call cry_calcIonLoss(pc,ldech,dest,betar,bgr,gammar,tmax,plen,&
                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
        pc = pc - half*dest*Ldech ! Energy loss to ionization while in CH [GeV]
        ci_x  = ci_x  + (half*(s_length-Sdech))*xp
        y  = y  + (half*(s_length-Sdech))*yp

        call cry_calcIonLoss(pc,s_length-sdech,dest,betar,bgr,gammar,tmax,plen,&
                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
        call cry_moveAM(nam,s_length-sdech,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
                        ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
        ci_x = ci_x + (half*(s_length-Sdech))*xp
        y = y + (half*(s_length-Sdech))*yp
      else
        ci_iProc = proc_CH
        xpin  = XP
        ypin  = YP

        ! check if a nuclear interaction happen while in CH
        call cry_moveCH(nam,l_chan,ci_x,xp,yp,pc,ci_rcurv,rcrit,ci_rho,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csect,&
                        ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_eUm,ci_collnt,ci_iProc)
        if(ci_iProc /= proc_CH) then
          ! if an nuclear interaction happened, move until the middle with initial xp,yp then
          ! propagate until the "crystal exit" with the new xp,yp accordingly with the rest
          ! of the code in "thin lens approx"
          ci_x = ci_x + (half*L_chan)*xpin
          y = y + (half*L_chan)*ypin
          ci_x = ci_x + (half*L_chan)*XP
          y = y + (half*L_chan)*YP

          call cry_calcIonLoss(pc,length,dest,betar,bgr,gammar,tmax,plen,&
                                ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          pc = pc - dest*length ! energy loss to ionization [GeV]
        else
          Dxp = tdefl + (half*ran_gauss(zero))*xpcrit ! Change angle[rad]
        
          xp  = Dxp
          ci_x = ci_x + L_chan*(sin_mb(half*Dxp)) ! Trajectory at channeling exit
          y   = y + s_length * yp

          call cry_calcIonLoss(pc,length,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          pc = pc - (half*dest)*length ! energy loss to ionization [GeV]
        end if
      end if

    else ! Option 2: VR

      ! good for channeling but don't channel (1-2)
      ci_iProc = proc_VR

      xp = xp + (0.45_fPrec*(xp_rel/xpcrit + one))*Ang_avr
      ci_x  = ci_x  + (half*s_length)*xp
      y  = y  + (half*s_length)*yp

      call cry_calcIonLoss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
                           ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
      call cry_moveAM(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
                      ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)

      ci_x = ci_x + (half*s_length)*xp
      y = y + (half*s_length)*yp

    end if

  else ! case 3-2: no good for channeling. check if the  can VR

    Lrefl = xp_rel*r ! Distance of refl. point [m]
    Srefl = sin_mb(xp_rel/two + ci_miscut)*Lrefl

    if(Lrefl > zero .and. Lrefl < L_chan) then ! VR point inside

      ! 2 options: volume capture and volume reflection

      if(coll_rand() > Vcapt .or. ZN == zero) then ! Option 1: VR

        ci_iProc = proc_VR
        ci_x  = ci_x + xp*Srefl
        y     = y + yp*Srefl
        Dxp   = Ang_avr
        xp    = xp + Dxp + Ang_rms*ran_gauss(zero)
        ci_x  = ci_x  + (half*xp)*(s_length - Srefl)
        y     = y  + (half*yp)*(s_length - Srefl)

        call cry_calcIonLoss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
        call cry_moveAM(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
                        ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
        ci_x = ci_x + (half*xp)*(s_length - Srefl)
        y = y + (half*yp)*(s_length - Srefl)

      else ! Option 2: VC

        ci_x = ci_x + xp*Srefl
        y = y + yp*Srefl

        TLdech2 = (const_dech/c1e1)*pc*(one-one/ratio)**2          ! Updated typical dechanneling length(m)
        Ldech   = TLdech2*(sqrt(c1m2 - log_mb(coll_rand())) - c1m1)**2 ! Updated DC length
        tdech   = Ldech/ci_rcurv
        Sdech   = Ldech*cos_mb(xp + half*tdech)

        if(Ldech < Length-Lrefl) then

          ci_iProc = proc_DC
          Dxp   = Ldech/ci_rcurv + (half*ran_gauss(zero))*xpcrit
          ci_x  = ci_x + Ldech*(sin_mb(half*Dxp+xp)) ! Trajectory at channeling exit
          y     = y + Sdech*yp
          xp    =  Dxp
          Red_S = (s_length - Srefl) - Sdech
          ci_x  = ci_x + (half*xp)*Red_S
          y     = y + (half*yp)*Red_S

          call cry_calcIonLoss(pc,srefl,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          pc = pc - dest*Srefl ! "added" energy loss before capture

          call cry_calcIonLoss(pc,sdech,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          pc = pc - (half*dest)*Sdech ! "added" energy loss while captured

          call cry_calcIonLoss(pc,red_s,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          call cry_moveAM(nam,red_s,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
                          ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
          ci_x = ci_x + (half*xp)*Red_S
          y = y + (half*yp)*Red_S

        else

          ci_iProc   = proc_VC
          Rlength = Length - Lrefl
          tchan   = Rlength/ci_rcurv
          Red_S   = Rlength*cos_mb(xp + half*tchan)

          call cry_calcIonLoss(pc,lrefl,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          pc   = pc - dest*Lrefl ! "added" energy loss before capture
          xpin = xp
          ypin = yp

          ! Check if a nuclear interaction happen while in ch
          call cry_moveCH(nam,rlength,ci_x,xp,yp,pc,ci_rcurv,rcrit,ci_rho,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csect,&
                          ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_eUm,ci_collnt,ci_iProc)
                          
          if(ci_iProc /= proc_VC) then
            ! if an nuclear interaction happened, move until the middle with initial xp,yp then propagate until
            ! the "crystal exit" with the new xp,yp aciordingly with the rest of the code in "thin lens approx"
            ci_x = ci_x + (half*Rlength)*xpin
            y = y + (half*Rlength)*ypin
            ci_x = ci_x + (half*Rlength)*XP
            y = y + (half*Rlength)*YP

            call cry_calcIonLoss(pc,rlength,dest,betar,bgr,gammar,tmax,plen,&
                                 ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
            pc = pc - dest*Rlength
          else
            Dxp = (Length-Lrefl)/ci_rcurv
            ci_x = ci_x + sin_mb(half*Dxp+xp)*Rlength ! Trajectory at channeling exit
            y   = y + red_S*yp
            xp  = tdefl + (half*ran_gauss(zero))*xpcrit ! [mrad]

            call cry_calcIonLoss(pc,rlength,dest,betar,bgr,gammar,tmax,plen,&
                                 ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
            pc = pc - (half*dest)*Rlength  ! "added" energy loss once captured
          end if
        end if
      end if

    else

      ! Case 3-3: move in amorphous substance (big input angles)
      ! Modified for transition vram daniele
      if(xp_rel > tdefl-ci_miscut + two*xpcrit .or. xp_rel < -xpcrit) then
        ci_iProc = proc_AM
        ci_x  = ci_x + (half*s_length)*xp
        y     = y + (half*s_length)*yp
        if(zn > zero) then
          call cry_calcIonLoss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          call cry_moveAM(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
                          ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
        end if
        ci_x = ci_x + (half*s_length)*xp
        y = y + (half*s_length)*yp
      else
        Pvr = (xp_rel-(tdefl-ci_miscut))/(two*xpcrit)
        if(coll_rand() > Pvr) then
          ci_iProc = proc_TRVR
          ci_x  = ci_x + xp*Srefl
          y     = y + yp*Srefl

          Dxp = (((-three*Ang_rms)*xp_rel)/(two*xpcrit) + Ang_avr) + ((three*Ang_rms)*(tdefl-ci_miscut))/(two*xpcrit)
          xp  = xp + Dxp
          ci_x = ci_x + (half*xp)*(s_length-Srefl)
          y   = y + (half*yp)*(s_length-Srefl)

          call cry_calcIonLoss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          call cry_moveAM(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
                          ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
          ci_x = ci_x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
        else
          ci_iProc = proc_TRAM
          ci_x = ci_x + xp*Srefl
          y = y + yp*Srefl
          Dxp = ((((-one*(13.6_fPrec/pc))*sqrt(s_length/ci_dlri))*c1m3)*xp_rel)/(two*xpcrit) + &
            (((13.6_fPrec/pc)*sqrt(s_length/ci_dlri))*c1m3)*(one+(tdefl-ci_miscut)/(two*xpcrit))
          xp = xp+Dxp
          ci_x  = ci_x + (half*xp)*(s_length-Srefl)
          y  = y + (half*yp)*(s_length-Srefl)

          call cry_calcIonLoss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
                              ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
          call cry_moveAM(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
                          ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
          ci_x = ci_x + (half*xp)*(s_length - Srefl)
          y = y + (half*yp)*(s_length - Srefl)
        end if
      end if
    end if
  end if

end subroutine cry_interact

! ================================================================================================ !
!  Subroutine for the calculazion of the energy loss by ionisation
! ================================================================================================ !
subroutine cry_calcIonLoss(pc,dz,EnLo,cc_betar,cc_bgr,cc_gammar,cc_tmax,cc_plen,cc_exenergy,cc_zatom,&
                           cc_rho,cc_anuc,cc_dlri,cc_dlyi,cc_ai,cc_eUm,cc_collnt)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  ! use coll_materials, only : zatom, exenergy, rho, anuc
  use mathlib_bouncer
  use physical_constants

  ! integer,          intent(in)  :: is
  real(kind=fPrec), intent(in)  :: pc
  real(kind=fPrec), intent(in)  :: dz
  real(kind=fPrec), intent(out) :: EnLo
  real(kind=fPrec), intent(in)  :: cc_betar
  real(kind=fPrec), intent(in)  :: cc_bgr
  real(kind=fPrec), intent(in)  :: cc_gammar
  real(kind=fPrec), intent(in)  :: cc_tmax
  real(kind=fPrec), intent(in)  :: cc_plen

  real(kind=fPrec), intent(in)  :: cc_exenergy
  real(kind=fPrec), intent(in)  :: cc_zatom
  real(kind=fPrec), intent(in)  :: cc_rho
  real(kind=fPrec), intent(in)  :: cc_anuc
  
  real(kind=fPrec), intent(in)  :: cc_dlri
  real(kind=fPrec), intent(in)  :: cc_dlyi
  real(kind=fPrec), intent(in)  :: cc_ai
  real(kind=fPrec), intent(in)  :: cc_eUm
  real(kind=fPrec), intent(in)  :: cc_collnt

  real(kind=fPrec) thl,tt,cs_tail,prob_tail
  real(kind=fPrec), parameter :: k = 0.307075_fPrec ! Constant in front bethe-bloch [mev g^-1 cm^2]

  thl       = (((((four*k)*cc_zatom)*dz)*c1e2)*cc_rho)/(cc_anuc*cc_betar**2) ! [MeV]
  EnLo      = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ( &
    half*log_mb(((((two*pmae)*cc_bgr)*cc_bgr)*cc_tmax)/(c1e6*cc_exenergy**2)) - &
    cc_betar**2 - log_mb(cc_plen/(cc_exenergy*c1e3)) - log_mb(cc_bgr) + half    &
  )
  EnLo      = ((EnLo*cc_rho)*c1m1)*dz ! [GeV]
  Tt        = (EnLo*c1e3)+thl          ! [MeV]

  cs_tail   = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ((half*((one/Tt)-(one/cc_tmax))) - &
    (log_mb(cc_tmax/Tt)*(cc_betar**2)/(two*cc_tmax)) + ((cc_tmax-Tt)/((four*(cc_gammar**2))*(pmap**2))))
  prob_tail = ((cs_tail*cc_rho)*dz)*c1e2

  if(coll_rand() < prob_tail) then
    EnLo = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ( &
      half*log_mb((two*pmae*cc_bgr*cc_bgr*cc_tmax)/(c1e6*cc_exenergy**2)) -      &
      cc_betar**2 - log_mb(cc_plen/(cc_exenergy*c1e3)) - log_mb(cc_bgr) + half + &
      cc_tMax**2/(eight*(cc_gammar**2)*(pmap**2)) &
    )
    EnLo = (EnLo*cc_rho)*c1m1 ! [GeV/m]
  else
    EnLo = EnLo/dz ! [GeV/m]
  end if

end subroutine cry_calcIonLoss

! ================================================================================================ !
!  Subroutine for the movement in the amorphous
! ================================================================================================ !
subroutine cry_moveAM(nam,dz,dei,dly,dlr,xp,yp,pc,cm_anuc,cm_zatom,cm_emr,cm_hcut,cm_bnref,&
                    cm_csref0,cm_csref1,cm_csref4,cm_csref5,cm_collnt,cm_iProc)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  ! use coll_materials, only : anuc, hcut, bnref, csref
  use mathlib_bouncer
  use physical_constants

  ! integer,          intent(in)    :: is
  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(in)    :: dei
  real(kind=fPrec), intent(in)    :: dly
  real(kind=fPrec), intent(in)    :: dlr
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc

  real(kind=fPrec), intent(in)    :: cm_csref0
  real(kind=fPrec), intent(in)    :: cm_csref1
  real(kind=fPrec), intent(in)    :: cm_csref4
  real(kind=fPrec), intent(in)    :: cm_csref5

  real(kind=fPrec), intent(in)    :: cm_anuc
  real(kind=fPrec), intent(in)    :: cm_zatom
  real(kind=fPrec), intent(in)    :: cm_emr
  real(kind=fPrec), intent(in)    :: cm_hcut
  real(kind=fPrec), intent(in)    :: cm_bnref
  real(kind=fPrec), intent(in)    :: cm_collnt

  integer,          intent(inout) :: cm_iProc

  integer i,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn,cs(0:5),cprob(0:5),freep,zlm,xp_in,yp_in,xm2,xln15s,tz,tx,tlow, &
    thigh,teta,pptot,ppsd,ppel,pc_in,kymcs,kxmcs,ecmsq,dya,bsd,bpp,aran

  xp_in = xp
  yp_in = yp
  pc_in = pc

  ! New treatment of scattering routine based on standard sixtrack routine
  ! useful calculations for cross-section and event topology calculation
  ecmsq  = ((two*pmap)*c1m3)*pc
  xln15s = log_mb(0.15_fPrec*ecmsq)

  ! New models, see Claudia's thesis
  pptot = (0.041084_fPrec - 0.0023302_fPrec*log_mb(ecmsq)) + 0.00031514_fPrec*log_mb(ecmsq)**2
  ppel  = (11.7_fPrec - 1.59_fPrec*log_mb(ecmsq) + 0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ppsd  = (4.3_fPrec + 0.3_fPrec*log_mb(ecmsq))/c1e3
  bpp   = 7.156_fPrec + 1.439_fPrec*log_mb(sqrt(ecmsq))

  ! Distribution for Ruth. scatt.
  tlow      = tlcut_cry
  thigh     = cm_hcut
  emr_curr_cry = cm_emr
  zatom_curr_cry = cm_zatom
  call funlxp(cry_ruth,cgen_cry(1,1),tlow,thigh)

  ! Cross-section calculation
  ! freep: number of nucleons involved in single scattering
  freep = freeco_cry * cm_anuc**(one/three)

  ! Compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3) = freep*ppel
  cs(4) = freep*ppsd

  ! Correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0) = cm_csref0 + freep*(pptot - pptref_cry)
  bn    = (cm_bnref*cs(0))/cm_csref0

  ! Also correct inel-CS
  cs(1) = (cm_csref1*cs(0))/cm_csref0

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = ((cs(0) - cs(1)) - cs(3)) - cs(4)
  cs(5) = cm_csref5

  ! Now add Coulomb
  cs(0) = cs(0) + cs(5)

  ! Calculate cumulative probability
  cprob(:) = zero
  cprob(5) = one
  do i=1,4
    cprob(i) = cprob(i-1) + cs(i)/cs(0)
  end do

  ! Multiple Coulomb Scattering
  xp  = xp*c1e3
  yp  = yp*c1e3
  pc  = pc - dei*dz ! Energy lost because of ionization process[GeV]

  dya   = (13.6_fPrec/pc)*sqrt(dz/dlr) ! RMS of coloumb scattering MCS (mrad)
  kxmcs = dya*ran_gauss(zero)
  kymcs = dya*ran_gauss(zero)

  xp = xp+kxmcs
  yp = yp+kymcs

  if(nam == 0) return ! Turn on/off nuclear interactions

  ! Can nuclear interaction happen?
  zlm = -cm_collnt*log_mb(coll_rand())

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = coll_rand()
    i=1
10  if(aran > cprob(i)) then
      i = i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction
    t = 0 ! default value to cover ichoix=1
    select case(ichoix)
    case(1) ! Deep inelastic, impinging p disappeared
      cm_iProc = proc_absorbed

    case(2) ! p-n elastic
      cm_iProc = proc_pne
      t     = -log_mb(coll_rand())/bn

    case(3) ! p-p elastic
      cm_iProc = proc_ppe
      t     = -log_mb(coll_rand())/bpp

    case(4) ! Single diffractive
      cm_iProc = proc_diff
      xm2   = exp_mb(coll_rand()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = 7.0_fPrec*bpp/12.0_fPrec
      end if
      t = -log_mb(coll_rand())/bsd

    case(5)
      cm_iProc      = proc_ruth
      length_cry = 1
      call funlux(cgen_cry(1,1),xran_cry,length_cry)
      t = xran_cry(1)

    end select

    ! Calculate the related kick
    if(ichoix == 4) then
      teta = sqrt(t)/pc_in ! DIFF has changed PC
    else
      teta = sqrt(t)/pc
    end if

    tx = (teta*ran_gauss(zero))*c1e3
    tz = (teta*ran_gauss(zero))*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz
  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine cry_moveAM

! ================================================================================================ !
!  Subroutine for check if a nuclear interaction happen while in channeling
! ================================================================================================ !
subroutine cry_moveCH(nam,dz,ch_x,xp,yp,pc,r,rc,ch_rho,ch_anuc,ch_zatom,ch_emr,ch_hcut,ch_bnref,ch_csect,&
                      ch_csref0,ch_csref1,ch_csref4,ch_csref5,ch_eUm,ch_collnt,ch_iProc)

  use crcoall
  use mod_ranlux
  use mod_funlux
  use floatPrecision
  use coll_common, only : coll_debug
  use mathlib_bouncer
  use physical_constants

  integer,          intent(in)    :: nam
  real(kind=fPrec), intent(in)    :: dz
  real(kind=fPrec), intent(inout) :: ch_x
  real(kind=fPrec), intent(inout) :: xp
  real(kind=fPrec), intent(inout) :: yp
  real(kind=fPrec), intent(inout) :: pc
  real(kind=fPrec), intent(in)    :: r
  real(kind=fPrec), intent(in)    :: rc

  real(kind=fPrec), intent(in)    :: ch_rho
  real(kind=fPrec), intent(in)    :: ch_anuc
  real(kind=fPrec), intent(in)    :: ch_zatom
  real(kind=fPrec), intent(in)    :: ch_emr
  real(kind=fPrec), intent(in)    :: ch_hcut
  real(kind=fPrec), intent(in)    :: ch_bnref

  real(kind=fPrec), intent(in)    :: ch_csref0
  real(kind=fPrec), intent(in)    :: ch_csref1
  real(kind=fPrec), intent(in)    :: ch_csref4
  real(kind=fPrec), intent(in)    :: ch_csref5

  real(kind=fPrec), intent(in)    :: ch_csect
  real(kind=fPrec), intent(in)    :: ch_eUm
  real(kind=fPrec), intent(in)    :: ch_collnt

  integer,          intent(inout) :: ch_iProc

  integer i,np,length_cry,ichoix
  real(kind=fPrec) t,xran_cry(1),bn,cs(0:5),cprob(0:5),freep,zlm,xp_in,yp_in,xminU,xm2,xln15s,x_min,&
    x_max,x_i,Umin,Ueff,tz,tx,tlow,thigh,teta,rho_min,rho_max,pv,pptot,ppsd,ppel,PC_in,nuc_cl_l,&
    N_am,Et,ecmsq,Ec,csref_inel_rsc,csref_tot_rsc,bsd,bpp,aran,avrrho

  xp_in = xp
  yp_in = yp
  pc_in = pc

  ! New treatment of scattering routine based on standard sixtrack routine

  ! Useful calculations for cross-section and event topology calculation
  ecmsq  = ((two*pmap)*c1m3)*pc
  xln15s = log_mb(0.15_fPrec*ecmsq)

  ! New models, see Claudia's thesis
  pptot = (0.041084_fPrec - 0.0023302_fPrec*log_mb(ecmsq)) + 0.00031514_fPrec*log_mb(ecmsq)**2
  ppel  = (11.7_fPrec - 1.59_fPrec*log_mb(ecmsq) + 0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ppsd  = (4.3_fPrec + 0.3_fPrec*log_mb(ecmsq))/c1e3
  bpp   = 7.156_fPrec + 1.439_fPrec*log_mb(sqrt(ecmsq))

  ! Distribution for Ruth. scatt.
  tlow      = tlcut_cry
  thigh     = ch_hcut
  emr_curr_cry = ch_emr
  zatom_curr_cry = ch_zatom
  call funlxp(cry_ruth,cgen_cry(1,1),tlow,thigh)

  ! Rescale the total and inelastic cross-section accordigly to the average density seen
  x_i = ch_x
  np  = int(x_i/dp)    ! Calculate in which crystalline plane the particle enters
  x_i = x_i - Np*dP    ! Rescale the incoming x at the left crystalline plane
  x_i = x_i - (dP/two) ! Rescale the incoming x in the middle of crystalline planes

  pv   = pc**2/sqrt(pc**2 + (pmap*c1m3)**2)*c1e9          ! Calculate pv=P/E
  Ueff = ch_eUm*((two*x_i)/dp)*((two*x_i)/dp) + pv*x_i/r ! Calculate effective potential
  Et   = (pv*xp**2)/two + Ueff                            ! Calculate transverse energy
  Ec   = (ch_eUm*(one-rc/r))*(one-rc/r)                  ! Calculate critical energy in bent crystals

  ! To avoid negative Et
  xminU = ((-dp**2*pc)*c1e9)/(eight*ch_eUm*r)
  Umin  = abs((ch_eUm*((two*xminU)/dp))*((two*xminU)/dP) + pv*xminU/R)
  Et    = Et + Umin
  Ec    = Ec + Umin

  ! Calculate min e max of the trajectory between crystalline planes
  x_min = (-(dP/two)*Rc)/R - (dP/two)*sqrt(Et/Ec)
  x_Max = (-(dP/two)*Rc)/R + (dP/two)*sqrt(Et/Ec)

  ! Change ref. frame and go back with 0 on the crystalline plane on the left
  x_min = x_min - dp/two
  x_max = x_max - dp/two

  ! Calculate the "normal density" in m^-3
  N_am  = ((ch_rho*6.022e23_fPrec)*c1e6)/ch_anuc

  ! Calculate atomic density at min and max of the trajectory oscillation
  rho_max = ((N_am*dp)/two)*(erf(x_max/sqrt(two*u1**2)) - erf((dP-x_Max)/sqrt(two*u1**2)))
  rho_min = ((N_am*dP)/two)*(erf(x_min/sqrt(two*u1**2)) - erf((dP-x_min)/sqrt(two*u1**2)))

  ! "zero-approximation" of average nuclear density seen along the trajectory
  avrrho  = (rho_max - rho_min)/(x_max - x_min)
  avrrho  = (two*avrrho)/N_am

  csref_tot_rsc  = ch_csref0*avrrho ! Rescaled total ref cs
  csref_inel_rsc = ch_csref1*avrrho ! Rescaled inelastic ref cs

  ! Cross-section calculation
  freep = freeco_cry * ch_anuc**(one/three)

  ! compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
  cs(3) = freep*ppel
  cs(4) = freep*ppsd

  ! correct TOT-CSec for energy dependence of qel
  ! TOT CS is here without a Coulomb contribution
  cs(0) = csref_tot_rsc + freep*(pptot - pptref_cry)

  ! Also correct inel-CS
  if(csref_tot_rsc == zero) then
    cs(1) = zero
  else
    cs(1) = (csref_inel_rsc*cs(0))/csref_tot_rsc
  end if

  ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
  cs(2) = ((cs(0) - cs(1)) - cs(3)) - cs(4)
  cs(5) = ch_csref5

  ! Now add Coulomb
  cs(0) = cs(0) + cs(5)

  ! Calculate cumulative probability
  cprob(:) = zero
  cprob(5) = one
  if(cs(0) == zero) then
    do i=1,4
      cprob(i) = cprob(i-1)
    end do
  else
    do i=1,4
      cprob(i) = cprob(i-1) + cs(i)/cs(0)
    end do
  end if

  ! Multiple Coulomb Scattering
  xp = xp*c1e3
  yp = yp*c1e3

  ! Turn on/off nuclear interactions
  if(nam == 0) return

  ! Can nuclear interaction happen?
  ! Rescaled nuclear collision length
  if(avrrho == zero) then
    nuc_cl_l = c1e6
  else
    nuc_cl_l = ch_collnt/avrrho
  end if
  zlm = -nuc_cl_l*log_mb(coll_rand())

  ! write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

  if(zlm < dz) then
    ! Choose nuclear interaction
    aran = coll_rand()
    i=1
10  if(aran > cprob(i)) then
      i=i+1
      goto 10
    end if
    ichoix = i

    ! Do the interaction
    t = 0 ! default value to cover ichoix=1
    select case(ichoix)
    case(1) ! deep inelastic, impinging p disappeared
      ch_iProc = proc_ch_absorbed

    case(2) ! p-n elastic
      ch_iProc = proc_ch_pne
      bn    = (ch_bnref*cs(0))/csref_tot_rsc
      t     = -log_mb(coll_rand())/bn

    case(3) ! p-p elastic
      ch_iProc = proc_ch_ppe
      t     = -log_mb(coll_rand())/bpp

    case(4) ! Single diffractive
      ch_iProc = proc_ch_diff
      xm2   = exp_mb(coll_rand()*xln15s)
      pc    = pc*(one - xm2/ecmsq)
      if(xm2 < two) then
        bsd = two*bpp
      else if(xm2 >= two .and. xm2 <= five) then
        bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
      else if(xm2 > five) then
        bsd = (seven*bpp)/12.0_fPrec
      end if
      t = -log_mb(coll_rand())/bsd

    case(5)
      ch_iProc      = proc_ch_ruth
      length_cry = 1
      call funlux(cgen_cry(1,1),xran_cry,length_cry)
      t = xran_cry(1)

    end select

    ! Calculate the related kick -----------
    if(ichoix == 4) then
      teta = sqrt(t)/pc_in ! DIFF has changed PC!!!
    else
      teta = sqrt(t)/pc
    end if

    tx = (teta*ran_gauss(zero))*c1e3
    tz = (teta*ran_gauss(zero))*c1e3

    ! Change p angle
    xp = xp + tx
    yp = yp + tz

  end if

  xp = xp/c1e3
  yp = yp/c1e3

end subroutine cry_moveCH

! ================================================================================================ !
! Definition of rutherford scattering formula
! ================================================================================================ !
real(kind=fPrec) function cry_ruth(t_cry)

  use floatPrecision
  ! use coll_materials
  use mathlib_bouncer

  real(kind=fPrec), intent(in) :: t_cry
  real(kind=fPrec), parameter  :: cnorm  = 2.607e-4_fPrec
  real(kind=fPrec), parameter  :: cnform = 0.8561e3_fPrec

  cry_ruth = (cnorm*exp_mb(((-one*t_cry)*cnform)*emr_curr_cry**2))*(zatom_curr_cry/t_cry)**2

end function cry_ruth

end module coll_crystal
