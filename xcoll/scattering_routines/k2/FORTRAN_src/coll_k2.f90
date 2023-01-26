

! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use floatPrecision
  use numerical_constants

  implicit none

  !integer,          private, save :: mat   ! Current material
  !integer,          private, save :: mcurr ! Current material, used for Rutherford scattering integration
  real(kind=fPrec), private, save :: zatom_curr
  real(kind=fPrec), private, save :: emr_curr
  real(kind=fPrec), private, save :: zlm
  real(kind=fPrec), private, save :: zlm1
  real(kind=fPrec), private, save :: p0
  real(kind=fPrec), private, save :: length
  real(kind=fPrec), private, save :: xln15s
  real(kind=fPrec), private, save :: pptot
  real(kind=fPrec), private, save :: ppsd
  real(kind=fPrec), private, save :: ppel
  real(kind=fPrec), private, save :: ecmsq
  real(kind=fPrec), private, save :: bpp
  real(kind=fPrec), private, save :: xpsd
  real(kind=fPrec), private, save :: zpsd
  real(kind=fPrec), private, save :: psd

  ! RB DM 2014 added variables for FLUKA output
  real(kind=fPrec), private, save :: xInt,xpInt,yInt,ypInt,sInt
  real(kind=fPrec), private, save :: x,xp,z,zp,dpop

contains

!>
!! "Merlin" scattering collimation configuration
!! This routine pre-calcuates some varibles for
!! the nuclear properties
!<
!subroutine k2coll_merlinInit
!
!  use coll_materials
!
!  integer i
!
!  ! Compute the electron densnity and plasma energy for each material
!  do i=1, nmat
!    edens(i) = k2coll_calcElectronDensity(zatom(i),rho(i),anuc(i))
!    pleng(i) = k2coll_calcPlasmaEnergy(edens(i))
!  end do
!
!end subroutine k2coll_merlinInit

!>
!! K2 scattering collimation configuration
!<
subroutine k2coll_init
!nothing currently
end subroutine k2coll_init

! ================================================================================================ !
!  Collimation K2 Routine
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  G. ROBERT-DEMOLAIZE, November 1st, 2004
!  Based on routines by JBJ. Changed by RA 2001
! ================================================================================================ !
subroutine k2coll_collimate(run_exenergy, run_anuc, run_zatom, run_emr, run_rho, &
                run_hcut, run_bnref, run_csref0, run_csref1, &
                run_csref5, run_radl, run_dlri, run_dlyi, run_eum, run_ai, &
                run_collnt, is_crystal, c_length, c_rotation, c_aperture, c_offset, c_tilt,  &
                x_in, xp_in, y_in, yp_in, p_in, s_in, enom, lhit, part_abs_local, &
                impact, indiv, lint, onesided, nhit_stage, j_slices, nabs_type, linside)

  use, intrinsic :: iso_fortran_env, only : int16
  use parpro
  use crcoall
  !use coll_db
  use coll_common
  use coll_crystal, only : cry_doCrystal
  use coll_materials
  use mod_common, only : iexact, napx, unit208
  use mod_common_main, only : partID, naa
  use mathlib_bouncer
  use mod_ranlux

  !integer,          intent(in)    :: matid        ! Material ID
  real(kind=fPrec), intent(in)    :: run_exenergy
  real(kind=fPrec), intent(in)    :: run_anuc
  real(kind=fPrec), intent(in)    :: run_zatom
  real(kind=fPrec), intent(in)    :: run_emr
  real(kind=fPrec), intent(in)    :: run_rho
  real(kind=fPrec), intent(in)    :: run_hcut
  real(kind=fPrec), intent(in)    :: run_bnref
  real(kind=fPrec), intent(in)    :: run_csref0
  real(kind=fPrec), intent(in)    :: run_csref1
  real(kind=fPrec), intent(in)    :: run_csref5
  real(kind=fPrec), intent(in)    :: run_radl
  real(kind=fPrec), intent(in)    :: run_dlri
  real(kind=fPrec), intent(in)    :: run_dlyi
  real(kind=fPrec), intent(in)    :: run_eum
  real(kind=fPrec), intent(in)    :: run_ai
  real(kind=fPrec), intent(in)    :: run_collnt
  logical,          intent(in)    :: is_crystal

  real(kind=fPrec), intent(in)    :: c_length     ! Collimator length in m
  real(kind=fPrec), intent(in)    :: c_rotation   ! Collimator rotation angle vs vertical in radians
  real(kind=fPrec), intent(in)    :: c_aperture   ! Collimator aperture in m
  real(kind=fPrec), intent(in)    :: c_offset     ! Collimator offset in m
  real(kind=fPrec), intent(inout) :: c_tilt(2)    ! Collimator tilt in radians

  real(kind=fPrec), intent(inout) :: x_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: xp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: y_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: yp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: p_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: s_in(npart)  ! Particle coordinate

  real(kind=fPrec), intent(in)    :: enom         ! Reference momentum in GeV
  logical,          intent(in)    :: onesided

  integer,          intent(inout) :: lhit(npart)
  integer,          intent(inout) :: part_abs_local(npart)
  integer,          intent(inout) :: nabs_type(npart)
  integer,          intent(inout) :: nhit_stage(npart)
  real(kind=fPrec), intent(inout) :: indiv(npart)
  real(kind=fPrec), intent(inout) :: lint(npart)
  real(kind=fPrec), intent(inout) :: impact(npart)
  logical,          intent(inout) :: linside(napx)

  logical isImp
  integer j,nabs,nhit,j_slices
  real(kind=fPrec) keeps,fracab,drift_length,mirror,tiltangle
  real(kind=fPrec) x00,z00,p,sp,s,s_impact
  real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,s_flk,zpj
  real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump
  real(kind=fPrec) cRot,sRot,cRRot,sRRot
  real(kind=fPrec) xIn,xpIn,yIn,ypIn,xOut,xpOut,yOut,ypOut,sImp,sOut
  real(kind=fPrec) x_in0,xp_in0

  ! ien0,ien1: particle energy entering/leaving the collimator
  ! energy in MeV
  real(kind=fPrec)    :: ien0, ien1
  integer(kind=int16) :: nnuc0,nnuc1

  ! Initilaisation
  !mat = matid
  length = c_length
  p0     = enom

  ! Initialise scattering processes
  call k2coll_scatin(p0, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, run_csref0, run_csref1, run_csref5)

  nhit   = 0
  nabs   = 0
  fracab = zero
  mirror = one

  ! Compute rotation factors for collimator rotation
  cRot   = cos_mb(c_rotation)
  sRot   = sin_mb(c_rotation)
  cRRot  = cos_mb(-c_rotation)
  sRRot  = sin_mb(-c_rotation)

  !Set energy and nucleon change variables as with the coupling
  nnuc0 = 0
  ien0  = zero
  nnuc1 = 0
  ien1  = zero

  do j=1,napx

    if(part_abs_local(j) /= 0) then
      ! Don't do scattering process for particles already absorbed
      cycle
    end if

    impact(j) = -one
    lint(j)   = -one
    indiv(j)  = -one

    x      = x_in(j)
    xp     = xp_in(j)
    xp_in0 = xp_in(j)
    z      = y_in(j)
    zp     = yp_in(j)
    p      = p_in(j)
    sp     = zero
    dpop   = (p - p0)/p0
    x_flk  = zero
    y_flk  = zero
    xp_flk = zero
    yp_flk = zero

    ! Transform particle coordinates to get into collimator coordinate  system
    ! First do rotation into collimator frame
    x  =  x_in(j)*cRot + sRot*y_in(j)
    z  =  y_in(j)*cRot - sRot*x_in(j)
    xp = xp_in(j)*cRot + sRot*yp_in(j)
    zp = yp_in(j)*cRot - sRot*xp_in(j)

    ! For one-sided collimators consider only positive X. For negative X jump to the next particle
    if(onesided .and. x < zero) then
      cycle
    end if

! Log input energy + nucleons as per the FLUKA coupling
    nnuc0   = nnuc0 + naa(j)
    ien0    = ien0 + rcp(j) * c1e3


    ! Now mirror at the horizontal axis for negative X offset
    if(x < zero) then
      mirror    = -one
      tiltangle = -one*c_tilt(2)
    else
      mirror    = one
      tiltangle = c_tilt(1)
    end if
    x  = mirror*x
    xp = mirror*xp

    ! Shift with opening and offset
    x = (x - c_aperture/two) - mirror*c_offset

    ! Include collimator tilt
    if(tiltangle > zero) then
      xp = xp - tiltangle
    end if
    if(tiltangle < zero) then
      x  = x + sin_mb(tiltangle) * c_length
      xp = xp - tiltangle
    end if

    ! CRY Only: x_in0 has to be assigned after the change of reference frame
    x_in0 = x

    ! After finishing the coordinate transformation, or the coordinate manipulations in case of pencil beams,
    ! save the initial coordinates of the impacting particles
    xIn  = x
    xpIn = xp
    yIn  = z
    ypIn = zp

    ! particle passing above the jaw are discarded => take new event
    ! entering by the face, shorten the length (zlm) and keep track of
    ! entrance longitudinal coordinate (keeps) for histograms

    ! The definition is that the collimator jaw is at x>=0.

    ! 1) Check whether particle hits the collimator
    isImp = .false.
    s     = zero
    keeps = zero
    zlm   = -one*length

    if(is_crystal) then ! This is a crystal collimator

      call cry_doCrystal(j,x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit,&
                         part_abs_local,impact,indiv,c_length, run_exenergy, run_anuc, &
                         run_zatom, run_emr, run_rho, run_hcut, run_bnref, run_csref0, &
                         run_csref1, run_csref5, run_dlri, run_dlyi, run_eum, run_ai, run_collnt)

      if(nabs /= 0) then
        part_abs_local(j)  = 1
        lint(j)                = zlm
      end if

      sImp  = (s - c_length) + sImp
      sOut  = s
      xOut  = x
      xpOut = xp
      yOut  = z
      ypOut = zp

    else ! "Normal" collimator

      if(x >= zero) then
        ! Particle hits collimator and we assume interaction length ZLM equal
        ! to collimator length (what if it would leave collimator after
        ! small length due to angle???)
        zlm       = length
        impact(j) = x
        indiv(j)  = xp
      else if(xp <= zero) then
        ! Particle does not hit collimator. Interaction length ZLM is zero.
        zlm = zero
      else
        ! Calculate s-coordinate of interaction point
        s = (-one*x)/xp
        if(s <= zero) then
          write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
          !call prror
        end if
        if(s < length) then
          zlm       = length - s
          impact(j) = zero
          indiv(j)  = xp
        else
          zlm = zero
        end if
      end if

      ! First do the drift part
      ! DRIFT PART
      drift_length = length - zlm
      if(drift_length > zero) then
        if(iexact) then
          zpj = sqrt(one-xp**2-zp**2)
          x   = x  + drift_length*(xp/zpj)
          z   = z  + drift_length*(zp/zpj)
          sp  = sp + drift_length
        else
          x  = x  + xp* drift_length
          z  = z  + zp * drift_length
          sp = sp + drift_length
        end if
      end if

      ! Now do the scattering part
      if(zlm > zero) then
        if(.not.linside(j)) then
          ! first time particle hits collimator: entering jaw
          linside(j) = .true.
          if(dowrite_impact) then
            if(tiltangle > zero) then
              x_Dump = (x + c_aperture/two + tiltangle*sp)*mirror + c_offset
            else
              x_Dump = (x + c_aperture/two + tiltangle*(sp - c_length))*mirror + c_offset
            end if
            xpDump = (xp + tiltangle)*mirror
            y_Dump = z
            ypDump = zp
            s_Dump = sp+real(j_slices-1,fPrec)*c_length
          end if
        end if

        s_impact = sp
        nhit = nhit + 1
        call k2coll_jaw(s,nabs,partID(j), run_exenergy, run_anuc, run_zatom, run_rho, run_radl)

        nabs_type(j) = nabs
        lhit(j)  = 1

        isImp = .true.
        sImp  = s_impact+(real(j_slices,fPrec)-one)*c_length
        sOut  = (s+sp)+(real(j_slices,fPrec)-one)*c_length
        xOut  = x
        xpOut = xp
        yOut  = z
        ypOut = zp

        ! Writeout should be done for both inelastic and single diffractive. doing all transformations
        ! in x_flk and making the set to 99.99 mm conditional for nabs=1
        if(dowrite_impact .or. nabs == 1 .or. nabs == 4) then
          ! Transform back to lab system for writeout.
          ! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

          x_flk  = xInt
          xp_flk = xpInt

          if(tiltangle > zero) then
            x_flk  = x_flk  + tiltangle*(sInt+sp)
            xp_flk = xp_flk + tiltangle
          else if(tiltangle < zero) then
            xp_flk = xp_flk + tiltangle
            x_flk  = x_flk  - sin_mb(tiltangle) * (length-(sInt+sp))
          end if

          x_flk  = (x_flk + c_aperture/two) + mirror*c_offset
          x_flk  = mirror*x_flk
          xp_flk = mirror*xp_flk
          y_flk  = (  yInt*cRRot -  x_flk*sRRot)*c1e3
          yp_flk = ( ypInt*cRRot - xp_flk*sRRot)*c1e3
          x_flk  = ( x_flk*cRRot +   yInt*sRRot)*c1e3
          xp_flk = (xp_flk*cRRot +  ypInt*sRRot)*c1e3
          s_flk  = (sInt+sp)+(real(j_slices,fPrec)-one)*c_length

          ! Finally, the actual coordinate change to 99 mm
          if(nabs == 1) then
            fracab  = fracab + 1
            x       = 99.99e-3_fPrec
            z       = 99.99e-3_fPrec
            lint(j) = zlm
            part_abs_local(j)  = 1
          end if
        end if
      end if ! Collimator jaw interaction

      if(nabs /= 1 .and. zlm > zero) then
        ! Do the rest drift, if particle left collimator early
        drift_length = (length-(s+sp))
        if(drift_length > c1m15) then
          linside(j) = .false.
          if(dowrite_impact) then
            if(tiltangle > zero) then
              x_Dump = (x + c_aperture/two + tiltangle*(s+sp))*mirror + c_offset
            else
              x_Dump = (x + c_aperture/two + tiltangle*(s+sp-c_length))*mirror + c_offset
            end if
            xpDump = (xp+tiltangle)*mirror
            y_Dump = z
            ypDump = zp
            s_Dump = s+sp+real(j_slices-1,fPrec)*c_length

          end if
          if(iexact) then
            zpj = sqrt(one-xp**2-zp**2)
            x   = x  + drift_length*(xp/zpj)
            z   = z  + drift_length*(zp/zpj)
            sp  = sp + drift_length
          else
            x  = x  + xp * drift_length
            z  = z  + zp * drift_length
            sp = sp + drift_length
          end if
        end if
        lint(j) = zlm - drift_length
      end if

    end if ! Collimator isCrystal

    ! Transform back to particle coordinates with opening and offset
    if(x < 99.0e-3_fPrec) then
      ! Include collimator tilt
      if(tiltangle > zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle < zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
        x  = x  - sin_mb(tiltangle) * c_length
      end if

      ! Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x   = (x + c_aperture/two) + mirror*c_offset

      ! Now mirror at the horizontal axis for negative X offset
      x  = mirror * x
      xp = mirror * xp

      ! Last do rotation into collimator frame
      x_in(j)  =  x*cRRot +  z*sRRot
      y_in(j)  =  z*cRRot -  x*sRRot
      xp_in(j) = xp*cRRot + zp*sRRot
      yp_in(j) = zp*cRRot - xp*sRRot

! Log output energy + nucleons as per the FLUKA coupling
! Do not log dead particles
      nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
      ien1        = ien1  + rcp(j) * c1e3                   ! outcoming energy

      if(is_crystal) then
        p_in(j) = p
        s_in(j) = s_in(j) + s
      else
        p_in(j) = (one + dpop) * p0
        s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length
      end if
    else
      x_in(j) = x
      y_in(j) = z
    end if

  end do ! End of loop over all particles

end subroutine k2coll_collimate

!>
!! k2coll_scatin(plab)
!! Configure the K2 scattering routine cross sections
!<
subroutine k2coll_scatin(plab, sc_anuc, sc_zatom, sc_emr, sc_rho, sc_hcut, sc_bnref, sc_csref0, sc_csref1, sc_csref5)

  use mod_funlux
  use coll_common
  use coll_materials
  use mathlib_bouncer
  use physical_constants
  use mod_units
  use crcoall

  real(kind=fPrec), intent(in) :: plab
  real(kind=fPrec), intent(in) :: sc_anuc
  real(kind=fPrec), intent(in) :: sc_zatom
  real(kind=fPrec), intent(in) :: sc_emr
  real(kind=fPrec), intent(in) :: sc_rho
  real(kind=fPrec), intent(in) :: sc_hcut
  real(kind=fPrec), intent(in) :: sc_bnref
  real(kind=fPrec), intent(in) :: sc_csref0
  real(kind=fPrec), intent(in) :: sc_csref1
  real(kind=fPrec), intent(in) :: sc_csref5

  real(kind=fPrec), parameter :: tlcut = 0.0009982_fPrec
  integer i

  integer csUnit
  character(len=23), parameter :: cs_fileName = "MaterialInformation.txt"
  logical csErr

  ecmsq = (two*(pmap*c1m3)) * plab
  xln15s = log_mb(0.15_fPrec*ecmsq)
  ! Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
  pptot = 0.041084_fPrec-0.0023302_fPrec*log_mb(ecmsq)+0.00031514_fPrec*log_mb(ecmsq)**2
  ! Claudia used the fit from TOTEM for ppel (in barn)
  ppel = (11.7_fPrec-1.59_fPrec*log_mb(ecmsq)+0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ! Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
  ppsd = (4.3_fPrec+0.3_fPrec*log_mb(ecmsq))/c1e3


  ! Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
  bpp = 7.156_fPrec + 1.439_fPrec*log_mb(sqrt(ecmsq))

  ! Unmeasured tungsten data, computed with lead data and power laws
  !bnref(4) = bnref(5)*(anuc(4)/anuc(5))**(two/three)
  !emr(4)   = emr(5)  *(anuc(4)/anuc(5))**(one/three)

  ! Compute cross-sections (CS) and probabilities + Interaction length
  ! Last two material treated below statement number 100
    zatom_curr = sc_zatom
    emr_curr = sc_emr
    ! Prepare for Rutherford differential distribution
    call funlxp(k2coll_ruth, cgen(1), tlcut, sc_hcut)

    ! freep: number of nucleons involved in single scattering
    freep = freeco * sc_anuc**(one/three)

    ! compute pp and pn el+single diff contributions to cross-section
    ! (both added : quasi-elastic or qel later)
    csect(3) = freep * ppel
    csect(4) = freep * ppsd

    ! correct TOT-CSec for energy dependence of qel
    ! TOT CS is here without a Coulomb contribution
    csect(0) = sc_csref0 + freep * (pptot - pptref)
    bn       = (sc_bnref * csect(0)) / sc_csref0

    ! also correct inel-CS
    csect(1) = (sc_csref1 * csect(0)) / sc_csref0

    ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect(2) = ((csect(0) - csect(1)) - csect(3)) - csect(4)
    csect(5) = sc_csref5

    ! Now add Coulomb
    csect(0) = csect(0) + csect(5)

    ! Interaction length in meter
    xintl = (c1m2*sc_anuc)/(((fnavo * sc_rho)*csect(0))*1e-24_fPrec)

    ! Filling CProb with cumulated normalised Cross-sections
    do i=1,4
      cprob(i) = cprob(i-1) + csect(i)/csect(0)
    end do

  ! Last two materials for 'vaccum' (nmat-1) and 'full black' (nmat)
  !cprob(1,nmat-1) = one
  !cprob(1,nmat)   = one
  !xintl(nmat-1)   = c1e12
  !xintl(nmat)     = zero

end subroutine k2coll_scatin

!>
!! jaw(s,nabs,icoll,iturn,ipart)
!!     RB: adding as input arguments to jaw variables icoll,iturn,ipart
!!         these are only used for the writeout of particle histories
!!
!! Input:   ZLM is interaction length
!!          MAT is choice of material
!!
!! Output:  nabs = 1   Particle is absorped
!!          nabs = 4   Single-diffractive scattering
!!          dpop       Adjusted for momentum loss (dE/dx)
!!          s          Exit longitudinal position
!!
!! Physics:  If monte carlo interaction length greater than input
!!           interaction length, then use input interaction length
!!           Is that justified???
!<
subroutine k2coll_jaw(s, nabs, ipart, j_exenergy, j_anuc, j_zatom, j_rho, j_radl)

  use mod_ranlux
  use coll_common
  use coll_materials
  use mathlib_bouncer

  real(kind=fPrec), intent(inout) :: s
  integer,          intent(inout) :: nabs
  integer,          intent(in)    :: ipart
  real(kind=fPrec), intent(in)    :: j_exenergy
  real(kind=fPrec), intent(in)    :: j_anuc
  real(kind=fPrec), intent(in)    :: j_zatom
  real(kind=fPrec), intent(in)    :: j_rho
  real(kind=fPrec), intent(in)    :: j_radl

  real(kind=fPrec) m_dpodx,p,rlen,t,dxp,dzp,p1,zpBef,xpBef,pBef
  integer inter,nabs_tmp

  ! Note that the input parameter is dpop. Here the momentum p is constructed out of this input.
  p    = p0*(one+dpop)
  nabs = 0
  nabs_tmp = nabs


  ! Initialize the interaction length to input interaction length
  rlen = zlm

  ! Do a step for a point-like interaction.
  ! Get monte-carlo interaction length.
10 continue
  zlm1     = (-one*xintl)*log_mb(coll_rand())
  nabs_tmp = 0  ! type of interaction reset before following scattering process
  xpBef    = xp ! save angles and momentum before scattering
  zpBef    = zp
  pBef     = p

  ! If the monte-carlo interaction length is longer than the
  ! remaining collimator length, then put it to the remaining
  ! length, do multiple coulomb scattering and return.
  ! LAST STEP IN ITERATION LOOP
  if(zlm1 > rlen) then
    zlm1 = rlen
    call k2coll_mcs(s, j_radl)
    s = (zlm-rlen)+s
    call k2coll_calcIonLoss(p,rlen,m_dpodx, j_exenergy, j_anuc, j_zatom, j_rho)  ! DM routine to include tail
    p = p-m_dpodx*s

    dpop = (p-p0)/p0
    return
  end if

  ! Otherwise do multi-coulomb scattering.
  ! REGULAR STEP IN ITERATION LOOP
  call k2coll_mcs(s, j_radl)

  ! Check if particle is outside of collimator (X.LT.0) after
  ! MCS. If yes, calculate output longitudinal position (s),
  ! reduce momentum (output as dpop) and return.
  ! PARTICLE LEFT COLLIMATOR BEFORE ITS END.
  if(x <= zero) then
    s = (zlm-rlen)+s

    call k2coll_calcIonLoss(p,rlen,m_dpodx, j_exenergy, j_anuc, j_zatom, j_rho)
    p = p-m_dpodx*s
    dpop = (p-p0)/p0

    return
  end if

  ! Check whether particle is absorbed. If yes, calculate output
  ! longitudinal position (s), reduce momentum (output as dpop)
  ! and return.
  ! PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

  inter    = k2coll_ichoix()
  nabs     = inter
  nabs_tmp = nabs

  ! RB, DM: save coordinates before interaction for writeout to FLUKA_impacts.dat
  xInt  = x
  xpInt = xp
  yInt  = z
  ypInt = zp
  sInt  = (zlm-rlen)+zlm1

  if(inter == 1) then
    s = (zlm-rlen)+zlm1

    call k2coll_calcIonLoss(p,rlen,m_dpodx, j_exenergy, j_anuc, j_zatom, j_rho)
    p = p-m_dpodx*s

    dpop = (p-p0)/p0

    return
  end if

  ! Now treat the other types of interaction, as determined by ICHOIX:

  ! Nuclear-Elastic:          inter = 2
  ! pp Elastic:               inter = 3
  ! Single-Diffractive:       inter = 4    (changes momentum p)
  ! Coulomb:                  inter = 5

  ! As the single-diffractive interaction changes the momentum, save input momentum in p1.
  p1 = p

  ! Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.
  t = k2coll_gettran(inter,p)

  ! Tetat calculates from the rms transverse momentum transfer in
  ! monte-carlo fashion the angle changes for x and z planes. The
  ! angle change is proportional to SQRT(t) and 1/p, as expected.
  call k2coll_tetat(t,p,dxp,dzp)

  ! Apply angle changes
  xp = xp+dxp
  zp = zp+dzp

  ! Treat single-diffractive scattering.
  if(inter == 4) then

    ! added update for s
    s    = (zlm-rlen)+zlm1
    xpsd = dxp
    zpsd = dzp
    psd  = p1

    ! Add this code to get the momentum transfer also in the calling routine
    dpop = (p-p0)/p0
  end if

  ! Calculate the remaining interaction length and close the iteration loop.
  rlen = rlen-zlm1
  goto 10

end subroutine k2coll_jaw

!>
!! mcs(s)
!!++  Input:   zlm1   Monte-carlo interaction length
!!
!!++  Output:  s      Longitudinal position
!!++           p0     Reference momentum
!!++           dpop   Relative momentum offset
!!
!!     collimator: x>0 and y<zlm1
!<
subroutine k2coll_mcs(s, mc_radl)

  use coll_materials

  real(kind=fPrec), intent(inout) :: s
  real(kind=fPrec), intent(in)    :: mc_radl

  real(kind=fPrec) theta,rlen0,rlen,ae,be!,radl_mat,rad_len

  real(kind=fPrec), parameter :: h   = 0.001_fPrec
  real(kind=fPrec), parameter :: dh  = 0.0001_fPrec
  real(kind=fPrec), parameter :: bn0 = 0.4330127019_fPrec

  !radl_mat = radl(mat)
  theta    = 13.6e-3_fPrec/(p0*(one+dpop))
  !rad_len  = radl(mat)

  x     = (x/theta)/mc_radl
  xp    = xp/theta
  z     = (z/theta)/mc_radl
  zp    = zp/theta
  rlen0 = zlm1/mc_radl
  rlen  = rlen0

10 continue
  ae = bn0*x
  be = bn0*xp

  call k2coll_soln3(ae,be,dh,rlen,s)
  if(s < h) s = h

  call k2coll_scamcs(x,xp,s)
  if(x <= zero) then
    s = (rlen0-rlen)+s
    goto 20
  end if
  if(s+dh >= rlen) then
    s = rlen0
    goto 20
  end if
  rlen = rlen-s
  goto 10

20 continue
  call k2coll_scamcs(z,zp,s)
  s  = s*mc_radl
  x  = (x*theta)*mc_radl
  xp = xp*theta
  z  = (z*theta)*mc_radl
  zp = zp*theta

end subroutine k2coll_mcs

!>
!! k2coll_calcIonLoss(IS,PC,DZ,EnLo)
!! subroutine for the calculazion of the energy loss by ionization
!! Either mean energy loss from Bethe-Bloch, or higher energy loss, according to finite probability from cross section
!! Written by DM for crystals, introduced in main code by RB
!! Updated and improved for numerical stability by VKBO
!<
subroutine k2coll_calcIonLoss(PC, DZ, EnLo, ci_exenergy, ci_anuc, ci_zatom, ci_rho)

  use mod_ranlux
  use coll_materials
  use mathlib_bouncer
  use physical_constants

  !integer,          intent(in)  :: IS   ! IS material ID
  real(kind=fPrec), intent(in)  :: PC   ! PC momentum in GeV
  real(kind=fPrec), intent(in)  :: DZ   ! DZ length traversed in material (meters)
  real(kind=fPrec), intent(out) :: EnLo ! EnLo energy loss in GeV/meter
  real(kind=fPrec),    intent(in)    :: ci_exenergy
  real(kind=fPrec),    intent(in)    :: ci_anuc
  real(kind=fPrec),    intent(in)    :: ci_zatom
  real(kind=fPrec),    intent(in)    :: ci_rho

  real(kind=fPrec) exEn,thl,Tt,cs_tail,prob_tail,enr,mom,betar,gammar,bgr,kine,Tmax,plen
  real(kind=fPrec), parameter :: k = 0.307075_fPrec ! Constant in front of Bethe-Bloch [MeV g^-1 cm^2]

  mom    = PC*c1e3                     ! [GeV/c] -> [MeV/c]
  enr    = (mom*mom + pmap*pmap)**half ! [MeV]
  gammar = enr/pmap
  betar  = mom/enr
  bgr    = betar*gammar
  kine   = ((two*pmae)*bgr)*bgr

  ! Mean excitation energy
  exEn = ci_exenergy*c1e3 ! [MeV]

  ! Tmax is max energy loss from kinematics
  Tmax = kine/(one + (two*gammar)*(pmae/pmap) + (pmae/pmap)**2) ! [MeV]

  ! Plasma energy - see PDG 2010 table 27.1
  plen = (((ci_rho*ci_zatom)/ci_anuc)**half)*28.816e-6_fPrec ! [MeV]

  ! Calculate threshold energy
  ! Above this threshold, the cross section for high energy loss is calculated and then
  ! a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
  ! below threshold, only the standard Bethe-Bloch is used (all particles get average energy loss)

  ! thl is 2*width of Landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
  thl = ((((four*(k*ci_zatom))*DZ)*c1e2)*ci_rho)/(ci_anuc*betar**2) ! [MeV]

  ! Bethe-Bloch mean energy loss
  EnLo = ((k*ci_zatom)/(ci_anuc*betar**2)) * ( &
    half*log_mb((kine*Tmax)/(exEn*exEn)) - betar**2 - log_mb(plen/exEn) - log_mb(bgr) + half &
  )
  EnLo = ((EnLo*ci_rho)*c1m1)*DZ ! [GeV]

  ! Threshold Tt is Bethe-Bloch + 2*width of Landau distribution
  Tt = EnLo*c1e3 + thl ! [MeV]

  ! Cross section - see Alfredo's presentation for derivation
  cs_tail = ((k*ci_zatom)/(ci_anuc*betar**2)) * ( &
    half*((one/Tt)-(one/Tmax)) - (log_mb(Tmax/Tt)*betar**2)/(two*Tmax) + (Tmax-Tt)/((four*gammar**2)*pmap**2) &
  )

  ! Probability of being in tail: cross section * density * path length
  prob_tail = ((cs_tail*ci_rho)*DZ)*c1e2

  ! Determine based on random number if tail energy loss occurs.
  if(coll_rand() < prob_tail) then
    EnLo = ((k*ci_zatom)/(ci_anuc*betar**2)) * ( &
      half*log_mb((kine*Tmax)/(exEn*exEn)) - betar**2 - log_mb(plen/exEn) - log_mb(bgr) + &
      half + TMax**2/((eight*gammar**2)*pmap**2) &
    )
    EnLo = (EnLo*ci_rho)*c1m1 ! [GeV/m]
  else
    ! If tail energy loss does not occur, just use the standard Bethe-Bloch
    EnLo = EnLo/DZ  ! [GeV/m]
  end if

end subroutine k2coll_calcIonLoss

!>
!! k2coll_tetat(t,p,tx,tz)
!! Generate sine and cosine of an angle uniform in [0,2pi](see RPP)
!<
subroutine k2coll_tetat(t, p, tx, tz)

  use mod_ranlux

  real(kind=fPrec), intent(in)  :: t
  real(kind=fPrec), intent(in)  :: p
  real(kind=fPrec), intent(out) :: tx
  real(kind=fPrec), intent(out) :: tz

  real(kind=fPrec) va,vb,va2,vb2,r2,teta

  teta = sqrt(t)/p

10 continue
  va  = two*coll_rand() - one
  vb  = coll_rand()
  va2 = va**2
  vb2 = vb**2
  r2  = va2 + vb2
  if(r2 > one) goto 10
  tx  = (teta*((two*va)*vb))/r2
  tz  = (teta*(va2 - vb2))/r2

end subroutine k2coll_tetat

!>
!! k2coll_soln3(a,b,dh,smax,s)
!<
subroutine k2coll_soln3(a, b, dh, smax, s)

  real(kind=fPrec), intent(in)    :: a
  real(kind=fPrec), intent(in)    :: b
  real(kind=fPrec), intent(in)    :: dh
  real(kind=fPrec), intent(in)    :: smax
  real(kind=fPrec), intent(inout) :: s

  real(kind=fPrec) c

  if(b == zero) then
    s = a**0.6666666666666667_fPrec
  ! s = a**(two/three)
    if(s > smax) s = smax
    return
  end if

  if(a == zero) then
    if(b > zero) then
      s = b**2
    else
      s = zero
    end if
    if(s > smax) s=smax
    return
  end if

  if(b > zero) then
    if(smax**3 <= (a + b*smax)**2) then
      s = smax
      return
    else
      s = smax*half
      call k2coll_iterat(a,b,dh,s)
    end if
  else
    c = (-one*a)/b
    if(smax < c) then
      if(smax**3 <= (a + b*smax)**2) then
        s = smax
        return
      else
        s = smax*half
        call k2coll_iterat(a,b,dh,s)
      end if
    else
      s = c*half
      call k2coll_iterat(a,b,dh,s)
    end if
  end if

end subroutine k2coll_soln3

!>
!! k2coll_scamcs(xx,xxp,s)
!<
subroutine k2coll_scamcs(xx, xxp, s)

  use mathlib_bouncer
  use mod_ranlux

  real(kind=fPrec), intent(inout) :: xx
  real(kind=fPrec), intent(inout) :: xxp
  real(kind=fPrec), intent(in)    :: s

  real(kind=fPrec) v1,v2,r2,a,z1,z2,ss,x0,xp0,sss

  x0  = xx
  xp0 = xxp

10 continue
  v1 = two*coll_rand() - one
  v2 = two*coll_rand() - one
  r2 = v1**2 + v2**2
  if(r2 >= one) goto 10

  a   = sqrt((-two*log_mb(r2))/r2)
  z1  = v1*a
  z2  = v2*a
  ss  = sqrt(s)
  sss = one + 0.038_fPrec*log_mb(s)
  xx  = x0  + s*(xp0 + ((half*ss)*sss)*(z2 + z1*0.577350269_fPrec))
  xxp = xp0 + (ss*z2)*sss

end subroutine k2coll_scamcs

subroutine k2coll_iterat(a, b, dh, s)

  real(kind=fPrec), intent(in)    :: a
  real(kind=fPrec), intent(in)    :: b
  real(kind=fPrec), intent(in)    :: dh
  real(kind=fPrec), intent(inout) :: s

  real(kind=fPrec) ds

  ds = s

10 continue
  ds = ds*half

  if(s**3 < (a+b*s)**2) then
    s = s+ds
  else
    s = s-ds
  end if

  if(ds < dh) then
    return
  else
    goto 10
  end if

end subroutine k2coll_iterat

!>
!! k2coll_ruth(t)
!! Calculate the rutherford scattering cross section
!<
real(kind=fPrec) function k2coll_ruth(t)

  use mathlib_bouncer
  use coll_materials

  real(kind=fPrec), intent(in) :: t

  ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug
  real(kind=fPrec), parameter :: cnorm  = 2.607e-5_fPrec
  real(kind=fPrec), parameter :: cnform = 0.8561e3_fPrec

  k2coll_ruth = (cnorm*exp_mb(((-one*t)*cnform)*emr_curr**2)) * (zatom_curr/t)**2

end function k2coll_ruth

!>
!! k2coll_gettran(inter,xmat,p)
!! This function determines: GETTRAN - rms transverse momentum transfer
!! Note: For single-diffractive scattering the vector p of momentum
!! is modified (energy loss is applied)
!<
real(kind=fPrec) function k2coll_gettran(inter, p)

  use mathlib_bouncer
  use mod_ranlux
  use mod_funlux
  use coll_materials

  integer,          intent(in)    :: inter
  !integer,          intent(in)    :: xmat
  real(kind=fPrec), intent(inout) :: p

  real(kind=fPrec) xm2,bsd,xran(1)

  ! Neither if-statements below have an else, so defaulting function return to zero.
  k2coll_gettran = zero

  select case(inter)
  case(2) ! Nuclear Elastic
    k2coll_gettran = (-one*log_mb(coll_rand()))/bn
  case(3) ! pp Elastic
    k2coll_gettran = (-one*log_mb(coll_rand()))/bpp
  case(4) ! Single Diffractive
    xm2 = exp_mb(coll_rand() * xln15s)
    p   = p * (one - xm2/ecmsq)
    if(xm2 < two) then
      bsd = two * bpp
    else if(xm2 >= two .and. xm2 <= five) then
      bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*bpp)/36.0_fPrec
    else
      bsd = (seven*bpp)/12.0_fPrec
    end if
    k2coll_gettran = (-one*log_mb(coll_rand()))/bsd
  case(5) ! Coulomb
    call funlux(cgen(1), xran, 1)
    k2coll_gettran = xran(1)
  end select

end function k2coll_gettran

!>
!! k2coll_ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
integer function k2coll_ichoix()

  use mod_ranlux
  use coll_materials

  !integer, intent(in) :: ma
  integer i
  real(kind=fPrec) aran

  aran = coll_rand()
  i    = 1
10 continue
  if(aran > cprob(i)) then
    i = i+1
    goto 10
  end if

  k2coll_ichoix = i

end function k2coll_ichoix

!>
!! k2coll_calcElectronDensity(AtomicNumber, Density, AtomicMass)
!! Function to calculate the electron density in a material
!! Should give the number per cubic meter
!<
real(kind=fPrec) function k2coll_calcElectronDensity(AtomicNumber, Density, AtomicMass)

  real(kind=fPrec) , intent(in) :: AtomicNumber
  real(kind=fPrec) , intent(in) :: Density
  real(kind=fPrec) , intent(in) :: AtomicMass

  real(kind=fPrec), parameter :: Avogadro = 6.022140857e23_fPrec
  real(kind=fPrec) PartA, PartB

  PartA = (AtomicNumber*Avogadro) * Density
  PartB = AtomicMass * c1m6 ! 1e-6 factor converts to n/m^-3
  k2coll_calcElectronDensity = PartA/PartB

end function k2coll_calcElectronDensity

!>
!! k2coll_calcPlasmaEnergy(ElectronDensity)
!! Function to calculate the plasma energy in a material
!! CalculatePlasmaEnergy = (PlanckConstantBar * sqrt((ElectronDensity *(ElectronCharge**2)) / &
!!& (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
!<
real(kind=fPrec) function k2coll_calcPlasmaEnergy(ElectronDensity)

  use physical_constants

  real(kind=fPrec), intent(in) :: ElectronDensity

  real(kind=fPrec) sqrtAB,PartA

  ! Values from the 2016 PDG
  real(kind=fPrec), parameter :: PlanckConstantBar = 1.054571800e-34_fPrec
  real(kind=fPrec), parameter :: ElectronCharge = echarge 
  real(kind=fPrec), parameter :: ElectronCharge2 = ElectronCharge*ElectronCharge
  real(kind=fPrec), parameter :: ElectronMass = 9.10938356e-31_fPrec
  real(kind=fPrec), parameter :: SpeedOfLight2 = clight*clight
  real(kind=fPrec), parameter :: FreeSpacePermeability = 16.0e-7_fPrec*atan(one)
  real(kind=fPrec), parameter :: FSPC2 = FreeSpacePermeability*SpeedOfLight2
  real(kind=fPrec), parameter :: FreeSpacePermittivity = one/FSPC2
  real(kind=fPrec), parameter :: PartB = ElectronMass * FreeSpacePermittivity

  PartA  = ElectronDensity * ElectronCharge2
  sqrtAB = sqrt(PartA/PartB)

  k2coll_calcPlasmaEnergy = ((PlanckConstantBar*sqrtAB)/ElectronCharge)*c1m9

end function k2coll_calcPlasmaEnergy


end module coll_k2
