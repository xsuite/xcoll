

! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use floatPrecision
  use numerical_constants

  implicit none

  real(kind=fPrec), save :: zatom_curr ! Current zatom, used for Rutherford scattering integration
  real(kind=fPrec), save :: emr_curr ! Current emr, used for Rutherford scattering integration

contains

subroutine k2coll_init
!nothing currently
end subroutine k2coll_init

! ================================================================================================ !
!  Collimation K2 Routine
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  G. ROBERT-DEMOLAIZE, November 1st, 2004
!  Based on routines by JBJ. Changed by RA 2001
! ================================================================================================ !

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
! subroutine k2coll_jaw(s, nabs, j_exenergy, j_anuc, j_zatom, j_rho, j_radl, &
!                     j_cprob, j_xintl, j_bn, j_cgen, j_ecmsq, j_xln15s, j_bpp, j_zlm, &
!                     j_p0, j_x, j_xp, j_z, j_zp, j_dpop, j_xInt, j_xpInt, j_yInt, j_ypInt, j_sInt)

!   use mod_ranlux, only : coll_rand
!   use coll_common
!   !
!   use mathlib_bouncer

!   real(kind=fPrec), intent(inout) :: s
!   integer,          intent(inout) :: nabs
!   real(kind=fPrec), intent(inout) :: j_exenergy
!   real(kind=fPrec), intent(in)    :: j_anuc
!   real(kind=fPrec), intent(in)    :: j_zatom
!   real(kind=fPrec), intent(in)    :: j_rho
!   real(kind=fPrec), intent(in)    :: j_radl
!   real(kind=fPrec), intent(in)    :: j_cprob(0:5)
!   real(kind=fPrec), intent(in)    :: j_xintl
!   real(kind=fPrec), intent(in)    :: j_bn
!   real(kind=fPrec), intent(in)    :: j_cgen(200)
!   real(kind=fPrec), intent(in)    :: j_ecmsq
!   real(kind=fPrec), intent(in)    :: j_xln15s
!   real(kind=fPrec), intent(in)    :: j_bpp
!   real(kind=fPrec), intent(in)    :: j_zlm
!   real(kind=fPrec), intent(in)    :: j_p0

!   real(kind=fPrec), intent(inout) :: j_xInt
!   real(kind=fPrec), intent(inout) :: j_xpInt
!   real(kind=fPrec), intent(inout) :: j_yInt
!   real(kind=fPrec), intent(inout) :: j_ypInt
!   real(kind=fPrec), intent(inout) :: j_sInt

!   real(kind=fPrec), intent(inout) :: j_x
!   real(kind=fPrec), intent(inout) :: j_xp
!   real(kind=fPrec), intent(inout) :: j_z
!   real(kind=fPrec), intent(inout) :: j_zp
!   real(kind=fPrec), intent(inout) :: j_dpop

!   real(kind=fPrec) m_dpodx,p,rlen,t,dxp,dzp,p1,zpBef,xpBef,pBef,j_zlm1,xpsd,zpsd,psd
!   integer inter,nabs_tmp

!   ! Note that the input parameter is dpop. Here the momentum p is constructed out of this input.
!   p    = j_p0*(one+j_dpop)
!   nabs = 0
!   nabs_tmp = nabs
    

!   ! Initialize the interaction length to input interaction length
!   rlen = j_zlm

!   ! Do a step for a point-like interaction.
!   ! Get monte-carlo interaction length.
! 10 continue
!   j_zlm1     = (-one*j_xintl)*log_mb(coll_rand())
!   nabs_tmp = 0  ! type of interaction reset before following scattering process
!   xpBef    = j_xp ! save angles and momentum before scattering
!   zpBef    = j_zp
!   pBef     = p

!   ! If the monte-carlo interaction length is longer than the
!   ! remaining collimator length, then put it to the remaining
!   ! length, do multiple coulomb scattering and return.
!   ! LAST STEP IN ITERATION LOOP
!   if(j_zlm1 > rlen) then
!     j_zlm1 = rlen
!     call k2coll_mcs(s,j_radl,j_zlm1,j_p0,j_x,j_xp,j_z,j_zp,j_dpop)
!     s = (j_zlm-rlen)+s
!     call k2coll_calcIonLoss(p,rlen,j_exenergy,j_anuc,j_zatom,j_rho,m_dpodx)  ! DM routine to include tail
!     p = p-m_dpodx*s

!     j_dpop = (p-j_p0)/j_p0
!     return
!   end if
!     ! Otherwise do multi-coulomb scattering.
!   ! REGULAR STEP IN ITERATION LOOP
  ! call k2coll_mcs(s,j_radl,j_zlm1,j_p0,j_x,j_xp,j_z,j_zp,j_dpop)
  !   ! Check if particle is outside of collimator (X.LT.0) after
  ! ! MCS. If yes, calculate output longitudinal position (s),
  ! ! reduce momentum (output as dpop) and return.
  ! ! PARTICLE LEFT COLLIMATOR BEFORE ITS END.
  ! if(j_x <= zero) then
  !   s = (j_zlm-rlen)+s
  !       call k2coll_calcIonLoss(p,rlen,j_exenergy,j_anuc,j_zatom,j_rho,m_dpodx)
  !   p = p-m_dpodx*s
  !   j_dpop = (p-j_p0)/j_p0
  !       return
  ! end if

  ! ! Check whether particle is absorbed. If yes, calculate output
  ! ! longitudinal position (s), reduce momentum (output as dpop)
  ! ! and return.
  ! ! PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.
  !   inter    = k2coll_ichoix(j_cprob)
  !   nabs     = inter
  ! nabs_tmp = nabs

  ! ! RB, DM: save coordinates before interaction for writeout to FLUKA_impacts.dat
  ! j_xInt  = j_x
  ! j_xpInt = j_xp
  ! j_yInt  = j_z
  ! j_ypInt = j_zp
  ! j_sInt  = (j_zlm-rlen)+j_zlm1

  ! if(inter == 1) then
  !   s = (j_zlm-rlen)+j_zlm1
  !       call k2coll_calcIonLoss(p,rlen,j_exenergy,j_anuc,j_zatom,j_rho,m_dpodx)
  !       p = p-m_dpodx*s

  !   j_dpop = (p-j_p0)/j_p0

  !   return
  ! end if

  ! ! Now treat the other types of interaction, as determined by ICHOIX:

  ! ! Nuclear-Elastic:          inter = 2
  ! ! pp Elastic:               inter = 3
  ! ! Single-Diffractive:       inter = 4    (changes momentum p)
  ! ! Coulomb:                  inter = 5

  ! ! As the single-diffractive interaction changes the momentum, save input momentum in p1.
  ! p1 = p
  ! ! Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.
  ! t = k2coll_gettran(inter,p,j_bn,j_cgen, j_ecmsq, j_xln15s, j_bpp)

  ! ! Tetat calculates from the rms transverse momentum transfer in
  ! ! monte-carlo fashion the angle changes for x and z planes. The
  ! ! angle change is proportional to SQRT(t) and 1/p, as expected.
  ! call k2coll_tetat(t,p,dxp,dzp)
  ! ! Apply angle changes
  ! j_xp = j_xp+dxp
  ! j_zp = j_zp+dzp

  ! ! Treat single-diffractive scattering.
  ! if(inter == 4) then

  !   ! added update for s
  !   s    = (j_zlm-rlen)+j_zlm1
  !   xpsd = dxp
  !   zpsd = dzp
  !   psd  = p1

  !   ! Add this code to get the momentum transfer also in the calling routine
  !   j_dpop = (p-j_p0)/j_p0
  ! end if

  ! ! Calculate the remaining interaction length and close the iteration loop.
  ! rlen = rlen-j_zlm1
  ! goto 10

! end subroutine k2coll_jaw

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
subroutine k2coll_mcs(s, mc_radl, mc_zlm1, mc_p0, mc_x, mc_xp, mc_z, mc_zp, mc_dpop)

  !use coll_materials

  real(kind=fPrec), intent(inout) :: s
  real(kind=fPrec), intent(in)    :: mc_radl
  real(kind=fPrec), intent(in)    :: mc_zlm1
  real(kind=fPrec), intent(in)    :: mc_p0

  real(kind=fPrec), intent(inout) :: mc_x
  real(kind=fPrec), intent(inout) :: mc_xp
  real(kind=fPrec), intent(inout) :: mc_z
  real(kind=fPrec), intent(inout) :: mc_zp
  real(kind=fPrec), intent(inout) :: mc_dpop

  real(kind=fPrec) theta,rlen0,rlen,ae,be,rad_len

  real(kind=fPrec), parameter :: h   = 0.001_fPrec
  real(kind=fPrec), parameter :: dh  = 0.0001_fPrec
  real(kind=fPrec), parameter :: bn0 = 0.4330127019_fPrec

  ! radl_mat = mc_radl
  theta    = 13.6e-3_fPrec/(mc_p0*(one+mc_dpop)) ! dpop   = (p - p0)/p0
  rad_len  = mc_radl


  mc_x     = (mc_x/theta)/mc_radl
  mc_xp    = mc_xp/theta
  mc_z     = (mc_z/theta)/mc_radl
  mc_zp    = mc_zp/theta
  rlen0 = mc_zlm1/mc_radl
  rlen  = rlen0

10 continue
  ae = bn0*mc_x
  be = bn0*mc_xp

    call k2coll_soln3(ae,be,dh,rlen,s)
    if(s < h) s = h
    call k2coll_scamcs(mc_x,mc_xp,s)
    if(mc_x <= zero) then
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
    call k2coll_scamcs(mc_z,mc_zp,s)
    s  = s*mc_radl
  mc_x  = (mc_x*theta)*mc_radl
  mc_xp = mc_xp*theta
  mc_z  = (mc_z*theta)*mc_radl
  mc_zp = mc_zp*theta

end subroutine k2coll_mcs

!>
!! k2coll_calcIonLoss(IS,PC,DZ,il_anuc,il_zatom,il_rho,EnLo)
!! subroutine for the calculazion of the energy loss by ionization
!! Either mean energy loss from Bethe-Bloch, or higher energy loss, according to finite probability from cross section
!! Written by DM for crystals, introduced in main code by RB
!! Updated and improved for numerical stability by VKBO
!<
! subroutine k2coll_calcIonLoss(PC, DZ, il_exenergy, il_anuc, il_zatom, il_rho, EnLo)

!   use mod_ranlux
!   !use coll_materials
!   use mathlib_bouncer
!   use physical_constants

!   !integer,          intent(in)  :: IS           ! IS material ID
!   real(kind=fPrec), intent(in)  :: PC           ! PC momentum in GeV
!   real(kind=fPrec), intent(in)  :: DZ           ! DZ length traversed in material (meters)
!   real(kind=fPrec), intent(inout)  :: il_exenergy  ! il_exenergy
!   real(kind=fPrec), intent(in)  :: il_anuc      ! il_anuc 
!   real(kind=fPrec), intent(in)  :: il_zatom     ! il_zatom
!   real(kind=fPrec), intent(in)  :: il_rho       ! il_rho
!   real(kind=fPrec), intent(out) :: EnLo         ! EnLo energy loss in GeV/meter

!   real(kind=fPrec) exEn,thl,Tt,cs_tail,prob_tail,enr,mom,betar,gammar,bgr,kine,Tmax,plen
!   real(kind=fPrec), parameter :: k = 0.307075_fPrec ! Constant in front of Bethe-Bloch [MeV g^-1 cm^2]

!   mom    = PC*c1e3                     ! [GeV/c] -> [MeV/c]
!   enr    = (mom*mom + pmap*pmap)**half ! [MeV]
!   gammar = enr/pmap
!   betar  = mom/enr
!   bgr    = betar*gammar
!   kine   = ((two*pmae)*bgr)*bgr

!   ! Mean excitation energy
!   exEn = il_exenergy*c1e3 ! [MeV]

!   ! Tmax is max energy loss from kinematics
!   Tmax = kine/(one + (two*gammar)*(pmae/pmap) + (pmae/pmap)**2) ! [MeV]

!   ! Plasma energy - see PDG 2010 table 27.1
!   plen = (((il_rho*il_zatom)/il_anuc)**half)*28.816e-6_fPrec ! [MeV]

!   ! Calculate threshold energy
!   ! Above this threshold, the cross section for high energy loss is calculated and then
!   ! a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
!   ! below threshold, only the standard Bethe-Bloch is used (all particles get average energy loss)

!   ! thl is 2*width of Landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
!   thl = ((((four*(k*il_zatom))*DZ)*c1e2)*il_rho)/(il_anuc*betar**2) ! [MeV]

!   ! Bethe-Bloch mean energy loss
!   EnLo = ((k*il_zatom)/(il_anuc*betar**2)) * ( &
!     half*log_mb((kine*Tmax)/(exEn*exEn)) - betar**2 - log_mb(plen/exEn) - log_mb(bgr) + half &
!   )
!   EnLo = ((EnLo*il_rho)*c1m1)*DZ ! [GeV]

!   ! Threshold Tt is Bethe-Bloch + 2*width of Landau distribution
!   Tt = EnLo*c1e3 + thl ! [MeV]

!   ! Cross section - see Alfredo's presentation for derivation
!   cs_tail = ((k*il_zatom)/(il_anuc*betar**2)) * ( &
!     half*((one/Tt)-(one/Tmax)) - (log_mb(Tmax/Tt)*betar**2)/(two*Tmax) + (Tmax-Tt)/((four*gammar**2)*pmap**2) &
!   )

!   ! Probability of being in tail: cross section * density * path length
!   prob_tail = ((cs_tail*il_rho)*DZ)*c1e2

!   ! Determine based on random number if tail energy loss occurs.
!   if(coll_rand() < prob_tail) then
!     EnLo = ((k*il_zatom)/(il_anuc*betar**2)) * ( &
!       half*log_mb((kine*Tmax)/(exEn*exEn)) - betar**2 - log_mb(plen/exEn) - log_mb(bgr) + &
!       half + TMax**2/((eight*gammar**2)*pmap**2) &
!     )
!     EnLo = (EnLo*il_rho)*c1m1 ! [GeV/m]
!   else
!     ! If tail energy loss does not occur, just use the standard Bethe-Bloch
!     EnLo = EnLo/DZ  ! [GeV/m]
!   end if


! end subroutine k2coll_calcIonLoss

!>
!! k2coll_tetat(t,p,tx,tz)
!! Generate sine and cosine of an angle uniform in [0,2pi](see RPP)
! !<
! subroutine k2coll_tetat(t, p, tx, tz)

!   use mod_ranlux, only : coll_rand

!   real(kind=fPrec), intent(in)  :: t
!   real(kind=fPrec), intent(in)  :: p
!   real(kind=fPrec), intent(out) :: tx
!   real(kind=fPrec), intent(out) :: tz

!   real(kind=fPrec) va,vb,va2,vb2,r2,teta

!   teta = sqrt(t)/p
! 10 continue
!   va  = two*coll_rand() - one
!   vb  = coll_rand()
!   va2 = va**2
!   vb2 = vb**2
!   r2  = va2 + vb2
!   if(r2 > one) goto 10
!   tx  = (teta*((two*va)*vb))/r2
!   tz  = (teta*(va2 - vb2))/r2

! end subroutine k2coll_tetat

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
! HACK coll_zatom is used as global; it cannot be passed as a function argument
! because this function is passed into funlxp which builds the random distribution
! However, the latter expects a function with one argument
real(kind=fPrec) function k2coll_ruth(t)

  use mathlib_bouncer
  !use coll_materials

  real(kind=fPrec), intent(in) :: t
  !real(kind=fPrec), intent(in) :: ru_emr

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
! XMAT AND MAT ARE THE SAME, EQUAL TO MAIN MATID
! real(kind=fPrec) function k2coll_gettran(inter, p, tt_bn, tt_cgen, tt_ecmsq, tt_xln15s, tt_bpp)

!   use mathlib_bouncer
!   use mod_ranlux, only : coll_rand
!   use mod_funlux, only : funlux  
!   !use coll_materials

!   integer,          intent(in)    :: inter
!   real(kind=fPrec), intent(inout) :: p
!   real(kind=fPrec), intent(in)    :: tt_bn
!   real(kind=fPrec), intent(in)    :: tt_cgen(200)
!   real(kind=fPrec), intent(in)    :: tt_xln15s
!   real(kind=fPrec), intent(in)    :: tt_ecmsq
!   real(kind=fPrec), intent(in)    :: tt_bpp

!   real(kind=fPrec) xm2,bsd,xran(1)

!   ! Neither if-statements below have an else, so defaulting function return to zero.
!   k2coll_gettran = zero


!   select case(inter)
!   case(2) ! Nuclear Elastic
!     k2coll_gettran = (-one*log_mb(coll_rand()))/tt_bn
!   case(3) ! pp Elastic
!     k2coll_gettran = (-one*log_mb(coll_rand()))/tt_bpp
!   case(4) ! Single Diffractive
!     xm2 = exp_mb(coll_rand() * tt_xln15s)
!     p   = p * (one - xm2/tt_ecmsq)
!     if(xm2 < two) then
!       bsd = two * tt_bpp
!     else if(xm2 >= two .and. xm2 <= five) then
!       bsd = ((106.0_fPrec - 17.0_fPrec*xm2)*tt_bpp)/36.0_fPrec
!     else
!       bsd = (seven*tt_bpp)/12.0_fPrec
!     end if
!     k2coll_gettran = (-one*log_mb(coll_rand()))/bsd
!   case(5) ! Coulomb
!     call funlux(tt_cgen(1), xran, 1)
!     k2coll_gettran = xran(1)
!   end select

! end function k2coll_gettran

!>
!! k2coll_ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
! integer function k2coll_ichoix(ich_cprob)

!   use mod_ranlux, only : coll_rand 
!   !use coll_materials
!   real(kind=fPrec), intent(in) :: ich_cprob(0:5)      ! Cprob to choose an interaction in iChoix

!   !integer, intent(in) :: ma
!   integer i
!   real(kind=fPrec) aran

!   aran = coll_rand()
!   i    = 1
! 10 continue
!   if(aran > ich_cprob(i)) then
!     i = i+1
!     goto 10
!   end if

!   k2coll_ichoix = i

! end function k2coll_ichoix


end module coll_k2
