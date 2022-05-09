!subroutine pyk2_init(n_alloc, random_generator_seed)
subroutine pyk2_init(random_generator_seed)
  
  !use floatPrecision
  !use numerical_constants
  ! use crcoall    NODIG ??
  !use mod_alloc ,        only : alloc      !to allocate partID etc
  use mod_ranlux ,       only : rluxgo     ! for ranlux init
  use coll_common ,      only : rnd_seed !, coll_expandArrays
  !use coll_materials ! for collmat_init
  !use coll_k2        ! for scattering

  implicit none

  ! integer, intent(in)          :: n_alloc
  integer, intent(in)          :: random_generator_seed

  ! Set default values for collimator materials
 ! call collmat_init

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

  ! call coll_expandArrays(n_alloc)
end subroutine


subroutine initialise_random(random_generator_seed, cgen, zatom, emr, hcut)
  ! use parpro ,           only : npart
  ! use mod_common ,       only : napx !, aa0
  !use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init
  use mod_funlux
  use coll_k2,           only : k2coll_ruth, zatom_curr, emr_curr

  implicit none

  integer, intent(in)         :: random_generator_seed
  real(kind=8), intent(inout) :: cgen(200)
  real(kind=8), intent(in)    :: zatom
  real(kind=8), intent(in)    :: emr
  real(kind=8), intent(in)    :: hcut

  real(kind=8), parameter     :: tlcut = 0.0009982
  !integer j

  ! ####################
  ! ## initialisation ##
  ! ####################
  !character(len=:),    allocatable   :: numpart
  !numpart="20000"
  !read(numpart,*) napx
  ! npart=num_particles

  if(random_generator_seed .ge. 0) then
        call rluxgo(3, random_generator_seed, 0, 0)
  end if

  ! npart = num_particles
  
  ! napx=npart  ! this decreases after absorptions!

  ! Prepare for Rutherford differential distribution
  !mcurr = mat ! HACK> mcurr is global, and coll_zatom too which is used inside k2coll_ruth
  zatom_curr = zatom
  emr_curr = emr
  call funlxp(k2coll_ruth, cgen(1), tlcut, hcut)

end subroutine

! subroutine crystal_collimation()


! subroutine pyk2_run( &
!                     val_part_hit, &
!                     val_part_abs, &
!                     val_part_impact, &
!                     val_part_indiv, &
!                     val_part_linteract, &
!                     val_nabs_type, &
!                     val_linside, &
!                     run_exenergy, &
!                     run_anuc, &
!                     run_zatom, &
!                     run_emr, &
!                     run_rho, &
!                     run_hcut, &
!                     run_bnref, &
!                     run_csref0, &
!                     run_csref1, &
!                     run_csref4, &
!                     run_csref5, &
!                     run_radl, &
!                     run_dlri, & 
!                     run_dlyi, &
!                     run_eUm, &
!                     run_ai, &
!                     run_collnt, &
!                     run_cprob, &
!                     run_xintl, &
!                     run_bn, &
!                     run_ecmsq, &
!                     run_xln15s, &
!                     run_bpp, &
!                     run_cgen, &
!                     is_crystal, &
!                     c_length, &
!                     length, &
!                     p0, &
!                     nhit, &
!                     nabs, &
!                     fracab, &
!                     isImp, &
!                     s, &
!                     zlm, &
!                     x, &
!                     xp, &
!                     xp_in0, &
!                     z, &
!                     zp, &
!                     p, &
!                     sp, &
!                     dpop, &
!                     x_in0)

!   use coll_k2        ! for scattering

!   implicit none


!   ! ############################
!   ! ## variables declarations ##
!   ! ############################

!   integer(kind=4)  , intent(inout) :: val_part_hit
!   integer(kind=4)  , intent(inout) :: val_part_abs
!   real(kind=8) , intent(inout) :: val_part_impact
!   real(kind=8) , intent(inout) :: val_part_indiv
!   real(kind=8) , intent(inout) :: val_part_linteract
!   integer(kind=4)  , intent(inout) :: val_nabs_type
!   logical(kind=4)  , intent(inout) :: val_linside

!   real(kind=8)     , intent(inout) :: run_exenergy
!   real(kind=8)     , intent(in) :: run_anuc
!   real(kind=8)     , intent(in) :: run_zatom
!   real(kind=8)     , intent(in) :: run_emr
!   real(kind=8)     , intent(in) :: run_rho
!   real(kind=8)     , intent(in) :: run_hcut
!   real(kind=8)     , intent(in) :: run_bnref

!   real(kind=8)     , intent(in) :: run_csref0
!   real(kind=8)     , intent(in) :: run_csref1
!   real(kind=8)     , intent(in) :: run_csref4
!   real(kind=8)     , intent(in) :: run_csref5

!   real(kind=8)     , intent(in) :: run_radl
!   real(kind=8)     , intent(in) :: run_dlri
!   real(kind=8)     , intent(in) :: run_dlyi
!   real(kind=8)     , intent(in) :: run_eUm
!   real(kind=8)     , intent(in) :: run_ai
!   real(kind=8)     , intent(in) :: run_collnt
!   real(kind=8)     , intent(in) :: run_cprob(0:5)
!   real(kind=8)     , intent(in) :: run_xintl
!   real(kind=8)     , intent(inout) :: run_bn
!   real(kind=8)     , intent(in) :: run_ecmsq
!   real(kind=8)     , intent(in) :: run_xln15s
!   real(kind=8)     , intent(in) :: run_bpp
!   real(kind=8)     , intent(in) :: run_cgen(200)

!   logical(kind=4)  , intent(in) :: is_crystal
!   real(kind=8) ,    intent(in) :: c_length
!   real(kind=8),  intent(inout) :: length
!   real(kind=8),  intent(inout) :: p0

!   integer,          intent(inout) :: nhit
!   integer,          intent(inout) :: nabs
!   real(kind=8), intent(inout) :: fracab

!   logical(kind=4), intent(inout) :: isImp
!   real(kind=8),    intent(inout) :: s
!   real(kind=8),    intent(inout) :: zlm
!   real(kind=8),    intent(inout) :: sp

!   real(kind=8),    intent(inout) :: x
!   real(kind=8),    intent(inout) :: xp
!   real(kind=8),    intent(inout) :: xp_in0
!   real(kind=8),    intent(inout) :: z
!   real(kind=8),    intent(inout) :: zp
!   real(kind=8),    intent(inout) :: p
!   real(kind=8),    intent(inout) :: dpop  
!   real(kind=8),    intent(inout) :: x_in0


!   ! needs to be passed from cry_startElement
!   integer cry_proc, cry_proc_prev, cry_proc_tmp
!   cry_proc = -1
!   cry_proc_prev = -1
!   cry_proc_tmp = -1

!   call k2coll_collimate( &
!      run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, &
!      run_csref0, run_csref1, run_csref4, run_csref5, run_radl, run_dlri, &
!      run_dlyi, run_eUm, run_ai, run_collnt, run_cprob, run_xintl, run_bn, &
!      run_ecmsq, run_xln15s, run_bpp, run_cgen, is_crystal, &
!      c_length, val_part_hit, val_part_abs, &
!      val_part_impact, val_part_indiv, val_part_linteract, &
!      val_nabs_type, val_linside, length, p0, nhit, &
!      nabs, fracab, &
!     cry_proc, cry_proc_prev, cry_proc_tmp, &
!     isImp, s, zlm, sp, &
!     x, xp, xp_in0, z, zp, p, dpop, x_in0)

! end subroutine



subroutine pyk2_crystal( &
  val_part_hit, &
  val_part_abs, &
  val_part_impact, &
  val_part_indiv, &
  val_part_linteract, &
  run_exenergy, &
  run_anuc, &
  run_zatom, &
  run_emr, &
  run_rho, &
  run_hcut, &
  run_bnref, &
  run_csref0, &
  run_csref1, &
  run_csref4, &
  run_csref5, &
  run_dlri, & 
  run_dlyi, &
  run_eUm, &
  run_ai, &
  run_collnt, &
  run_bn, &
  c_length, &
  nhit, &
  nabs, &
  isImp, &
  s, &
  zlm, &
  x, &
  xp, &
  xp_in0, &
  z, &
  zp, &
  p, &
  x_in0)

use coll_k2        ! for scattering
use coll_crystal, only : cry_doCrystal

implicit none


! ############################
! ## variables declarations ##
! ############################

integer(kind=4)  , intent(inout) :: val_part_hit
integer(kind=4)  , intent(inout) :: val_part_abs
real(kind=8) , intent(inout) :: val_part_impact
real(kind=8) , intent(inout) :: val_part_indiv
real(kind=8) , intent(inout) :: val_part_linteract

real(kind=8)     , intent(inout) :: run_exenergy
real(kind=8)     , intent(in) :: run_anuc
real(kind=8)     , intent(in) :: run_zatom
real(kind=8)     , intent(in) :: run_emr
real(kind=8)     , intent(in) :: run_rho
real(kind=8)     , intent(in) :: run_hcut
real(kind=8)     , intent(in) :: run_bnref

real(kind=8)     , intent(in) :: run_csref0
real(kind=8)     , intent(in) :: run_csref1
real(kind=8)     , intent(in) :: run_csref4
real(kind=8)     , intent(in) :: run_csref5

real(kind=8)     , intent(in) :: run_dlri
real(kind=8)     , intent(in) :: run_dlyi
real(kind=8)     , intent(in) :: run_eUm
real(kind=8)     , intent(in) :: run_ai
real(kind=8)     , intent(in) :: run_collnt
real(kind=8)     , intent(inout) :: run_bn
real(kind=8) ,    intent(in) :: c_length

integer,          intent(inout) :: nhit
integer,          intent(inout) :: nabs
logical(kind=4), intent(inout) :: isImp
real(kind=8),    intent(inout) :: s
real(kind=8),    intent(inout) :: zlm

real(kind=8),    intent(inout) :: x
real(kind=8),    intent(inout) :: xp
real(kind=8),    intent(inout) :: xp_in0
real(kind=8),    intent(inout) :: z
real(kind=8),    intent(inout) :: zp
real(kind=8),    intent(inout) :: p
real(kind=8),    intent(inout) :: x_in0

real(kind=fPrec) xOut,xpOut,yOut,ypOut,sImp,sOut

! needs to be passed from cry_startElement
integer cry_proc, cry_proc_prev, cry_proc_tmp
cry_proc = -1
cry_proc_prev = -1
cry_proc_tmp = -1

call cry_doCrystal(x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,val_part_hit,&
val_part_abs,val_part_impact,val_part_indiv,c_length,run_exenergy,run_anuc,run_zatom,run_emr,run_rho,&
run_hcut,run_bnref,run_csref0,run_csref1,run_csref4,run_csref5,run_dlri,run_dlyi,&
run_eUm,run_ai,run_collnt,run_bn, cry_proc, cry_proc_prev, cry_proc_tmp)

if(nabs /= 0) then
val_part_abs = 1
val_part_linteract = zlm
end if

sImp  = (s - c_length) + sImp
sOut  = s
xOut  = x
xpOut = xp
yOut  = z
ypOut = zp

end subroutine




subroutine pyk2_jaw( &
  val_part_hit, &
  val_part_abs, &
  val_part_impact, &
  val_part_indiv, &
  val_part_linteract, &
  val_nabs_type, &
  val_linside, &
  run_exenergy, &
  run_anuc, &
  run_zatom, &
  run_rho, &
  run_radl, &
  run_cprob, &
  run_xintl, &
  run_bn, &
  run_ecmsq, &
  run_xln15s, &
  run_bpp, &
  run_cgen, &
  length, &
  p0, &
  nhit, &
  nabs, &
  fracab, &
  isImp, &
  s, &
  zlm, &
  x, &
  xp, &
  z, &
  zp, &
  sp, &
  dpop)

use coll_k2        ! for scattering

implicit none


! ############################
! ## variables declarations ##
! ############################

integer(kind=4)  , intent(inout) :: val_part_hit
integer(kind=4)  , intent(inout) :: val_part_abs
real(kind=8) , intent(inout) :: val_part_impact
real(kind=8) , intent(inout) :: val_part_indiv
real(kind=8) , intent(inout) :: val_part_linteract
integer(kind=4)  , intent(inout) :: val_nabs_type
logical(kind=4)  , intent(inout) :: val_linside

real(kind=8)     , intent(inout) :: run_exenergy
real(kind=8)     , intent(in) :: run_anuc
real(kind=8)     , intent(in) :: run_zatom
real(kind=8)     , intent(in) :: run_rho

real(kind=8)     , intent(in) :: run_radl
real(kind=8)     , intent(in) :: run_cprob(0:5)
real(kind=8)     , intent(in) :: run_xintl
real(kind=8)     , intent(inout) :: run_bn
real(kind=8)     , intent(in) :: run_ecmsq
real(kind=8)     , intent(in) :: run_xln15s
real(kind=8)     , intent(in) :: run_bpp
real(kind=8)     , intent(in) :: run_cgen(200)

real(kind=8),  intent(inout) :: length
real(kind=8),  intent(inout) :: p0

integer,          intent(inout) :: nhit
integer,          intent(inout) :: nabs
real(kind=8), intent(inout) :: fracab

logical(kind=4), intent(inout) :: isImp
real(kind=8),    intent(inout) :: s
real(kind=8),    intent(inout) :: zlm
real(kind=8),    intent(inout) :: sp

real(kind=8),    intent(inout) :: x
real(kind=8),    intent(inout) :: xp
real(kind=8),    intent(inout) :: z
real(kind=8),    intent(inout) :: zp
real(kind=8),    intent(inout) :: dpop  

real(kind=fPrec) s_impact,drift_length
real(kind=fPrec) xOut,xpOut,yOut,ypOut,sImp,sOut
real(kind=fPrec) xInt,xpInt,yInt,ypInt,sInt

if(x >= zero) then
  ! Particle hits collimator and we assume interaction length ZLM equal
  ! to collimator length (what if it would leave collimator after
  ! small length due to angle???)
  zlm  = length
  val_part_impact = x
  val_part_indiv  = xp
else if(xp <= zero) then
  ! Particle does not hit collimator. Interaction length ZLM is zero.
  zlm = zero
else
  ! Calculate s-coordinate of interaction point
  s = (-one*x)/xp
  if(s <= zero) then
    ! write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
    !call prror
  end if
  if(s < length) then
    zlm       = length - s
    val_part_impact = zero
    val_part_indiv  = xp
  else
    zlm = zero
  end if
end if

! First do the drift part
! DRIFT PART
drift_length = length - zlm
if(drift_length > zero) then
  x  = x  + xp* drift_length
  z  = z  + zp * drift_length
  sp = sp + drift_length
end if

! Now do the scattering part
if(zlm > zero) then
  if(.not.val_linside) then
    ! first time particle hits collimator: entering jaw
    val_linside = .true.
  end if

  s_impact = sp
  nhit = nhit + 1
  call k2coll_jaw(s,nabs,run_exenergy,run_anuc,run_zatom,run_rho,run_radl,&
                    run_cprob,run_xintl,run_bn,run_cgen,run_ecmsq,run_xln15s,run_bpp,zlm,p0,&
                    x,xp,z,zp,dpop,xInt,xpInt,yInt,ypInt,sInt)
  val_nabs_type = nabs
  val_part_hit  = 1

  isImp = .true.
  sImp  = s_impact
  sOut  = (s+sp)
  xOut  = x
  xpOut = xp
  yOut  = z
  ypOut = zp

  ! Writeout should be done for both inelastic and single diffractive. doing all transformations
  ! in x_flk and making the set to 99.99 mm conditional for nabs=1
  if(nabs == 1 .or. nabs == 4) then
    ! Transform back to lab system for writeout.
    ! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

    ! Finally, the actual coordinate change to 99 mm
    if(nabs == 1) then
      fracab  = fracab + 1
      x       = 99.99e-3_fPrec
      z       = 99.99e-3_fPrec
      val_part_linteract = zlm
      val_part_abs = 1
    end if
  end if
end if ! Collimator jaw interaction

if(nabs /= 1 .and. zlm > zero) then
  ! Do the rest drift, if particle left collimator early
  drift_length = (length-(s+sp))
  if(drift_length > c1m15) then
    val_linside = .false.
    x  = x  + xp * drift_length
    z  = z  + zp * drift_length
    sp = sp + drift_length
  end if
  val_part_linteract = zlm - drift_length
end if

end subroutine

