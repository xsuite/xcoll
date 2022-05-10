!subroutine pyk2_init(n_alloc, random_generator_seed)
subroutine pyk2_init(random_generator_seed)
  use mod_ranlux ,       only : rluxgo     ! for ranlux init
  use coll_common ,      only : rnd_seed !, coll_expandArrays
  
  implicit none

  integer, intent(in)          :: random_generator_seed

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)
end subroutine


subroutine initialise_random(random_generator_seed, cgen, zatom, emr, hcut)
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

  if(random_generator_seed .ge. 0) then
        call rluxgo(3, random_generator_seed, 0, 0)
  end if

  zatom_curr = zatom
  emr_curr = emr
  call funlxp(k2coll_ruth, cgen(1), tlcut, hcut)

end subroutine


subroutine pyk2_crystal( &
  val_part_hit, &
  val_part_abs, &
  val_part_impact, &
  val_part_indiv, &
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

real(kind=fPrec) sImp

! needs to be passed from cry_startElement
integer cry_proc, cry_proc_prev, cry_proc_tmp
cry_proc = -1
cry_proc_prev = -1
cry_proc_tmp = -1

call cry_doCrystal(x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,val_part_hit,&
val_part_abs,val_part_impact,val_part_indiv,c_length,run_exenergy,run_anuc,run_zatom,run_emr,run_rho,&
run_hcut,run_bnref,run_csref0,run_csref1,run_csref4,run_csref5,run_dlri,run_dlyi,&
run_eUm,run_ai,run_collnt,run_bn, cry_proc, cry_proc_prev, cry_proc_tmp)

end subroutine



subroutine pyk2_jaw( &
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
  p0, &
  nabs, &
  s, &
  zlm, &
  x, &
  xp, &
  z, &
  zp, &
  dpop)

use coll_k2     ! for scattering

implicit none


! ############################
! ## variables declarations ##
! ############################

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
real(kind=8),  intent(inout) :: p0
integer,          intent(inout) :: nabs
real(kind=8),    intent(inout) :: s
real(kind=8),    intent(inout) :: zlm
real(kind=8),    intent(inout) :: x
real(kind=8),    intent(inout) :: xp
real(kind=8),    intent(inout) :: z
real(kind=8),    intent(inout) :: zp
real(kind=8),    intent(inout) :: dpop  

real(kind=8) xInt,xpInt,yInt,ypInt,sInt

call k2coll_jaw(s,nabs,run_exenergy,run_anuc,run_zatom,run_rho,run_radl,&
                    run_cprob,run_xintl,run_bn,run_cgen,run_ecmsq,run_xln15s,run_bpp,zlm,p0,&
                    x,xp,z,zp,dpop,xInt,xpInt,yInt,ypInt,sInt)
end subroutine

