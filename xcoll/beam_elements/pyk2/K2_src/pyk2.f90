subroutine pyk2_init(n_alloc, random_generator_seed)
  use floatPrecision
  use numerical_constants
  ! use crcoall    NODIG ??
!  use parpro ,           only : npart
!  use mod_alloc ,        only : alloc      ! to allocate partID etc
 ! use mod_common ,       only : iexact, napx, unit208, aa0
 ! use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed !, rcx, rcxp, rcy, rcyp, rcp, rcs, &
!                                coll_expandArrays
  use coll_materials ! for collmat_init
!  use coll_k2        ! for scattering
  use coll_crystal, only: cry_init

  implicit none

  integer, intent(in)          :: n_alloc
  integer, intent(in)          :: random_generator_seed

  ! Set default values for collimator materials
  call collmat_init
  call cry_init

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

!  call coll_expandArrays(n_alloc)
!  call alloc(naa, n_alloc, aa0, "naa")
!  call alloc(partID, n_alloc, 0, "partID")
!  call alloc(parentID, n_alloc, 0, "parentID")
!  call alloc(pairID, 2, n_alloc, 0, "pairID")

end subroutine

subroutine pyk2_startcry(c_length, new_length, c_rotation, cryTilt, cryBend, cryThick, cryXDim, &
                         cryYDim, cryOrient, cryMiscut)
  use coll_crystal, only: cry_startElement

  real(kind=8), intent(in)    :: c_length     ! Collimator length in m
  real(kind=8), intent(inout) :: new_length
  real(kind=8), intent(in)    :: c_rotation   ! Collimator rotation angle vs vertical in radians
  real(kind=8), intent(in)    :: cryTilt
  real(kind=8), intent(in)    :: cryBend
  real(kind=8), intent(in)    :: cryThick
  real(kind=8), intent(in)    :: cryXDim
  real(kind=8), intent(in)    :: cryYDim
  real(kind=8), intent(in)    :: cryOrient
  real(kind=8), intent(in)    :: cryMiscut
  
  call cry_startElement(c_length, new_length, c_rotation, cryTilt, cryBend, cryThick, cryXDim, &
                        cryYDim, cryOrient, cryMiscut)
end subroutine


subroutine pyk2_docrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs, &
  lhit,part_abs,impact,indiv,c_length,exenergy,zatom,rho,anuc,dlri,dlyi,ai,eUm,collnt,&
  hcut,csref0,csref1,csref4,csref5,nmat,bnref,csect)
  
  use coll_crystal, only: cry_doCrystal
  use parpro
  ! use coll_common, only : cry_proc, cry_proc_prev, cry_proc_tmp
  use mathlib_bouncer

  real(kind=8), intent(inout) :: x,xp
  real(kind=8), intent(inout) :: z,zp
  real(kind=8), intent(inout) :: s,p
  real(kind=8), intent(inout) :: x0,xp0
  real(kind=8), intent(inout) :: zlm,s_imp
  integer,      intent(inout) :: nhit,nabs
  integer,      intent(inout) :: lhit
  integer,      intent(inout) :: part_abs
  real(kind=8), intent(inout) :: impact
  real(kind=8), intent(inout) :: indiv
  real(kind=8), intent(in)    :: c_length
  logical,      intent(inout) :: isImp

  real(kind=8), intent(in)    :: exenergy
  real(kind=8), intent(in)    :: zatom
  real(kind=8), intent(in)    :: rho
  real(kind=8), intent(in)    :: anuc

  real(kind=8), intent(in)    :: dlri
  real(kind=8), intent(in)    :: dlyi
  real(kind=8), intent(in)    :: ai
  real(kind=8), intent(in)    :: eUm
  real(kind=8), intent(in)    :: collnt

  real(kind=8), intent(in)    :: hcut
  real(kind=8), intent(in)    :: csref0
  real(kind=8), intent(in)    :: csref1
  real(kind=8), intent(in)    :: csref4
  real(kind=8), intent(in)    :: csref5
  integer,      intent(in)    :: nmat
  real(kind=8), intent(in)    :: bnref
  real(kind=8), intent(in)    :: csect
  
  call cry_doCrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs,lhit,part_abs,impact,indiv,&
                    c_length,exenergy,zatom,rho,anuc,dlri,dlyi,ai,eUm,collnt,&
                    hcut,bnref,csref0,csref1,csref4,csref5,nmat,csect)

end subroutine
