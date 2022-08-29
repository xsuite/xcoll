subroutine pyk2_init(random_generator_seed)
  use mod_ranlux ,       only : rluxgo     ! for ranlux init
  use coll_common ,      only : rnd_seed !, rcx, rcxp, rcy, rcyp, rcp, rcs, &

  implicit none

  integer, intent(in)          :: random_generator_seed

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

end subroutine

subroutine pyk2_startcry(c_length, new_length, cryTilt, cryBend, cryThick, cryXDim, &
                         cryYDim, cryOrient, cryMiscut)
  use coll_crystal, only: cry_startElement

  real(kind=8), intent(in)    :: c_length     ! Collimator length in m
  real(kind=8), intent(inout) :: new_length
  real(kind=8), intent(in)    :: cryTilt
  real(kind=8), intent(in)    :: cryBend
  real(kind=8), intent(in)    :: cryThick
  real(kind=8), intent(in)    :: cryXDim
  real(kind=8), intent(in)    :: cryYDim
  integer,      intent(in)    :: cryOrient
  real(kind=8), intent(in)    :: cryMiscut
  
  call cry_startElement(c_length, new_length, cryTilt, cryBend, cryThick, cryXDim, &
                        cryYDim, cryOrient, cryMiscut)
end subroutine


subroutine pyk2_docrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs, &
  lhit,part_abs,impact,indiv,c_length,exenergy,rho,anuc,zatom,emr,dlri,dlyi,ai,eUm,collnt,&
  hcut,csref0,csref1,csref4,csref5,bnref,csect)
  
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
  real(kind=8), intent(in)    :: rho
  real(kind=8), intent(in)    :: anuc
  real(kind=8), intent(in)    :: zatom
  real(kind=8), intent(in)    :: emr

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
  real(kind=8), intent(in)    :: bnref
  real(kind=8), intent(in)    :: csect

  call cry_doCrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs,lhit,part_abs,impact,indiv,&
                    c_length,exenergy,rho,anuc,zatom,emr,dlri,dlyi,ai,eUm,collnt,&
                    hcut,bnref,csref0,csref1,csref4,csref5,csect)

end subroutine
