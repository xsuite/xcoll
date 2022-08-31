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

! subroutine pyk2_startcry(c_length, new_length, cryTilt, cryBend, cryThick, cryXDim, &
!                          cryYDim, cryOrient, cryMiscut)
!   use coll_crystal, only: cry_startElement

!   real(kind=8), intent(in)    :: c_length     ! Collimator length in m
!   real(kind=8), intent(inout) :: new_length
!   real(kind=8), intent(in)    :: cryTilt
!   real(kind=8), intent(in)    :: cryBend
!   real(kind=8), intent(in)    :: cryThick
!   real(kind=8), intent(in)    :: cryXDim
!   real(kind=8), intent(in)    :: cryYDim
!   integer,      intent(in)    :: cryOrient
!   real(kind=8), intent(in)    :: cryMiscut
  
!   call cry_startElement(c_length, new_length, cryTilt, cryBend, cryThick, cryXDim, &
!                         cryYDim, cryOrient, cryMiscut)
! end subroutine


subroutine pyk2_docrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs, &
  lhit,part_abs,impact,indiv,c_length,exenergy,rho,anuc,zatom,emr,dlri,dlyi,ai,eUm,collnt,&
  hcut,csref0,csref1,csref4,csref5,bnref,csect,cry_tilt,c_rcurv,&
  c_alayer,c_xmax,c_ymax,c_orient,c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,&
  c_spTilt,c_cnTilt,c_snTilt,iProc,n_chan,n_VR,n_amorphous)
  
  use coll_crystal, only: cry_doCrystal
  use parpro
  ! use coll_common, only : cry_proc, cry_proc_prev, cry_proc_tmp
  use mathlib_bouncer

  implicit none

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

  real(kind=8), intent(in)    :: cry_tilt
  real(kind=8), intent(in)    :: c_rcurv
  real(kind=8), intent(in)    :: c_alayer
  real(kind=8), intent(in)    :: c_xmax
  real(kind=8), intent(in)    :: c_ymax
  integer,      intent(in)    :: c_orient
  real(kind=8), intent(in)    :: c_miscut
  real(kind=8), intent(in)    :: cry_bend
  real(kind=8), intent(in)    :: c_cBend
  real(kind=8), intent(in)    :: c_sBend
  real(kind=8), intent(in)    :: c_cpTilt
  real(kind=8), intent(in)    :: c_spTilt
  real(kind=8), intent(in)    :: c_cnTilt
  real(kind=8), intent(in)    :: c_snTilt
  integer,      intent(inout) :: iProc
  integer,      intent(inout) :: n_chan
  integer,      intent(inout) :: n_VR
  integer,      intent(inout) :: n_amorphous

  call cry_doCrystal(x,xp,z,zp,s,p,x0,xp0,zlm,s_imp,isImp,nhit,nabs,lhit,part_abs,impact,indiv,&
                    c_length,exenergy,rho,anuc,zatom,emr,dlri,dlyi,ai,eUm,collnt,&
                    hcut,bnref,csref0,csref1,csref4,csref5,csect,cry_tilt,c_rcurv,&
                    c_alayer,c_xmax,c_ymax,c_orient,c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,&
                    c_spTilt,c_cnTilt,c_snTilt,iProc,n_chan,n_VR,n_amorphous)

end subroutine


subroutine pyk2_cryinteract(ci_x,xp,y,yp,pc,length,s_P,x_P,ci_exenergy,ci_rho,ci_anuc,ci_zatom,ci_emr,&
                      ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt,ci_hcut,ci_csref0,ci_csref1,ci_csref4,&
                      ci_csref5,ci_bnref,ci_csect,cry_tilt,c_rcurv,c_alayer,c_xmax,c_ymax,c_orient,&
                      c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,c_spTilt,c_cnTilt,c_snTilt,iProc)
    use coll_crystal, only: cry_interact

    implicit none
  
    ! integer,          intent(in)    :: is  ! Material number
    real(kind=8), intent(inout) :: ci_x
    real(kind=8), intent(inout) :: xp
    real(kind=8), intent(inout) :: y
    real(kind=8), intent(inout) :: yp
    real(kind=8), intent(inout) :: pc
    real(kind=8), intent(in)    :: length
    real(kind=8), intent(in)    :: s_P
    real(kind=8), intent(in)    :: x_P
  
    real(kind=8), intent(in)    :: ci_exenergy
    real(kind=8), intent(in)    :: ci_rho
    real(kind=8), intent(in)    :: ci_anuc
    real(kind=8), intent(in)    :: ci_zatom
    real(kind=8), intent(in)    :: ci_emr
    
    real(kind=8), intent(in)    :: ci_dlri
    real(kind=8), intent(in)    :: ci_dlyi
    real(kind=8), intent(in)    :: ci_ai
    real(kind=8), intent(in)    :: ci_eUm
    real(kind=8), intent(in)    :: ci_collnt
    real(kind=8), intent(in)    :: ci_hcut
    real(kind=8), intent(in)    :: ci_csref0
    real(kind=8), intent(in)    :: ci_csref1
    real(kind=8), intent(in)    :: ci_csref4
    real(kind=8), intent(in)    :: ci_csref5
  
    real(kind=8), intent(in)    :: ci_bnref
    real(kind=8), intent(in)    :: ci_csect

    real(kind=8), intent(in)    :: cry_tilt
    real(kind=8), intent(in)    :: c_rcurv
    real(kind=8), intent(in)    :: c_alayer
    real(kind=8), intent(in)    :: c_xmax
    real(kind=8), intent(in)    :: c_ymax
    integer,      intent(in)    :: c_orient
    real(kind=8), intent(in)    :: c_miscut
    real(kind=8), intent(in)    :: cry_bend
    real(kind=8), intent(in)    :: c_cBend
    real(kind=8), intent(in)    :: c_sBend
    real(kind=8), intent(in)    :: c_cpTilt
    real(kind=8), intent(in)    :: c_spTilt
    real(kind=8), intent(in)    :: c_cnTilt
    real(kind=8), intent(in)    :: c_snTilt
    integer,      intent(inout) :: iProc


    call cry_interact(ci_x,xp,y,yp,pc,length,s_P,x_P,ci_exenergy,ci_rho,ci_anuc,ci_zatom,ci_emr,&
                      ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt,ci_hcut,ci_csref0,ci_csref1,ci_csref4,&
                      ci_csref5,ci_bnref,ci_csect,cry_tilt,c_rcurv,c_alayer,c_xmax,c_ymax,c_orient,&
                      c_miscut,cry_bend,c_cBend,c_sBend,c_cpTilt,c_spTilt,c_cnTilt,c_snTilt, iProc)

end subroutine pyk2_cryinteract