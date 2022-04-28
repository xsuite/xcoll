subroutine pyk2_init(n_alloc, random_generator_seed)
  use floatPrecision
  use numerical_constants
  ! use crcoall    NODIG ??
  use parpro ,           only : npart
  use mod_alloc ,        only : alloc      ! to allocate partID etc
  use mod_common ,       only : iexact, napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, &
                                coll_expandArrays
  !use coll_materials ! for collmat_init
  use coll_k2         ! for scattering

  implicit none

  integer, intent(in)          :: n_alloc
  integer, intent(in)          :: random_generator_seed

  ! Set default values for collimator materials
 ! call collmat_init

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

  call coll_expandArrays(n_alloc)
  call alloc(naa, n_alloc, aa0, "naa")
  call alloc(partID, n_alloc, 0, "partID")
  call alloc(parentID, n_alloc, 0, "parentID")
  call alloc(pairID, 2, n_alloc, 0, "pairID")

end subroutine








subroutine pyk2_start_run(num_particles, cgen, hcut, zatom, emr)
  use floatPrecision
  use numerical_constants

  use parpro ,           only : npart
  use mod_common ,       only : napx, unit208
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_funlux,        only : funlxp
  
  implicit none

  integer, intent(in)             :: num_particles
  real(kind=fPrec), intent(inout) :: cgen(200)
  real(kind=fPrec), intent(in)    :: zatom   !
  real(kind=fPrec), intent(in)    :: emr     !
  real(kind=fPrec), intent(in)    :: hcut    ! 
  
  npart=num_particles
  unit208=109 

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do
  
  napx=npart  ! this decreases after absorptions!

  
  real(kind=fPrec), parameter :: tlcut = 0.0009982_fPrec
  ! Prepare for Rutherford differential distribution
  !mcurr = mat ! HACK> mcurr is global, and coll_zatom too which is used inside k2coll_ruth
  zatom_curr = coll_zatom
  emr_curr = coll_emr

  call funlxp(k2coll_ruth, cgen(1), tlcut, coll_hcut)
  
end subroutine


subroutine pyk2_doCrystal(j,x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit,&
  part_abs_local,impact,indiv,c_length,run_exenergy,run_anuc,run_zatom,run_emr,run_rho,&
  run_hcut,run_bnref,run_csref0,run_csref1,run_csref4,run_csref5,run_dlri,run_dlyi,&
  run_eUm,run_ai,run_runnt,run_bn)

  ! ... all variable definitions go here ...
  integer,          intent(in)    :: j
  real(kind=fPrec), intent(inout) :: x,xp
  real(kind=fPrec), intent(inout) :: z,zp
  real(kind=fPrec), intent(inout) :: s,p
  real(kind=fPrec), intent(inout) :: x0,xp0
  real(kind=fPrec), intent(inout) :: zlm,s_imp
  integer,          intent(inout) :: nhit,nabs
  integer,          intent(inout) :: lhit(npart)
  integer,          intent(inout) :: part_abs(npart)
  real(kind=fPrec), intent(inout) :: impact(npart)
  real(kind=fPrec), intent(inout) :: indiv(npart)
  real(kind=fPrec), intent(in)    :: c_length
  real(kind=fPrec), intent(inout) :: run_exenergy
  real(kind=fPrec), intent(in)    :: run_anuc
  real(kind=fPrec), intent(in)    :: run_zatom
  real(kind=fPrec), intent(in)    :: run_emr
  real(kind=fPrec), intent(in)    :: run_rho
  real(kind=fPrec), intent(in)    :: run_hcut
  real(kind=fPrec), intent(in)    :: run_bnref
  real(kind=fPrec), intent(in)    :: run_csref0
  real(kind=fPrec), intent(in)    :: run_csref1
  real(kind=fPrec), intent(in)    :: run_csref4
  real(kind=fPrec), intent(in)    :: run_csref5
  real(kind=fPrec), intent(in)    :: run_dlri
  real(kind=fPrec), intent(in)    :: run_dlyi
  real(kind=fPrec), intent(in)    :: run_eUm
  real(kind=fPrec), intent(in)    :: run_ai
  real(kind=fPrec), intent(in)    :: run_collnt
  real(kind=fPrec), intent(inout) :: run_bn

  logical,          intent(inout) :: isImp

  use module coll_crystal

  do j=1,npart
    rcx(j) = x_particles(j)
    rcxp(j) = xp_particles(j)
    rcy(j) = y_particles(j)
    rcyp(j) = yp_particles(j)
    rcs(j) = s_particles(j)
    rcp(j) = p_particles(j)
  end do

  call cry_doCrystal(j,rcx,rcxp,rcz,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit,&
  part_abs_local,impact,indiv,c_length,coll_exenergy,coll_anuc,coll_zatom,coll_emr,coll_rho,&
  coll_hcut,coll_bnref,coll_csref0,coll_csref1,coll_csref4,coll_csref5,coll_dlri,coll_dlyi,&
  coll_eUm,coll_ai,coll_collnt,coll_bn)

  do j=1,npart
    x_particles(j) = rcx(j)
    xp_particles(j) = rcxp(j)
    y_particles(j) = rcy(j)
    yp_particles(j) = rcyp(j)
    s_particles(j) = rcs(j)
    p_particles(j) = rcp(j)
 end do
end subroutine


subroutine pyk2_jaw(s, nabs, ipart, j_exenergy, j_anuc, j_zatom, j_rho, j_radl, &
  j_cprob, j_xintl, j_bn, j_cgen, j_ecmsq, j_xln15s, j_bpp)

end subroutine pyk2_jaw







subroutine pyk2_run(num_particles, &
                    x_particles, &
                    xp_particles, &
                    y_particles, &
                    yp_particles, &
                    s_particles, &
                    p_particles, &
                    part_hit, &
                    part_abs, &
                    part_impact, &
                    part_indiv, &
                    part_linteract, &
                    nhit_stage, &
                    nabs_type, &
                    linside, &
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
                    run_radl, &
                    run_dlri, & 
                    run_dlyi, &
                    run_eUm, &
                    run_ai, &
                    run_collnt, &
                    run_cprob, &
                    run_xintl, &
                    run_bn, &
                    run_ecmsq, &
                    run_xln15s, &
                    run_bpp, &
                    is_crystal, &
                    c_length, &
                    c_rotation, &
                    c_aperture, &
                    c_offset, &
                    c_tilt, &
                    c_enom, &
                    onesided, &
                    random_generator_seed)

  use floatPrecision
  use numerical_constants

  use parpro ,           only : npart
  use mod_alloc ,        only : alloc      ! to allocate partID etc
  use mod_common ,       only : iexact, napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, coll_expandArrays
  ! //use coll_materials ! for collmat_init
  use coll_k2        ! for scattering

  implicit none


  ! ############################
  ! ## variables declarations ##
  ! ############################

  integer, intent(in)          :: num_particles
  real(kind=8), intent(inout)  :: x_particles(num_particles)
  real(kind=8), intent(inout)  :: xp_particles(num_particles)
  real(kind=8), intent(inout)  :: y_particles(num_particles)
  real(kind=8), intent(inout)  :: yp_particles(num_particles)
  real(kind=8), intent(inout)  :: s_particles(num_particles)
  real(kind=8), intent(inout)  :: p_particles(num_particles)

  integer(kind=4)  , intent(inout) :: part_hit(num_particles)
  integer(kind=4)  , intent(inout) :: part_abs(num_particles)
  real(kind=8) , intent(inout) :: part_impact(num_particles)
  real(kind=8) , intent(inout) :: part_indiv(num_particles)
  real(kind=8) , intent(inout) :: part_linteract(num_particles)
  integer(kind=4)  , intent(inout) :: nhit_stage(num_particles)
  integer(kind=4)  , intent(inout) :: nabs_type(num_particles)
  logical(kind=4)  , intent(inout) :: linside(num_particles)

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

  real(kind=8)     , intent(in) :: run_radl
  real(kind=8)     , intent(in) :: run_dlri
  real(kind=8)     , intent(in) :: run_dlyi
  real(kind=8)     , intent(in) :: run_eUm
  real(kind=8)     , intent(in) :: run_ai
  real(kind=8)     , intent(in) :: run_collnt
  real(kind=8)     , intent(in) :: run_cprob(0:5)
  real(kind=8)     , intent(in) :: run_xintl
  real(kind=8)     , intent(inout) :: run_bn
  real(kind=8)     , intent(in) :: run_ecmsq
  real(kind=8)     , intent(in) :: run_xln15s
  real(kind=8)     , intent(in) :: run_bpp

  logical(kind=4)  , intent(in) :: is_crystal
  real(kind=8) ,    intent(in) :: c_length
  real(kind=8) ,    intent(in) :: c_rotation
  real(kind=8) ,    intent(in) :: c_aperture
  real(kind=8) ,    intent(in) :: c_offset
  real(kind=8) , intent(inout) :: c_tilt(2)
  real(kind=8) ,    intent(in) :: c_enom
  logical(kind=4) ,  intent(in):: onesided
  integer, intent(in)          :: random_generator_seed

  integer j



  ! ####################
  ! ## initialisation ##
  ! ####################
  !character(len=:),    allocatable   :: numpart
  !numpart="20000"
  !read(numpart,*) napx
  npart=num_particles

  if(random_generator_seed .ge. 0) then
        call rluxgo(3, random_generator_seed, 0, 0)
  end if

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do
  
  napx=npart  ! this decreases after absorptions!
  unit208=109 

  do j=1,npart
    rcx(j) = x_particles(j)
    rcxp(j) = xp_particles(j)
    rcy(j) = y_particles(j)
    rcyp(j) = yp_particles(j)
    rcs(j) = s_particles(j)
    rcp(j) = p_particles(j)
  end do

  call k2coll_collimate( &
     run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, &
     run_csref0, run_csref1, run_csref4, run_csref5, run_radl, run_dlri, &
     run_dlyi, run_eUm, run_ai, run_collnt, run_cprob, run_xintl, run_bn, &
     run_ecmsq, run_xln15s, run_bpp, is_crystal, &
     c_length, c_rotation, c_aperture, c_offset, c_tilt, &
     rcx, rcxp, rcy, rcyp, rcp, rcs, &
     c_enom, part_hit, part_abs, &
     part_impact, part_indiv, part_linteract, &
     onesided, nhit_stage, 1, nabs_type, linside)

  do j=1,npart
     x_particles(j) = rcx(j)
     xp_particles(j) = rcxp(j)
     y_particles(j) = rcy(j)
     yp_particles(j) = rcyp(j)
     s_particles(j) = rcs(j)
     p_particles(j) = rcp(j)
  end do
end subroutine 
