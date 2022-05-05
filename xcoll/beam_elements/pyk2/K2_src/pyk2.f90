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


subroutine pyk2_run( &
                    x_in, &
                    xp_in, &
                    y_in, &
                    yp_in, &
                    s_in, &
                    p_in, &
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
                    run_cgen, &
                    is_crystal, &
                    c_length, &
                    c_aperture, &
                    c_offset, &
                    c_tilt, &
                    onesided, &
                    length, &
                    p0, &
                    nhit, &
                    nabs, &
                    fracab, &
                    mirror, &
                    cRot, &
                    sRot, &
                    cRRot, &
                    sRRot, &
                    nnuc0, &
                    nnuc1, &
                    ien0, &
                    ien1, &
                    isImp, &
                    s, &
                    keeps, &
                    zlm, &
                    x, &
                    xp, &
                    xp_in0, &
                    z, &
                    zp, &
                    p, &
                    sp, &
                    dpop, &
                    x_flk, &
                    y_flk, &
                    xp_flk, &
                    yp_flk, &
                    x_in0, &
                    xIn, &
                    xpIn, &
                    yIn, &
                    ypIn, &
                    tiltangle)

  ! use parpro ,           only : npart
  ! use coll_common ,      only : rcx, rcxp, rcy, rcyp, rcp, rcs, coll_expandArrays
  ! //use coll_materials ! for collmat_init
  use coll_k2        ! for scattering

  implicit none


  ! ############################
  ! ## variables declarations ##
  ! ############################

  ! integer, intent(in)          :: npart
  real(kind=8), intent(inout)  :: x_in
  real(kind=8), intent(inout)  :: xp_in
  real(kind=8), intent(inout)  :: y_in
  real(kind=8), intent(inout)  :: yp_in
  real(kind=8), intent(inout)  :: s_in
  real(kind=8), intent(inout)  :: p_in

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
  real(kind=8)     , intent(in) :: run_cgen(200)

  logical(kind=4)  , intent(in) :: is_crystal
  real(kind=8) ,    intent(in) :: c_length
  real(kind=8) ,    intent(in) :: c_aperture
  real(kind=8) ,    intent(in) :: c_offset
  real(kind=8) , intent(inout) :: c_tilt(2)
  logical(kind=4) ,  intent(in):: onesided
  real(kind=8),  intent(inout) :: length
  real(kind=8),  intent(inout) :: p0

  integer,          intent(inout) :: nhit
  integer,          intent(inout) :: nabs
  integer(kind=8),          intent(inout) :: nnuc0
  integer(kind=8),          intent(inout) :: nnuc1
  real(kind=8), intent(inout) :: ien0
  real(kind=8), intent(inout) :: ien1
  real(kind=8), intent(inout) :: fracab
  real(kind=8), intent(inout) :: mirror
  real(kind=8), intent(inout) :: cRot
  real(kind=8), intent(inout) :: sRot
  real(kind=8), intent(inout) :: cRRot
  real(kind=8), intent(inout) :: sRRot

  logical(kind=4), intent(inout) :: isImp
  real(kind=8),    intent(inout) :: s
  real(kind=8),    intent(in) :: keeps
  real(kind=8),    intent(inout) :: zlm
  real(kind=8),    intent(inout) :: sp
  real(kind=8),    intent(inout) :: x_flk
  real(kind=8),    intent(inout) :: y_flk
  real(kind=8),    intent(inout) :: xp_flk
  real(kind=8),    intent(inout) :: yp_flk

  real(kind=8),    intent(inout) :: x
  real(kind=8),    intent(inout) :: xp
  real(kind=8),    intent(inout) :: xp_in0
  real(kind=8),    intent(inout) :: z
  real(kind=8),    intent(inout) :: zp
  real(kind=8),    intent(inout) :: p
  real(kind=8),    intent(inout) :: dpop  
  real(kind=8),    intent(inout) :: x_in0
  real(kind=8),    intent(inout) :: xIn
  real(kind=8),    intent(inout) :: xpIn
  real(kind=8),    intent(inout) :: yIn
  real(kind=8),    intent(inout) :: ypIn
  real(kind=8),    intent(in)    :: tiltangle


  ! needs to be passed from cry_startElement
  integer cry_proc, cry_proc_prev, cry_proc_tmp
  cry_proc = -1
  cry_proc_prev = -1
  cry_proc_tmp = -1

  call k2coll_collimate( &
     run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, &
     run_csref0, run_csref1, run_csref4, run_csref5, run_radl, run_dlri, &
     run_dlyi, run_eUm, run_ai, run_collnt, run_cprob, run_xintl, run_bn, &
     run_ecmsq, run_xln15s, run_bpp, run_cgen, is_crystal, &
     c_length, c_aperture, c_offset, c_tilt, &
     x_in, xp_in, y_in, yp_in, p_in, s_in, &
     val_part_hit, val_part_abs, &
     val_part_impact, val_part_indiv, val_part_linteract, &
     onesided, 1, val_nabs_type, val_linside, length, p0, nhit, &
     nabs, fracab, mirror, cRot, sRot, cRRot, sRRot, nnuc0, &
    nnuc1, ien0, ien1, cry_proc, cry_proc_prev, cry_proc_tmp, &
    isImp, s, keeps, zlm, sp, x_flk, y_flk, xp_flk, yp_flk, &
    x, xp, xp_in0, z, zp, p, dpop, x_in0, xIn, xpIn, yIn, ypIn, tiltangle)

end subroutine 

