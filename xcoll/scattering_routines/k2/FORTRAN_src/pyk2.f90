subroutine pyk2_init(n_alloc, colldb_input_fname, random_generator_seed, num_coll, &
                     betax, betay, alphax, alphay, orbx, orby, orbxp, orbyp, gamma, emit)
  use floatPrecision
  use numerical_constants
  use parpro
  use mod_alloc
  use mod_common
  use mod_common_main
  use mod_units
  use mod_common_track
  use mod_ranlux
  use collimation
  use coll_common
  use coll_materials ! for collmat_init
  use coll_db        ! for cdb_readCollDB
  use coll_k2        ! for scattering
  use coll_crystal   ! for crystal scattering


  implicit none

  integer, intent(in)          :: n_alloc
  character(64), intent(in)   :: colldb_input_fname
  integer, intent(in)          :: random_generator_seed
  integer, intent(in)          :: num_coll
  real(kind=8), intent(in)     :: betax(num_coll)
  real(kind=8), intent(in)     :: betay(num_coll)
  real(kind=8), intent(in)     :: alphax(num_coll)
  real(kind=8), intent(in)     :: alphay(num_coll)
  real(kind=8), intent(in)     :: orbx(num_coll)
  real(kind=8), intent(in)     :: orby(num_coll)
  real(kind=8), intent(in)     :: orbxp(num_coll)
  real(kind=8), intent(in)     :: orbyp(num_coll)
  real(kind=8), intent(in)     :: gamma
  real(kind=8), intent(in)     :: emit   ! emittance in um
  integer j

  gammar = 1./gamma
  do j=1,num_coll
    ic(j) = j
    tbetax(j)  = betax(j)
    tbetay(j)  = betay(j)
    talphax(j) = alphax(j)
    talphay(j) = alphay(j)
    torbx(j)   = orbx(j)
    torby(j)   = orby(j)
    torbxp(j)  = orbxp(j)
    torbyp(j)  = orbyp(j)
  end do

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do

  emitnx0_dist = emit
  emitny0_dist = emit
  emitnx0_collgap = emit
  emitny0_collgap = emit
  call collmat_init
  cdb_fileName=colldb_input_fname
  call cdb_readCollDB

  ! Then do any implementation specific initial loading
  call k2coll_init
  if(coll_hasCrystal) then
    call cry_init
  end if
  ! Open the edep file
  if(unit208 == -1) then
    call f_requestUnit(fort208,unit208)
    call f_open(unit=unit208,file=fort208,formatted=.true.,mode="w",status="replace")
  end if

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

  ipencil = -1
  call coll_expandArrays(n_alloc)
  call alloc(naa, n_alloc, aa0, "naa")
  call alloc(partID, n_alloc, 0, "partID")
  call alloc(parentID, n_alloc, 0, "parentID")
  call alloc(pairID, 2, n_alloc, 0, "pairID")

end subroutine



subroutine pyk2_track(num_particles, x_particles, xp_particles, &
                      y_particles, yp_particles, p_particles, e_particles, &
                      ix, hit, absorbed)

  use floatPrecision
  use numerical_constants
  use parpro
  use mod_common
  use mod_common_main
  use coll_common
  use collimation
  use coll_k2

  implicit none

  integer, intent(in)          :: num_particles
  real(kind=8), intent(inout)  :: x_particles(num_particles)
  real(kind=8), intent(inout)  :: xp_particles(num_particles)
  real(kind=8), intent(inout)  :: y_particles(num_particles)
  real(kind=8), intent(inout)  :: yp_particles(num_particles)
  real(kind=8), intent(inout)  :: p_particles(num_particles)
  real(kind=8), intent(inout)  :: e_particles(num_particles)
  integer,      intent(in)     :: ix
  integer(kind=4), intent(inout) :: hit(num_particles)
  integer(kind=4), intent(inout) :: absorbed(num_particles)

  real(kind=8) stracki
  integer j

  npart = num_particles
  c_ix = ix

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do

  do j=1,npart
    xv1(j) = x_particles(j)
    yv1(j) = xp_particles(j)
    xv2(j) = y_particles(j)
    yv2(j) = yp_particles(j)
    ejv(j) = e_particles(j)
    ejfv(j) = p_particles(j)
  end do

  stracki = 0.
  call coll_doCollimator(stracki)

  do j=1,npart
    x_particles(j)  = xv1(j)
    xp_particles(j) = yv1(j)
    y_particles(j)  = xv2(j)
    yp_particles(j) = yv2(j)
    p_particles(j)  = ejv(j)
    hit(j)          = part_hit_pos(j)
    absorbed(j)     = part_abs_pos(j)
  end do

end subroutine

