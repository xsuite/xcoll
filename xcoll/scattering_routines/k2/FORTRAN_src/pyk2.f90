subroutine pyk2_init(n_alloc, colldb_input_fname, random_generator_seed, num_coll, &
                     alphax, alphay, betax, betay, orbx, orby, orbxp, orbyp, &
                     nemitt_x, nemitt_y, e_ref, p_ref, m_ref, beta_ref, gamma_ref)
  use floatPrecision
  use numerical_constants
  use parpro
  use mod_alloc
  use mod_common
  use mod_commond2
  use mod_common_main
  use mod_common_track
  use mod_units
  use mod_ranlux
  use collimation
  use coll_common
  use coll_materials ! for collmat_init
  use coll_db        ! for cdb_readCollDB
  use coll_k2        ! for scattering
  use coll_crystal   ! for crystal scattering

  implicit none

  integer(kind=4), intent(in)  :: n_alloc
  character(64),   intent(in)  :: colldb_input_fname
  integer(kind=4), intent(in)  :: random_generator_seed
  integer(kind=4), intent(in)  :: num_coll
  real(kind=8),    intent(in)  :: alphax(num_coll)
  real(kind=8),    intent(in)  :: alphay(num_coll)
  real(kind=8),    intent(in)  :: betax(num_coll)
  real(kind=8),    intent(in)  :: betay(num_coll)
  real(kind=8),    intent(in)  :: orbx(num_coll)
  real(kind=8),    intent(in)  :: orby(num_coll)
  real(kind=8),    intent(in)  :: orbxp(num_coll)
  real(kind=8),    intent(in)  :: orbyp(num_coll)
  real(kind=8),    intent(in)  :: nemitt_x   ! emittance in m
  real(kind=8),    intent(in)  :: nemitt_y   ! emittance in m
  real(kind=8),    intent(in)  :: e_ref
  real(kind=8),    intent(in)  :: p_ref
  real(kind=8),    intent(in)  :: m_ref
  real(kind=8),    intent(in)  :: beta_ref
  real(kind=8),    intent(in)  :: gamma_ref
  integer(kind=4) j

  call f_initUnits
  units_beQuiet = .false. ! Allow mod_units to write to lout now
  call mod_common_expand_arrays(num_coll,1,num_coll,n_alloc,1)
  call mod_commont_expand_arrays(num_coll,n_alloc)
  call mod_commonmn_expand_arrays(num_coll,n_alloc)
  call mod_commond2_expand_arrays(num_coll)
  call collimation_expand_arrays(n_alloc,num_coll)
  call cdb_expand_arrays(num_coll)

  nucm0 = m_ref / c1e6
  e0 = e_ref / c1e6
  e0f = p_ref / c1e6
  beta0 = beta_ref
  gamma0 = gamma_ref
  gammar = 1./gamma0

  do j=1,num_coll
    ic(j) = j
    talphax(j) = alphax(j)
    talphay(j) = alphay(j)
    tbetax(j)  = betax(j)
    tbetay(j)  = betay(j)
    torbx(j)   = orbx(j)  * c1e3   ! closed orbit in mm in SixTrack
    torby(j)   = orby(j)  * c1e3
    torbxp(j)  = orbxp(j) * c1e3
    torbyp(j)  = orbyp(j) * c1e3
  end do

  !call alloc(naa, n_alloc, aa0, "naa")
  !call alloc(partID, n_alloc, 0, "partID")
  !call alloc(parentID, n_alloc, 0, "parentID")
  !call alloc(pairID, 2, n_alloc, 0, "pairID")
  !do j=1,n_alloc
  !  naa(j)      = 1
  !  partID(j)   = j
  !  parentID(j) = j
  !  pairID(1,j) = (j+1)/2    ! The pairID of particle j
  !  pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  !end do

  emitnx0_dist    = nemitt_x * c1e6  ! emittance in um in K2
  emitny0_dist    = nemitt_y * c1e6
  emitnx0_collgap = nemitt_x * c1e6
  emitny0_collgap = nemitt_y * c1e6
  call collmat_init
  cdb_fileName = colldb_input_fname
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

end subroutine



subroutine pyk2_track(num_particles, x_particles, xp_particles, y_particles, &
                      yp_particles, e_particles, p_particles, delta_particles, &
                      rvv_particles, rpp_particles, idx, hit, absorbed)

  use floatPrecision
  use numerical_constants
  use parpro
  use mod_common
  use mod_common_main
  use coll_common
  use collimation
  use coll_k2

  implicit none

  integer(kind=4), intent(in)    :: num_particles
  real(kind=8),    intent(inout) :: x_particles(num_particles)
  real(kind=8),    intent(inout) :: xp_particles(num_particles)
  real(kind=8),    intent(inout) :: y_particles(num_particles)
  real(kind=8),    intent(inout) :: yp_particles(num_particles)
  real(kind=8),    intent(inout) :: e_particles(num_particles)
  real(kind=8),    intent(in)    :: p_particles(num_particles)
  real(kind=8),    intent(in)    :: delta_particles(num_particles)
  real(kind=8),    intent(in)    :: rvv_particles(num_particles)
  real(kind=8),    intent(in)    :: rpp_particles(num_particles)
  integer(kind=4), intent(in)    :: idx
  integer(kind=4), intent(inout) :: hit(num_particles)
  integer(kind=4), intent(inout) :: absorbed(num_particles)

  real(kind=8) stracki
  integer(kind=4) j

  npart = num_particles
  napx  = npart
  c_ix  = idx

  !do j=1,npart
  !  naa(j) = aa0
  !  partID(j)   = j
  !  parentID(j) = j
  !  pairID(1,j) = (j+1)/2    ! The pairID of particle j
  !  pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  !end do

  do j=1,npart
    xv1(j)    = x_particles(j)  * c1e3
    yv1(j)    = xp_particles(j) * c1e3
    xv2(j)    = y_particles(j)  * c1e3
    yv2(j)    = yp_particles(j) * c1e3
    ejv(j)    = e_particles(j)  / c1e6
    ejfv(j)   = p_particles(j)  / c1e6
    dpsv(j)   = delta_particles(j)
    rvv(j)    = 1./rvv_particles(j)
    oidpsv(j) = rpp_particles(j)
    nucm(j)   = nucm0
    mtc(j)    = 1
  end do

  stracki = 0.
  call coll_doCollimator(stracki)

  do j=1,npart
    x_particles(j)  = xv1(j) / c1e3
    xp_particles(j) = yv1(j) / c1e3
    y_particles(j)  = xv2(j) / c1e3
    yp_particles(j) = yv2(j) / c1e3
    e_particles(j)  = ejv(j) * c1e6
    hit(j)          = part_hit_pos(j)
    absorbed(j)     = part_abs_pos(j)
  end do

end subroutine

