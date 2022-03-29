program main
  use floatPrecision
  use numerical_constants
  ! use crcoall    NODIG ??
  use parpro ,           only : npart
  use mod_alloc ,        only : alloc      ! to allocate partID etc
  use mod_common ,       only : iexact, napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, coll_expandArrays
  use coll_materials ! for collmat_init
  use coll_db        ! for cdb_readCollDB
  use coll_k2        ! for scattering

  use files  ! for testing

  implicit none


  ! ####################
  ! ## test variables ##
  ! ####################
  character(len=80)  :: filename
  integer j


  ! ########################
  ! ## function variables ##
  ! ########################
  integer(kind=4)  :: icoll
  integer(kind=4)  :: iturn
  integer(kind=4)  :: ie
  real(kind=fPrec) :: c_length
  real(kind=fPrec) :: c_rotation
  real(kind=fPrec) :: c_aperture
  real(kind=fPrec) :: c_offset
  real(kind=fPrec) :: c_tilt(2)
  !integer, save    :: rnd_seed   = 0        ! also defined in coll_common
  !real(kind=fPrec) :: rcx(20000)            ! also defined in coll_common
  !real(kind=fPrec) :: rcxp(20000)           ! also defined in coll_common
  !real(kind=fPrec) :: rcy(20000)            ! also defined in coll_common
  !real(kind=fPrec) :: rcyp(20000)           ! also defined in coll_common
  !real(kind=fPrec) :: rcp(20000)            ! also defined in coll_common
  !real(kind=fPrec) :: rcs(20000)            ! also defined in coll_common
  real(kind=fPrec) :: c_enom
  integer(kind=4)  :: part_hit_pos(20000)
  integer(kind=4)  :: part_hit_turn(20000)
  integer(kind=4)  :: part_abs_pos(20000)
  integer(kind=4)  :: part_abs_turn(20000)
  real(kind=fPrec) :: part_impact(20000)
  real(kind=fPrec) :: part_indiv(20000)
  real(kind=fPrec) :: part_linteract(20000)
  logical(kind=4)  :: onesided
  integer(kind=4)  :: nhit_stage(20000)
  integer(kind=4)  :: nabs_type(20000)
  logical(kind=4)  :: linside(20000)


  ! ####################
  ! ## initialisation ##
  ! ####################
  !character(len=:),    allocatable   :: numpart
  !numpart="20000"
  !read(numpart,*) napx
  npart=20000
  
  call alloc(naa, npart, aa0, "naa")
  call alloc(partID, npart, 0, "partID")
  call alloc(parentID, npart, 0, "parentID")
  call alloc(pairID, 2, npart, 0, "pairID")
  do j=1,npart
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do
  ! call alloc(part_abs_turn, npart, 0, "part_abs_turn")
  ! call alloc(rcx,  npart, zero, "rcx")
  ! call alloc(rcxp, npart, zero, "rcxp")
  ! call alloc(rcy,  npart, zero, "rcy")
  ! call alloc(rcyp, npart, zero, "rcyp")
  ! call alloc(rcp,  npart, zero, "rcp")
  ! call alloc(rcs,  npart, zero, "rcs")
  ! call alloc(cry_proc, npart, -1, "cry_proc")
  ! call alloc(cry_proc_prev, npart, -1, "cry_proc_prev")
  ! call alloc(cry_proc_tmp, npart, -1, "cry_proc_tmp")
  ! call alloc(rcx0,           npart,  zero, "rcx0")
  ! call alloc(rcxp0,          npart,  zero, "rcxp0")
  ! call alloc(rcy0,           npart,  zero, "rcy0")
  ! call alloc(rcyp0,          npart,  zero, "rcyp0")
  ! call alloc(rcp0,           npart,  zero, "rcp0")
  ! call alloc(part_hit_pos,   npart,  0,    "part_hit_pos")
  ! call alloc(part_hit_turn,  npart,  0,    "part_hit_turn")
  ! call alloc(part_abs_pos,   npart,  0,    "part_abs_pos")
  ! call alloc(part_select,    npart,  1,    "part_select")
  ! call alloc(nabs_type,      npart,  0,    "nabs_type")
  ! call alloc(nhit_stage,     npart,  0,    "nhit_stage")
  ! call alloc(part_impact,    npart,  zero, "part_impact")
  ! call alloc(part_indiv,     npart, -c1m6, "part_indiv")
  ! call alloc(part_linteract, npart,  zero, "part_linteract")




  
  napx=npart  ! this decreases after absorptions!
  unit208=109
  rnd_seed = 7569

  icoll = 31
  iturn = 1
  ie =1
  c_length = 0.59999999999999998
  c_rotation = 0
  c_aperture = 0.0025711021962573095
  c_offset = 0
  c_tilt = (0, 0)
  c_enom = 7000000
  onesided = .FALSE.
  linside(:) = .FALSE.


  call coll_expandArrays(npart)
  filename="rcx.dump"
  call realfromfile(filename,rcx)
  filename="rcxp.dump"
  call realfromfile(filename,rcxp)
  filename="rcy.dump"
  call realfromfile(filename,rcy)
  filename="rcyp.dump"
  call realfromfile(filename,rcyp)
  filename="rcp.dump"
  call realfromfile(filename,rcp)
  filename="rcs.dump"
  call realfromfile(filename,rcs)
  filename="part_hit_pos.dump"
  call intfromfile(filename,part_hit_pos)
  filename="part_hit_turn.dump"
  call intfromfile(filename,part_hit_turn)
  filename="part_abs_pos.dump"
  call intfromfile(filename,part_abs_pos)
  filename="part_abs_turn.dump"
  call intfromfile(filename,part_abs_turn)
  filename="part_impact.dump"
  call realfromfile(filename,part_impact)
  filename="part_indiv.dump"
  call realfromfile(filename,part_indiv)
  filename="part_linteract.dump"
  call realfromfile(filename,part_linteract)
  filename="nhit_stage.dump"
  call intfromfile(filename,nhit_stage)
  filename="nabs_type.dump"
  call intfromfile(filename,nabs_type)

! print *, c1m3
! print *, rcx
! print *, rcyp

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

  ! Set default values for collimator materials
  call collmat_init
  cdb_fileName="CollDB-RunIII.dat"
  call cdb_readCollDB


  call k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, &
     c_tilt, rcx, rcxp, rcy, rcyp, rcp, rcs, c_enom*c1m3, part_hit_pos, part_hit_turn, &
     part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract,             &
     onesided, nhit_stage, 1, nabs_type, linside)

  filename="rcx.dump_after"
  call realtofile(rcx,filename)
  filename="rcxp.dump_after"
  call realtofile(rcxp,filename)
  filename="rcy.dump_after"
  call realtofile(rcy,filename)
  filename="rcyp.dump_after"
  call realtofile(rcyp,filename)
  filename="rcp.dump_after"
  call realtofile(rcp,filename)
  filename="rcs.dump_after"
  call realtofile(rcs,filename)
  filename="part_hit_pos.dump_after"
  call inttofile(part_hit_pos,filename)
  filename="part_hit_turn.dump_after"
  call inttofile(part_hit_turn,filename)
  filename="part_abs_pos.dump_after"
  call inttofile(part_abs_pos,filename)
  filename="part_abs_turn.dump_after"
  call inttofile(part_abs_turn,filename)
  filename="part_impact.dump_after"
  call realtofile(part_impact,filename)
  filename="part_indiv.dump_after"
  call realtofile(part_indiv,filename)
  filename="part_linteract.dump_after"
  call realtofile(part_linteract,filename)
  filename="nhit_stage.dump_after"
  call inttofile(nhit_stage,filename)
  filename="nabs_type.dump_after"
  call inttofile(nabs_type,filename)

end program main

