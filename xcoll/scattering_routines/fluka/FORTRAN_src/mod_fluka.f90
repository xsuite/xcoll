module mod_fluka

  use floatPrecision
  use numerical_constants
  use mod_alloc

  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 18-01-2016
  ! fortran 90 module for coupling SixTrack to FLUKA
  ! NOTA BENE:
  !    napx  (SixTrack) -> npart     (mod_fluka)
  !    npart (SixTrack) -> max_npart (mod_fluka)

  implicit none
  private

  public :: fluka_mod_init
  public :: fluka_mod_expand_arrays
  public :: fluka_mod_end

  public :: fluka_connect
  public :: fluka_end

  public :: fluka_send_receive
  public :: fluka_send
  public :: fluka_receive
  public :: fluka_shuffleLostParticles
  public :: fluka_set_synch_part
  public :: fluka_init_max_uid
  public :: fluka_is_running

  public :: fluka_close

  ! HION variables that are only used for FLUKA
  ! ien0,ien1: ion energy entering/leaving the collimator
  real(kind=8),    public :: ien0, ien1
  integer(kind=4), public :: nnuc0,nnuc1

  ! FlukaIO Connection parameters
  character(len = 255), public  :: fluka_host
  integer, public :: fluka_port
  character(len = 255), parameter :: fluka_net_nfo_file = 'network.nfo'

  ! FlukaIO interface
  external ntinit, ntconnect, ntend
  external ntsendp,     &
           ntsendeob,   &
           ntsendeoc,   &
           ntsendipt,   &
           ntrecv,      &
           ntwait,      &
           ntsendnpart, &
           ntsendbrhono
  external ntrtimeout, ntwtimeout

  integer(kind=4) :: ntconnect,   &
                         ntsendp,     &
                         ntsendeob,   &
                         ntsendeoc,   &
                         ntsendipt,   &
                         ntrecv,      &
                         ntwait,      &
                         ntsendnpart, &
                         ntsendbrhono,&
                         ntend
  integer(kind=4) :: ntrtimeout, ntwtimeout

  ! FlukaIO Message types
  integer(kind=1), parameter :: FLUKA_PART = 1, &
                                FLUKA_EOB  = 2, &
                                FLUKA_EOC  = 3, &
                                FLUKA_CONF = 4, &
                                FLUKA_IPT  = 5, &
                                FLUKA_HSK  = 6, &
                                FLUKA_NPART= 7, &
                                FLUKA_BRHO = 8
  ! connection ID
  integer(kind=4) :: fluka_cid

  ! FLUKA input block
  logical, public :: fluka_enable    = .false.         ! enable coupling
  logical, public :: fluka_connected = .false.         ! fluka is connected
  integer, public :: fluka_debug_level = 1             ! write debug messages

  ! fluka insertions
  logical, public :: fluka_inside = .false.                    ! Are we in a fluka insertion?
  integer(kind=4), public, allocatable :: fluka_type(:)        ! type of insertion (one per SINGLE ELEMENT)
  integer(kind=4), public, allocatable :: fluka_geo_index(:)   ! index of insertion (one per SINGLE ELEMENT)
  real(kind=8), public, allocatable :: fluka_synch_length(:)   ! length of insertion [m] (one per SINGLE ELEMENT)
  ! current fluka insertion:
  integer :: fluka_i=-1              ! lattice entry (FLUKA_ELEMENT/FLUKA_EXIT)
  integer :: fluka_ix=-1             ! single element entry (FLUKA_ELEMENT/FLUKA_EXIT)
  integer :: fluka_nturn=-1          ! turn
  integer :: fluka_last_sent_mess=-1 ! last sent message
  integer :: fluka_last_rcvd_mess=-1 ! last received message
  ! recognised insertion types
  integer(kind=4), parameter, public :: FLUKA_NONE    = 0, & ! no insertion
                                            FLUKA_ELEMENT = 1, & ! insertion covers only the present SINGLE ELEMENT
                                            FLUKA_ENTRY   = 2, & ! SINGLE ELEMENT marking the start of the insertion
                                            FLUKA_EXIT    = 3    ! SINGLE ELEMENT marking the end   of the insertion
  ! ancillary tracking values
  integer(kind=4), public :: fluka_max_npart                          ! Maximum number of particles (array size)
  integer,          public, allocatable    :: pids(:)         ! Particle ID moved from hisixtrack, to be harmonised

  ! Useful values
  integer :: fluka_nsent     ! Temporary count of sent particles
  integer :: fluka_nrecv     ! Temporary count of received particles
  real(kind=8), public :: fluka_clight ! [m/s]

  ! Reference particle
  real(kind=8), public :: fluka_e0     ! [GeV]
  real(kind=8), public :: fluka_pc0    ! [GeV/c]
  real(kind=8), public :: fluka_mass0  ! [GeV/c2]
  real(kind=8), public :: fluka_brho0  ! [Tm]
  integer(kind=4),          public :: fluka_chrg0  ! []
  integer(kind=4),          public :: fluka_a0     ! nucelon number (hisix)
  integer(kind=4),          public :: fluka_z0     ! nuclear charge

  save

contains

  !----------------------------------------------------------------------------
  ! set the module up
  subroutine fluka_mod_init(npart, nele, clight)

    use mod_units
    use mod_common, only : fort208, unit208
    implicit none

    ! interface variables
    integer :: npart, nele
    real(kind=8) :: clight

    ! temporary variables
    integer :: j

    if(fluka_debug_level > 1) then
       write(lout,*) "fluka_mod_init npart=", npart, ", nele=", nele, ", clight=", clight
       flush(lout)
    end if

    fluka_max_npart = npart
    fluka_clight    = clight

    call alloc(pids,               npart, 0, "pids")
    call alloc(fluka_type,         nele, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele, zero, 'fluka_synch_length')

    if(unit208 == -1) then
      call f_requestUnit(fort208,unit208)
      call f_open(unit=unit208,file=fort208,formatted=.true.,mode="w")
    end if

  end subroutine fluka_mod_init

  subroutine fluka_mod_expand_arrays(npart_new, nele_new)

    use parpro, only : npart

    implicit none

    integer :: npart_new, nele_new, j

    call alloc(pids,               npart_new, 0, "pids")
    call alloc(fluka_type,         nele_new, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele_new, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele_new, zero, 'fluka_synch_length')

    fluka_max_npart = npart_new

  end subroutine fluka_mod_expand_arrays

  !----------------------------------------------------------------------------
  ! un-set the module
  subroutine fluka_mod_end()
    implicit none
    call dealloc(pids,"pids")
    call dealloc(fluka_type,'fluka_type')
    call dealloc(fluka_geo_index,'fluka_geo_index')
    call dealloc(fluka_synch_length,'fluka_synch_length')
  end subroutine fluka_mod_end

  !----------------------------------------------------------------------------
  ! acquire info for network communication
  subroutine fluka_read_config(net_nfo_file, host, port)

    use mod_units

    implicit none

    ! interface variables
    character(len=255) :: net_nfo_file
    character(len=255) :: host
    integer :: port
    integer :: net_nfo_unit
    integer :: ios

    call f_requestUnit(net_nfo_file, net_nfo_unit)
    call f_open(net_nfo_unit, file=net_nfo_file, formatted=.true., mode="rw", status='old')
    read(unit=net_nfo_unit, fmt=*, iostat=ios) host
    if(ios .ne. 0) then
      write(lout,'(A)') 'FLUKA> ERROR Could not read the host name from network.nfo'
      ERROR STOP 'ENDED WITH ERROR.'
    end if

    read(unit=net_nfo_unit, fmt=*, iostat=ios) port
    if(ios .ne. 0) then
      write(lout,'(A)') 'FLUKA> ERROR Could not read the port number from network.nfo'
      write(lout,'(A)') 'FLUKA>       Is the FLUKA server running and has it had time to write the port number?'
      ERROR STOP 'ENDED WITH ERROR.'
    end if

    if(fluka_debug_level > 1) then
       write(lout,*) "FLUKA>   host=", host, " port=", port
       flush(lout)
    end if

    call f_close(net_nfo_unit)

  end subroutine fluka_read_config

  !----------------------------------------------------------------------------
  ! start communication with fluka
  integer function fluka_connect(timeout_sec)
    implicit none

    integer(kind=4) :: timeout_sec
    integer :: rtimeout
    integer :: wtimeout

    call fluka_read_config(fluka_net_nfo_file, fluka_host, fluka_port)
    call ntinit()
    fluka_cid = ntconnect(fluka_host, fluka_port)
    fluka_connect = fluka_cid

    rtimeout = ntrtimeout(fluka_cid, timeout_sec)
    wtimeout = ntwtimeout(fluka_cid, timeout_sec)

  end function fluka_connect

  !----------------------------------------------------------------------------
  ! close communication with fluka
  subroutine fluka_end()
    implicit none

    ! Finish connection
    integer(kind=4) :: n

    ! Fluka I/O parameters
    integer(kind=4)         :: flid, flgen
    real(kind=8)  :: flwgt, flx, fly, flz, flxp, flyp, flzp, flpc, flm, flt
    integer(kind=4)         :: flaa, flzz
    integer(kind=2)         :: flq
    integer(kind=1)          :: mtype

    integer(kind=4)         :: flpdgid
    real(kind=8)            :: flsx, flsy, flsz

    if(fluka_debug_level > 0) then
        write(lout,'(A)') 'FlukaIO> Sending End of Computation signal'
        flush(lout)
    end if

    ! Send end of computation
    n = ntsendeoc(fluka_cid)
    if(n.lt.0) then
      write(lout,'(A,i0,A)') "FlukaIO> ERROR ", n, " - Error sending End of Computation"
      flush(lout)
      return
    end if

    ! Wait end of comp
    n = ntwait(fluka_cid, mtype, &
          flid, flgen, flwgt, flx, fly, flz, flxp, flyp, flzp, &
          flm, flpc, flt, flpdgid, flq, flsx, flsy, flsz)
    if(n.lt.0) then
      write(lout,'(A,i0,A)') "FlukaIO> ERROR ", n, " - Server timed out while waiting End of Computation"
      flush(lout)
      return
    end if
    if(mtype.ne.FLUKA_EOC) then
      write(lout,*) "FlukaIO> WARNING Received unexpected message at shutdown"
      flush(lout)
    end if

    ! At this point both ends agreed to disconnect

    ! Close connection
    n = ntend(fluka_cid)

    return
  end subroutine fluka_end

  !----------------------------------------------------------------------------
  ! send and receive particles from Fluka
  integer function fluka_send_receive(turn, ipt, el, npart, max_part, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                   partID, parentID, partWeight, spinx, spiny, spinz)
    implicit none

    ! Parameters
    integer(kind=4) :: turn
    integer(kind=2) :: ipt
    integer           :: npart
    integer           :: max_part
    real(kind=8)  :: el

    real(kind=8) :: xv1(max_part)
    real(kind=8) :: yv1(max_part)
    real(kind=8) :: xv2(max_part)
    real(kind=8) :: yv2(max_part)
    real(kind=8) :: s(max_part)
    real(kind=8) :: etot(max_part)

    real(kind=8) :: mass(max_part)
    integer(kind=4) :: aa(max_part)
    integer(kind=4) :: zz(max_part)
    integer(kind=2) :: qq(max_part)
    integer(kind=4) :: pdg_id(max_part)
    integer(kind=4) :: partID(max_part)
    integer(kind=4) :: parentID(max_part)
    real(kind=8) :: partWeight(max_part)
    real(kind=8) :: spinx(max_part)
    real(kind=8) :: spiny(max_part)
    real(kind=8) :: spinz(max_part)

    fluka_send_receive = fluka_send(turn, ipt, el, npart, max_part, &
         xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
         partID, parentID, partWeight, spinx, spiny, spinz)
    if(fluka_debug_level > 1) then
        write(lout,*) "    Sent all particles (return code ", fluka_send_receive, ")"
        flush(lout)
    end if
    if(fluka_send_receive.lt.0) return

    fluka_send_receive = fluka_receive(turn, ipt, el, npart, max_part, &
         xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
         partID, parentID, partWeight, spinx, spiny, spinz)
    if(fluka_debug_level > 1) then
        write(lout,*) "    Received all particles (return code ", fluka_send_receive, ")"
        flush(lout)
    end if
  end function fluka_send_receive

  !----------------------------------------------------------------------------
  ! just send particles to Fluka
  integer function fluka_send(turn, ipt, el, npart, max_part, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                   partID, parentID, partWeight, spinx, spiny, spinz)
    implicit none

    ! Interface variables
    integer(kind=4) :: turn
    integer(kind=2) :: ipt
    integer           :: npart
    integer           :: max_part
    real(kind=8)  :: el

    real(kind=8) :: xv1(max_part)
    real(kind=8) :: yv1(max_part)
    real(kind=8) :: xv2(max_part)
    real(kind=8) :: yv2(max_part)
    real(kind=8) :: s(max_part)
    real(kind=8) :: etot(max_part)

    real(kind=8) :: mass(max_part)
    integer(kind=4) :: aa(max_part)
    integer(kind=4) :: zz(max_part)
    integer(kind=2) :: qq(max_part)
    integer(kind=4) :: pdg_id(max_part)
    integer(kind=4) :: partID(max_part)
    integer(kind=4) :: parentID(max_part)
    real(kind=8) :: partWeight(max_part)
    real(kind=8) :: spinx(max_part)
    real(kind=8) :: spiny(max_part)
    real(kind=8) :: spinz(max_part)

    ! Fluka I/O parameters
    integer(kind=4) :: flid, flgen
    real(kind=8)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=4) :: flaa, flzz
    integer(kind=2) :: flq
    integer(kind=1)  :: mtype

    integer(kind=4) :: flpdgid
    real(kind=8)    :: flsx, flsy, flsz

    ! Auxiliary variables
    integer :: j
    integer(kind=4) :: n

    fluka_send = 0
    fluka_last_rcvd_mess = -1

    n = ntsendipt(fluka_cid, turn, ipt)
    if(n.lt.0) then
        write(lout,*) "FlukaIO> ERROR (return code ", n, ") - Error sending Insertion Point"
        fluka_cid = -1
        fluka_send = n
        return
    end if
    fluka_last_sent_mess=FLUKA_IPT
    if(fluka_debug_level > 1) then
        write(lout,*) "FLUKA> Sent insertion point (return code ", n, ")"
        flush(lout)
    end if

    fluka_nsent = 0
    fluka_nrecv = 0
    mtype = 0

    do j=1, npart

        flid  = partID(j)
        flgen = parentID(j)
        flwgt = partWeight(j)

        flx   = xv1(j) * c1m1  ! from [mm] to [cm]
        fly   = xv2(j) * c1m1  ! from [mm] to [cm]
        flz   = zero

        flxp  = yv1(j) * c1m3 ! from [1.0E-03] to [1.0]
        flyp  = yv2(j) * c1m3 ! from [1.0E-03] to [1.0]
        ! director cosines:
        ! full transformation:
        flzp = sqrt(one - flxp**2 - flyp**2)
        ! flzp  = sqrt( one / ( flxp**2 + flyp**2 + one ) )
        !      ! taylor expansion, for consistency with drifts in SixTrack:
        !      flzp  = 1d0 / ( 1d0 + ( flxp**2+flyp**2 )/2d0 )
        !flxp  = flxp * flzp
        !flyp  = flyp * flzp

        ! total energy:
        flet  = etot(j) * c1m3 ! from [MeV] to [GeV]
        ! longitudinal phase:
        flt   = -s(j) * c1m3 / ( (fluka_pc0/fluka_e0)*fluka_clight ) ! from [mm] to [s]


        ! Ion properties (PH for hiSix)
        flm   = mass(j) * c1m3      ! unit is [GeV]
        flaa  = aa(j)
        flzz  = zz(j)

        flpdgid = pdg_id(j)
        flq  = qq(j)
        flsx = spinx(j)
        flsy = spiny(j)
        flsz = spinz(j)

        if(fluka_debug_level > 2) then
            write(lout,*) "    Sending particle: ", flid, flgen, flwgt, flx, fly, flz, flxp, flyp, flzp, flm, flet, flt
            write(lout,*) flsx, flsy, flsz, flaa, flzz, flq, flpdgid
            flush(lout)
        end if

        ! Send particle         TODO: it seems flet should be pc!!
        n = ntsendp(fluka_cid, &
            flid, flgen, flwgt, &
            flx, fly, flz, &
            flxp, flyp, flzp, &
            flm, flet, flt, &
            flpdgid, flq, flsx, flsy, flsz)
        !      write(lout,*) "Sent particle: ", n, fluka_nsent, FLUKA_PART

        if(n.lt.0) then
            write(lout,*) "FlukaIO> ERROR (return code ", n, ") - Error sending Particle"
            fluka_cid = -1
            fluka_send = -1
            return
        end if

        fluka_nsent = fluka_nsent + 1
        fluka_last_sent_mess=FLUKA_PART

    end do

    ! Send end of batch
    n = ntsendeob(fluka_cid)

    if(n.lt.0) then
        write(lout,*) "FlukaIO> ERROR (return code ", n, ") - Error sending End of Batch"
        fluka_cid = -1
        fluka_send = -1
        return
    end if
    fluka_last_sent_mess=FLUKA_EOB

  end function fluka_send

  !----------------------------------------------------------------------------
  ! just receive particles from Fluka
  ! The call from fluka.s90 is:
  ! fluka_receive( nturn, fluka_geo_index(ix), eltot, napx, xv1(:), yv1(:), xv2(:), yv2(:), sigmv, ejv, naa(:), nzz(:), nucm(:))
  ! When the above arrays are made allocatable, the below variables will need updating - see mod_commonmn and mod_hions
  integer function fluka_receive(turn, ipt, el, napx, max_part, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                                 partID, parentID, partWeight, spinx, spiny, spinz)

    use parpro
    use mod_pdgid
    use mod_common_main, only : MaximumPartID

    implicit none

    ! Interface variables
    integer(kind=4) :: turn
    integer(kind=2) :: ipt
    integer           :: napx
    integer           :: max_part
    real(kind=8)  :: el

    real(kind=8) :: xv1(max_part)
    real(kind=8) :: yv1(max_part)
    real(kind=8) :: xv2(max_part)
    real(kind=8) :: yv2(max_part)
    real(kind=8) :: s(max_part)
    real(kind=8) :: etot(max_part)

    real(kind=8) :: mass(max_part)
    integer(kind=4) :: aa(max_part)
    integer(kind=4) :: zz(max_part)
    integer(kind=2) :: qq(max_part)
    integer(kind=4) :: pdg_id(max_part)
    integer(kind=4) :: partID(max_part)
    integer(kind=4) :: parentID(max_part)
    real(kind=8) :: partWeight(max_part)
    real(kind=8) :: spinx(max_part)
    real(kind=8) :: spiny(max_part)
    real(kind=8) :: spinz(max_part)

    ! Fluka I/O parameters
    integer(kind=4) :: flid, flgen
    real(kind=8)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=4) :: flaa, flzz
    integer(kind=2) :: flq
    integer(kind=1)  :: mtype

    integer(kind=4)         :: flpdgid
    real(kind=8)            :: flsx, flsy, flsz

    ! Auxiliary variables
    integer(kind=4) :: n, j

    fluka_receive = 0
    fluka_last_sent_mess = -1

    fluka_nrecv = 0
    mtype = 0

    ! assign default values
    do j = 1, npart
        partID(j) = j
        parentID(j) = j

        partWeight(j) = one

        xv1 (j) = zero
        xv2 (j) = zero
        yv1 (j) = zero
        yv2 (j) = zero
        etot(j) = zero
        s   (j) = zero
        ! hisix: we should also parse m0,A0,Z0
        aa  (j) = 1
        zz  (j) = 1
        mass(j) = zero
        qq  (j) = 1
        pdg_id(j) = 0
        spinx = zero
        spiny = zero
        spinz = zero
    end do

    ! Wait until end of turn (Synchronize)
    do while(mtype.ne.FLUKA_EOB)
        n = ntwait(fluka_cid, mtype, &
              flid, flgen, flwgt, &
              flx, fly, flz, &
              flxp, flyp, flzp, &
              flm, flet, flt, &
              flpdgid, flq, flsx, flsy, flsz)
        if(n.lt.0) then
            write(lout,*) "FlukaIO> ERROR (return code ", n ,") - Server timed out while waiting for message"
            fluka_cid = -1
            fluka_receive = n
            return
        end if

        if(mtype.eq.FLUKA_PART) then

            fluka_nrecv = fluka_nrecv + 1
            fluka_last_rcvd_mess = FLUKA_PART

            if(fluka_nrecv .gt. npart) then
                write(lout,*) 'FlukaIO> ERROR - HIT END OF PARTICLE ARRAY'
                return
            end if

            call CalculateAZ(flpdgid, flaa, flzz)

            if(fluka_debug_level > 2) then
                write(lout,*) "    Received particle: ", flid, flgen, flwgt, flx, fly, flz, flxp, flyp, flzp, flm, flet, flt
                write(lout,*) flsx, flsy, flsz, flaa, flzz, flq, flpdgid
                flush(lout)
            end if

            partID(fluka_nrecv)      = flid
            parentID(fluka_nrecv)    = flgen
            if (partID(fluka_nrecv).gt.MaximumPartID) then
                MaximumPartID = partID(fluka_nrecv)
            end if

            partWeight(fluka_nrecv)  = flwgt
            xv1(fluka_nrecv)         = flx * c1e1   ! from [cm]  to [mm]
            xv2(fluka_nrecv)         = fly * c1e1   ! from [cm]  to [mm]
            yv1(fluka_nrecv)         = flxp * c1e3 ! from director cosine to x' [1.0E-03]
            yv2(fluka_nrecv)         = flyp * c1e3 ! from director cosine to x' [1.0E-03]
            etot(fluka_nrecv)         = flet * c1e3  ! from [GeV] to [MeV]
            s(fluka_nrecv)            = ( el - (fluka_pc0/fluka_e0)*(flt*fluka_clight) ) * c1e3 ! from [s] to [mm]
            aa(fluka_nrecv)           = flaa          !PH for hiSix
            zz(fluka_nrecv)           = flzz          !PH for hiSix
            mass(fluka_nrecv)         = flm  * c1e3  ! from [GeV] to [MeV]         !PH for hiSix
            qq(fluka_nrecv)           = flq
            pdg_id(fluka_nrecv)       = flpdgid
            spinx(fluka_nrecv)        = flsx
            spiny(fluka_nrecv)        = flsy
            spinz(fluka_nrecv)        = flsz

    !       The conversion is now done inside the coupling server
    !       call GetPDGid_fromFLUKA(-2, pdg_id(fluka_nrecv), flaa, flzz)
        end if

      !Finished waiting end of turn
    end do

    napx = fluka_nrecv
    fluka_last_rcvd_mess = FLUKA_EOB

    if(fluka_debug_level > 1) then
        write(lout,*) "FlukaIO> turn = ", turn, " ipt = ", ipt, " sent = ", fluka_nsent, &
                      " received = ", fluka_nrecv, " max_uid = ", MaximumPartID
        flush(lout)
    end if

  end function fluka_receive

  !----------------------------------------------------------------------------
  ! compact ancillary tracking arrays
  subroutine fluka_shuffleLostParticles(tnapx, j)

    integer, intent(in) :: tnapx
    integer, intent(in) :: j

    if(fluka_debug_level > 1) then
      write(lout, *) 'FLUKA> fluka_shuffleLostParticles called with napx (lnapx for SixTrack) = ', tnapx, ', j = ', j
      flush(lout)
    end if

  end subroutine fluka_shuffleLostParticles

  !----------------------------------------------------------------------------
  ! set reference particle properties (mainly for longitudinal dynamics)
  integer function fluka_set_synch_part( e0, pc0, mass0, a0, z0, q0 )
    implicit none

    ! interface variables
    real(kind=8) :: e0, pc0, mass0
    integer(kind=4) :: a0, z0, q0

    ! Auxiliary variables
    integer(kind=4) :: n

    fluka_set_synch_part = 0

    fluka_e0    = e0    *c1m3 ! from  [MeV]    to [GeV]
    fluka_pc0   = pc0   *c1m3 ! from  [MeV/c]  to [GeV/c]
    fluka_mass0 = mass0 *c1m3 ! from  [MeV/c2] to [GeV/c2]
    fluka_a0 = a0
    fluka_z0 = z0
    fluka_chrg0 = q0

    ! update magnetic rigidity, unless division by clight and 10^-9
    fluka_brho0 = fluka_pc0 / real(fluka_chrg0,real64)

    if(fluka_debug_level > 1) then
        write(lout,*) '    - total en        [GeV]:',fluka_e0
        write(lout,*) '    - momentum      [GeV/c]:',fluka_pc0
        write(lout,*) '    - mass         [GeV/c2]:',fluka_mass0
        write(lout,*) '    - A mass number      []:',fluka_a0
        write(lout,*) '    - Z number           []:',fluka_z0
        write(lout,*) '    - charge state      [e]:',fluka_chrg0
        write(lout,*) '    - rigidity     [Tm/0.3]:', fluka_brho0
        write(lout,*) '      in proper units  [Tm]:', fluka_brho0 / ( fluka_clight*c1m9 )
        flush(lout)
    end if

    ! inform Fluka about the new magnetic rigidity
    ! This step cannot be skipped!!!
    ! The fluka server is waiting for this message, so that BRHONO can be set.
    n = ntsendbrhono(fluka_cid, fluka_brho0)
    if (n .lt. 0) then
      fluka_set_synch_part = -1
      return
    end if

  end function fluka_set_synch_part

  !----------------------------------------------------------------------------
  ! set max ID
  integer function fluka_init_max_uid( npart )

    use mod_common_main, only : MaximumPartID

    implicit none

    ! interface variables
    integer(kind=4) :: npart

    ! Auxiliary variables
    integer(kind=4) :: n

    fluka_init_max_uid = 0

    MaximumPartID = npart

    n = ntsendnpart(fluka_cid, npart)
    if (n .lt. 0) then
      fluka_init_max_uid = -1
      return
    end if

  end function fluka_init_max_uid


  !----------------------------------------------------------------------------
  ! check if fluka is running, ie if it created the
  integer function fluka_is_running()
    implicit none

    ! temporary variables
    logical :: lexist

    fluka_is_running = 0
    inquire( file=fluka_net_nfo_file, exist=lexist)

    if (.not.lexist) then
       write(lout,*) 'FLUKA> ERROR File containing network infos ', fluka_net_nfo_file, ' does not exist!!'
       flush(lout)
       fluka_is_running = -1
    endif

  end function fluka_is_running

!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     clean closure of communication with fluka and un-set mod_fluka
!     inserted in main code by the 'fluka' compilation flag
subroutine fluka_close

     use crcoall

     implicit none

     integer fluka_con

     if(fluka_debug_level > 0) then
         write(lout,'(A)') 'FLUKA> Call to fluka_close'
         flush(lout)
     end if
     if (fluka_enable) then
       if (fluka_last_sent_mess==FLUKA_EOB .and. fluka_last_rcvd_mess.eq.-1) then
         ! FLUKA is still crunching stuff
         if(fluka_debug_level > 0) then
             write(lout,'(A)') 'FLUKA> fluka_close called while particles are still on the Fluka side'
             write(lout,'(A)') 'FLUKA>    dummy wait to have a clean closure'
             flush(lout)
         end if
       end if
       fluka_con = fluka_is_running()
       if(fluka_con.eq.0) then
         if( .not. fluka_connected ) then
!              temporarily connect to fluka, to properly terminate the run
           fluka_con = fluka_connect(3600)
           if(fluka_con.lt.0) then
!                no hope to properly close the run
             write(lout,'(A)') 'FLUKA> ERROR Unable to connect to fluka while closing the simulation:'
             write(lout,'(A)') 'FLUKA>       please, manually kill all its instances'
             flush(lout)
             goto 1982
           endif
           if(fluka_debug_level > 0) then
               write(lout,'(A)') 'FLUKA> Successfully connected to Fluka server (only temporarily)'
               flush(lout)
           end if
         end if
         call fluka_end
       end if
     end if
1982 call fluka_mod_end
     flush(lout)
end subroutine fluka_close

end module mod_fluka

