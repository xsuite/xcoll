subroutine pyfluka_init(n_alloc)

    ! Initialize FLUKA
    ! Same initialization as was done in SixTrack main_cr.f90
    !use crcoall
    !use parpro
    !use mod_settings
    !use mod_common
    !use mod_common_main
    !use mod_commons
    !use mod_common_track

    use mod_fluka
    !, only : fluka_enable, fluka_mod_init
    use physical_constants, only : clight

    implicit none
    integer, intent(in)    :: n_alloc

    call fluka_mod_init(n_alloc, 500, clight)
    fluka_enable = .true.

end subroutine


subroutine pyfluka_connect()
    use crcoall
    use mod_fluka
    !, only : fluka_connect, fluka_connected

    implicit none

    integer fluka_con

    ! start connection to FLUKA and initialise max ID
    if(fluka_enable) then
       fluka_con = fluka_is_running()
       if(fluka_con == -1) then
          write(lerr,"(a)") "FLUKA> ERROR Fluka is expected to run but it is NOT actually the case"
          write(fluka_log_unit,*) "# Fluka is expected to run but it is NOT actually the case"
          call prror
       end if
       write(lout,"(a)") "FLUKA> Initializing FlukaIO interface ..."
       write(fluka_log_unit,*) "# Initializing FlukaIO interface ..."
       fluka_con = fluka_connect()
       if(fluka_con == -1) then
          write(lerr,"(a)") "FLUKA> ERROR Cannot connect to Fluka server"
          write(fluka_log_unit,*) "# Error connecting to Fluka server"
          call prror
       endif
       write(lout,"(a)") "FLUKA> Successfully connected to Fluka server"
       write(fluka_log_unit,*) "# Successfully connected to Fluka server"
       fluka_connected = .true.
    endif

end subroutine


subroutine pyfluka_close()
    use mod_fluka
    !, only : fluka_close

    implicit none

    call fluka_close

end subroutine
      

subroutine pyfluka_set_n_alloc(n_alloc)
    use crcoall
    use mod_fluka
    !, only : fluka_init_max_uid, fluka_enable

    implicit none
    integer, intent(in)    :: n_alloc
    integer fluka_con

    ! P.Garcia Ortega, A.Mereghetti and V.Vlachoudis, for the FLUKA Team
    ! last modified: 26-08-2014
    ! send n_alloc to fluka
    if(fluka_enable) then
       write(lout,"(a,i0)") "FLUKA> Sending n_alloc = ",n_alloc
       write(fluka_log_unit,*) "# Sending n_alloc: ", n_alloc
       fluka_con = fluka_init_max_uid( n_alloc )

       if(fluka_con < 0) then
          write(lerr,"(a,i0,a,i0,a)") "FLUKA> ERROR ", fluka_con, ": Failed to send n_alloc ",n_alloc," to fluka "
          write(fluka_log_unit, *) "# failed to send n_alloc to fluka ",n_alloc
          call prror
       end if

       write(lout,"(a)") "FLUKA> Sending n_alloc successful"
       write(fluka_log_unit,*) "# Sending n_alloc successful;"
       flush(lout)
       flush(fluka_log_unit)
    end if

end subroutine


subroutine pyfluka_set_synch_part(part_e0, part_pc0, part_mass0, part_a0, part_z0, part_q0)
    use crcoall
    use mod_common
    !use floatPrecision
    use mod_fluka

    !use, intrinsic :: ISO_FORTRAN_ENV, only : int32

    implicit none
    integer fluka_con

    real(kind=8),        intent(in) :: part_e0, part_pc0, part_mass0
    integer(kind=int32), intent(in) :: part_a0, part_z0, part_q0

    ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
    ! last modified: 18-01-2016
    ! initialise energy/momentum/rest mass of reference particle in mod_fluka
    !     and synch magnetic rigidity with Fluka (for the time being, consider
    !     only protons);
    if(fluka_enable) then
       write(lout,"(a)") "FLUKA> Updating the reference particle"
       write(fluka_log_unit,*) "# Updating ref particle"
       flush(lout)
       flush(fluka_log_unit)

       !e0 = part_e0
       !e0f = part_pc0
       !nucm0 = part_mass0
       !aa0 = part_a0
       !zz0 = part_z0
       !qq0 = part_q0

       !write(lout,*) 'e0 = ', e0
       !write(lout,*) 'e0f = ', e0f
       !write(lout,*) 'nucm0 = ', nucm0
       !write(lout,*) 'aa0 = ', aa0
       !write(lout,*) 'zz0 = ', zz0
       !write(lout,*) 'qq0 = ', qq0
       !flush(lout)

       !fluka_con = fluka_set_synch_part( e0, e0f, nucm0, aa0, zz0, qq0)
       fluka_con = fluka_set_synch_part(part_e0, part_pc0, part_mass0, part_a0, part_z0, part_q0)

       if(fluka_con < 0) then
          write(lerr,"(a,i0,a)") "FLUKA> ERROR ", fluka_con, ": Failed to update the reference particle"
          write(fluka_log_unit,*) "# failed to update ref particle"
          call prror
       end if

       write(lout,"(a)") "FLUKA> Updating the reference particle successful"
       write(fluka_log_unit,*) "# Updating ref particle successful;"
       flush(lout)
       flush(fluka_log_unit)
    end if

end subroutine


subroutine track_fluka(turn, fluka_id, length, part_p0c, part_e0, alive_part, max_part, x_part, xp_part, y_part, yp_part, &
                       zeta_part, e_part, m_part, q_part, A_part, Z_part, pdgid_part, part_id, parent_id, &
                       part_weight, spin_x_part, spin_y_part, spin_z_part)

    use floatPrecision
    use numerical_constants, only : zero, one, c1e3, c1m3
    use crcoall
    use parpro
    use mod_common
    use mod_common_track
    use mod_common_main
    use mod_fluka

    implicit none

    integer(kind=int32), intent(in)    :: turn
    integer(kind=int32), intent(in)    :: fluka_id
    real(kind=8),        intent(in)    :: length
    real(kind=8),        intent(in)    :: part_p0c
    real(kind=8),        intent(in)    :: part_e0
    integer,             intent(in)    :: alive_part           ! napx
    integer,             intent(in)    :: max_part             ! npart
    real(kind=8),        intent(inout) :: x_part(max_part)     ! [mm]    xv1
    real(kind=8),        intent(inout) :: xp_part(max_part)    ! [1e-3]  yv1
    real(kind=8),        intent(inout) :: y_part(max_part)     ! [mm]    xv2
    real(kind=8),        intent(inout) :: yp_part(max_part)    ! [1e-3]  yv2
    real(kind=8),        intent(inout) :: zeta_part(max_part)  ! [mm]    sigmv
    real(kind=8),        intent(inout) :: e_part(max_part)     ! [MeV]   ejv   (ejfv is momentum, dpsv is delta, oidpsv is 1/(1+d))
    real(kind=8),        intent(inout) :: m_part(max_part)     ! [MeV]   nucm
    integer(kind=int32), intent(inout) :: q_part(max_part)     !         nqq     Charge
    integer(kind=int32), intent(inout) :: A_part(max_part)     !         naa     Ion atomic mass
    integer(kind=int32), intent(inout) :: Z_part(max_part)     !         nzz     Ion atomic number
    integer(kind=int32), intent(inout) :: pdgid_part(max_part) !         pdgid   Particle PDGid
    integer(kind=int32), intent(inout) :: part_id(max_part)
    integer(kind=int32), intent(inout) :: parent_id(max_part)
    real(kind=8),        intent(inout) :: part_weight(max_part)
    real(kind=8),        intent(inout) :: spin_x_part(max_part)  ! spin_x  ! x component of the particle spin
    real(kind=8),        intent(inout) :: spin_y_part(max_part)  ! spin_y  ! y component of the particle spin
    real(kind=8),        intent(inout) :: spin_z_part(max_part)  ! spin_z  ! z component of the particle spin

    integer ret

    npart = max_part
    napx = alive_part
    fluka_pc0 = part_p0c
    fluka_e0 = part_e0

    ret = fluka_send_receive(turn, fluka_id, length, alive_part, max_part, x_part, y_part, xp_part, yp_part, &
                           zeta_part, e_part, A_part, Z_part, m_part, q_part, pdgid_part, &
                           part_id, parent_id, part_weight, spin_x_part, spin_y_part, spin_z_part )
    napx = alive_part

    if (ret.lt.0) then
        write(fluka_log_unit,*) 'FLUKA> ERROR ', ret, ' in Fluka communication returned by fluka_send_receive...'
        write(fluka_log_unit,*) 'ENDED WITH ERROR.'
    end if

    return
end subroutine track_fluka


subroutine prror
  use crcoall
  use mod_fluka

  implicit none

  ! These should not go to lerr
  write(lout,"(a)") ""
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") "    +      ERROR DETECTED!      +"
  write(lout,"(a)") "    + RUN TERMINATED ABNORMALLY +"
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") ""

  call fluka_close

end subroutine prror
