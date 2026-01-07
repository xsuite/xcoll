subroutine pyfluka_init(n_alloc, debug_level, cwd_path)
    use mod_fluka
    !, only : fluka_enable, fluka_mod_init
    use physical_constants, only : clight
    use mod_units, only : units_path

    implicit none
    integer, intent(in)    :: n_alloc
    integer, intent(in)    :: debug_level
    character(len=255), intent(in) :: cwd_path

    ! NB: In SixTrack, npart was passed, not n_alloc.
    ! (Needed for e.g. avoiding to re-compile?)
    units_path = cwd_path
    fluka_debug_level  = debug_level

    call fluka_mod_init(n_alloc, 500, clight)
    fluka_enable = .true.
end subroutine


subroutine pyfluka_connect(timeout_sec)
    use crcoall
    use mod_fluka
    !, only : fluka_connect, fluka_connected

    implicit none

    integer fluka_con
    integer(kind=4), intent(in) :: timeout_sec

    ! start connection to FLUKA and initialise max ID
    if(fluka_enable) then
       fluka_con = fluka_is_running()
       if(fluka_con == -1) then
          write(lout,"(a)") "FLUKA> ERROR Fluka is expected to run but it is NOT actually the case"
          flush(lout)
          call fluka_close
       end if
       if(fluka_debug_level > 0) then
           write(unit_pyfluka,"(a)") "FLUKA> Initializing FlukaIO interface"
           flush(unit_pyfluka)
       end if
       fluka_con = fluka_connect(timeout_sec)
       if(fluka_con == -1) then
          write(lout,"(a)") "FLUKA> ERROR Cannot connect to Fluka server"
          flush(lout)
          call fluka_close
       endif
       if(fluka_debug_level > 0) then
          write(unit_pyfluka,"(a)") "FLUKA> Successfully connected to Fluka server"
          flush(unit_pyfluka)
       end if
       fluka_connected = .true.
    endif

end subroutine


subroutine pyfluka_close()
    use mod_fluka
    !, only : fluka_close

    implicit none

    call fluka_close

end subroutine


subroutine pyfluka_init_max_uid(npart)
    use crcoall
    use mod_fluka
    !, only : fluka_init_max_uid, fluka_enable

    implicit none
    integer, intent(in)    :: npart
    integer fluka_con

    ! P.Garcia Ortega, A.Mereghetti and V.Vlachoudis, for the FLUKA Team
    ! last modified: 26-08-2014
    ! send npart to fluka
    if(fluka_enable) then
       if(fluka_debug_level > 0) then
           write(unit_pyfluka,"(a,i0)") "FLUKA> Sending npart = ",npart
           flush(unit_pyfluka)
       end if
       ! IMPORTANT: The call to fluka_init_max_uid is absolutely needed.
       ! The FLUKA server looks (in order!) for the corresponding message.
       fluka_con = fluka_init_max_uid( npart )

       if(fluka_con < 0) then
          write(lout,"(a,i0,a,i0,a)") "FLUKA> ERROR ", fluka_con, ": Failed to send npart ",npart," to fluka "
          flush(lout)
          call fluka_close
       end if

       if(fluka_debug_level > 0) then
           write(unit_pyfluka,"(a)") "FLUKA> Sending npart successful"
           flush(unit_pyfluka)
       end if
       flush(lout)
    end if

end subroutine


subroutine pyfluka_set_synch_part(part_e0, part_pc0, part_mass0, part_a0, part_z0, part_q0)
    use crcoall
    use mod_common
    use mod_fluka

    implicit none
    integer fluka_con

    real(kind=8),    intent(in) :: part_e0, part_pc0, part_mass0
    integer(kind=4), intent(in) :: part_a0, part_z0, part_q0

    ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
    ! last modified: 18-01-2016
    ! initialise energy/momentum/rest mass of reference particle in mod_fluka
    !     and synch magnetic rigidity with Fluka (for the time being, consider
    !     only protons);
    if(fluka_enable) then
       if(fluka_debug_level > 0) then
           write(unit_pyfluka,"(a)") "FLUKA> Updating the reference particle"
           flush(unit_pyfluka)
       end if

       ! Optional: Let's also set the reference particle in mod_common, like was done in SixTrack.
       ! Default values for e0 and e0f are 0!
       ! Will only be needed when/if e0f etc become used in mod_fluka.f90.
       e0 = part_e0
       e0f = part_pc0
       nucm0 = part_mass0
       aa0 = part_a0
       zz0 = part_z0
       qq0 = part_q0

       ! IMPORTANT: The call to fluka_set_synch_part is absolutely needed.
       ! The FLUKA server looks (in order!) for the corresponding message.
       fluka_con = fluka_set_synch_part(part_e0, part_pc0, part_mass0, part_a0, part_z0, part_q0)

       if(fluka_con < 0) then
          write(lout,"(a,i0,a)") "FLUKA> ERROR ", fluka_con, ": Failed to update the reference particle"
          flush(lout)
          call fluka_close
       end if

       if(fluka_debug_level > 0) then
           write(unit_pyfluka,"(a)") "FLUKA> Updating the reference particle successful"
           flush(unit_pyfluka)
       end if
    end if

end subroutine


subroutine track_fluka(turn, fluka_id, length, alive_part, max_part, x_part, xp_part, y_part, yp_part, &
                       zeta_part, e_part, m_part, q_part, A_part, Z_part, pdg_id_part, part_id, parent_id, &
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

    integer(kind=4), intent(in)    :: turn
    integer(kind=2), intent(in)    :: fluka_id
    real(kind=8),    intent(in)    :: length
    integer,         intent(inout) :: alive_part           ! napx
    integer,         intent(in)    :: max_part             ! npart
    real(kind=8),    intent(inout) :: x_part(max_part)     ! [mm]    xv1
    real(kind=8),    intent(inout) :: xp_part(max_part)    ! [1e-3]  yv1
    real(kind=8),    intent(inout) :: y_part(max_part)     ! [mm]    xv2
    real(kind=8),    intent(inout) :: yp_part(max_part)    ! [1e-3]  yv2
    real(kind=8),    intent(inout) :: zeta_part(max_part)  ! [mm]    sigmv
    real(kind=8),    intent(inout) :: e_part(max_part)     ! [MeV]   ejv   (ejfv is momentum, dpsv is delta, oidpsv is 1/(1+d))
    real(kind=8),    intent(inout) :: m_part(max_part)     ! [MeV]   nucm
    integer(kind=2), intent(inout) :: q_part(max_part)     !         nqq     Charge
    integer(kind=4), intent(inout) :: A_part(max_part)     !         naa     Ion atomic mass
    integer(kind=4), intent(inout) :: Z_part(max_part)     !         nzz     Ion atomic number
    integer(kind=4), intent(inout) :: pdg_id_part(max_part) !         pdgid   Particle PDGid
    integer(kind=4), intent(inout) :: part_id(max_part)
    integer(kind=4), intent(inout) :: parent_id(max_part)
    real(kind=8),    intent(inout) :: part_weight(max_part)
    real(kind=8),    intent(inout) :: spin_x_part(max_part)  ! spin_x  ! x component of the particle spin
    real(kind=8),    intent(inout) :: spin_y_part(max_part)  ! spin_y  ! y component of the particle spin
    real(kind=8),    intent(inout) :: spin_z_part(max_part)  ! spin_z  ! z component of the particle spin

    integer ret

    npart = max_part
    napx = alive_part

    ret = fluka_send_receive(turn, fluka_id, length, alive_part, max_part, x_part, y_part, xp_part, yp_part, &
                           zeta_part, e_part, A_part, Z_part, m_part, q_part, pdg_id_part, &
                           part_id, parent_id, part_weight, spin_x_part, spin_y_part, spin_z_part )
    napx = alive_part

    if (ret.lt.0) then
        write(lout,*) 'FLUKA> ERROR ', ret, ' in Fluka communication returned by fluka_send_receive...'
        write(lout,*) 'ENDED WITH ERROR.'
        flush(lout)
    end if

    return
end subroutine track_fluka

