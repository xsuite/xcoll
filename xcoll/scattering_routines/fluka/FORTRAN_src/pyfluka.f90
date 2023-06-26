subroutine pyfluka_init(n_alloc)
    use mod_fluka
    !, only : fluka_enable, fluka_mod_init
    use physical_constants, only : clight

    implicit none
    integer, intent(in)    :: n_alloc

    call fluka_mod_init(n_alloc, 500, clight)
    fluka_enable = .true.

end subroutine


subroutine pyfluka_connect()
    use mod_fluka
    !, only : fluka_connect, fluka_connected

    implicit none

    integer fluka_con

    fluka_con = fluka_connect()
    if(fluka_con == -1) then
      PRINT *, "ERROR Cannot connect to Fluka server"
    endif
    PRINT *, "Successfully connected to Fluka server"
    fluka_connected = .true.

end subroutine


subroutine pyfluka_close()
    use mod_fluka
    !, only : fluka_close

    implicit none

    call fluka_close

end subroutine


subroutine pyfluka_set_n_alloc(n_alloc)
    use mod_fluka
    !, only : fluka_init_max_uid, fluka_enable

    implicit none
    integer, intent(in)    :: n_alloc
    integer fluka_con

    if(fluka_enable) then
        PRINT *, "Changing n_alloc in FLUKA"
        fluka_con = fluka_init_max_uid( n_alloc )
        if(fluka_con < 0) then
            PRINT *, "Failed to change n_alloc in FLUKA"
        end if
        PRINT *, "Succesfully changed n_alloc in FLUKA"
    end if
end subroutine
