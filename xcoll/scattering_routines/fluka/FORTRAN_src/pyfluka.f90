subroutine pyfluka_init()
  use mod_fluka, only : fluka_enable
  implicit none


  fluka_enable = .true.
  PRINT *, "lol"

end subroutine

