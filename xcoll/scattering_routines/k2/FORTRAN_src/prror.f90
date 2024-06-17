subroutine prror

  use crcoall
#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  ! These should not go to lerr
  write(lout,"(a)") ""
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") "    +      ERROR DETECTED!      +"
  write(lout,"(a)") "    + RUN TERMINATED ABNORMALLY +"
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") ""

#ifdef FLUKA
  call fluka_close
#endif

#ifdef CR
  call abend("ERROR")
#else
  call closeUnits
  stop 1
#endif

end subroutine prror


subroutine closeUnits

  use crcoall
  use mod_units,   only : units_maxUnit, f_close
  use, intrinsic :: iso_fortran_env

  implicit none

  integer i
  logical isOpen

  ! Then iterate through 1 to units_maxUnit
  do i=1, units_maxUnit
    ! Do not close the following units:
    if(i == output_unit .or. i == input_unit .or. i == error_unit .or. i == lerr .or. i == lout .or. i == crlog) cycle
    inquire(unit=i, opened=isOpen)
    if(isOpen) call f_close(i)
  end do

end subroutine closeUnits
