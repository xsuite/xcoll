

! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use floatPrecision
  use numerical_constants

  implicit none

  real(kind=fPrec), save :: zatom_curr ! Current zatom, used for Rutherford scattering integration
  real(kind=fPrec), save :: emr_curr ! Current emr, used for Rutherford scattering integration

contains

!>
!! k2coll_ruth(t)
!! Calculate the rutherford scattering cross section
!<
! HACK coll_zatom is used as global; it cannot be passed as a function argument
! because this function is passed into funlxp which builds the random distribution
! However, the latter expects a function with one argument
real(kind=fPrec) function k2coll_ruth(t)

  use mathlib_bouncer
  !use coll_materials

  real(kind=fPrec), intent(in) :: t
  !real(kind=fPrec), intent(in) :: ru_emr

  ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug
  real(kind=fPrec), parameter :: cnorm  = 2.607e-5_fPrec
  real(kind=fPrec), parameter :: cnform = 0.8561e3_fPrec

  k2coll_ruth = (cnorm*exp_mb(((-one*t)*cnform)*emr_curr**2)) * (zatom_curr/t)**2

end function k2coll_ruth

end module coll_k2
