

! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  implicit none

  real(kind=8), save :: zatom_curr ! Current zatom, used for Rutherford scattering integration
  real(kind=8), save :: emr_curr ! Current emr, used for Rutherford scattering integration
  real(kind=8), save :: cnorm_curr ! Current emr, used for Rutherford scattering integration
  !integer,      save :: rnd_seed

contains

!>
!! k2coll_ruth(t)
!! Calculate the rutherford scattering cross section
!<
! HACK coll_zatom is used as global; it cannot be passed as a function argument
! because this function is passed into funlxp which builds the random distribution
! However, the latter expects a function with one argument
real(kind=8) function k2coll_ruth(t)

  real(kind=8), intent(in) :: t
  real(kind=8), parameter :: cnform = 0.8561e3

  k2coll_ruth = (cnorm_curr*exp(((-1*t)*cnform)*emr_curr**2)) * (zatom_curr/t)**2

end function k2coll_ruth

end module coll_k2
