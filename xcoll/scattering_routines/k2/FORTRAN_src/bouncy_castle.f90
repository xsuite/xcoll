

! =================================================================================================
!  Mathlib Bouncer
!  Last modified: 2018-10-23
! =================================================================================================
module mathlib_bouncer

  use floatPrecision

  implicit none
  private

  public :: sin_mb, asin_mb, sinh_mb, cos_mb, acos_mb, cosh_mb
  public :: tan_mb, atan_mb, atan2_mb, exp_mb, log_mb, log10_mb, isnan_mb

  ! For linking with 1


  private :: acos_rn, asin_rn, atan2_rn

  ! Interface definitions for the other functions in 1
  interface

    real(kind=c_double) pure function exp_rn(arg) bind(C,name="exp_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function exp_rn

    real(kind=c_double) pure function log_rn(arg) bind(C,name="log_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log_rn

    real(kind=c_double) pure function log10_rn(arg) bind(C,name="log10_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function log10_rn

    real(kind=c_double) pure function atan_rn(arg) bind(C,name="atan_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function atan_rn

    real(kind=c_double) pure function tan_rn(arg) bind(C,name="tan_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function tan_rn

    real(kind=c_double) pure function sin_rn(arg) bind(C,name="sin_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sin_rn

    real(kind=c_double) pure function cos_rn(arg) bind(C,name="cos_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cos_rn

    real(kind=c_double) pure function sinh_rn(arg) bind(C,name="sinh_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function sinh_rn

    real(kind=c_double) pure function cosh_rn(arg) bind(C,name="cosh_rn")
      use, intrinsic :: iso_c_binding, only : c_double
      real(kind=c_double), intent(in), value :: arg
    end function cosh_rn




  end interface

contains

  ! ========================================================================== !
  !  Definition of the MathlibBouncer (_mb) callable functions
  ! ========================================================================== !

  logical pure elemental function isnan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    isnan_mb = arg /= arg
  end function isnan_mb

  real(kind=fPrec) pure function sin_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    sin_mb=sin_rn(arg)
  end function sin_mb

  real(kind=fPrec) pure function asin_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    asin_mb=asin_rn(arg)
  end function asin_mb

  real(kind=fPrec) pure function sinh_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    sinh_mb=sinh_rn(arg)
  end function sinh_mb

  real(kind=fPrec) pure function cos_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    cos_mb=cos_rn(arg)
  end function cos_mb

  real(kind=fPrec) pure function acos_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    acos_mb=acos_rn(arg)
  end function acos_mb

  real(kind=fPrec) pure function cosh_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    cosh_mb=cosh_rn(arg)
  end function cosh_mb

  real(kind=fPrec) pure function tan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    tan_mb=tan_rn(arg)
  end function tan_mb

  real(kind=fPrec) pure function atan_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    atan_mb=atan_rn(arg)
  end function atan_mb

  real(kind=fPrec) pure function atan2_mb(y,x)
    real(kind=fPrec), intent(in) :: y,x
    atan2_mb=atan2_rn(y,x)
  end function atan2_mb

  real(kind=fPrec) pure function exp_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    exp_mb=exp_rn(arg)
  end function exp_mb

  real(kind=fPrec) pure function log_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    log_mb=log_rn(arg)
  end function log_mb

  real(kind=fPrec) pure function log10_mb(arg)
    real(kind=fPrec), intent(in) :: arg
    log10_mb=log10_rn(arg)
  end function log10_mb


  ! ========================================================================== !
  !  ROUND NEAR SPEACIAL FUNCTIONS
  ! ========================================================================== !

  real(kind=real64) pure function acos_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      acos_rn=x
    elseif(abs(x) == 0.0d0) then
      acos_rn=pi2
    else
      ! Try using (1-x)*(1+x) in case x is very small or close to 1. Write a test program?
      acos_rn=atan_rn(sqrt((1.0d0-x)*(1.0d0+x))/x)
      if(x < 0d0) then
        acos_rn=pi+acos_rn
      end if
    end if
  end function acos_rn

  real(kind=real64) pure function asin_rn(x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x
    if(x /= x) then ! Check if NaN
      asin_rn=x
      return
    end if
    if(abs(x) == 1.0d0) then
      asin_rn=sign(pi2,x)
    else
      ! Try using (1-x)*(1+x) in case x is very small or close to 1. Write a test program?
      asin_rn=atan_rn(x/sqrt((1.0d0-x)*(1.0d0+x)))
    end if
  end function asin_rn

  real(kind=real64) pure function atan2_rn(y,x)
    use, intrinsic :: iso_fortran_env, only : real64
    use numerical_constants, only : pi, pi2
    real(kind=real64), intent(in) :: x,y
    if(x == 0d0) then
      if(y == 0d0) then
        ! Let the internal atan2 handle this case according to ieee
        atan2_rn=atan2(y,x)
      else
        atan2_rn=sign(pi2,y)
      end if
    else
      if(y == 0d0) then
        if(x > 0d0) then
          ! Let the internal atan2 handle this case according to ieee
          atan2_rn=atan2(y,x)
        else
          atan2_rn=pi
        end if
      else
        atan2_rn=atan_rn(y/x)
        if(x < 0d0) then
          atan2_rn=sign(pi,y)+atan2_rn
        end if
      end if
    end if
  end function atan2_rn

  ! ========================================================================== !
  !  ROUND UP SPEACIAL FUNCTIONS
  ! ========================================================================== !


  ! ========================================================================== !
  !  ROUND DOWN SPEACIAL FUNCTIONS
  ! ========================================================================== !


  ! ========================================================================== !
  !  ROUND ZERO SPEACIAL FUNCTIONS
  ! ========================================================================== !



end module mathlib_bouncer
