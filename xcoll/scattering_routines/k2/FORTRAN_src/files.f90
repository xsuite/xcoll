MODULE files
  USE floatPrecision

CONTAINS

SUBROUTINE realfromfile(filename,x)
  IMPLICIT NONE
  CHARACTER(len=80), INTENT(in)  :: filename  ! input
  REAL(kind=fPrec), INTENT(out)  :: x(20000)  ! output

  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="read")
  DO  i = 1, nrows
    READ(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE realfromfile

SUBROUTINE intfromfile(filename,x)
  IMPLICIT NONE
  CHARACTER(len=80), INTENT(in) :: filename   ! input
  INTEGER(kind=4), INTENT(out)  :: x(20000)   ! output
  

  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="read")
  DO  i = 1, nrows
    READ(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE intfromfile

SUBROUTINE boolfromfile(filename,x)
  IMPLICIT NONE
  CHARACTER(len=80), INTENT(in) :: filename   ! input
  LOGICAL(kind=4), INTENT(out)  :: x(20000)   ! output

  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="read")
  DO  i = 1, nrows
    READ(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE boolfromfile

SUBROUTINE realtofile(x,filename)
  IMPLICIT NONE
  REAL(kind=fPrec), INTENT(in)   :: x(20000)  ! input
  CHARACTER(len=80), INTENT(in)  :: filename  ! input
  
  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="write")
  DO  i = 1, nrows
    WRITE(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE realtofile

SUBROUTINE inttofile(x,filename)
  IMPLICIT NONE
  INTEGER(kind=4), INTENT(in)    :: x(20000)  ! input
  CHARACTER(len=80), INTENT(in)  :: filename  ! input
  
  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="write")
  DO  i = 1, nrows
    WRITE(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE inttofile

SUBROUTINE booltofile(x,filename)
  IMPLICIT NONE
  LOGICAL(kind=4), INTENT(in)    :: x(20000)  ! input
  CHARACTER(len=80), INTENT(in)  :: filename  ! input
  
  INTEGER :: nrows=20000
  INTEGER :: i

  OPEN(UNIT=11, FILE=filename, ACTION="write")
  DO  i = 1, nrows
    WRITE(11,*)  x(i)
  END DO
  CLOSE(UNIT=11)

END SUBROUTINE booltofile

END MODULE files