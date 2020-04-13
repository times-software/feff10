!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xlogx.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COMPLEX*16 FUNCTION xLogx(x)
  DOUBLE PRECISION,INTENT(IN) :: x
  INTEGER sign
  DOUBLE PRECISION,PARAMETER :: tol = 1.d-10

  IF(ABS(x).gt.tol) THEN
     xLogx = x*LOG(x)
  ELSEIF(x.gt.0.d0) THEN
     ! linear interpolation to zero.
     IF(x.gt.0) THEN
        sign = 1
     ELSE
        sign = -1
     END IF
     xLogx = x*LOG(sign*tol)
  END IF

END FUNCTION xLogx
     
