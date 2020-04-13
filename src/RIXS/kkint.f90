COMPLEX*16 FUNCTION KKInt(a, b, x0, x1, gam, x)
  IMPLICIT NONE
  COMPLEX*16 a, b
  REAL(8) x0, x1, gam, Pi, x
  INTEGER ix
  COMPLEX*16, PARAMETER :: coni = (0.d0,1.d0)
  Pi = 3.141592653d0
  
! Integrate f(x')/(x' - x + i gam) from x0 to x1 where
! f(x') = a*x' + b.
! x1 should be > x0
  KKInt = 0.d0
  IF(x0.GE.x1) THEN
     print*,'x0= ',x0
     print*,'x1= ',x1
     PRINT*, 'Error in KKInt: x0 > x1.'
     STOP
  END IF
  IF((x.NE.x0).AND.(x.NE.x1)) THEN
     KKInt = a*( x1 - x0 + (gam+coni*x)*( ATAN(gam/(x-x0)) - ATAN(gam/(x-x1)) &
          & + 0.5d0*coni*LOG( (gam**2 + (x-x0)**2)/(gam**2 + (x-x1)**2) ) - Pi ) )
     IF((x.LT.x0).OR.(x.GT.x1)) THEN
        KKInt = KKInt + a*Pi*(gam + coni*x)
     END IF
     KKInt = (KKInt + b*Log( (x1-x+coni*gam)/(x0-x+coni*gam) ))
  ELSEIF(x.EQ.x0) THEN
     KKInt = a*( x1 - x0 - (gam + coni*x0)*(ATAN(gam/(x0-x1)) + &
          & 0.5d0*coni*LOG((gam**2 + (x0-x1)**2)/gam**2) + 0.5d0*Pi)) + &
          & b*LOG(gam/(gam - coni*(x1-x0)))
  ELSE
     KKInt = a*( x1 - x0 - (gam + coni*x1)*(ATAN(gam/(x0-x1)) + &
          & 0.5d0*coni*LOG(gam**2/(gam**2 + (x0-x1)**2)) + 0.5d0*Pi) ) + &
          & b*LOG((gam - coni*(x1-x0))/gam)
  END IF
  !KKInt = - DBLE(KKInt) - coni*DIMAG(KKInt)
END FUNCTION KKInt
