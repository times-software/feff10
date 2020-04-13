PROGRAM test
  COMPLEX*16 A(100,100), A0
  REAL(8) x(100), y(100), x0, y0, dx, dy
  INTEGER i1, i2

  DO i1 = 1, 100
     x(i1) = i1**2
     y(i1) = i1**2
  END DO
  DO i1 = 1, 100
     DO i2 = 1, 100
        A(i1,i2) = SIN(1.d0/1000.d0*x(i1))*SIN(1.d0/1000.d0*y(i2))
        WRITE(14,'(3f20.10)') x(i1), y(i2), DBLE(A(i1,i2))
     END DO
     WRITE(14,*)
  END DO
  
  dx = (x(100) - x(1))/99.d0
  dy = (y(100) - y(1))/99.d0
  DO i1 = 1, 100
     x0 = x(1) + (i1 - 1)*dx
     DO i2 = 1, 100
        y0 = y(1) + (i2 - 1)*dy        
        CALL BLInterp2D(x,y,A,100,100,x0,y0,A0)
        WRITE(15,'(3f20.10)') x0, y0, DBLE(A0)
     END DO
     WRITE(15,*)
  END DO
END PROGRAM test
