SUBROUTINE BLInterp2D(x,y,A,nx,ny,lda,x0,y0,A0)
  INTEGER nx,ny,lda,ix, iy, imx, imy, ipx, ipy
  REAL(8) x(nx),y(ny),x0,y0, dx,dy,dxdy,tol
  COMPLEX*16 A(lda,lda), A0
  tol = 1.d-5
  IF((x0.LT.x(1)-tol).OR.(x0.GT.(x(nx)+tol))) THEN
     PRINT*, 'ERROR in Linterp2D: x0 out of range.'
     PRINT*, 'x0 =', x0
     STOP
  ELSE IF((y0.LT.y(1)-tol).OR.(y0.GT.y(ny)+tol)) THEN
     PRINT*, 'ERROR in Linterp2D: y0 out of range.'
     PRINT*, 'y0 =', y0
     STOP
  END IF

  ! Find im1, ip1 s.t. x(im1) < x0 < x(ip1)
  imx = -1
  ipx = nx*2
  DO ix = 1, nx
     IF(x(ix).GE.x0) THEN
        imx = ix - 1
        ipx = ix
        EXIT
     END IF
  END DO
  ! Same for y0
  imy = -1
  ipy = ny*2
  DO iy = 1, ny
     IF(y(iy).GE.y0) THEN
        imy = iy - 1
        ipy = iy
        EXIT
     END IF
  END DO
  ! Check for end points.
  IF(imx.LT.1) THEN
     imx = 1
     ipx = 2
  END IF
  IF(ipx.GT.nx) THEN
     ipx = nx
     imx = nx -1
  END IF
  IF(imy.LT.1) THEN
     imy = 1
     ipy = 2
  END IF
  IF(ipy.GT.ny) THEN
     ipy = ny
     imy = ny - 1
  END IF
  dx = x(ipx) - x(imx)
  dy = y(ipy) - y(imy)
  dxdy = dx*dy
  ! Linear interpolation in 2-D
  A0 = A(imx,imy)*(x(ipx) - x0)*(y(ipy) - y0) + &
       & A(ipx,imy)*(x0 - x(imx))*(y(ipy) - y0) + &
       & A(imx,ipy)*(x(ipx) - x0)*(y0 - y(imy)) + &
       & A(ipx,ipy)*(x0 - x(imx))*(y0 - y(imy))
  IF(dxdy.EQ.0.d0) THEN
     PRINT*, "ERROR: Duplicate points (x,y) in Linterp2D"
     STOP
  END IF
  A0 = A0/dxdy
  RETURN
END SUBROUTINE BLInterp2D
