!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: conv.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine conv(omega,xsec,ne1,vicorr)
!     multiply xsec by theta(omega-efermi) and
!     convolute xsec(omega) with  xloss/((omega-omega0)**2+xloss**2)/pi
!     the result is xsec0(omega0)

      use dimsmod, only: nex
	  use constants
      implicit double precision (a-h, o-z)
      dimension  omega(nex)
      complex*16 xsec(nex), xsec0(nex), xsecdx

      complex*16 conv1
      external conv1

      do 100 ie = 1,ne1
         xsec0(ie) = 0.0d0
         omega0 = omega(ie)
!        Add one more point to correct for the finite grid
!        at large energies. Use linear interpolation.
         dx = max( omega(ne1) - omega(ne1-1), 50*vicorr)
         xlast = omega(ne1)+dx
         dx = dx / ( omega(ne1) - omega(ne1-1))
         xsecdx = xsec(ne1)+ (xsec(ne1)-xsec(ne1-1)) * dx

!        first interval
         do 50  i = 1, ne1-1
            xsec0(ie) = xsec0(ie) +                                     &
     &      conv1(omega(i),omega(i+1),xsec(i),xsec(i+1),omega0,vicorr)
  50     continue
!        last interval
         xsec0(ie) = xsec0(ie) +                                        &
     &   conv1(omega(ne1),xlast,xsec(ne1),xsecdx,omega0,vicorr)
         xsec0(ie) = xsec0(ie) /real(pi)
  100 continue
      do 200 ie = 1, ne1
  200 xsec(ie) = xsec0(ie)

      return
      end

      complex*16 function conv1(x1,x2,y1,y2,x0,xloss)
!     convolution of function 1/(omega-omega0-i*xloss)/pi
!     makes linear interpolation for function between x1,x2 and
!     takes advantage that the integral can be taken analytically.
      implicit double precision (a-h, o-z)
      complex*16  y1, y2, t, coni,dum, a, b
      parameter (coni = (0.0,1.0))

      d = (x2-x1) / 2.0
      a = dble(y2-y1) / 2.0
      b = dble(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = dimag (dum)

      d = (x2-x1) / 2.0
      a = dimag(y2-y1) / 2.0
      b = dimag(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = conv1 + coni* dimag( dum)

      return
      end
