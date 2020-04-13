!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: polint.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine polint( xa, ya, n, x, y, dy)
!     draws a polynimial P(x) of order (n-1) through n points.
!     returns y = P(x) and dy - estimate of the error
!     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) stop 'failure in polint'    !pause to stop  KJ 6-2012
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
