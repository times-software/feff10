!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: besjh.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine besjh (x, lbmax, jl, hl)

!-----------------------------------------------------------------------
!
!     purpose:  to calculate the spherical bessel functions jl and hl
!               for l = 0 to lbmax (no offset)
!
!     arguments:
!       x = argument of jl and nl
!       lbmax
!       jl = jl bessel function (abramowitz conventions)
!       hl = hl^+ bessel function (messiah conventions) for Im x >=0
!       hl = hl^- bessel function (messiah conventions) for Im x < 0
!       jl and hl must be dimensioned 
!            complex*16 jl(0:lbmax), hl(0:lbmax), 
!
!     notes:  jl and hl should be calculated at least to 10 place
!             accuracy for the range 0<x<100 according to spot
!             checks with tables
!
!     error messages written with PRINT statement.
!
!     first coded by r. c. albers on 14 dec 82
!
!     version 3
!
!     last modified: 27 jan 83 by r. c. albers
!     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
!     modified again, siz, June 1992
!     rewritten for jl and hl by a.l. ankudinov feb 2000
!
!-----------------------------------------------------------------------
      use dimsmod, only: ltot
      implicit double precision (a-h, o-z)

      complex*16 x
      complex*16 jl(0:lbmax), nl(ltot+2)
      complex*16 hl(0:lbmax)
      complex*16 cjl(ltot+2), sjl(ltot+2)

      complex*16 xjl,xnl,asx,acx, epx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      !parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)
      parameter (xcut = 1.d0, xcut1 = 17.d0, xcut2 = 5.01d0)
      complex*16 coni
      parameter (coni=(0,1))

      if (dble(x) .lt. 0)  stop 'Re(x) is .lt. zero in besjh'

      lmax = min(lbmax, ltot+1)
      lmaxp1 = lmax + 1

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
!        case Re(x) < 1, just use series expansion
         do 10 ll = 0,lmax
            ifl = 0
            call bjnser (x,ll,xjl,xnl,ifl)
            jl(ll) = xjl
            hl(ll) = -xnl + coni*xjl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

!        case 1 <= Re(x) < 7.5

         call bjnser (x,lmax,xjl,xnl,1)
         jl(lmax) = xjl

         call bjnser (x,lmax-1,xjl,xnl,1)
         jl(lmax-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
!           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
!           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

!        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(l) = tlxp3 * jl(l+1) / x  -  jl(l+2)
   60    continue

         do 65 il = 1, lmaxp1
            l = il - 1
            hl(l) = -nl(il) + coni*jl(l)
   65    continue

      else
!        case Re(x) > 7.5
!        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
!        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
!        scrambled (mod 4) here, since n is integer.  These are hard-
!        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5                  &
     &            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6            &
     &            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7        &
     &            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5             &
     &             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6         &
     &             + 25740*xi4 - 55*xi2
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = DBLE(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         if (dimag(x).ge. 0.d0) then
           epx = exp(coni*x)
         else 
           epx = exp(-coni*x)
         endif
         do 110 ll = 0,lmax
            lp1 = ll + 1
            jl(ll) = asx*sjl(lp1)+acx*cjl(lp1)
            if (dimag(x).ge. 0.d0) then
              hl(ll) = (sjl(lp1)+coni*cjl(lp1)) * epx
            else
              hl(ll) = (sjl(lp1)-coni*cjl(lp1)) * epx
            endif
  110    continue
      endif

      return
      end
