!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bjnser.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine bjnser (x, l, jl, nl, ifl)

!-----------------------------------------------------------------------
!
!     subroutine: bjnser (x,l,jl,nl,ifl)
!
!     purpose:  to calculate the spherical bessel functions jl and nl
!
!     arguments:
!       x = argument of jl and nl
!       l = l value calculated (no offset)
!       jl = jl bessel function (abramowitz conventions)
!       nl = nl bessel function (abramowitz yl conventions)
!       ifl = 0 return both jl and nl
!             1 return jl only
!             2 return nl only
!
!     notes:  jl and nl are calculated by a series
!             expansion according to 10.1.2 and 10.1.3
!             in abramowitz and stegun (ninth printing),
!             page 437
!
!             error msgs written with PRINT statements.
!
!     first coded by r. c. albers on 26 jan 83
!
!     version 2
!
!     last modified: 27 jan 83 by r. c. albers
!
!-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      character*512 slog

      parameter (niter = 160, tol = 1.e-15)

      if (l .lt. 0) then
         call wlog(' l .lt. 0 in bjnser')
         stop 'bjnser 1'
      endif
      if (dble(x).lt. 0.) then
         write(slog,30) x
         call wlog(slog)
   30    format (' x = ', 1p, 2e14.6, ' is .le. 0 in bjnser')
         stop 'bjnser 2'
      endif

      lp1 = l+1
      u = x**2 / 2

!     make djl = 1 * 3 * 5 * ... * (2*l+1),
!          dnl = 1 * 3 * 5 * ... * (2*l-1)
      djl = 1
      fac = -1
      do 50 il = 1, lp1
         fac = fac + 2
         djl = fac * djl
   50 continue
      dnl = djl / (2*l+1)


      if (ifl .eq. 2)   goto 90
!     make jl
!     pj is term in { } in 10.1.2, del is last factor in the series
!     convergence test is (last factor)/(total term) <= tol
      pj = 1
      nf = 1
      nfac = 2*l + 3
      den = nfac
      sgn = -1
      ux = u
      do 60 il = 1, niter
         del = sgn*ux / den
         pj = pj + del
         trel = abs (del / pj)
         if (trel .le. tol)  goto 80
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
   60 continue
      stop  'jl does not converge in bjnser'
   80 jl = pj * (x**l) / djl

   90 if (ifl.eq.1) return
!     make nl
!     pn is term in { } in 10.1.3, del is last factor in the series
!     convergence test is (last factor)/(total term) <= tol
      pn = 1
      nf = 1
      nfac = 1 - 2*l
      den = nfac
      sgn = -1
      ux = u
      do 100  il = 1, niter
         del = sgn * ux / den
         pn = pn + del
         trel = abs (del / pn)
         if (trel .le. tol) goto 120
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
  100 continue
      stop  'nl does not converge in bjnser'
  120 nl = -pn * dnl / (x**lp1)

      return
      end
