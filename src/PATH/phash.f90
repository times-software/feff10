!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phash.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine phash (npat, ipat, rx, ry, rz, dhash)
!     hashes a path into double precision real dhash

      use dimsmod, only: npatx, natx
      double precision dhash
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      double precision xx

      parameter (iscale = 1000)
      parameter (factor = 16.12345678)
      parameter (facto2 = 8.57654321)

!     Hashing scheme: Assume about 15 significant digits in a double 
!     precision number.  This is 53 bit mantissa and 11 bits for sign 
!     and exponent, vax g_floating and probably most other machines.
!     With max of 9 legs, 47**9 = 1.12e15, so with a number less than 
!     47, we can use all these digits, scaling each leg's data by 
!     47**(j-1).  Actually, since our numbers can go up to about 10,000,
!     we should keep total number < 1.0e11, 17**9 = 1.18e11, which means
!     a factor a bit less than 17.  Choose 16.12345678, a non-integer,
!     to help avoid hash collisions.

!     iscale and 'int' below are to strip off trailing digits, which
!     may contain roundoff errors

      dhash = 0
      do 210  j = 1, npat
         xx = factor**(j-1)
         dhash = dhash + xx * (nint(rx(j)*iscale) +                     &
     &               nint(ry(j)*iscale)*0.894375 +                      &
     &               nint(rz(j)*iscale)*0.573498)
  210 continue
      do 220  j = 1, npat
         xx = facto2**(j-1)
         dhash = dhash + xx * iscale * ipot(ipat(j))
!        dhash = dhash + xx * ipot(ipat(j))
  220 continue
      dhash = dhash + npat * 40000000

      return
      end
