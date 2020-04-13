!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cpl0.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/11/20 00:09:46 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cpl0 (x, pl0, lmaxp1)
      implicit double precision (a-h, o-z)

!-----------------------------------------------------------------------
!
!     cpl0:  Calculate associated legendre polynomials p_l0(x)
!            by recursion.
!            Adapted from aslgndr.
!
!     first written: (25 june 86) by j. j. rehr
!
!     version 1 (25 june 86) (aslgndr)
!     version 2 (March, 1992) siz
!
!-----------------------------------------------------------------------

      dimension pl0 (lmaxp1)

      lmax = lmaxp1-1

!     calculate legendre polynomials p_l0(x) up to l=lmax
      pl0(1) = 1
	  if(lmaxp1.lt.2) return !KJ avoid exceeding array bounds when lmaxp1=0 e.g. for FPRIME
      pl0(2) = x
	  if(lmaxp1.lt.3) return !KJ same   11-2012
      do 10  il = 2, lmax
         l = il-1
         pl0(il+1) = ( (2*l+1)*x*pl0(il) - l*pl0(l) ) / il
   10 continue

      return
      end
