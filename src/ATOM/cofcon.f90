!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cofcon.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cofcon (a,b,p,q)
!     acceleration of the convergence in the iterative process
!     b is the part of final iteration n is a function of the error (p)
!     (p) at iteration n and the error (q) at the iteration n-1.
!     if the product p*q is positive  b is increased by 0.1
!                        zero b is unchanged
!                        negative b is decreased by 0.1
!     b is between 0.1 and 0.9
!                a = 1. - b
!     ** at the end makes q=p
!
      implicit double precision (a-h,o-z)

      if (p*q)  11,31,21
 11   if (b .ge. 0.2) b = b - 0.1
      go to 31

 21   if (b .le. 0.8) b = b + 0.1

 31   a = 1.0 - b
      q=p
      return
      end
