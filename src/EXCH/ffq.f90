!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ffq.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function ffq (q, ef, xk, wp, alph)
      implicit double precision (a-h,o-z)

!     input:  q, wp, alph, ef, xk
!             q is dimensionless, normalized to fermi momentum
!             xk is momentum in invBohrs
!     output: ffq only

      wq = sqrt (wp**2 + alph*q**2 + q**4)
      ffq = (wp+wq)/(q**2) + alph/(2*wp)
      ffq = ((ef*wp) / (4*xk))  * log(ffq)

      return
      end
