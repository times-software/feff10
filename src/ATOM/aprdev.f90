!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: aprdev.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function aprdev (a,b,l)
!     the result of this function is the coefficient for the term of 
!     power (l-1) for the product of two polynomes, whose coefficients
!     are in rows a and b 
 
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
 
      aprdev=0.0d00
      do 11 m=1,l
 11      aprdev=aprdev+a(m)*b(l+1-m)
      return
      end
