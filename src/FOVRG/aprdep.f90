!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: aprdep.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function aprdep (a,b,l)
!     need to be in library for ATOM and PHASE; renamed aprdev
!     the result of this function is the coefficient for the term of 
!     power (l-1) for the product of two polynomes, whose coefficients
!     are in rows a and b 
 
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
 
      aprdep=0.0d00
      do 11 m=1,l
 11      aprdep=aprdep+a(m)*b(l+1-m)
      return
      end
