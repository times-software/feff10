!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strap.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine strap (x, y, n, sum)

!     Trapeziodal integration of y(x), result in sum
!     SINGLE PRECISION
!     modified by ala to handle cases for E<Efermi
!     sum only positive numbers

      dimension x(n), y(n)

      sum = y(1) * abs(x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * abs(x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * abs(x(n) - x(n-1))
      sum = sum/2

      return
      end
