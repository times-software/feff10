!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cubic.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cubic (xk0, wp, alph, rad, qplus, qminus)

!     input:  xk0, wp, alph
!     output: rad, qplus, qminus

      implicit double precision (a-h, o-z)
      complex*16 s1,s13
      parameter (three = 3)
      parameter (third = 1/three)

!     this subroutine finds the roots of the equation
!     4xk0 * q^3  +  (alph-4xk0^2) * q^2  +  wp^2 = 0
!     see abramowitz and stegun pg 17 for formulae.

      a2 = (alph / (4*xk0**2)  -  1) * xk0
      a0 = wp**2 / (4*xk0)
      a1 = 0
      q = a1/3 - a2**2/9
      r = (a1*a2 - 3*a0)/6  -  a2**3/27
      rad = q**3 + r**2
      if (rad .gt. 0) then
         qplus = 0
         qminus = 0
         return
      endif

      s13 = dcmplx (r, sqrt(-rad))
      s1 = s13 ** third
      qz1 = 2*s1 - a2/3
!     qz2 = -(s1 + sqrt(three)*dimag(s1) + a2/3)
      qz3 = -(s1 - sqrt(three)*dimag(s1) + a2/3)
      qplus = qz1
      qminus = qz3

      return
      end
