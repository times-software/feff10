!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: exjlnl.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine exjlnl (z, l, jl, nl)

!     purpose:  to calculate the spherical bessel functions jl and nl
!               for l = 0 to 6  using exact analytic expression
!
!     arguments:
!       z = argument of jl and nl
!       l = integer order of spherical bessel function
!       jl = jl bessel function (abramowitz conventions)
!       nl = nl bessel function (abramowitz yl conventions)
!            Note that this nl = abramowitz yl.
!
!       analytic expressions from abramowitz 10.1.11 and 10.1.12
!       recurrence relation to get analytic j4,n4  eqns 10.1.19-22 ala

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

!     Exact formulae unstable for very small z, so use series
!     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3)  then
         call bjnser (z, l, jl, nl, 0)
      else
!        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3/z**3 - 1/z)*sinz - 3*cosz/z**2
            nl = (-3/z**3 + 1/z)*cosz - 3*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15/z**4 - 6/z**2)*sinz + (-15/z**3 + 1/z)*cosz
            nl = (-15/z**4 + 6/z**2)*cosz + (-15/z**3 + 1/z)*sinz

         elseif (l .eq. 4)  then
            jl = ( 105/z**5 - 45/z**3 + 1/z )*sinz +                    &
     &                ( -105/z**4 + 10/z**2 )*cosz
            nl = (-105/z**5 + 45/z**3 - 1/z )*cosz +                    &
     &                ( -105/z**4 + 10/z**2 )*sinz

         elseif (l .eq. 5)  then
            jl = ( 945/z**6 - 420/z**4 + 15/z**2 )*sinz +               &
     &              ( -945/z**5 + 105/z**3 - 1/z )*cosz
            nl = (-945/z**6 + 420/z**4 - 15/z**2 )*cosz +               &
     &              ( -945/z**5 + 105/z**3 - 1/z )*sinz

         elseif (l .eq. 6)  then
!    Modified the expression to 1260/z**4 previous value was wrong
!    Author: comorado date:10.14.2014
            jl = ( 10395/z**7 - 4725/z**5 + 210/z**3 - 1/z )*sinz +     &
     &                 ( -10395/z**6 + 1260/z**4 - 21/z**2 )*cosz
            nl = (-10395/z**7 + 4725/z**5 - 210/z**3 + 1/z )*cosz +     &
     &                 ( -10395/z**6 + 1260/z**4 - 21/z**2 )*sinz
!    Author: comorado
!    Added on 10.14.2014 for compton calculation of higher l_fms max
!    Used Abramowitz 10.1.13 expression for jl(z) + Mathematica 9 to do
!    Integral analytically
         elseif (l .eq. 7)  then
          jl = ( 135135/z**8 - 62370/z**6 + 3150/z**4 - 28/z**2 )*sinz +  &
     &             ( -135135/z**7 + 17325/z**5 - 378/z**3 + 1/z )*cosz

          nl = (-135135/z**8 + 62370/z**6 - 3150/z**4 + 28/z**2 )*cosz +  &
     &             ( -135135/z**7 + 17325/z**5 - 378/z**3 + 1/z )*sinz

         elseif (l .eq. 8)  then
          jl = ( 2027025/z**9 - 945945/z**7 + 51975/z**5 - 630/z**3 + 1/z )*sinz +  &
     &                ( -2027025/z**8 + 270270/z**6 - 6930/z**4 + 36/z**2 )*cosz
          nl = (-2027025/z**9 + 945945/z**7 - 51975/z**5 + 630/z**3 - 1/z )*cosz + &
     &                  (-2027025/z**8 + 270270/z**6 - 6930/z**4 + 36/z**2)*sinz

         elseif (l .eq. 9)  then
          jl = ( 34459425/z**10 - 16216200/z**8 + 945945/z**6 - 13860/z**4 + 45/z**2 )*sinz +  &
     &                ( -34459425/z**9 + 4729725/z**7 - 135135/z**5 + 990/z**3 - 1/z )*cosz
          nl = (-34459425/z**10 + 16216200/z**8 - 945945/z**6 + 13860/z**4 - 45/z**2 )*cosz + &
     &                ( -34459425/z**9 + 4729725/z**7 - 135135/z**5 + 990/z**3 - 1/z )*sinz

         else
            stop 'exjlnl, l out of range'
         endif
      endif

      return
      end
