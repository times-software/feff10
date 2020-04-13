!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: exconv.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/02/03 07:17:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine exconv (omega, nk, efermi, s02, erelax, wp, xmu)
      parameter (nfinex = 601)
!     convolution of xmu(e) with excitation spectrum, which is modeled
!     by: f(e) = s02*delta(e) + theta(e)*exp(-e/ed)*x1 + fp(e)
!     plasmon contribution modeled by fp(e) = theta(e-wp)*exp(-e/wp)*x2
!     normalization factors x1, x2 and distribution parameter ed are
!     found from conditions: 1) integral d(e)*f(e) = 1 
!     2) integral d(e)*fp(e) = wwp  0<=wwp<1  s02+wwp<=1
!     3) integral d(e)*e*f(e) = erelax
!     Input:
!       omega - energy grid (e)
!       nk    - number of points in energy grid
!       efermi- fermi level
!       s02   - overlap with fully relaxed orbitals
!       erelax- relaxation energy 
!       wp    - plasmon frequency
!       xmu   - original absorption coefficient
!     Output
!       xmu  - result of convolution, rewritten at the end. 
!     This subroutine uses the fact, that if convolution is made for
!     e(i), then for e(i+1), the convolution integral with exp(-e/ed)
!     for e<e(i) is simply scaled by exp((e(i)-e(i+1)) / ed). This makes
!     this convolution fast.
!     written by ala. december 1995
!     
      implicit double precision (a-h,o-z)
      dimension  omega(nk), xmu(nk)
!     work space
      dimension  slope(nfinex), dmu(nfinex), xmup(nfinex)

      if (s02 .ge. 0.999) return
      if (wp .le. 0.0) wp = 0.00001
      if (nk .gt. nfinex) then
         call par_stop('check nfinex in subroutine exconv')
      endif
!     change weight for plasmon excitation here
      wwp = 0.00
!     sm1 - weight for shake up (off) processes
      sm1 = 1.0 - s02 - wwp
      edp = wp
      ed = (erelax - wwp * (wp + edp)) / sm1
      i0 = locat (efermi, nk, omega)
      do 10 i = 1, i0
         slope(i) = 0.0
         dmu(i) = 0.0
 10   continue
      do 20 i = i0, nk - 1
 20     slope(i) = ed * (xmu(i+1) - xmu(i)) / (omega(i+1) - omega(i))
      call terp (omega, xmu, nk, 1, efermi, xmuf)

!     start induction
      xmult = exp ((efermi - omega(i0+1)) / ed)
      dmu(i0+1) = xmu(i0+1) - slope(i0) - xmult * (xmuf - slope(i0))
      do 50 i = i0 + 1, nk - 1
         xmult = exp ((omega(i) - omega(i+1)) / ed)
         dmu(i+1) = xmu(i+1) - slope(i) + xmult*(dmu(i)-xmu(i)+slope(i))
 50   continue
      do 55 i = 1, nk
 55   xmup(i) = s02 * xmu(i) + sm1 * dmu(i)

!     do convolution with plasmon pole
      do 60 i = i0, nk - 1
 60     slope(i) = slope(i) / ed * edp
      xmult = exp ((efermi - omega(i0+1)) / edp)
      dmu(i0+1) = xmu(i0+1) - slope(i0) - xmult * (xmuf - slope(i0))
      do 70 i = i0 + 1, nk - 1
         xmult = exp ((omega(i) - omega(i+1)) / edp)
         dmu(i+1) = xmu(i+1) - slope(i) + xmult*(dmu(i)-xmu(i)+slope(i))
 70   continue

      do 90 i = 1, nk
        en = omega(i) - wp
        j0 = locat(en, nk, omega)
        if (en .gt. efermi) then
           xmult = exp ((omega(j0) - en) / edp)
           dif = xmu(j0) - slope(j0)
           xmup(i) = xmup(i) + wwp * (xmult * (dmu(j0) - dif) + dif +   &
     &                            slope(j0) * (en - omega(j0)) / edp)
        endif
 90   continue

      do 200 i = 1, nk
 200  xmu(i) = xmup(i)

      return
      end
