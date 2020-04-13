!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xmult.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xmult (k, kp, ls, lb, xm1, xm2)

      use constants
      implicit double precision (a-h, o-z)
      complex*16 xm1, xm2, alslb
!     xm1, xm2 both either real or pure imaginary
      integer a, ap

!     see Grant eq. 6.30. calculate the factors 
!     <k|alpha*A( l, L)|k'> = (-)**(j-m) * 3j( j L j'; -m p m')*R_k,k'
!     R_k,k'(l,L) = \int dr (xm1*P_k*Q_k'+ xm2*Q_k*P_k') * j_l(wr)

!     set the factor in front of bessel function (eq.6.26)
      if (ls+1.eq.lb) then
!        e.g. dipole and quadrupole transition
         aa = (2*lb-1) * (lb+1) / 2.d0
         alslb = coni**ls * sqrt(aa)
      elseif (ls-1.eq.lb) then
!        e.g. cross dipole-octupole
         aa = (2*lb+3) * lb / 2.d0
         alslb = coni**ls * sqrt(aa)
      elseif (ls.eq.lb) then
!        e.g. magnetic dipole
         alslb = coni**ls * (2*lb+1) /sqrt(2.d0)
      else
         alslb = 0
      endif

!     set all angular momenta
      j2 = 2*abs(k) -1
      a = 1
      if (k.gt.0) a=-1
      jp2 = 2*abs(kp) -1
      ap = 1
      if (kp.gt.0) ap=-1

!     calculate xm1 (beta=1 in eq.6.30)
!     check out 2 Kronecker symbols
      lam = (j2-a) / 2
      lamp = (jp2+ap) / 2
      if ( 2*lam.eq.j2-a .and. 2*lamp.eq.jp2+ap) then
         call ninej (lam, lamp, ls, j2,jp2, lb, aa)
         xm1 = alslb * aa * cwig3j(lam, ls, lamp, 0, 0, 1) * (-1)**lam  &
     &        * sqrt(6.d0*(j2+1)*(jp2+1)*(2*lb+1)*(2*lam+1)*(2*lamp+1) )
         xm1 = xm1 * coni
      else
         xm1 = 0
      endif

!     calculate xm2 (beta=-1 in eq.6.30)
!     check out 2 Kronecker symbols
      lam = (j2+a) / 2
      lamp = (jp2-ap) / 2
      if ( 2*lam.eq.j2+a .and. 2*lamp.eq.jp2-ap) then
         call ninej (lam, lamp, ls, j2,jp2, lb, aa)
         xm2 = alslb * aa * cwig3j(lam, ls, lamp, 0, 0, 1) * (-1)**lam  &
     &       * sqrt(6.d0*(j2+1)*(jp2+1)*(2*lb+1)*(2*lam+1)*(2*lamp+1) )
!        factor -1 due to complex conjugation of i*Q_k
         xm2 = - coni * xm2
      else
         xm2 = 0
      endif

      return
      end

      subroutine ninej (lam, lamp, ls, j2,jp2, lb, aa)
      implicit double precision (a-h, o-z)
!     calculate 9j-symbol in 6.30 of Grant using eq. C.41 in Messiah

      if (ls.gt.lb) then
        aa = - (ls+lb+1)* sixj(1,2,2*lb,ls+lb,2*ls) *                   &
     &       sixj(2*lb, ls+lb, 2*lamp, jp2, j2) *                       &
     &       sixj(ls+lb,2*ls, 2*lam, j2, 2*lamp)
      elseif (ls.lt.lb) then
        aa = - (ls+lb+1)* sixj(1,2,2*lb,ls+lb,2*ls) *                   &
     &       sixj(ls+lb, 2*lb, jp2, 2*lamp, j2) *                       &
     &       sixj(2*ls, ls+lb, j2, 2*lam, 2*lamp)
      else
!       ls=lb (magnetic dipole)
        aa = -(2*ls+2) * sixj(1,2,2*lb,2*lb+1,2*lb) *                   &
     &       sixj(2*lb, 2*lb+1, 2*lamp, jp2, j2) *                      &
     &       sixj(2*lb, 2*lb+1, j2, 2*lam, 2*lamp)
        aa = aa -(2*ls) * sixj(1,2,2*lb,2*lb-1,2*lb) *                  &
     &       sixj(2*lb-1, 2*lb, jp2, 2*lamp, j2) *                      &
     &       sixj(2*lb-1, 2*lb, 2*lam, j2, 2*lamp)
      endif

      return
      end

      double precision function sixj(j1,j2,j3,j4,j5)
      implicit double precision (a-h, o-z)
!     calculate 6j symbols in eq. c.38, c39 of Messiah
!     all input angular momenta are multiplied by 2 and
!     j2 should be equal to j1+1
      integer g2

      aa = 0
      if (j2.eq.j1+1) then
        if (j4.eq.j3+1) then
!         eq.c.38
          g2 = j5 - 1
          if (g2.ge.abs(j1-j3) .and. g2.le.j1+j3) then
            aa = (1.d0 + (g2+j1-j3)/2.d0) * (1.d0 +(g2-j1+j3)/2.d0) /   &
     &           (j1+1) /(j1+2)/(j3+1)/(j3+2)
            aa = sqrt(aa) * (-1)**(nint(1+(g2+j1+j3)/2.d0))
          endif
        elseif(j3.eq.j4+1) then
!         eq.c.39
          g2 = j5
          if (g2.ge.abs(j1-j4) .and. g2.le.j1+j4) then
            aa = (1.d0 - (g2-j1-j4)/2.d0) * (2.d0 +(g2+j1+j4)/2.d0) /   &
     &           (j1+1) /(j1+2)/(j4+1)/(j4+2)
            aa = sqrt(aa) * (-1)**(nint(1+(g2+j1+j4)/2.d0))
          endif
        endif
      endif
      sixj = aa

      return
      end

