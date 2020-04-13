!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: snlm.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine snlm (lmaxp1, mmaxp1)

      use dimsmod, only: ltot, mtot
      use nlm
      implicit double precision(a-h,o-z)

!     Set nlm, legendre normalization factors, xnlm in common /nlm/
!     Calculates legendre norm factors
!     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)

!     flg(i) = i! * afac**i, set in factst
      dimension flg(0:210)

      call factst (afac, flg)

!     initialize xnlm explicitly
      do 5  il = 1, ltot+1
      do 5  im = 1, mtot+1
         xnlm(il,im) = 0
    5 continue

      do 10  il = 1, lmaxp1
         mmxp1 = min (mmaxp1, il)
         do 10  im = 1, mmxp1
            l = il-1
            m = im-1
            cnlm = (2*l+1) * flg(l-m) / flg(l+m)
            cnlm = sqrt(cnlm) * afac**m
            xnlm(il,im) = cnlm
   10 continue

      return
      end
      subroutine factst (afac, flg)
      implicit double precision (a-h,o-z)

!     FACTorial SeT, flg(i) = i! * afac**i
      dimension flg(0:210)

!     afac = 1/64 works with double precision on a VAX
      afac = 1./64.

      flg(0) = 1
      flg(1) = afac

      do 10  i = 2, 210
   10 flg(i) = flg(i-1) * i * afac

      return
      end
