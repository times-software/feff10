!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sclmz.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sclmz (rho, lmaxp1, mmaxp1, ileg)

      use dimsmod, only: legtot
	  use constants
      use clmz
      implicit double precision (a-h, o-z)

!     Set CLM(Z) for current leg.
!     Makes clm(z) (eq B11).  Fills array clmi in /clmz/ for ileg,
!     elements clm(0,0) -> clm(lmax+1,mmax+1).
!     If mmaxp1 > lmaxp1, fills m only to lmaxp1.

!     calculates energy dependent factors
!     c(il,im) = c_l^(m)z**m/m! = c_lm    by recursion
!     c_l+1,m = c_l-1,m-(2l+1)z(c_l,m-c_l,m-1, l ne m
!     c_m,m = (-z)**m (2m)!/(2**m m!) with z = 1/i rho
!
!     To test pw approx, set z = 0

      complex*16 rho(legtot)
      complex*16 z, cmm

      cmm = 1
      z = -coni / rho(ileg)

      clmi(1,1,ileg) = (1,0)
      clmi(2,1,ileg) = clmi(1,1,ileg) - z

      lmax = lmaxp1-1

      do 10  il = 2, lmax
         clmi(il+1,1,ileg) =                                            &
     &           clmi(il-1,1,ileg) - z*(2*il-1)*clmi(il,1,ileg)
   10 continue
      mmxp1 = min (mmaxp1, lmaxp1)
      do 20  im = 2, mmxp1
         m = im-1
         imp1 = im+1
         cmm = -cmm * (2*m-1) * z
         clmi(im,im,ileg) = cmm
         clmi(imp1,im,ileg) = cmm * (2*m+1) * (1-im*z)
         do 20  il = imp1, lmax
            l = il-1
            clmi(il+1,im,ileg) = clmi(l,im,ileg) -                      &
     &          (2*l+1) * z * (clmi(il,im,ileg) + clmi(il,m,ileg))
   20 continue

      return
      end
