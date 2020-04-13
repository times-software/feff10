!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mmtrxi.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mmtrxi ( rkk, lam1x, bmati, ie, ileg, ilegp, lind)
!     calculates matrix M in Rehr,Albers paper.
!     in polarization case
      use dimsmod, only: mtot, nex, ltot, ntot
	  use constants
      use nlm
      use lambda
      use clmz
      use fmatrx
      use rotmat
      use pdata
      implicit double precision (a-h, o-z)

!     all commons except for /fmat/ are inputs

!     inputs:
!       lam1x:  limits on lambda and lambda'
!       ie:  energy grid points
!       ileg, ilegp: leg and leg'
!
!     Inputs from common:
!        phases, use ph(ie,...,ilegp), and lmax(ie,ilegp)
!        lambda arrays
!        rotation matrix for ilegp
!        clmz for ileg and ilegp
!        path data, eta(ilegp) and ipot(ilegp)
!        xnlm array
!
!     Output:  fmati(...,ilegp) in common /fmatrx/ is set for
!              current energy point.

!     calculate scattering amplitude matrices
!     f(lam,lam') = sum_l tl gam(l,m,n)dri(l,m,m',ileg)gamt(l,m',n')
!                 *cexp(-i*m*eta),  eta = gamma+alpha'
!     lam lt lam1x, lam' lt lam2x such that m(lam) lt l0, n(lam) lt l0
!     gam = (-)**m c_l,n+m*xnlm, gamt = (2l+1)*c_ln/xnlm,
!     gamtl = gamt*tl


      complex*16 cam, camt, tltl, bmati
      dimension bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8), lind(8)
      complex*16  rkk(nex,8)
      complex*16 gam(ltot+1,mtot+1,ntot+1),                             &
     &           gamtl(ltot+1,mtot+1,ntot+1)

!     calculate factors gam and gamtl
 
!     set limits for orbital momentum
      lmn = ltot
      lmx = 0
      do 10 k1 = 1,8
        if (lind(k1).gt.lmx) lmx = lind(k1)
        if (lind(k1).lt.lmn .and. lind(k1).ge.0) lmn = lind(k1)
  10  continue
      iln = lmn + 1
      ilx = lmx + 1
 
      do 30  il = iln, ilx
         tltl = 2*il - 1
         do 20  lam = 1, lam1x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 20
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 50  lam2 = 1, lam1x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = abs(m2) + 1
            fmati(lam1,lam2,ilegp) = 0.0d0

            do 40 k1 = 1, 8
            do 40 k2 = 1, 8
               l1 = lind(k1) + 1
               l2 = lind(k2) + 1
               if (l1.gt.0.and.l2.gt.0 .and. iam1.le.l1.and.iam2.le.l2) &
     &           fmati(lam1,lam2,ilegp) = fmati(lam1,lam2,ilegp) +      &
     &           bmati(m1,k1, m2,k2) * rkk(ie,k1) * rkk(ie,k2) *        &
     &           gam( l1, iam1, in1) * gamtl( l2, iam2, in2)
   40       continue
            fmati(lam1,lam2,ilegp) = fmati(lam1,lam2,ilegp) *           &
     &      exp(-coni*eta(ileg)*m1)
   50    continue
   60 continue

      return
      end
