!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmtrxi.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fmtrxi (lam1x, lam2x, ie, ileg, ilegp)

!     all commons except for /fmat/ are inputs

!     inputs:
!       lam1x, lam2x:  limits on lambda and lambda'
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
      use dimsmod !, only: ltot, mtot, lamtot, ntot, legtot, nphx=>nphu, nex
	  use constants
      use nlm
      use lambda
      use clmz
      use fmatrx
      use rotmat
      use pdata
      implicit double precision (a-h, o-z)

      complex*16 cam, camt, cterm, tltl
      complex*16 gam(ltot+1,mtot+1,ntot+1),                             &
     &           gamtl(ltot+1,mtot+1,ntot+1), tl

!     calculate factors gam and gamtl
      iln = 1
      ilx = lmax(ie,ipot(ilegp)) + 1
      do 30  il = iln, ilx
         ll = il - 1
!        do j-average t-matrix
         tl = (exp(2*coni*ph(ie,-ll,ipot(ilegp))) - 1) / (2*coni)
         tltl = (exp(2*coni*ph(ie,ll,ipot(ilegp))) - 1) / (2*coni)
         tltl = tl * (ll+1) + tltl * ll
         lam12x = max (lam1x, lam2x)
         do 20  lam = 1, lam12x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 10
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
   10       if (lam .gt. lam2x) goto 20
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 60  lam2 = 1, lam2x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = iabs(m2) + 1
            imn1 = iam1 + in1 - 1
            cterm = 0
            ilmin = max (iam1, iam2, imn1, in2, iln)
            do 40  il = ilmin, ilx
!              skip terms with mu > l (NB il=l+1, so mu=il is mu>l)
               if (abs(m1).ge.il .or. abs(m2).ge.il)  goto 40
               m1d = m1 + mtot+1
               m2d = m2 + mtot+1

               cterm = cterm + gam(il,iam1,in1)*gamtl(il,iam2,in2)      &
     &                         *dri(il,m1d,m2d,ilegp)

   40       continue
            if (eta(ileg) .ne. 0.0) then
               m1 = mlam(lam1)
               cterm = cterm * exp(-coni*eta(ileg)*m1)
            endif
!           Above was org coding, change to use eta(ilegp) as test
!           based on algebra check.  July 20, 1992, siz&jjr
!           Changed back with redifinition of eta(see rdpath.f)
!           which is more convinient in polarization case.
!           August 8,1993, ala.
!           if (eta(ilegp) .ne. 0.0) then
!              m1 = mlam(lam1)
!              cterm = cterm * exp(-coni*eta(ilegp)*m1)
!           endif
            fmati(lam1,lam2,ilegp) = cterm
   60 continue

!     test of fmati(lam,lam',ileg)
!     plot fmat(lam,lam') = csqrt((z/2)**(m1-m2))*fmat

      return
      end
