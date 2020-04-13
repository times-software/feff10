!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mmtrxijas.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/03/30 04:50:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine mmtrxi ( rkk, lam1x, bmati, ie, ileg, ilegp, lind)
      subroutine mmtrxijas (kfinmax,indmax,l2lp,rkk,lam1x,jinit, &
           hbmatl,hbmatr, &
           ie,ileg,ilegp,lind,fmatl,fmatr,ldecmx,ldecs, &
           lgfmatr,lgfmatl)
!     calculates matrix M in Rehr,Albers paper. 
!     in polarization case
!     ala J. A. Soininen

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

      use dimsmod, only: mtot, lamtot, ltot, ntot
	  use constants
	  use global_inp,only: nq,qweights=>qw 
      use nlm
      use lambda
      use clmz
      use fmatrx
      use rotmat
      use pdata
      implicit double precision (a-h, o-z)

      complex*16 cam, camt, tltl, bmati
!      dimension bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8), 
      integer jinit,mj
      integer lind(kfinmax)
      complex*16  rkk(nex,nq,kfinmax)
      complex*16  hbmatl(-jinit:jinit,-mtot:mtot,nq,kfinmax)
      complex*16  hbmatr(-jinit:jinit,-mtot:mtot,nq,kfinmax)
      complex*16  fmatl(-jinit:jinit,nq,lamtot),fmatr(-jinit:jinit,nq,lamtot)

      integer ldecmx,ldecs
      complex*16  lgfmatl(-jinit:jinit,nq,0:ldecs,lamtot), lgfmatr(-jinit:jinit,nq,0:ldecs,lamtot)
      complex*16 itol
      complex*16 gam(ltot+1,mtot+1,ntot+1), gamtl(ltot+1,mtot+1,ntot+1)
      integer l2lp,ll,iq

!     calculate factors gam and gamtl
 
!     set limits for orbital momentum
      lmn = ltot
      lmx = 0
      do 10 k1 = 1,indmax
        if (lind(k1).gt.lmx) lmx = lind(k1)
        if (lind(k1).lt.lmn .and. lind(k1).ge.0) lmn = lind(k1)
  10  continue
      lmx=min(lmx,l2lp)
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


            fmatr(:,:,:)=dcmplx(0.0d0,0.0d0)
            fmatl(:,:,:)=dcmplx(0.0d0,0.0d0)
                  lgfmatr(:,:,:,:)=dcmplx(0.0d0,0.0d0)
                  lgfmatl(:,:,:,:)=dcmplx(0.0d0,0.0d0)


      do lam1=1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do k1=1,indmax
            l1=lind(k1)
            ll=l1
            l1=l1+1
            if ( l1.gt.0 .and. iam1.le.l1 .and. (l1-1).le.l2lp) then  
			   do iq=1,nq
               do mj=-jinit,jinit,2
                  fmatr(mj,iq,lam1)=fmatr(mj,iq,lam1)+rkk(ie,iq,k1)*hbmatr(mj,m1,iq,k1)*gamtl(l1,iam1,in1)  !KJ Aleksi has a bug - he has iq in 2nd position instead of 3rd, which is WRONG
                  fmatl(mj,iq,lam1)=fmatl(mj,iq,lam1)+rkk(ie,iq,k1)*hbmatl(mj,m1,iq,k1)*gam(l1,iam1,in1)
                  if (ll.le.ldecmx) then
                     lgfmatl(mj,iq,ll,lam1)=lgfmatl(mj,iq,ll,lam1)+rkk(ie,iq,k1)*hbmatl(mj,m1,iq,k1)*gam(l1,iam1,in1)
                     lgfmatr(mj,iq,ll,lam1)=lgfmatr(mj,iq,ll,lam1)+rkk(ie,iq,k1)*hbmatr(mj,m1,iq,k1)*gamtl(l1,iam1,in1)
                  end if
               enddo
			   enddo
            end if

         end do
         if (ldecmx.ge.0) then 
            do ll=0,ldecmx
			   do iq=1,nq
               do mj=-jinit,jinit,2
                  lgfmatl(mj,iq,ll,lam1)=lgfmatl(mj,iq,ll,lam1)*exp(-coni*eta(ileg)*m1)*qweights(iq)
               enddo
			   enddo
            end do
         end if
		 do iq=1,nq
         do mj=-jinit,jinit,2
            fmatl(mj,iq,lam1)=fmatl(mj,iq,lam1)*exp(-coni*eta(ileg)*m1)*qweights(iq)
         enddo
		 enddo
      end do



      return
      end
