!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mmtrxijas0.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/03/30 04:50:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine mmtrxi ( rkk, lam1x, bmati, ie, ileg, ilegp, lind)
      subroutine mmtrxijas0 (kfinmax,indmax,l2lp,rkk,lam1x,jinit, &
           hbmatrs,ie,ileg,ilegp,lind,fmats,ldecmx,ldecs,lgfmats)
!     calculates matrix M in Rehr,Albers paper. 
!     in polarization case
!     ala J. A. Soininen
      use dimsmod, only: mtot, lamtot, ntot, ltot
	  use constants
	  use global_inp,only: nq,qweights=>qw 
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
      integer jinit,mj,iq
      integer lind(kfinmax)
      complex*16  rkk(nex,nq,kfinmax)
      complex*16  hbmatrs(-jinit:jinit,0:1,-mtot:mtot,-mtot:mtot,kfinmax)
      complex*16  fmats(-jinit:jinit,0:1,lamtot,lamtot)

      integer ldecmx,ldecs
      complex*16  lgfmats(-jinit:jinit,0:1,0:ldecs,lamtot,lamtot)
      complex*16 itol
      complex*16 gam(ltot+1,mtot+1,ntot+1),gamtl(ltot+1,mtot+1,ntot+1),qsum(kfinmax)
      integer l2lp,ll
      integer lam2,lam1,is2
      integer m1,in1,iam1
      integer m2,in2,iam2
      integer lj
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
      
	  
!    do the sum over the momentum transfers here

       do k1=1,indmax
          qsum(k1)=dcmplx(0.0d0,0.0d0)
          do iq=1,nq
             qsum(k1)=qsum(k1) +rkk(ie,iq,k1)*rkk(ie,iq,k1)*qweights(iq)
          end do
       end do
       
	  

	  fmats(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	  lgfmats(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)


      do lam1=1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do lam2=1,lam1x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = abs(m2) + 1
            do k1=1,indmax
               l1=lind(k1)
               ll=l1
               lj=lind(k1)
               l1=l1+1
               if ( l1.gt.0 .and. iam1.le.l1 .and. (l1-1).le.l2lp) then
                  do is2=0,1
                     do mj=-jinit,jinit,2
                        fmats(mj,is2,lam2,lam1) &
                             =fmats(mj,is2,lam2,lam1)+qsum(k1)*hbmatrs(mj,is2,m2,m1,k1)* &
                             gamtl(l1,iam2,in2)*gam(l1,iam1,in1)/dble(2*lj+1)

                        if (ll.le.ldecmx) then
                           lgfmats(mj,is2,ll,lam2,lam1) &
                                =lgfmats(mj,is2,ll,lam2,lam1)+qsum(k1)* &
                                hbmatrs(mj,is2,m2,m1,k1)*gamtl(l1,iam2,in2)*gam(l1,iam1,in1)/dble(2*lj+1)
                        end if
                     end do
                  end do
               end if
            end do
         end do !lam2
         if (ldecmx.ge.0) then 
            do lam2=1,lam1x
               do ll=0,ldecmx
                  do is2=0,1
                     do mj=-jinit,jinit,2
                        lgfmats(mj,is2,ll,lam2,lam1)=lgfmats(mj,is2,ll,lam2,lam1)*exp(-coni*eta(ileg)*m1)
                     end do
                  end do
               end do
            end do
         end if
      end do ! lam 1


      return
      end
