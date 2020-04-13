!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mmtrjas.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2011/12/11 02:39:17 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mmtrjas( hbmatr,hbmatl, ipol, ispin, ptz, nsp,pha,qbeta) !KJ ,jinit,kfinmax,indmax) !KJ 9-09,jmax,ljmax,indmax,lind)
!     calculates the part of matrix M which does not depend on energy
!     point.( see Rehr and Albers paper)
!     for path expansion always neglect spin-flip processes
!     to simplify calculations; (bmati does not have spin indices)
      use dimsmod, only: mtot, nspx=>nspu
	  use constants
      use nrixs_inp
	  use global_inp,only: nq
      use rotmat
      use pdata

      implicit none !KJdouble precision (a-h, o-z)

!     all commons are inputs
!     Inputs from common:
!        kinit: quantum number kappa for initial orbital
!        rotation matrix for ilegp
!        path data, eta(ilegp) and ipot(ilegp)
!        mtot,l0
!        polarization data : ipol, ispin, ptz
!     Output:  bmati(...) 
      integer,intent(in) :: ipol,nsp !KJ ,jinit,kfinmax,indmax
	  integer minit,is,ispin,lxx,k1,mu1,mu1d,m1,m1d,ims1,is1,mj1, mlim

      complex*16 ptz(-1:1, -1:1)

      complex*16 hbmatr(-jinit:jinit,nq,-mtot:mtot,kfinmax)
      complex*16 hbmatl(-jinit:jinit,nq,-mtot:mtot,kfinmax)
      complex*16 bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8)
      real*8 hbmat(0:1,kfinmax,-jinit:jinit)
	  complex pha(nq),pham
	  real*8 qbeta(nq)
	  complex*16 rotdri(2*mtot+1),pheta,phetam,rotlm1m2
	  real*8 rot 
	  real*8,external :: rotwig

      complex*16 itolj,itolg

!KJ 9-09      dimension kiind(kfinmax), lgind(kfinmax)!,ljind(kfinmax) JK - 09/2009
!KJ      integer lind(kfinmax)
      integer l1,l1p,iq,m2,m2d
      integer mj
      logical ltrace


      hbmatr=dcmplx(0.0d0,0.0d0)
	  hbmatl=dcmplx(0.0d0,0.0d0)
      hbmat=dcmplx(0.d0,0.d0)
      do minit=-jinit,jinit,2
         call bcoefjas(kinit,minit,ltot,hbmat(0,1,minit))
      end do
	  
	  is = 0
      if (ispin.eq.1) is = nspx - 1

!     ilinit = initial orb. momentum + 1. 

      lxx=MIN(mtot,MAXVAL(lind))
      if (ipol.eq.0) then 
         write(6,*) "ipol=0 in mmtrjas"
      end if
!
!     lefthand side matrix
      pheta=exp(-coni*eta(0))
	  
      do k1=1,indmax
         l1=lind(k1)
         mlim=MIN(lxx,l1) ! JK - mlim introduced to limit m1, m2, since m1d can't be less than 1.
         l1p=l1+1
         itolj=dcmplx(0.0d0,1.0d0)**ljind(k1)
         itolg=dcmplx(0.0d0,1.0d0)**l1
         do mu1=-lxx,lxx
            mu1d = mu1+mtot+1
			do iq=1,nq
			!here we apply the rotation for the q-vector
            do m1=-mlim,mlim ! JK - used to be -l1, l1
               m1d = m1 + mtot+1
			   !if(m1d.gt.0) then !KJ 12-2011 added "if" since rotdri(0) doesn't exist ...
			   rotdri(m1d)=dcmplx(0.0d0,0.0d0)
			   do m2=-mlim, mlim ! JK Used to be -l1,l1
                      m2d=m2 + mtot+1
                      pham=pha(iq)**m2
                      phetam=pheta**m2
                      rot=rotwig(qbeta(iq),l1,m2,m1,1)
                      rotlm1m2=rot*pham
                      rotdri(m1d)=rotdri(m1d)+phetam*dri(l1p,mu1d,m2d,nsc+2)*rotlm1m2
			   end do
			   !endif
			end do
!              rotating is done			   
			do m1=-mlim, mlim !JK - Used to be -l1,l1
			   m1d=m1+mtot+1
               do ims1=1,nsp
                  is1=ims1-1
                  mj1=2*m1+(2*is1-1)
                  if (abs(mj1).le.jinit) then
                     hbmatl(mj1,iq,mu1,k1)=hbmatl(mj1,iq,mu1,k1)+hbmat(is1,k1,mj1)* rotdri(m1d)*dconjg(itolg)*dconjg(itolj)
                  end if
               end do    
            end do
         end do
		 enddo
		 
      end do
!
!     righthand side matrix
!
      pheta=exp(-coni*eta(nsc+2))

      do k1=1,indmax
         l1=MIN(lind(k1),mtot)
         l1p=l1+1
         itolj=dcmplx(0.0d0,1.0d0)**ljind(k1)
         itolg=dcmplx(0.0d0,1.0d0)**l1
         do mu1=-lxx,lxx
            mu1d = mu1+mtot+1
			do iq=1,nq
               !Here we apply the rotation for the q-vector        
                do m1=-mlim, mlim !-l1,l1 
                   m1d = m1 + mtot+1
!			       if(m1d.gt.0) then !KJ 12-2011 added "if" since rotdri(0) doesn't exist ...				   
                   rotdri(m1d)=dcmplx(0.0d0,0.0d0)
                   do m2=-mlim, mlim !-l1,l1
                      m2d=m2 + mtot+1
                      pham=pha(iq)**m2
                      phetam=pheta**m2
                      rot=rotwig(qbeta(iq),l1,m2,m1,1)
                      rotlm1m2=conjg(rot*pham)
                      rotdri(m1d)=rotdri(m1d)+phetam*dri(l1p,m2d,mu1d,nleg)*rotlm1m2
                   end do
!				   endif
                end do			

            do m1=-mlim, mlim !-l1,l1
               m1d = m1 + mtot+1
               do ims1=1,nsp
                  is1=ims1-1
                  mj1=2*m1+(2*is1-1)
                  if (abs(mj1).le.jinit) then
					hbmatr(mj1,iq,mu1,k1)=hbmatr(mj1,iq,mu1,k1)+hbmat(is1,k1,mj1)* &
					rotdri(m1d)*(itolg)*(itolj)
                  end if
               end do    
            enddo
			enddo
         end do
      end do


      return
      end


