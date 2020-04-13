!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mmtrjas0.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mmtrjas0( hbmatrs, ipol, ispin, ptz, nsp,lx1)
!KJ   subroutine mmtrjas0( hbmatrs, ipol, ispin, le2, angks, ptz, &
!KJ        nsp,jinit,jmax,kfinmax,ljmax,lx1,indmax,lind)
!     calculates the part of matrix M which does not depend on energy
!     point.( see Rehr and Albers paper)
!     for path expansion always neglect spin-flip processes
!     to simplify calculations; (bmati does not have spin indices)
      use dimsmod, only: mtot, nspx=>nspu
	  use constants
	  use nrixs_inp,kiind=>kind
      use rotmat
      use pdata

      implicit none

!     all commons are inputs
!     Inputs from common:
!        kinit: quantum number kappa for initial orbital
!        rotation matrix for ilegp
!        path data, eta(ilegp) and ipot(ilegp)
!        mtot,l0
!        polarization data : ipol, ispin, ptz
!     Output:  bmati(...) 

      complex*16 ptz(-1:1, -1:1)

      complex*16 bmati, bmat,hbmatr,hbmatrs,hbmatl
      
      dimension hbmatrs(-jinit:jinit,0:1,-mtot:mtot,-mtot:mtot,kfinmax)
      dimension bmati(-mtot:mtot, 1:8, -mtot:mtot, 1:8)
      real*8 hbmat(0:1,kfinmax,-jinit:jinit)

      complex*16 itolj,itolg

      integer jind(kfinmax),indmat(kfinmax)
      integer l1,l1p,l2
      integer mj
      logical ltrace

      integer ainit,lj,lt,ji
      integer jfin,jfinmax,jfinmin
      integer jpartest
      integer kfin,indt,ind
      integer afin
      integer ms2,jj,ims1,is1,ims2,is2,mj1,mj2,ipol,nsp,lfin
      integer lxx,lx1,k1
      integer is,ispin
      integer mu1,mu1d,m1,m1d
      integer mu2,mu2d,m2,m2d
      integer jspin,ient
      real*8 res3j1,res3j2
      real*8, external :: cwig3j

      
	  hbmatrs=dcmplx(0.0d0,0.d0)
      jspin=1
      ient=2

      ji=2*abs(kinit)-1
      if (kinit.gt.0) then
         ainit=-1
      else
         ainit=1
      end if
    
      ind=0

      do lj=0,ljmax
         lt=2*lj
!     triangle
         jfinmax=min(lt+ji,jmax)
         jfinmin=max(abs(lt-ji),1)
         do jfin=jfinmin,jfinmax,2
!     use parity as an way to check/get afin  
            jpartest=ji+jfin+lt
            if (mod(jpartest,4).eq.0) then
               afin=-ainit
            else
               afin=ainit
            end if
            if (afin.gt.0) then 
               lfin=jfin-1
            else
               lfin=jfin+1
            end if
            
            indt=0
            if (lfin.le.2*lx1) then 
               ind=ind+1
               indt=indt+1
               if (ind.gt.kfinmax) then
                  write(6,*) "Unexpected error in bcoefjas"
                  stop
               end if
               indmat(ind)=indt
               jind(ind)=jfin

            end if
               
         end do                 ! jfin
      end do                    ! lj

      is = 0
      if (ispin.eq.1) is = nspx - 1

!
!     The combination matrix
!
!
      do k1=1,indmax
         l1=lind(k1)
         jfin=jind(k1)
         l1p=l1+1
         itolj=dcmplx(0.0d0,1.0d0)**ljind(k1)
         itolg=dcmplx(0.0d0,1.0d0)**l1
         do mu1=-mtot,mtot
            mu1d = mu1+mtot+1
            do mu2=-mtot,mtot
               mu2d = mu2+mtot+1
               do m1=-l1,l1
                  m1d = m1 + mtot+1
                  do ims1=1,nsp
                     is1=ims1-1
                     mj1=2*m1+(2*is1-1)
                     if (abs(mj1).le.jfin ) then
                        res3j1=cwig3j(jspin,jfin,2*l1,2*is1-1,-mj1,ient)
                        jj=l1+mj1-jspin
                        jj=jj/2
                        if (mod(jj,2).ne.0) then
                           res3j1=-res3j1
                        end if
                        res3j1=res3j1
                        do ims2=1,nsp
                           is2=ims2-1
                           m2=(mj1-(2*is2-1))/2
                           if (m2.le.l1) then 
                           m2d=m2+mtot+1
						   res3j2=cwig3j(jspin,jfin,2*l1,2*is2-1,-mj1,ient)
                           jj=l1+mj1-jspin
                           jj=jj/2
                           if (mod(jj,2).ne.0) then
                              res3j2=-res3j2
                           end if
                           res3j2=res3j2
                           
						   hbmatrs(mj1,ims2-1,mu2,mu1,k1)=hbmatrs(mj1,ims2-1,mu2,mu1,k1)+ &
                             res3j1*res3j2*exp(-coni*eta(0)*m1)*dri(l1p,mu1d,m1d,nsc+2)* &
                          exp(-coni*eta(nsc+2)*m2)*dri(l1p,m2d,mu2d,nleg)

			              end if
                        end do
                     end if    
                  end do
               end do
            end do
         end do
      end do  

      return
      end


