!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getholeorb0.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getholeorb0(dx,dxnew,ihole,jnew,iz,xion,iunf,dgc,dpc,dgcx0,dpcx0)
      use constants
      use dimsmod, only: nphx=>nphu, nrptx
	  implicit none

!KJ      dimension dgc(251,30,0:nphx), dpc(251,30,0:nphx)
      double precision dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1) !KJ 7-09 definitions have changed ...
      double precision dgcx0(nrptx), dpcx0(nrptx)
      double precision gc0(nrptx), pc0(nrptx)
      double precision xion(0:nphx)
      integer iz(0:nphx)
      integer nqn(30), nk(30), iorb(-4:3)
	  double precision xnel(30), xnval(30), xmag(30)
      double precision xorg(nrptx), xnew(nrptx)
      double precision xxx(nrptx), rrr(nrptx)
      real*8, parameter :: xx00 = 8.80d0
	  integer ihole,iholep,iunf,norb,norbco
	  integer i,imax,j,jmax,jnew
	  double precision delta,dx,dxnew
 
!     first have to get iholep i.e. index of the state ihole 
!     (apparently usually == ihole but not always)
!
      call getorb (iz(0), ihole, xion(0), iunf, norb, norbco, iorb, iholep, nqn, nk, xnel, xnval, xmag, 0) !KJ 2-2011 added iph=0

!     copy the wf
      do i=1,251
         gc0(i)=dgc(i,iholep,0)
         pc0(i)=dpc(i,iholep,0)

      end do

!     now interpolate
      imax = 0
      do 100  i = 251, 1, -1
         if ( abs(gc0(i)) .ge. 1.0d-11 .or. abs(pc0(i)) .ge. 1.0d-11 )  then
            imax = i
            goto 16
         endif
  100 continue
      call wlog(' Should never see this line from sub fixdsp')
   16 continue
!     jmax is the first point where both dpc and dgc are zero in the original grid
      jmax = imax + 1
      if (jmax.gt.251) jmax = 251

      delta = dx
      do j=1,251
         xxx(j) = -xx00 + (j-1)*delta
         rrr(j) = exp (-xx00 + (j-1)*delta)
      end do
      do 10  j = 1, jmax
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      do j=1,nrptx
         xxx(j) = -xx00 + (j-1)*delta
         rrr(j) = exp (-xx00 + (j-1)*delta)
      end do

      do 20  j = 1, jnew
         xnew(j) = xxx(j)
 20   continue
      
!     interpolate to new grid using x, only inside of rmax
      do 30  j = 1, jnew
         call terp (xorg, gc0,  jmax, 3, xnew(j), dgcx0(j))
         call terp (xorg, pc0,  jmax, 3, xnew(j), dpcx0(j))
   30 continue

!     and zero the arrays past rmax
      dgcx0(jnew+1:nrptx)=0.d0
	  dpcx0(jnew+1:nrptx)=0.d0

      return 
      end
