!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mkrmu.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mkrmu(xmu,xmu0,rmu,wpts,npts)
! This subroutine does a Cramers-Kronig transform on the array xmu
! and returns an array rmu which is the real part of the analytic
! function whose imaginary part is xmu.  This is needed to get the proper
! phase shift for a convolution with a real function.
      implicit none
      integer npts,i,j
      double precision xmu(npts),xmu0(npts),rmu(npts),wpts(npts)
      double precision dw,pi
      pi=dacos(-1.d0)
      do j=1,npts
        rmu(j)=0.d0
        do i=1,npts
          if (i.eq.1) then
            dw=wpts(2)-wpts(1)
          elseif (i.eq.npts) then
            dw=wpts(npts)-wpts(npts-1)
          else
            dw=(wpts(i+1)-wpts(i-1))/2.d0
          endif
          if (i.ne.j) then
            rmu(j)=rmu(j)+dw*(xmu(i)-xmu0(i))/(wpts(i)-wpts(j))
          endif
        enddo
        rmu(j)=rmu(j)/pi
      enddo
      rmu(20)=(rmu(20)+rmu(21))/2.d0
      rmu(21)=rmu(20)
 500  format(1x,5(e12.5,1x))
      return
      end
