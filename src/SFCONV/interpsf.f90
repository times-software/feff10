!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: interpsf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine interpsf(npts,epts,wpts,spectf,cspec)
! Interpolates the spectral function calculated on a minimal grid
! to a uniform grid that can be handled by the convolution subroutine.
! input: npts - number of grid points the spectral function will
!               be interpolated on to.
!        epts - energy values of the minimal grid
!        wpts - energy values of the uniform grid
!        spectf - the spectral function on the minimal grid
!        cspec - the spectral function on the uniform grid
      implicit none
      integer npts,nsfpts,i,j
!      parameter (nsfpts=80)
      parameter (nsfpts=110)
      double precision spectf(8,nsfpts),cspec(npts),epts(nsfpts),       &
     &                 wpts(npts),wmin,wmax,dw,sfhi,sflo,delw
      double precision pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      common /convsf/ pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      wmin=epts(1)
      wmax=epts(nsfpts)
      dw=(wmax-wmin)/(npts-1)
      wpts(1)=wmin
      do i=2,npts
        wpts(i)=wmin+dw*(i-1)
      enddo
      cspec(1)=spectf(2,1)+spectf(5,1)-2.d0*spectf(4,1)
      cspec(npts)=spectf(2,nsfpts)+spectf(5,nsfpts)                     &
     &            -2.d0*spectf(4,nsfpts)
      do i=2,npts-1
        do j=2,nsfpts
          if (wpts(i).ge.epts(j-1).and.wpts(i).lt.epts(j)) then
            delw=wpts(i)-epts(j-1)
            sfhi=spectf(2,j)+spectf(5,j)-2.d0*spectf(4,j)
            sflo=spectf(2,j-1)+spectf(5,j-1)-2.d0*spectf(4,j-1)
            cspec(i)=sflo+(sfhi-sflo)*delw/(epts(j)-epts(j-1))
            goto 10
          endif
        enddo
 10     continue
      enddo
      return
      end
