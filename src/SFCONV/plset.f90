!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: plset.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
! Sets calculation parameters for a given pole in epsilon^{-1}
!        ipl - pole to set parameters for
!        nplmax - array dimensioning (maximum number of poles)
!        plengy - energy of poles (array)
!        plwt - weight of poles (array)
!        ompl - energy of selected pole
!        wt - weight of selected pole
!        adisp - dispersion parameter for selected pole
!        conc - electron density that would give plasma frequency
!               equal to pole energy
!        rs - Wigner Seitz radius that would give plasma frequency
!               equal to pole energy
!        ef - Fermi energy that would give plasma frequency
!               equal to pole energy
!        qf - Fermi momentum that would give plasma frequency
!               equal to pole energy
      implicit none
      integer ipl,nplmax
      double precision plengy(nplmax),plwt(nplmax),plbrd(nplmax)
      double precision ompl,wt,brd
!      double precision pi,adisp,conc,rs,ef,qf
      ompl=plengy(ipl)
      wt=plwt(ipl)
      brd=plbrd(ipl)
!      adisp=(0.75d0*pi/ompl**2)**(2.d0/3.d0)/3.d0
!      rs=(3.d0/ompl**2)**(1.d0/3.d0)
!      conc=3.d0/(4.d0*pi*(rs**3))
!      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
!      ef=qf*qf/2.d0
      return
      end

      subroutine ppset(rs,pi,qf,ef,omp)
! Sets electron gas parameters for a given Wigner-Seitz radius rs.
!        rs - Wigner Seitz radius that would give plasma frequency
!               equal to pole energy
!        omp - plasma frequency
!        ef - Fermi energy that would give plasma frequency
!               equal to pole energy
!        qf - Fermi momentum that would give plasma frequency
!               equal to pole energy
      implicit none
      double precision rs,qf,ef,conc,omp,pi
      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
      ef=qf*qf/2.d0
      conc=3.d0/(4.d0*pi*(rs**3))
      omp=dsqrt(4.d0*pi*conc)
      return
      end

