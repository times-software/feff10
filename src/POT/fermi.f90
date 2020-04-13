!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fermi.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermi (rhoint, vint, xmu, rs, xf)

      use constants
      IMPLICIT NONE
      !implicit double precision (a-h, o-z)
      
      DOUBLE PRECISION rhoint, vint, xmu, rs, xf
      DOUBLE PRECISION den
!     calculate fermi level of the system (mu) according to formula
!     mu=vcoulomb(interstitial)+vxc(interstitial)+kf(interstitial)^2
!     formula  2.13 in lee and beni, phys. rev. b15,2862(1977)

!     note that vint includes both coulomb and ground state
!     exchange-correlation potentials

!     den is the interstitial density
!     rs is the density parameter
!     xf is the interstital fermi momentum
!     xmu is the fermi level in hartrees

      den = rhoint / (4*pi)
      !rs = (3 / (4*pi*den)) ** third
      rs = (3 / rhoint) ** third
      xf = fa / rs
      xmu = vint + xf**2 / 2

      return
      end
