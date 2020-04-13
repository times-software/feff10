!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sigte3.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigte3 (iz1,iz2, sig2, alphat, thetad, reff, sig1,sig3)
!     single scattering only.

!     input: sig2
!     iz1, iz2 are iz at central atom and neighbor
!     alphat coeef of thermal expansion at high T
!     reff

!     output: sig1 sig3
      implicit double precision (a-h, o-z)
      real reff

!     con=hbar**2/kB*amu)*10**20   in ang**2 units
!     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
!     kB = 1.380 6581 d-23
      parameter (con = 48.508459393094)
      parameter (hbar = 1.054572666e-34)
      parameter (amu = 1.660540e-27)
      parameter (xkb = 1.3806581e-23)

      ami=atwtd(iz1)*amu
      amj=atwtd(iz2)*amu

!     reduced mass
      xmu = 1 / (1/ami + 1/amj)
!     Einstein frequency
      omega = (2 * xkb * thetad) / (3 * hbar)
      xks = xmu * omega**2
      xk3 = xks**2 * reff * alphat / (3 * xkb)
      sig02 = hbar * omega / xks
      sig1 = -3 * (xk3 / xks) * sig2
      sig3 = 2 - (4.0/3.0) * (sig02 / sig2)**2
      sig3 = sig3 * (sig1 * sig2)

      return
      end
