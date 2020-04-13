!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: quinn.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine quinn (x, rs, wp, ef, ei)

      use constants
      implicit double precision (a-h, o-z)

!     input  x, rs, wp, ef
!     output ei

!***********************************************************************
!
!     quinn: calculates low energy gamma (approx. proportional to e**2)
!             formula taken from john j. quinn, phys. rev. 126,
!             1453 (1962); equation (7).
!             a cut-off is set up at quinn's cutoff + ef = ekc; it is a
!             rounded inverted step function (a fermi function)
!             theta = 1/( 1 + exp((e-ekc)/gam)) )
!             where the rounding factor gam is set to be about 0.3 ekc.
!     modified by j. rehr (oct 1991) based on coding of r. albers
!     subroutines quinn.f and quinnc.f
!
!     variables:
!        x  = p/pf
!        rs = ws density parameter
!        ei = imaginary self energy
!        pfqryd = quinn's prefactor in atomic-rydberg units
!        wkc = quinn's plasmon threshold
!
!***********************************************************************

      parameter (alphaq = 1/ fa)

!     calculate quinn prefactor in atomin Hartree units
      pisqrt = sqrt(pi)
      pfq = pisqrt / (32 * (alphaq*rs)**1.5)
      temp1 = atan (sqrt (pi / (alphaq*rs)))
      temp2 = sqrt(alphaq*rs/pi) / (1 + alphaq*rs/pi)
      pfq = pfq * (temp1 + temp2)

!     calculate quinn cutoff
!     wkc = quinn's plasmon threshold
!     wkc is cut-off of quinn, pr126, 1453, 1962, eq. (11)
!     in formulae below wp=omegap/ef
      wkc = (sqrt(1+wp) - 1)**2
      wkc = (1 + (6./5.) * wkc / wp**2) * wp * ef

!     we add fermi energy to get correct energy for
!     plasma excitations to turn on
      ekc = wkc + ef

!     calculate gamma
!     gamryd = 2 * (pfqryd/x) * (x**2-1)**2
      gam = (pfq/x) * (x**2-1)**2

!     put in fermi function cutoff
      eabs = ef * x**2
      arg = (eabs-ekc) / (0.3*ekc)
      f = 0
      if (arg .lt. 80)  f = 1 / (1 + exp(arg))

      ei = -gam * f / 2

      return
      end
