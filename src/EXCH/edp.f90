!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: edp.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!
!     this subroutine calculates the ' energy dependent
!     exchange-correlation potential' (or 'dirac- hara potential')
!     ref.: paper by s.h.chou, j.j.rehr, e.a.stern, e.r.davidson (1986)
!
!     inputs:    rs in a.u.
!                xk momentum in a.u.
!     outputs:   vr --- dirac potential (Hartrees)
!     written by j. mustre 8/31/87
!**********************************************************************

      subroutine edp (rs, xk, vr)
      use par
	  use constants
      implicit double precision (a-h, o-z)

      vr = 0.0d0
      if (rs .le. 100.0) then
!       p = sqrt (k^2 + kf^2) is the local momentum, and x = p / kf
!       Reference formula 23 in Role of Inelastic effects in EXAFS
!       by Rehr and Chou. EXAFS1 conference editted by Bianconi.
!       x is local momentum in units of fermi momentum

        xf = fa / rs
        x = xk / xf
        x = x + 1.0e-5
!       set to fermi level if below fermi level
        if (x .lt. 1.00001) x = 1.00001
        c = abs( (1+x) / (1-x) )
        c = log(c)
        vr = - (xf/pi) * (1 + c * (1-x**2) / (2*x))
      endif

      return
      end
