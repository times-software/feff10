!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mcritk.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mcritk (npat, ipat, ri, beta, indbet,                  &
     &      ipot, nncrit, fbetac, xlamc, ckspc, xout, xcalcx)

      use dimsmod, only: npatx, nphx=>nphu, natx
	  use constants
	  dimension ipat(npatx)
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension xlamc(necrit)

!c    xcalcx is max xcalc encountered so far.  Set to -1 to reset it --
!c    otherwise it gets passed in and out as mcritk gets called.
!     calculation of xcalcx changed by ala. It is calculated only
!     a first call, i.e. for the NN SS path, and is not recalculated

!     We may want path in heap so that other paths built from this
!     path will be considered, but do not want this path to be
!     written out for itself.  Decide that now and save the flag
!     in the heap, so we won't have to re-calculate the mpprm
!     path parameters later.

!     Do not want it for output if last atom is central atom,
!     use xout = -1 as flag for undefined, don't keep it.
      if (ipat(npat) .eq. 0)  then
         xout = -1
         return
      endif

!     Make xout, output inportance factor.  This is sum over p of
!     (product of f(beta)/rho for the scatterers) * 
!                                 (cos(beta0)/rho(npat+1).
!     Compare this to xoutx, max xout encountered so far.
!     Use mean free path factor, exp(-rtot/xlam)
!     Multiply by 100 so we can think in percent.

      xcalc = 0
      rtot = 0
      do 410  i = 1, npat+1
         rtot = rtot + ri(i)
  410 continue
      do 460  icrit = 1, nncrit
         rho = ri(npat+1) * ckspc(icrit)
!        when beta(0)=90 degrees, get zero, so fudge with cos=.2
         x = max (abs(beta(npat+1)), 0.3) / rho
         do 420  iat = 1, npat
            rho = ri(iat) * ckspc(icrit)
            ipot0 = ipot(ipat(iat))
            x = x * fbetac(indbet(iat),ipot0,icrit) / rho
  420    continue
         x = x * exp (-rtot/xlamc(icrit))
         xcalc = xcalc + x
  460 continue
      if (xcalcx.le.0)  xcalcx = xcalc
      xout = 100 * xcalc / xcalcx
      return
      end
