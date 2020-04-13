!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: outcrt.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/12/11 23:20:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine outcrt (npat, ipat, ckspc,                             &
     &    nncrit, fbetac, xlamc, ne, ik0, cksp,                         &
     &    fbeta, xlam, ipot,                                            &
     &    xport, xheap, xheapr,                                         &
     &    xout, xcalcx)

!     This make pw importance factor for pathsd, also recalculates
!     pathfinder criteria for output.  Pathfinder recalculation
!     is hacked from ccrit, so be sure to update this if ccrit
!     is changed.

      use dimsmod, only: npatx, natx, nphx=>nphu, nex
	  use constants
      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      dimension xlamc(necrit), xlam(nex)

!     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension xporti(nex)
      parameter (eps = 1.0e-6)

!     Space for variables for time reversed path (used in xheapr
!     calculation below)
      dimension ipat0(npatx)
      dimension ri0(npatx+1), indbe0(npatx+1)

!     mrb is 'efficient' way to get only ri and beta
!     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

!     Make index into fbeta array (this is nearest cos(beta) grid point,
!     code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = tmp / 0.025
         del = tmp - n*0.025
         if (del .gt. 0.0125)  n = n+1
         if (beta(i) .lt. 0)  n = -n
         indbet(i) = n
  290 continue

!     Make pw importance factor by integrating over all points
!     above the edge
!     Path importance factor is integral d|p| of
!        (product of f(beta)/rho for the scatterers) * cos(beta0)/rho0
!     Include mean free path factor, exp(-rtot/xlam)
      rtot = 0
      do 510  i = 1, npat+1
         rtot = rtot + ri(i)
  510 continue

      do 560  ie = ik0, ne
         rho = ri(npat+1) * cksp(ie)
         crit = max (abs(beta(npat+1)), 0.3) / rho
         do 520  iat = 1, npat
            rho = ri(iat) * cksp(ie)
            ipot0 = ipot(ipat(iat))
            crit = crit * fbeta(indbet(iat),ipot0,ie) / rho
  520    continue
         crit = crit * exp (-rtot/xlam(ie))
         xporti(ie) =  abs(crit)
  560 continue


!     integrate from ik0 to ne
      nmax = ne - ik0 + 1
      call strap (cksp(ik0), xporti(ik0), nmax, xport)

!     Stuff for  output.
!     Heap crit thing (see ccrit and mcrith for comments)
!     If a path got time reversed, its xheap may be smaller than
!     it was before it got time-reversed.  So calculate it both
!     ways.
!     xheap for path, xheapr for time-reversed path

      xheap  = -1
      xheapr = -1
      call mcrith (npat, ipat, ri, indbet,                              &
     &             ipot, nncrit, fbetac, ckspc, xheap)

!     Prepare arrays for time reversed path and make xheapr
!     See timrev.f for details on indexing here.

      nleg = npat+1
!     ri
      do 200  i = 1, nleg
         ri0(i) = ri(nleg+1-i)
  200 continue
!     indbet  and ipat
      indbe0(nleg) = indbet(nleg)
      do 210  i = 1, nleg-1
         indbe0(i) = indbet(nleg-i)
         ipat0(i) = ipat(nleg-i)
  210 continue

      call mcrith (npat, ipat0, ri0, indbe0,                            &
     &             ipot, nncrit, fbetac, ckspc, xheapr)

!     Keep crit thing (see mcritk for comments)
      call mcritk (npat, ipat, ri, beta, indbet,                        &
     &             ipot, nncrit, fbetac, xlamc, ckspc, xout, xcalcx)

      return
      end
