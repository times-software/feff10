!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mcrith.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mcrith (npat, ipat, ri, indbet,                        &
     &                   ipot, nncrit, fbetac, ckspc, xheap)

      use dimsmod, only: npatx, nphx=>nphu, natx
	  use constants
      dimension ipat(npatx)
      dimension ri(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)

!     Decide if we want the path added to the heap.

      if (ipat(npat) .eq. 0 .or. npat.le.2)  then
!        Partial path is used for xheap, not defined for ss and
!        triangles.  Special case: central atom added to end of path 
!        necessary for complete tree, but not a real path, again,
!        xheap not defined.  Return -1 as not-defined flag.
         xheap = -1
      else
!        Calculate xheap and see if we want to add path to heap.
!        Factor for comparison is sum over nncrit of
!        f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
!        Compare this to sum(1/p), multiply by 100 so we can think 
!        in percent.  Allow for degeneracy when setting crit.
         xheap = 0
         spinv = 0
         do 340  icrit = 1, nncrit
            x = ckspc(icrit) ** (-(npat-1)) * ri(npat-1)
            do 320  i = 1, npat-2
               ipot0 = ipot(ipat(i))
               x = x * fbetac(indbet(i),ipot0,icrit) / ri(i)
  320       continue
            spinv = spinv + 1/ckspc(icrit)
            xheap = xheap + x
  340    continue
         xheap = 100 * xheap / spinv

!        Factor for comparison is sum over nncrit of
!        New xheap:
!        Full chi is
! f(beta1)*f(beta2)*..*f(beta npat)cos(beta0)/(rho1*rho2*..*rho nleg).
! Some of this stuff may change when the path is modified --
! we can't use rho nleg or nleg-1, beta0, beta(npat) or beta(npat-1).
! We DO want to normalize wrt first ss path, f(pi)/(rho nn)**2.
!
! So save f(pi)/(rho nn)**2, 
! calculate 
! f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
! divide nn ss term by stuff we left out -- beta(npat), beta(npat-1),
! cos(beta0), rho nleg, rho nleg-1.
!
! Sum this over nncrit and try it out.
!
!        Sum over nncrit of
!        1/(rho1+rho2+..+rho npat-1).
!        reff = 0
!        do 350  i = 1, npat-1
!           reff = reff + ri(i)
! 350    continue
!        xss = 0
!        do 360  icrit = 1, nncrit
!           rho = ckspc(icrit) * reff
!           xss = xss + 1/rho
! 360    continue
!        xheap = 100 * xheap / xss
      endif

      return
      end
