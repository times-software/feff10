!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ccrit.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ccrit (npat, ipat, ckspc,                              &
     &    fbetac, xlamc, rmax, pcrith, pcritk, nncrit, ipot,            &
     &    rpath, lheap, lkeep, xcalcx, iclus)

!     lheap to add to heap, lkeep if keep path at output.
!     NB, if lheap is false, lkeep is not used (since path
!     won't be in the heap).
      use dimsmod, only: natx, nphx=>nphu, npatx
	  use constants
      logical lheap, lkeep
      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension xlamc(necrit), iclus(0:natx)

!     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)


!     mrb is efficient way to get only ri and beta
!     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

      rpath = 0
      do 300  i = 1, npat+1
         rpath = rpath + ri(i)
  300 continue

!     If we can decide only on rpath, do it here...
      if (rpath .gt. rmax)  then
         lheap = .false.
         lkeep = .false.
         return
      endif

!     If last atom central atom, do put in heap, don't use it
!     as an actual path at output
      if (ipat(npat).eq.0)  then
         lheap = .true.
         lkeep = .false.
         return
      endif

!     Make index into fbetac array (this is nearest cos(beta) grid 
!     point, code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = tmp / 0.025
         del = tmp - n*0.025
         if (del .gt. 0.0125)  n = n+1
         if (beta(i) .lt. 0)  n = -n
         indbet(i) = n
  290 continue

!     Decide if we want the path added to the heap if necessary.
!     (Not necessary if no pcrith in use.)
      if (pcrith .gt. 0)  then

         call mcrith (npat, ipat, ri, indbet,                           &
     &                ipot, nncrit, fbetac, ckspc, xheap)

!        xheap = -1 if not defined for this path (too few legs, etc.)
         if (xheap .ge. 0  .and.  xheap .lt. pcrith)  then
!           Do not want path in heap
            lheap = .false.
            lkeep = .false.
            return
         endif
      endif
!     Keep this path in the heap
      lheap = .true.

!     We may want path in heap so that other paths built from this
!     path will be considered, but do not want this path to be
!     written out for itself.  Decide that now and save the flag
!     in the heap, so we won't have to re-calculate the mpprm
!     path parameters later.

!     Skip calc if pcritk < 0
      if (pcritk .le. 0)  then
         lkeep = .true.
         goto 999
      endif

!     Make xout, output inportance factor.
      call mcritk (npat, ipat, ri, beta, indbet,                        &
     &             ipot, nncrit, fbetac, xlamc, ckspc, xout, xcalcx)

!     See if path wanted for output
!     Do not want it if last atom is central atom (xout = -1) or
!     if xout is too small
      lkeep = .false.
      if (xout .ge. pcritk)  lkeep = .true.

!     If path is entirely inside a cluster do not keep it
  999 nclus=0
      do 700 i=1,npat
  700 nclus=nclus+iclus(ipat(i))
      if (nclus.eq.0) lkeep = .false.

      return
      end
