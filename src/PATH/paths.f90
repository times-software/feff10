!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: paths.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2012/11/15 01:16:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine paths (ckspc, fbetac, xlamc, pcritk, pcrith, critpw,   &
     &                  nncrit, rmax, nlegxx, rfms,                     &
     &                  nat, ratdp, iphat, ibounc)

!     finds multiple scattering paths
!     This is single precision, units are Angstroms.  BE CAREFUL!

!     pcrith is cut-off fraction used when building paths
!            (path criterion for heap)
!     pcritk is cut-off fraction used on output
!            (path criterion for keeping)

      use dimsmod, only: npx => npx_path, natx, nx, nphx=>nphu, npatx
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension xlamc(necrit)

!     This common in pathsd, mpprm
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)
      double precision ratdp(3,natx)
      integer  iphat(natx), ibounc(natx)

      dimension m(-1:natx,0:natx)
      dimension mindex(natx+1)
!     Used for packed integers
      dimension iout(3)

!     ok true if all paths to rmax found.  If heap full, npx exceeded,
!     etc., last general shell may be incomplete, set ok=.false.
      logical ok
!     is label nfound, etc, written yet?
      logical wlabel

!     Heap data structure:
!     index is the pointer to the element of the data structure.
!     Each element contains
!        r        total path length
!                 Note that r is sorted along with index -- this keeps
!                 the heap maintenance routines fast.
!        mi, mj   m matrix elements used to place last atom in this path
!        npat     number of atoms in this path
!        ipat(npatx) indices of atoms in this path
!     next is the index of the next data structure element available.
!     If an element is freed, npat is the index of the free element
!     to use after using current "next" element.

!     nx is max number in heap
!      integer    nx
!     parameter (nx = 10000)
!KJ      parameter (nx = 60000)
!KJ 2014: nx moved to dimsmod

!     r also used in making m matrix, must have nx >= natx+1
      integer   index(nx), np, n, ip, i, iat0, idum
	  
!      integer npx	  
!      parameter (npx = 100000)
!KJ      parameter (npx = 4000000)
!KJ 2014 npx moved to dimsmod

      dimension r(nx), mi(nx), mj(nx)
      dimension npat(nx)
      dimension ipat (npatx,nx)
!     Keep this path on output
      logical keep1(nx), kp1tmp
!     to remmember atoms outside rfms
      dimension iclus(0:natx) 

!     Used with ipack, so need ipat(8)
      dimension ipat0(8)

!     paths are typically about 10 or 20 Ang
      parameter (big = 1.0e3)

      character*80  title

!     Returned from criterion checker, false if path fails criterion
      logical keep
      character*512 slog

      external sdist

!     transform to single precision for pathfinder
!     also count in pathfinder starts with 0, not with 1
      do 15 iat = 1, nat
        j = iat-1
        do 10 i = 1, 3
          rat(i, j) = real (ratdp(i,iat))
 10     continue
        ipot(j) = iphat(iat)
        i1b(j) = ibounc(iat)
 15   continue
      nat = nat - 1

!     nlegxx is max number of legs user wants to consider.
!     nlegs = npat+1, so set npatxx = min (npatx, nlegxx-1)
      npatxx = min (npatx, nlegxx-1)
!     Input rmax is one-way distances
      rmax = rmax*2
!     ratx is distance to most distant atom, used to check rmax
      ratx = 0
      iat0 = -1
      i1b(0) = 0

!     find index for the central atom
      do 20 iat = 0, nat
         if (ipot(iat).eq.0 .and. iat0.lt.0) iat0=iat
   20 continue

!     iclus = 0 for atoms inside rfms cluster, 1 for atoms outside
      do 21 iat = 0,nat
        rtmp = sdist(rat(1,iat),rat(1,iat0))
        if (rtmp.gt.rfms) iclus(iat)=1
        iclus(iat)=0
        if (rtmp.gt.rfms) iclus(iat)=1
   21 continue

      if (iat0.ne. 0) then
!       permute atoms 0 and iat0
!       do not need to permute i1b, since all of them are 1 in
!       this case, except i1b(0) = 0, which we want to keep.
        do 25 j=1,3
          temp = rat(j,0)
          rat(j,0) = rat(j,iat0)
          rat(j,iat0) = temp
   25   continue
        idum = ipot(0)
        ipot(0) = ipot(iat0)
        ipot(iat0) = idum
        idum = iclus(0)
        iclus(0) = iclus(iat0)
        iclus(iat0) = idum
      endif

!KJ 12-2011 : This warning is completely useless.  ratx is set to 0 above, so the warning is ALWAYS printed.  It means nothing.
!KJ!     Warn user if rmax > dist to most distant atom
!KJ!     1.01 to avoid roundoff error, matches rdinp where rmax default set
!KJ      if (rmax/2 .gt. 1.01 * ratx)  then
!KJ         call wlog('   WARNING:  rmax > distance to most distant atom.')
!KJ         call wlog('             Some paths may be missing.')
!KJ         write(slog,22) rmax/2, ratx
!KJ         call wlog(slog)
!KJ   22    format('             rmax, ratx ', 1p, 2e13.5)
!KJ      endif

!     Count number of 1st bounce atoms (at least 1 required).
      n1b = 0
      do 30  i = 1, nat
         if (i1b(i) .gt. 0)  n1b = n1b + 1
   30 continue
      if (n1b .lt. 1)  call par_stop('At least one 1st bounce atom required.')

      if (rmax .ge. big)  call par_stop('RPATH rmax  error:  please reduce rmax to a reasonable value.')

!     Make title for this run, include carriage control because head (read above) includes carriage control.
      write(title,32)  rmax/2, pcritk, pcrith, critpw
   32 format('PATH  Rmax=', f6.3, ',  Keep_limit=', f5.2,', Heap_limit', f5.2,'  Pwcrit=', f5.2, '%')

      write(slog,34) rmax/2, pcritk, pcrith
      if(rmax.gt.0.1d0) call wlog(slog)
   34 format ('    Rmax', f8.4,'  keep and heap limits', 2f12.7)

      if(rmax.gt.0.1d0) call wlog('    Preparing neighbor table')
      
   36 format (1x, a)
!     prepare table telling distance from atom i to atom j and then
!     back to central atom
!     First bounce is m(-1,...), m(0,...) is bounces from central
!     atom that are not first bounces.
      do 60  i = -1, nat
         ir = i
         if (i .eq. -1)  ir = 0
         do 40  j = 0, nat
!           r begins with element 1 so sort routine later will work
            r(j+1) = sdist (rat(1,ir), rat(1,j))
            r(j+1) = r(j+1) + sdist (rat(1,j), rat(1,0))
!           we don't need m(i,i), since this will be = shortest
!           of the r(j), so just set it to something very big,
!           it will sort to the end of this row and it won't
!           bother us
            if (j .eq. ir)  r(j+1) = big
!           If we're doing first bounce, use only the allowed first
!           bounce paths.
            if (i .eq. -1)  then
               if (i1b(j) .le. 0)  r(j+1) = big
            endif
   40    continue

!        prepare row i of m table
!        m is a distance table ordered such that distance from
!               i to m(i,0) to 0 <
!               i to m(i,1) to 0 <
!               i    m(i,2)    0 <
!               :    :    :
!               i    m(i,nat)  0
!
!        That is, m(i,0) is index of atom that gives shortest path,
!                 m(i,1)                        next shortest path, etc.
!        Note that m(0,0) is shortest single bounce path.

!        Again, r and mindex go from 1 to nat+1, m goes from 0 to nat
         call sortir (nat+1, mindex, r)
         do 50  j = 0, nat
            m(i,j) = mindex(j+1)-1
   50    continue
   60 continue

!     label for nfound, heap size, etc written?
      wlabel = .false.
!     initialize heap data space "next" pointers
      do 70  i = 1, nx-1
         npat(i) = i+1
   70 continue
      npat(nx) = -1
!     initial condition:  make the first path
!     n    number in heap
!     nna  number skipped counter
!     nhx  number used in heap max, a counter
      n = 1
      nna = 0
      nhx = n
      nwrote = 0
      index(n) = 1
      ip = index(n)
      next = 2
      mi(ip) = -1
      mj(ip) = 0
      npat(ip) = 1
      ipat(npat(ip),1) = m(mi(ip),mj(ip))

!     Someday change keep and keep1 to lkeep and lheap to match
!     ccrit variable names.
!     Initialize keep criterion
      xcalcx = -1
      call ccrit (npat(ip), ipat(1,ip), ckspc, fbetac, xlamc, rmax, pcrith, pcritk, nncrit, ipot,            &
          r(n), keep, keep1(ip), xcalcx, iclus)

      open (file='paths.bin', unit=3, access='sequential', form='unformatted', status='unknown', iostat=ios)
      call chopen (ios, 'paths.bin', 'paths')
!     These strings are all char*80 and include carriage control
!     temporary fix for zero title lines: fix later
      nhead = 1
      write(3) nhead
      write(3) title
      write(3)  nat
      do i = 0, nat
         write(3) (rat(j,i),j=1,3), ipot(i), i1b(i)
      enddo

!     r is the heap, index is the pointer to the rest of the data
!     np is the number of paths found and saved
      np = 0
!     nbx  mpat max (Number of Bounces maX)
      nbx = 0

!     done if path at top of heap is longer than longest path we're interested in.
!     done if max number of paths we want have been found.
!     begin 'while not done' loop
      ok = .false.
  800 continue
         if (r(1) .gt. rmax  .or.  np .ge. npx .or. n.le.0)  then
!           n=0 means heap is empty
            if (n.le.0)  ok=.true.
			if (rmax .lt. 1.0) ok=.true.  !KJ 2-2012 assuming this means RPATH so small that we don't want any paths ...
            goto 2000
         endif

!        save element at top of heap in arrays labeled 0
!        dump to unit 3 (unformatted)
         ip = index(1)
         npat0 = npat(ip)
         do i = 1, npat0
            ipat0(i) = ipat(i,ip)
         enddo
         r0 = r(1)

!        Don't write out path if last atom is central atom, or if it doesn't meet pcritk
         if (ipat0(npat0).ne.0 .and. keep1(ip))  then
            np = np+1
!           pack integers
            call ipack (iout, npat0, ipat0)
            write(3)  r0, iout
            nwrote = nwrote+1
!           write status report to screen
            if (mod(np,1000) .eq. 0)  then
               if (.not. wlabel)  then
                  call wlog('    nfound  heapsize  maxheap  maxscatt   reff')
                  wlabel = .true.
               endif
               write(slog,132) np, n, nhx, nbx, r0/2
               call wlog(slog)
  132          format (4x, i10, i10, i10, i6, f11.4)
            endif
         endif

         if (np .ge. npx)  then
            write(slog,134) np
            call wlog(slog)
  134       format(i15, ' paths found.  (np > npx)')
            goto 2000
         endif

!        Make new path by replacing last atom in path from top of heap,
!        put this path on top of heap and buble it down.  If row is
!        finished, or new path is too long, don't add it, instead
!        move last path in heap to the top.
!        If working on row mi=-1 (first bounce atoms), don't
!        use them if not allowed 1st bounce atoms.
         mj(ip) = mj(ip) + 1
         if (mi(ip).eq.-1  .and.  i1b(m(mi(ip),mj(ip))).le.0)  then
!           not allowed first bounce atom
            r(1) = big
            keep = .false.
         elseif (mj(ip) .ge. nat)  then
!           we've finished a row of m matrix
            r(1) = big
            keep = .false.
         else
!           new path has same indices, etc.  Only need to replace last atom.
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc, fbetac, xlamc, rmax, pcrith, pcritk, nncrit,    &
                        ipot, r(1), keep, keep1(ip), xcalcx, iclus)
         endif

!        If r is bigger than rmax or keep=false, remove element from
!        heap by taking the last element in the heap and moving it to
!        the top.  Then bubble it down.  When removing an element
!        from the heap, be sure to save the newly freed up index.
!        r(1) and index(1) are new path, set above
         if (r(1).gt.rmax .and. keep)  then
            call wlog(' odd case rmax...')
         endif
         if (r(1).gt.rmax .or. .not.keep)  then
            index(1) = index(n)
            r(1) = r(n)
!           use npat as pointer to next free location
            npat(ip) = next
            next = ip
            n = n-1
!           nna is Number Not Added to heap
            nna = nna + 1
!           Maybe heap may be empty here, but that's alright
         endif
         if (npat(index(1)).gt.nbx .and. n.gt.0)  nbx = npat(index(1))

!        If heap is empty, don't call hdown.
         if (n.gt.0)  call hdown (r, index, n)

!        and make a new path by adding an atom onto the end of the path
!        we saved, put this at the end of the heap and bubble it up.
!        Do this only if it won't be too many bounces.
         if (npat0+1 .le. npatxx)  then
            ip = next
            if (ip .lt. 0)  then
!              call wlog('   Heap full')
               goto 2000
            endif
            next0 = npat(ip)
            do i = 1, npat0
               ipat(i,ip) = ipat0(i)
           enddo
            mi(ip) = ipat0(npat0)
            mj(ip) = 0
            npat(ip) = npat0+1
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc, fbetac, xlamc, rmax, pcrith, pcritk, nncrit,    &
                        ipot, rtmp, keep, kp1tmp, xcalcx, iclus)
            if (rtmp .gt. rmax  .and.  keep)  then
               call wlog(' odd case rmax and tmp...')
            endif
            if (rtmp .gt. rmax  .or.  .not.keep)  then
               npat(ip) = next0
               nna = nna+1
            else
!              add it to the heap
               next = next0
               n = n+1
               if (n .gt. nhx)  nhx = n
               index(n) = ip
               r(n) = rtmp
               keep1(ip) = kp1tmp
               if (npat(index(n)) .gt. nbx)  nbx = npat(index(n))
               call hup (r, index, n)
            endif
         endif

      goto 800 !add another path to the list
 2000 continue
!     end of 'while not done' loop


      if (.not. ok  .and. np.gt.0 )  then  !KJ added np condition 11-2012
         call wlog('   Internal path finder limit exceeded -- path list may be incomplete.')
      endif
      close (unit=3)
      write(slog,2010) np, nhx, nbx
 2010 format ('    Paths found', i9, 3x,'(maxheap, maxscatt', i8, i4, ')')
      if(rmax.gt.0.2d0) call wlog(slog)

!     restore the value of rmax
      rmax = rmax/2

      end
