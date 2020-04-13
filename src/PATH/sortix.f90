!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sortix.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sortir (n, index, r)

!     SORT by rearranges Indices, keys are Real numbers
!     Heap sort, following algorithm in Knuth using r as key
!     Knuth, The Art of Computer Programming,
!     Vol 3 / Sorting and Searching, pp 146-7
!     Array r is not modified, instead array index is returned
!     ordered so that r(index(1)) is smallest, etc.
!     rr is temporary r storage (Knuth's R), irr is index of stored r

      dimension r(n), index(n)

!     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
!     only 1 element is already sorted
      if (n .eq. 1)  return

!     H1: initialize
      l = n/2 + 1
      ir = n

!     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         rr = r(irr)
      else
         irr = index(ir)
         rr = r(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

!     H3: Prepare for sift-up
      j = l

!     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

!     H5: Find larger son of i
      if (r(index(j)) .lt. r(index(j+1)))  j = j+1

!     H6: Son larger than rr?
   60 continue
      if (rr .ge. r(index(j)))  goto 80

!     H7: Move son up
      index(i) = index(j)
      goto 40

!     H8: Store rr in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
      subroutine sortii (n, index, k)

!     SORT by rearranges Indices, keys are Integers
!     Heap sort, following algorithm in Knuth using r as key
!     Knuth, The Art of Computer Programming,
!     Vol 3 / Sorting and Searching, pp 146-7
!     Array r is not modified, instead array index is returned
!     ordered so that r(index(1)) is smallest, etc.
!     rr is temporary r storage (Knuth's R), irr is index of stored r

      dimension k(n)
      dimension index(n)

!     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
!     only 1 element is already sorted
      if (n .eq. 1)  return

!     H1: initialize
      l = n/2 + 1
      ir = n

!     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         kk = k(irr)
      else
         irr = index(ir)
         kk = k(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

!     H3: Prepare for sift-up
      j = l

!     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

!     H5: Find larger son of i
      if (k(index(j)) .lt. k(index(j+1)))  j = j+1

!     H6: Son larger than kk?
   60 continue
      if (kk .ge. k(index(j)))  goto 80

!     H7: Move son up
      index(i) = index(j)
      goto 40

!     H8: Store kk in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
      subroutine sortid (n, index, r)

!     SORT by rearranges Indices, keys are Double precision numbers
!     Heap sort, following algorithm in Knuth using r as key
!     Knuth, The Art of Computer Programming,
!     Vol 3 / Sorting and Searching, pp 146-7
!     Array r is not modified, instead array index is returned
!     ordered so that r(index(1)) is smallest, etc.
!     rr is temporary r storage (Knuth's R), irr is index of stored r

      implicit double precision (a-h, o-z)
      dimension r(n), index(n)

!     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
!     only 1 element is already sorted
      if (n .eq. 1)  return

!     H1: initialize
      l = n/2 + 1
      ir = n

!     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         rr = r(irr)
      else
         irr = index(ir)
         rr = r(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

!     H3: Prepare for sift-up
      j = l

!     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

!     H5: Find larger son of i
      if (r(index(j)) .lt. r(index(j+1)))  j = j+1

!     H6: Son larger than rr?
   60 continue
      if (rr .ge. r(index(j)))  goto 80

!     H7: Move son up
      index(i) = index(j)
      goto 40

!     H8: Store rr in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
