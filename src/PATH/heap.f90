!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: heap.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     These heap routines maintain a heap (array h) and an index
!     array (array ih) used to keep other data associated with the heap
!     elements.

      subroutine hup (h, ih, n)
!     heap is in order except for last element, which is new and must
!     be bubbled through to its proper location
!     new element is at i, j = index of parent
      integer  n,i,j
      integer  ih(n)
      dimension h(n)

      i = n

   10 j = i/2
!     if no parent, we're at the top of the heap, and done
      if (j .eq. 0)  return
      if (h(i) .lt. h(j))  then
         call swap (h(i), h(j))
         call iswap (ih(i), ih(j))
         i = j
         goto 10
      endif
      return
      end

      subroutine hdown (h, ih, n)
!     h is in order, except that 1st element has been replaced.
!     Bubble it down to its proper location.  New element is i,
!     children are j and k.

      integer  n,i,j,k
      integer  ih(n)
      dimension h(n)

      i = 1

   10 continue
      j = 2*i
      k = j + 1

!     if j > n, new element is at bottom, we're done
      if (j .gt. n)  return
!     handle case where new element has only one child
      if (k .gt. n)  k = j

      if (h(j) .gt. h(k))  j = k
!     j is now index of smallest of children

      if (h(i) .gt. h(j))  then
         call swap (h(i), h(j))
         call iswap (ih(i), ih(j))
         i = j
         goto 10
      endif

      return
      end

      subroutine swap (a, b)
      t = a
      a = b
      b = t
      return
      end

      subroutine iswap (i, j)
      integer  i,j,k
      k = i
      i = j
      j = k
      return
      end
