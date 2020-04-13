!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdhead.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdhead (io, nhead, head, lhead)
      implicit double precision (a-h, o-z)

!     Reads title line(s) from unit io.  Returns number of lines
!     read.  If more than nheadx lines, skips over them.  End-of-header
!     marker is a line of 1 blank, 71 '-'s.
!     lhead is length of each line w/o trailing blanks.
!     header lines returned will have 1st space on line blank for
!     carriage control

      character*80 head(nhead)
      dimension lhead(nhead)
      character*80  line
      integer istrln
      external istrln

      n = 0
      nheadx = nhead
      nhead = 0
   10 read(io,20)  line
   20    format(a)
         if (line(4:11) .eq. '--------')  goto 100
         n = n+1
         if (n .le. nheadx)  then
            head(n) = line
            lhead(n) = istrln(head(n))
            nhead = n
         endif
      goto 10
  100 continue
      return
      end
