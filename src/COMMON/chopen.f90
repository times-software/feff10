!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: chopen.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine chopen (ios, fname, mod)
!     Writes error msg and stops if error in ios flag from open
!     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog
      external istrln

!     open successful
      if (ios .le. 0)  return

!     error opening file, tell user and die.
      i = istrln(fname)
      j = istrln(mod)
      write(slog,100)  fname(1:i), mod(1:j)
      call wlog(slog)

  100 format (' Error opening file, ', a,                               &
     &        ' in module ', a)

      call wlog(' Fatal error')
      call par_stop('CHOPEN')
      end
