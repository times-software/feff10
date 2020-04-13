!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wlog.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wlog (string)
      use par
      character*(*) string
      integer istrln
      external istrln

!     This output routine is used to replace the PRINT statement
!     for output that "goes to the terminal", or to the log file.
!     If you use a window based system, you can modify this routine
!     to handle the running output elegantly.
!     Handle carriage control in the string you pass to wlog.
!
!     The log file is also written here, hard coded here.

!     The log file is unit 11.  The log file is opened in the
!     main program, program feff.

!     make sure not to write trailing blanks

   10 format (a)

!     Suppress output in sequential loops
      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print 10
         if (par_type .ne. 3) write(11,10)
      else
         print 10, string(1:il)
         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
!     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
