! Kevin Jorissen 2012.  Purpose: to pass exit codes to  an external program, e.g. the JFEFF Java GUI which launches FEFF9 modules and must be able to figure out if they succeed before launching the next module.
! There is no really reliable way to set exit codes in Fortran (compiler/platform dependencies; exit codes are only included in very latest Fortran standards (2008? 2010?).
! Therefore, we copy the WIEN2k approach.  Set error file in working directory at program launch.  Wipe it on successful termination.  The GUI can then check for the presence of a non-zero-size file:
! If such a file exists, the program did not exit cleanly, signifying a crash.  This approach is robust and does not depend on any Fortran programming to catch runtime exceptiosn, memory allocation problems, ...
! The downside is that the file will not contain much useful information.  However any runtime information will still be printed to the screen.

! There is already some error handling code present in FEFF, introduced by Josh; presumably only used in a few of the routines he contributed?  In any case, it is not set up to output to file rather than stdout/err,
! and I don't want to mess with it.  The simple code below is good enough for me.

module errorfile

implicit none
character*11,private :: ErrorFileName='.feff.error'
integer,private :: lun = 77


contains

subroutine OpenErrorfileAtLaunch(ModuleName)
   ! Open the errorfile and set a default message.  Call at the start of a module.
   character*(*),intent(in) :: ModuleName
   character*500 :: ErrorMessage
! Modified by FDV:
! Using shorter lines to help compile in Solaris Studio
   ErrorMessage= 'Starting FEFF9 module ' //ModuleName// &
                 '.  If this message is still here after the module' // &
                 ' finishes running, it must have crashed. The content ' // &
                 'of this file is wiped on successful termination.'
   call SetErrorfileMessage(ErrorMessage)
   return
end subroutine OpenErrorfileAtLaunch


subroutine WipeErrorfileAtFinish
  !Overwrite the error file with an empty, 0-byte file.  Call at the regular termination of a module.
  open(lun,file=ErrorFileName,status='replace',err=1000)
  close(lun)
  return
  1000 stop 'Unable to wipe errorfile in SetErrorfileMessage.  How ironic.'
end subroutine WipeErrorfileAtFinish



subroutine SetErrorfileMessage(ErrorMessage)
   !Write a user(programmer)-specified error message to the errorfile.  Useful for diagnostic purposes.
   character*(*),intent(in) :: ErrorMessage
   open(lun,file=ErrorFileName,status='unknown',err=1000)
   write(lun,*) ErrorMessage
   close(lun)
   return
   1000 stop 'Unable to open errorfile in SetErrorfileMessage'
end subroutine SetErrorfileMessage



end module errorfile
