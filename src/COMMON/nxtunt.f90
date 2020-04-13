!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: nxtunt.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!====================================================================
      integer function nxtunt(iunit)

!  this function returns the value of the next unopened unit number
!  equal to or larger than iunit.  it will return neither unit numbers
!  0, 5, or 6 nor a negative unit number
! $Id: nxtunt.f90,v 1.2 2010/02/23 23:52:06 hebhop Exp $
! $Log: nxtunt.f90,v $
! Revision 1.2  2010/02/23 23:52:06  hebhop
! Added header to all .f90 files to show version info, i.e. date, author, revision number of last modification.
!
! Revision 1.1  2007/12/08 08:06:16  hebhop
! Lots of changes.
!
! 1. Changed all filenames to *.f90
! 2. Changed Makefile to compile all files directly instead of making recursively.
! 3. Added DEP directory and dependency files.
! 4. Added Utilities to make dependency files. Note: these probably aren't robust and should be rewritten.
! 5. Added I/O modules and Error modules.
! 6. Changed a few routines to use I/O modules.
! 7. Split atomic calculation from SCF potentials calculation.
!    Atomic calculations now calculate things for all edges.
! 8. Added hack (soon to be option) to read muffin tin potentials from external file.
!
! Revision 1.1.1.1  2007/10/23 22:07:50  ytakimot
!
! FEFF9 Project
!
! Revision 1.1.1.1  2007/06/04 23:16:47  hebhop
! Initializing feff86. All references to so2conv have changed to sfconv
! Josh Kas
!
! Revision 1.1.1.1  2006/01/12 06:37:42  hebhop
! New version of feff. feff8.5 (Extension of feff8.4)
! Includes:
! 	1) All feff8.4 capabilities.
! 	2) Screened core hole (calculation of W).
! 	3) Multiple pole self energy calculation.
! 	4) Convolution with spectral function.
! New cards and options:
! 	1) NOHOLE 2      (screened hole)
! 	2) PLASMON ipl   (multiple pole self energy)
! 	3) SO2CONV       (convolve output with spectral function)
! 	4) SELF          (print on shell self energy as calculated by Luke)
! 	5) SFSE k0        (print off shell self energy Sigma(k0,e) )
!
! Revision 1.1.1.1  2000/02/11 02:23:58  alex
! Initialize feff82
!
! Revision 1.10  1999/04/02 21:32:47  newville
! cleaned up nxtunt (matt)
!
! Revision 1.9  1999/02/11 20:08:08  alex
! x39 version: dim.h + misc. small changes
!
! Revision 1.8  1998/12/29 23:59:07  alex
! feff8x35 version
!
! Revision 1.7  1998/11/19 03:23:11  alex
! feff8x32 version
!
! Revision 1.6  1998/10/26 14:11:16  ravel
! no comments beyond column 71
!
! Revision 1.5  1998/10/18 21:47:51  alex
! feff8x30 version implements Broyden algorithm for self-consistency
!
! Revision 1.4  1998/02/24 18:31:37  ravel
! I should really be more careful.  This is the last commitment done
!      cright.
!
! Revision 1.1.1.1  1997/04/27 20:18:03  ravel
! Initial import of xanes sources, version 0.37
!
! Revision 1.1  1996/06/23 16:05:02  bruce
! Initial revision
!

       integer iunit
       logical open

       nxtunt = max(1, iunit) - 1
 10    continue
       nxtunt = nxtunt + 1
       if ((nxtunt.eq.5).or.(nxtunt.eq.6)) nxtunt = 7
       inquire (unit=nxtunt, opened=open)
       if (open) go to 10
       return
!  end integer function nxtunt
       end

!====================================================================
