!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sfconv.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sub-program ffmod9
      program  ffmod9
!     subroutine ffmod9

!     Calculation of S_0^2 
!     written by Luke Campbell 2002
!     modified by Luke Campbell 2005 for new I/O structure

!     INPUT: s02.inp mod6.inp and any set of spectroscopy output files
!            (xmu.dat, chi.dat, chipNNNN.dat, feffNNNN.dat)
!     OUTPUT: specfunct.dat and the input spectroscopy files (overwritten)
	  use par
	  use sfconv_inp
	  use errorfile
      implicit none
	  integer ios


      call par_begin
      if (worker) go to 400
      call OpenErrorfileAtLaunch('sfconv')	
	  
!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='logsfconv.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logsfconv.dat', 'feff')

!     read  s02.inp
      call sfconv_read

      if (msfconv.eq.1) then 
         call wlog('Calculating S0^2 ...')
         call so2conv
         call wlog('Done with module: S0^2.'//char(13)//char(10))
      endif

  400 call par_barrier
      call par_end
      if(master)call WipeErrorfileAtFinish
!     sub-program ffmod9
      stop
!     return

      end
