!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ldos.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program  ffmod7
  ! L-projected density of states (LDOS) calculation
  !   written by a.ankudinov 1998
  !   modified by a.ankudinov 2001 for new I/O structure

  !   INPUT: ldos.inp geom.dat and pot.bin
  !   OUTPUT: files ldosNN.dat and others

  use DimsMod, only: init_dimensions
  use par
  use ldos_inp, only: mldos
  use errorfile
  implicit none
  integer ios
  real*8 wall_start, wall_end
  character*512 slog


  call par_begin
  if(master) call OpenErrorfileAtLaunch('ldos')
  call init_dimensions

  !     Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0

  ! open the log file
  if (master) then
     open (unit=11, file='logdos.dat', status='unknown', iostat=ios)
     call chopen (ios, 'logdos.dat', 'feff')
  else
     par_type = 2
  endif

  ! read ldos.inp, reciprocal.inp, geom.inp
  call reldos

  if (mldos.eq.1) then
     ! Do work
     call wlog(' Calculating LDOS ...')
     if (master.and.parallel_run)  write(slog,'(a,i5,a)') 'FEFF-MPI using ',numprocs,' parallel threads.'
     if (master.and.(.not.parallel_run))  write(slog,'(a)') 'FEFF-serial using 1 thread.'
     call wlog(slog)
     call ldos_driver
     call wlog(' Done with LDOS.')
  endif

  !--   Time at end of run
  call seconds(wall_end)
  if (master .and. mldos.eq.1) then
        if(parallel_run) write (6,'(a,f15.4,a,10x,a,e15.4,a)') &
          'total time ', wall_end - wall_start,'s', &
          '(communication time', wall_comm,'s)'
        call wlog('Done with module: LDOS.'//char(13)//char(10))
  endif
  if (master) close(11)
  
  call par_end
  if(master) call WipeErrorfileAtFinish

  stop
end program ffmod7




