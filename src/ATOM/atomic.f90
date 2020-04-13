!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: atomic.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/05/30 00:55:55 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program atomic_pot

  !     calculate  el. density and potential given atomic positions for cluster atoms or other similar information
  !     calculation can vary in complexity: self-consistency (on/off), spin dependency (on/off), etc..
  !       coded by a.l. ankudinov 2000, for modular code structure
  !       modified by a.l. ankudinov 2001, for new i/o structure
  !       modified again K. Jorissen 2009

  !     INPUT files: pot.inp, geom.dat
  !     OUTPUT file: apot.bin
  use DimsMod, only: init_dimensions
  use par
  use potential_inp,only: mpot
  use errorfile
  implicit none
  real*8 wall_start, wall_end
  integer ios

  call par_begin
  if(master)call OpenErrorfileAtLaunch('atomic')
  if(worker) goto 400
  call init_dimensions

  !     Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0

  !     open the log file, unit 11.  See subroutine wlog.
  if (master) then
     open (unit=11, file='log1.dat', status='unknown', iostat=ios)
     call chopen (ios, 'log1.dat', 'feff')
  else
     par_type = 2
  endif

  !     INPUT: read data in pot.inp  and geom.dat files
  call reapot(.false.)  !KJ 1-2012 don't launch kprep since we won't need it -> saves 90% of runtime ...

  if (mpot .eq. 1)  then
     call wlog('Calculating atomic potentials ...')
     call AtomicPotentials
     call wlog('Done with module: atomic potentials.'//char(13)//char(10))
  endif

  if (master) close (unit=11)

  !--   Time at end of run
  call seconds(wall_end)

400  call par_end

  if(master)call WipeErrorfileAtFinish
  stop
  
end program atomic_pot
