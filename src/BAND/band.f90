 program band


  use DimsMod, only: init_dimensions
  use par
  use band_inp
  use errorfile

  implicit none
  real*8 wall_start, wall_end
  integer :: ios


  call par_begin
  if(master)call OpenErrorfileAtLaunch('band')

  ! Initialize needed modules
  call init_dimensions

  ! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0


!     open the log file, unit 11.  See subroutine wlog.
      if (master) then
        open (unit=11, file='logband.dat', status='unknown', iostat=ios)
        call chopen (ios, 'logband.dat', 'feff')
      else
        par_type = 2
      endif

!     read  INPUT: files geom.dat, global.dat and mod3.inp
      call reafms
	  call band_read

      if (mband.eq.1)  then
		 if(master)call wlog('Calculating band structure ...')
         call bandtot    
         if(master)call wlog(' Done with module: band structure.'//char(13)//char(10))
      endif

  ! Time at end of run
  call seconds(wall_end)
  if (master .and. parallel_run) then
     write (6,*) 'total time     ', wall_end - wall_start
     write (6,*) 'communicate time', wall_comm
  endif

  call par_end
  if(master)call WipeErrorfileAtFinish
  stop
 end
