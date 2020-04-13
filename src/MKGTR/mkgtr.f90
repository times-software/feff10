!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mkgtr.f90,v $:
! $Revision: 1.13 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mkgtr

  !     full multiple scattering code (inversion of big matrix)
  !     written by a.ankudinov 2000 using earlier written subroutines
  !     coded by b.ravel
  !     modified by a.ankudinov 2001 for new matrix inversion algorithms
  !     and new I/O structure

  !     INPUT:  geom.inp, global.inp and mod3.inp
  !     OUTPUT:  fms.bin
  use DimsMod, only: init_dimensions
  use par
  use fms_inp,only: mfms
  use global_inp,only: do_nrixs
  use errorfile

  implicit none
  real*8 wall_start, wall_end
  integer :: ios

  call par_begin
  if(worker) goto 400
  call OpenErrorfileAtLaunch('mkgtr')	
  call init_dimensions
  ! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0

! Open the log file, unit 11.
     open (unit=11, file='log3.dat', status='unknown', access='append',iostat=ios)
     call chopen (ios, 'log3.dat', 'feff')
! Read  INPUT: files geom.dat, global.dat and mod3.inp
     call reafms
      
     if (mfms.eq.1 )  then
        call wlog("MKGTR: Tracing over Green's function ...")
		!Trace over G
		!KJ 7-09 feffq and regular feff different enough that I'll allow a separate subroutine :
		if (do_nrixs.eq.1) then  !NRIXS calculation
			call getgtrjas
		else  !regular FEFF calculation
			call getgtr
		endif
		call wlog('Done with module: MKGTR.'//char(13)//char(10))
     endif

! OUTPUT: gtr.bin file for the next modules
     close (unit=11)

! Time at end of run
  call seconds(wall_end)

400  call par_end
  if(master)call WipeErrorfileAtFinish
  stop

end program mkgtr
