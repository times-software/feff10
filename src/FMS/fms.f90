!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fms.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ffmod3
  ! Full multiple scattering code (inversion of big matrix)
  !
  ! INPUT:  geom.inp, global.inp and mod3.inp
  ! OUTPUT:  gg.bin
  !
  ! HISTORY:
  !   2008/11/03 jp.rinehimer - Modified for dynamic memory allocation
  !   2001 a.ankudinov  - New matrix inversion algorithms and new I/O structure
  !   2001 a.ankudinov  - Written using earlier subroutines coded by b.ravel

  use DimsMod, only: istatx, lx, nclusx, init_dimensions
  use par
  use fms_inp,only: mfms
  use stkets
  use rotx
  use lnlm
  use xstruc
  use errorfile
  use t3j
  use hubbard_inp

  implicit none
  real*8 wall_start, wall_end
  integer :: ios
  character*512 slog

  
  call par_begin
  if(master) call OpenErrorfileAtLaunch('fms')

  ! Initialize needed modules
  call init_dimensions
  call hubbard_init
  call init_stkets(istatx)
  call init_rotx(lx,nclusx)
  call init_lnlm(lx,nclusx)
  call init_xstruc(nclusx)
  call init_t3j(lx)
  

  ! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0
   
  ! Open the log file, unit 11.  See subroutine wlog.
  if (master) then
     open (unit=11, file='log3.dat', status='unknown', iostat=ios)
     call chopen (ios, 'log3.dat', 'feff')
  else
     par_type = 2
  endif


  ! Read  INPUT: files geom.dat, global.dat and mod3.inp
  call reafms

  if (mfms.eq.1 )  then
     if(master) call wlog("FMS calculation of full Green's function ...")
     ! Do fms inside rfms2 
         if (master.and.parallel_run)  write(slog,'(a,i5,a)') 'FEFF-MPI using ',numprocs,' parallel threads.'
         if (master.and.(.not.parallel_run))  write(slog,'(a,i5,a)') 'FEFF-serial using 1 thread.'
         call wlog(slog)
     call fmstot
  endif

  ! OUTPUT: gg.bin file for the next modules
  if (master) close (unit=11)

  ! Deallocate needed modules
  call kill_stkets
  call kill_rotx
  call kill_lnlm
  call kill_xstruc
  call kill_t3j

  ! Time at end of run
  call seconds(wall_end)
  if (master .and. mfms.eq.1) then
        if(parallel_run) write (6,'(a,f15.4,a,10x,a,e15.4,a)') &
         'total time ', wall_end - wall_start,'s', &
         '(communication time', wall_comm,'s)'
        call wlog('Done with module: FMS.'//char(13)//char(10))
  endif


  call par_end
  if(master) call WipeErrorfileAtFinish
  stop
end program ffmod3
