!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: pot.f90,v $:
! $Revision: 1.20 $
! $Author: hebhop $
! $Date: 2012/11/29 23:20:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program ffmod1

!     calculate  el. density and potential given atomic positions for
!     cluster atoms or other similar information
!     calculation can vary in complexity: self-consistency (on/off),
!     spin dependency (on/off), etc..
!       coded by a.l. ankudinov 2000, for modular code structure
!       modified by a.l. ankudinov 2001, for new i/o structure

!     INPUT files: mod1.inp, geom.inp
!     OUTPUT file: pot.bin

      use DimsMod, only: lx, nclusx, istatx, init_dimensions
	  use stkets !KJ
      use rotx
      use lnlm
      use xstruc
      use t3j  !KJ
      use par
	  use potential_inp,only: mpot
	  use errorfile
!     JJK - added use broydn to keep module in scope. 11/2012
      use broydn_workspace
      implicit none
	  integer ios
      real*8 wall_start, wall_end
      character*512 slog

      call par_begin
	  if(master)call OpenErrorfileAtLaunch('pot')
      ! Initialize needed modules
      call init_dimensions
      call init_stkets(istatx)
      call init_rotx(lx,nclusx)
      call init_lnlm(lx,nclusx)
      call init_xstruc(nclusx)
      call init_t3j(lx)

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


!     INPUT: read data in pot.inp  and geom.inp files
!     and transform it to atomic hartree units
      call reapot(.true.)
!     JJK - added call to initialize broydn vars. 11/2012
      call broydn_workspace_init
      if (mpot .eq. 1)  then
         if(master)call wlog('Calculating SCF potentials ...')
         if (master.and.parallel_run)  write(slog,'(a,i5,a)') 'FEFF-MPI using ',numprocs,' parallel threads.'
         if (master.and.(.not.parallel_run))  write(slog,'(a,i5,a)') 'FEFF-serial using 1 thread.'
         call wlog(slog)
         call pot
      endif

!     OUTPUT: subroutine pot writes main output file pot.bin
!     with information on potentials, necessary for other modules;
!     additional output files can be obtained using PRINT card

!--   Time at end of run
      call seconds(wall_end)
      if (master .and. mpot.eq.1) then
         if(parallel_run) write (6,'(a,f15.4,a,10x,a,e15.4,a)') &
           'total time ', wall_end - wall_start,'s', &
           '(communication time', wall_comm,'s)'
         call wlog('Done with module: potentials.'//char(13)//char(10))
      endif
      if(master) close (unit=11)

      call par_end

      if(master)call WipeErrorfileAtFinish
      stop

      end

