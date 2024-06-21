!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: screen.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     FFMOD8
!=======================================================================
!     sub-program screen
program  screenw2
  !     subroutine ffmod8
  !     Written by Y. Takimoto

  use DimsMod, only: nrptx, init_dimensions
  use constants
  use par  
  use potential_inp,only:nohole,potential_read
  use errorfile
  !use screenw_inp ! Add new input for dx to play with grid.
  implicit none

  double precision dx, x0, vr0, ri(nrptx)
  !     W, the result of calculation
  double precision vch(nrptx), wscrn(nrptx)
  integer i, ilast, ixc0, ios
  integer,parameter :: nrx = 2502  !251 !2501 !1251 

  call par_begin
  !if (worker) go to 400
  if(master) call OpenErrorfileAtLaunch('screen')
  call init_dimensions

  call potential_read
  !if (nohole.ne.2) goto 5

  if(master) call wlog('Calculating screened core-hole potential ...')
  if(master) open (unit=11, file='logscreen.dat', status='unknown', iostat=ios)
  if(master) call chopen (ios, 'logscreen.dat', 'feff')

  !=================================================================
  !     read  INPUT data files: geom.inc and ldos.inp.
  !=================================================================
  call rdgeom ! read geom.inc and ldos.inp
  call rdscreen

  !=================================================================
  !     Prepare, then start the calculation
  !=================================================================
  !     make radial grid
  !dx   = 0.05d0
  !dx   = 0.03d0
  !dx = 0.01d0
  dx   = 0.0075d0
  x0   = 8.8d0
  call setri(dx, x0, nrx, ri)

  vr0  = 0.0d0
  ixc0 = 0

  call prepw2(vr0, nrx, ri, dx, x0, ilast)

  if(master) close (unit=11)
  if(master) call wlog('Done with module: screened core-hole potential.'//char(13)//char(10))
5 continue

400 call par_barrier
  call par_end
  if(master)call WipeErrorfileAtFinish
  !     sub-pro exchange point
  stop
  !     return

end program screenw2
!=======================================================================
!     END FFMOD8
!=======================================================================
