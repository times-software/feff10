!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: screen.f90,v $:
! $Revision: 1.7 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program  crpa

  use DimsMod, only: nrptx, init_dimensions
  use constants
  use par  
  use potential_inp,only:nohole,potential_read
  use crpa_inp
  implicit none

  double precision dx, x0, vr0, ri(nrptx)
  !     W, the result of calculation
  double precision vch(nrptx), wscrn(nrptx)
  integer i, ilast, ixc0, ios
  integer,parameter :: nrx = 251

  call par_begin
  if (worker) go to 400
  call init_dimensions
  call crpa_read
  call potential_read

  if (CRPAI%do_CRPA .eq. 1) then
      open (unit=11, file='logscrn.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logscrn.dat', 'feff')
      !     read  INPUT data files: geom.inc and ldos.inp.
      call rdgeom
      call wlog(' Calculating Hubbard U.')
      !     make radial grid
      dx   = 0.05d0
      x0   = 8.8d0
      call setri(dx, x0, nrx, ri)
      vr0  = 0.0d0
      ixc0 = 0
      call prep(vr0, ixc0, nrx, ri, dx, x0, ilast, vch, wscrn,.TRUE.)
      call wlog(' Done with Hubbard U calculation.')
      !     open files for output
      open (unit=24, file='wscrn.dat', status='unknown', iostat=ios)
      do i = 1, ilast
         write(24, '(33e20.10)') ri(i), wscrn(i), vch(i)
      end do
      close (unit=24)
      close (unit=11)
endif

400 call par_barrier
  call par_end
  stop

end program crpa
