program compton
  !
  ! Compton Profile
  ! Written by B. A. Mattern 2011
  !
  ! INPUT: compton.inp pot.bin gg_slice.bin
  ! OUTPUT: rhozzp.dat jzzp.dat jpq.dat
  !
  use dimsmod, only: init_dimensions
  use atoms_inp
  use compton_inp
  use compton_mod
  use errorfile
  use par
  implicit none

  type(compton_grid) :: grid
  real*8 wall_start, wall_end
  integer ios
  character*512 slog

  call par_begin
  if(master)call OpenErrorfileAtLaunch('compton')
  call init_dimensions
  call compton_read
  if(run_compton_module .ne. 1) goto 400
  call compton_init

! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0


  ! Set up logging
  if (master) then
     open (unit=11, file='logcompton.dat', status='unknown', iostat=ios)
     call chopen (ios, 'logcompton.dat', 'feff')
  else
     par_type = 2
  endif

  if(master)call wlog('Calculating Compton scattering ...')
  if (master.and.parallel_run)  write(slog,'(a,i5,a)') 'FEFF-MPI using ',numprocs,' parallel threads.'
  if (master.and.(.not.parallel_run))  write(slog,'(a,i5,a)') 'FEFF-serial using 1 thread.'
  call wlog(slog)

  call compton_build_grid(grid) ! construct grid on which to calculate rho(r,r')

  if (do_rhozzp.and.master) then ! this is pretty quick, so not //ized
    call wlog("Calculating rho(z,z')")
    call calculate_rhozzp(grid)
  end if

  if (do_compton) then
    if(master)call wlog('Calculating Compton profile')
    call calculate_compton(grid)
  end if

  call seconds(wall_end)
  if (master) then
      if(parallel_run) write (6,'(a,f15.4,a,10x,a,e15.4,a)') &
         'total time ', wall_end - wall_start,'s', '(communication time', wall_comm,'s)'
      call wlog('Done with module: Compton scattering.'//char(13)//char(10))
  endif
  if (master) close(unit=11)
400  call par_end
  if(master)call WipeErrorfileAtFinish


contains

subroutine calculate_compton(grid)
  !
  ! Calculate J(pq)
  !
  ! First, rho(r,rp) is calculated on the supplied grid.
  ! This is integrated over the x-y plane to get J(z,zp)
  ! Finally, J(z,zp) is fourier transformed to obtain J(p_q)
  !
  ! J(z,zp,E) is saved out to jzzp.dat to allow quickly recalculating
  ! J(pq) on a different pq grid, or with a different apodization function.
  !
  use par
  use compton_mod
  use compton_inp
  implicit none

  type(compton_grid), intent(in) :: grid

  real*8, dimension(:,:), allocatable :: jzzp
  real*8 pq, J, smax_, phimax_, zmax_, zpmax_

  integer :: i, ns_, nphi_, nz_, nzp_, ne_
  logical :: jzzp_exists, jzzp_valid
  character*1 :: tmpstr

  allocate(jzzp(grid%nz, grid%nzp))

  ! Check if a cached calculation of J(z,z') exists in jzzp.dat, was calculated on the same
  ! grid as currently requested, and we don't specify "force_jzzp" in compton.inp
  inquire(file='jzzp.dat', exist=jzzp_exists)
  jzzp_valid = .false.
  if (jzzp_exists .and. .not. force_jzzp) then
    ! try loading it
    open(unit=20, file='jzzp.dat')
    read(20,*,err=800) tmpstr, ns_, nphi_, nz_, nzp_
    read(20,*,err=800) tmpstr, smax_, phimax_, zmax_, zpmax_
    read(20,*,err=800) jzzp

    ! ensure that saved grid is same as current one
    if (  (ns.eq.ns_).and.(nphi.eq.nphi_).and.(nz.eq.nz_).and.(nzp.eq.nzp_)    &
       &  .and.(abs(smax-smax_).lt.1e-6).and.(abs(phimax-phimax_).lt.1e-6)     &
       &  .and.(abs(zmax-zmax_).lt.1e-6).and.(abs(zpmax-zpmax_).lt.1e-6)) then
      if(master) call wlog('Reusing previously calculated j(z,z'')')
      jzzp_valid = .true.
    else
      if(master)then
	  print *, "smax", smax, smax_
	  print *, "phimax", phimax, phimax_
	  print *, "zmax", zmax, zmax_
	  print *, "zpmax", zpmax, zpmax_
      call wlog('Grid for saved j(z,z'',E) differs from current. Recalculating.')
      endif
    end if

    goto 810
    800 continue
    call wlog('Error reading in saved j(z,z'',E). Recalculating.')
    810 continue
    close(20)
  end if

  ! If we don't have a valid cached calculation of J(z,z'), then calculate it here (and cache)
  if (.not.jzzp_valid) then
    !if(master)call wlog('Calculating j(z,z'')')
    ! integrate over x,y plane
    call compton_jzzp(grid, jzzp)

    if (master) then
      ! save out jezzp for later
      call wlog('Saving j(z,z'')')
      open(unit=20, file='jzzp.dat')
      write(20, *) '#', ns, nphi, nz, nzp
      write(20, '(a)', advance='no') '#'
      write(20, '(20f20.13)') smax, phimax, zmax, zpmax
      write(20, *) jzzp
      close(20)
    end if
  end if

  ! Fourier transform jzzp to get the Compton profile J(p_q)
  if (master) then
    call wlog('Calculate j(pq)')
    ! perform fourier transform
    open(unit=20, file='compton.dat')
    call write_jpq_header(grid, 20)

    ! Loop over requested number of p_q points
    do i=1,npq
      if (npq .eq. 1) then
        pq = 0
      else
        pq = pqmax / (npq-1) * (i - 1)
      end if

      call jpq(grid, jzzp, pq, J)
      write (20,*) pq, J
    end do
    close(20)
  end if

end subroutine

subroutine write_jpq_header(grid, fileno)
  !
  ! Write out header information for compton.dat
  !
  use compton_mod
  use compton_inp
  implicit none
  type(compton_grid), intent(in) :: grid
  integer, intent(in) :: fileno

  write(fileno, *) '# Compton profile, J(pq)'
  write(fileno, *) '# ns:  ', grid%ns
  write(fileno, *) '# nphi:', grid%nphi
  write(fileno, *) '# nz:  ', grid%nz
  write(fileno, *) '# nzp: ', grid%nzp
  write(fileno, *) '# zpmax:', -grid%zp(1)
  write(fileno, *) '# temperature (eV):', temperature
  write(fileno, *) '#----------------------------'
  write(fileno, *) '# pq               J'
end subroutine

subroutine calculate_rhozzp(grid)
  !
  ! Calculate rho(z,zp) at fixed z=0.01 with z' varying from 0.01 to zpmax
  !
  ! This is primarily used for diagnostics and to determine where
  ! to place zpmax.
  !
  ! This could easily be made more flexible to caluclate other slices of
  ! rho(r,r').
  !
  !use dimsmod, only:
  use constants
  use atoms_inp
  use compton_mod
  use compton_inp
  use rhorrp_mod
  implicit none

  type(compton_grid) :: grid

  real*8, dimension(3) :: v, vp
  real*8 :: zp
  integer i, n
  real*8 rho

  open(unit=99, file='rhozzp.dat')
  n = 1000
  do i=1,n
    v(1) = 0.0
    v(2) = 0.0
    v(3) = 0.01

    vp(:) = 0.0
    vp(3) = .01 + grid%zp(grid%nzp)/(n-1) * (i-1)

    zp = vp(3)

    if (grid%rotate) then
      call rotate_in_place(grid%rotation_matrix, v)
      call rotate_in_place(grid%rotation_matrix, vp)
    end if

    call rhorrp(v, vp, rho)
    write(99,*) zp, rho
  end do
  close(99)
end subroutine

end program
