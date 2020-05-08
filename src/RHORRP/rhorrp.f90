program rhorrp_prog
  !
  ! Main program for calcalating real-space electronic densities
  !
  use rhorrp_mod
  use density_inp
  use constants
  use par
  use errorfile
  implicit none

  integer ios
  integer i, j

  call par_begin
  if(master)call OpenErrorfileAtLaunch('rhorrp')

  call density_inp_read
  if(ngrids.le.0) goto 400  ! This means empty density.inp because no RHORRP card in feff.inp

  ! open the log file, unit 11.  See subroutine wlog.
  if (master) then
    open (unit=11, file='logrhorrp.dat', status='unknown', iostat=ios)
    call chopen (ios, 'logrhorrp.dat', 'feff')
  else
    par_type = 2 ! supress wlog calls from workers
  endif
  if(master) call wlog("Calculating charge density output ...")

  call rhorrp_init
  !print *, "ngrids: ", ngrids
  !do i=1,ngrids
  !  print "(a, i3)", "grid ", i
  !  print "(a, i10)", "  ndims: ", grids(i)%ndims
  !  print "(a, (3f10.3))", "  origin: ", grids(i)%origin(:)
  !  do j=1,grids(i)%ndims
  !    print "(a, (3f10.3) i)", "  ", grids(i)%axes(:,j), grids(i)%npts(j)
  !  end do
  !end do

  call calculate_densities

  !call test_nearest_atom
  !call test_rhozzp
  !call test_rhoe
  call rhorrp_deinit


  if(master) call wlog("Done with module: charge density output (RHORRP)."//char(13)//char(10))
400 continue
  call par_end

  if(master) call WipeErrorfileAtFinish


contains

  subroutine calculate_densities
    !
    ! Calculate density for each grid specified in density.inp
    !
    integer i

    do i = 1,ngrids
      call calculate_density(grids(i))
    end do
  end subroutine

  function filename_is_binary(filename)
    !
    ! Check whether a string ends in ".bin"
    !
    character(*), intent(in) :: filename
    logical filename_is_binary
    integer :: idot
    character(4) :: extension

    idot = scan(filename, ".", .true.)
    if (idot > 0) then
      extension = filename(idot+1:)
      call lower(extension)
      filename_is_binary = extension == "bin"
    else
      filename_is_binary = .false.
    end if
  end function

  subroutine next_index(grid, idx)
    ! Increment multi-dimensional iterator index
    type(density_grid), intent(in) :: grid
    integer, intent(inout) :: idx(grid%ndims)

    integer i

    do i=1,grid%ndims
      if (idx(i) < grid%npts(i)) then
        idx(i) = idx(i) + 1
        exit
      else
        idx(i) = 1
      end if
    end do
  end subroutine

  subroutine point_at_index(grid, idx, pt)
    ! Calculate point in grid at index specified by idx
    type(density_grid), intent(in) :: grid
    integer, intent(in) :: idx(grid%ndims)
    double precision, intent(out) :: pt(3)

    integer i

    pt(:) = grid%origin(:)
    do i = 1, grid%ndims
      pt(:) = pt(:) + (grid%axes(:,i) / (grid%npts(i)-1)) * (idx(i) - 1)
    end do
  end subroutine

  subroutine calculate_density(grid)
    type(density_grid), intent(in) :: grid

    integer :: idx(grid%ndims), i, j, n
    integer :: totpts, istat
    double precision :: v(3)
    double precision, allocatable :: points(:,:), rho(:)
    character(8) :: spts
    integer :: proc_i1(0:numprocs-1), proc_i2(0:numprocs-1), i1, i2

    integer :: fd, pts_per_proc, extra_pts, process
    logical :: binary_output
    integer iat, iph
    double precision dv(3)
    logical PrintNearestAtomInfo


    ! Debug flag for printing info about nearest neighbour of to the current position of density 
    PrintNearestAtomInfo =.true.


    ! Total number of points to calculate
    totpts = product(grid%npts)

    ! Periodically log progress
    write (spts, "(i8)") totpts
    call wlog("Calculate density: "//trim(grid%filename)//" ("//trim(spts)//" total points)")

    allocate(points(3,totpts), rho(totpts))

    ! open output file early to ensure that we can write to it before spending time doing calculation
    fd = 22
    if (filename_is_binary(grid%filename)) then
      open(unit=fd, file=grid%filename, status='unknown', form='unformatted', iostat=istat)
      binary_output = .true.
    else
      open(unit=fd, file=grid%filename, status='unknown', iostat=istat)
      binary_output = .false.
    end if

    if (istat .ne. 0) then
      call wlog("Unable to write to file: " // trim(grid%filename) // ". Skipping.")
      return
    end if

    idx(:) = 1

    ! determine typical number of points per processor (final processor may have fewer)
    pts_per_proc = totpts / numprocs
    extra_pts = mod(totpts, numprocs)

    ! distribute points across processes (store first and last points for each process)
    ! all points between proc_i1(n) and proc_i2(n), inclusive, will be handled by process n
    i = 1
    do n=0,numprocs-1
      proc_i1(n) = i
      proc_i2(n) = i + pts_per_proc - 1
      if (n < extra_pts) proc_i2(n) = proc_i2(n) + 1
      i = proc_i2(n) + 1
    end do

    ! calculate density
    do i=1,totpts
      call point_at_index(grid, idx, points(:,i))
      if (proc_i1(this_process) <= i .and. i <= proc_i2(this_process)) then
        if (grid%core) then
          call atomic_density(points(:,i), 1, rho(i))
        else
          call rhorrp(points(:,i), points(:,i), rho(i))
        end if
      end if
      call next_index(grid, idx)
    end do

    call wlog("Calculation complete, transfer back to master.")
    ! transfer calculations back to master
    if (master) then
      do n=1,numprocs-1
        i1 = proc_i1(n)
        i2 = proc_i2(n)
        !print *, "recv from", n, i1, i2
        call par_recv_double(rho(i1:i2), i2-i1+1, n, 0)
      end do
    else
      i1 = proc_i1(this_process)
      i2 = proc_i2(this_process)
      !print *, "send from", this_process, i1, i2
      call par_send_double(rho(i1:i2), i2-i1+1, 0, 0)
    end if

    call par_barrier

    if (master) then
      ! convert units to angstroms for output
      do i=1,totpts
        points(:,i) = points(:,i)*bohr
        rho(i) = rho(i) / bohr**3
      end do

      call wlog("Saving density calculation.")

      ! now save to file

      if (binary_output) then
        ! ---- binary format ----
        ! For higher-dimensional / denser grids, we can save disk space by storing only the density values
        ! and the grid specification (which is sufficient to regenerate the grid points)
        write(fd) grid%ndims
        write(fd) grid%origin*bohr
        do i=1,grid%ndims
          write(fd) grid%axes(:,i)*bohr
          write(fd) grid%npts(i)
        end do
        write(fd) rho
        close(fd)
      else
        ! ---- ascii format ----
        ! To make plotting easier, we write out the coordinates for each point

        if (.NOT. PrintNearestAtomInfo) then
          do i=1,totpts
            write(fd, '(4(ES12.5,TR1))') (points(j,i), j=1,3), rho(i)
          end do
        else
          do i=1,totpts
            call nearest_atom(points(:,i)/bohr, iat, iph, dv, .false.)
            iat = iat - 1
            write(fd, '(7(ES12.5,TR1),TR1,I2,TR1,I1)') (points(j,i), j=1,3), rho(i), dv, iat, iph
          end do
        end if
        close(fd)
      end if
    end if

    deallocate(points, rho)
  end subroutine

!
! code below here is debugging code, and can be ignored / removed
!

subroutine test_nearest_atom
  double precision r
  double precision v(3), dv(3)
  integer iat, iph

  r = 0.0
  v(:) = 0
  do while (r < 10.0)
    v(3) = r
    call nearest_atom(v, iat, iph, dv, .true.)
    print *, r, iat, iph
    r = r + 0.5
  end do
end subroutine

end program

