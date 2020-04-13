module density_inp
  use constants
  implicit none

  public :: density_inp_read, density_inp_cleanup, density_grid
  private

  integer, parameter, public :: GRID_TYPE_BASIC = 1

  integer, parameter :: nwordsx = 50
  integer, parameter :: wordsz = 50

  type density_grid
    integer :: grid_type
    character(30) :: filename
    integer :: ndims
    double precision :: origin(3)
    integer, allocatable :: npts(:)
    double precision, allocatable :: axes(:,:)
    logical :: core
  end type

  integer, public :: ngrids
  type(density_grid), allocatable, public :: grids(:)
contains

  subroutine next_line(fd, nwords, words, ier)
    integer, intent(in) :: fd
    integer, intent(inout) :: nwords
    character(wordsz), intent(in) :: words(nwordsx)
    integer, intent(out) :: ier

    character(1024) line
    call rdcmt(fd,'#!*C',ier)
    IF(ier.EQ.0) read(fd,'(A)',iostat=ier) line
    if (ier .ne. 0) return

    call untab(line)
    nwords = nwordsx
    call bwords(line,nwords,words)
  end subroutine

  function is_command(word)
    character(wordsz) :: word
    logical :: is_command
    is_command = (word == "line" .or. &
                  word == "plane" .or. &
                  word == "volume")
  end function

  subroutine density_inp_read

    integer fd ! file descriptor
    character(wordsz) words(nwordsx)
    integer ier, nwords, iword, igrid, i, j
    type(density_grid)::grid
    logical ltest
    character :: slog(100)

    fd = 22

    ! count number of command lines to determine how many grids are present
    ngrids = 0
    open(fd, file='density.inp', status='unknown')
    do
      call next_line(fd, nwords, words, ier)
      if (ier == -1) exit
      if (ier .ne. 0) goto 99
      if ( is_command(words(1)) ) ngrids = ngrids + 1
    end do
    close(fd)

    ! now read in grids
    allocate(grids(ngrids))
    igrid = 1
    open(fd, file='density.inp', status='unknown')
    do
      call next_line(fd, nwords, words, ier)
      if (ier .eq. -1) goto 90 ! end of file
      if (ier .ne. 0) goto 99 ! error

      if (.not. is_command(words(1))) then
        !write (slog, "(a, a)") "Unknown density grid type ", words(1)
        !call wlog(slog)
        call wlog("Unknown density grid type: "// trim(words(1)))
        goto 99
      end if

      grid = grids(igrid)
      ! default grid type is basic
      grid%grid_type = GRID_TYPE_BASIC
      if (words(1) ==  "line") then
        grid%ndims = 1
      else if (words(1) == "plane") then
        grid%ndims = 2
      else if (words(1) == "volume") then
        grid%ndims = 3
      else
        write (slog, "(a, a)") "Unhandled density grid type ", words(1)
        call wlog(slog)
      end if

      grid%filename = words(2)

      iword = 3
      if (grid%grid_type == GRID_TYPE_BASIC) then
        ! read in reast of parameters
        allocate(grid%npts(grid%ndims))
        allocate(grid%axes(3,grid%ndims))

        ! read in origin
        read(words(3),*) grid%origin(1)
        read(words(4),*) grid%origin(2)
        read(words(5),*) grid%origin(3)

        grid%core = .false.
        if (nwords .ge. 6 .and. words(6) == "core") then
          grid%core = .true.
        end if

        ! now read in unit vectors and number of points
        do i = 1,grid%ndims
          call next_line(fd, nwords, words, ier)
          read(words(1),*) grid%axes(1,i)
          read(words(2),*) grid%axes(2,i)
          read(words(3),*) grid%axes(3,i)
          read(words(4),*) grid%npts(i)
        end do
      end if

      ! convert from Angstroms to a.u.
      grid%origin = grid%origin / bohr
      grid%axes = grid%axes / bohr

      grids(igrid) = grid
      igrid = igrid + 1
    end do

90  continue ! end of file
    close(fd)
    return

99  continue ! error in file
    close(fd)
    call par_stop("Error reading density.inp. Exiting.")
  end subroutine

  subroutine density_inp_cleanup
    integer igrid

    if (.not. allocated(grids)) return

    do igrid = 1, ngrids
      if (allocated(grids(igrid)%axes)) deallocate(grids(igrid)%axes)
      if (allocated(grids(igrid)%npts)) deallocate(grids(igrid)%npts)
    end do
    deallocate(grids)
  end subroutine
end module
