program mpi_hello

! include 'mpif.h'
  use mpi, only : mpi_init, mpi_comm_rank, mpi_comm_size, mpi_finalize
  use mpi, only : MPI_COMM_WORLD
! use mpi, only : MPI_IRECV

  implicit none

  integer ierr,my_rank,size
  character(len=40) :: nodename
  integer :: count0, count1, count_rate

  call mpi_init(ierr)

  call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,size,ierr)

! print rank and size to screen

  call system_clock(count0,count_rate)
! write(6,fmt='(2i18,a,a12)') count0, count_rate, '  ', 'time'
! if ( my_rank .eq. 0 ) then
!   write(6,fmt='(i18,a,a12)') count0/count_rate, '  ', 'time'
! end if
! write(6,fmt='(a,i5,a,i5)') &
!   'Hello World! I am rank ', my_rank, ' of size ', size
  call get_environment_variable('SLURMD_NODENAME',nodename)
  write(6,fmt='(3a,i6,i18)') &
    'NODE: ', trim(nodename), ' time ', my_rank, count0/count_rate

  call mpi_finalize(ierr)

end program mpi_hello
