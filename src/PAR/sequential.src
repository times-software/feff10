!  **************************************************
!  Parallel feff8 routines
!  Jim Sims
!  **************************************************

      subroutine par_begin
!  **************************************************
!  Initializations for parallel version(s)
!  **************************************************

      use par

!-- So cvd or dbx can attach to a running process
!     call sleep(30) 

      numprocs = 1
      my_rank = 0
      this_process = my_rank

      par_type = 0
      parallel_run = .false.
!-- The following variable will be used for IO that should only be
!-- done in one process.
      master = (my_rank .eq. 0)

      worker = (.not. master)
      if (worker) par_type = 1

      return
      end

      subroutine par_stop (string)
!  **************************************************
!  Abnormal termination of the parallel session
!  **************************************************
      use par
!     For abnormal exits 
!     If open, close unit = 11
!     Go to the barrier that workers are sitting at
!     Then everyone will call par_end and stop
      logical is_open
      character*(*) string

      inquire(unit=11,opened=is_open)
      if (is_open) then
        call wlog(string)
        close(unit=11)
      else if (string .ne. ' ') then
        print *,string
        print *,'Abnormal termination on processor ',this_process
      endif

      stop ' '
      end

      subroutine par_end
!  **************************************************
!  Terminate the parallel session
!  **************************************************
      return
      end

      subroutine par_barrier
!  **************************************************
!  Calls mpi_barrier
!  **************************************************
      return
      end

      subroutine par_send_int(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for integer arrays
!  **************************************************
      integer count,dest,tag
      integer buf(*)
      return
      end

      subroutine par_send_int_scalar(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for integer arrays
!  **************************************************
      integer count,dest,tag
      integer buf
      return
      end


      subroutine par_send_cmplx(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for complex arrays
!  **************************************************
      integer count,dest,tag
      complex buf(*)
      return
      end

      subroutine par_send_real(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for real arrays
!  **************************************************
      integer count,dest,tag
      real buf(*)
      return
      end

      subroutine par_send_dc(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for double_complex arrays
!  **************************************************
      integer count,dest,tag
      complex*16 buf(*)
      return
      end

      subroutine par_send_double(buf,count,dest,tag)
!  **************************************************
!  Call mpi_send for double arrays
!  **************************************************
      integer count,dest,tag
      double precision buf(*)
      return
      end

      subroutine par_recv_int(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for integer arrays
!  **************************************************
      integer count,source,tag
      integer buf(*)
      return
      end

      subroutine par_recv_int_scalar(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for integer arrays
!  **************************************************
      integer count,source,tag
      integer buf
      return
      end

      subroutine par_recv_cmplx(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for complex arrays
!  **************************************************
      integer count,source,tag
      complex buf(*)
      return
      end

      subroutine par_recv_real(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for real arrays
!  **************************************************
      integer count,source,tag
      real buf(*)
      return
      end

      subroutine par_recv_dc(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for double complex arrays
!  **************************************************
      integer count,source,tag
      complex*16 buf(*)
      return
      end

      subroutine par_recv_double(buf,count,source,tag)
!  **************************************************
!  Call mpi_recv for double arrays
!  **************************************************
      integer count,source,tag
      double precision buf(*)
      return
      end

      subroutine par_bcast_int(buf,count,source)
!  **************************************************
!  Call mpi_bcast for integer arrays
!  **************************************************
      integer count,source
      integer buf(*)
      return
      end

      subroutine par_bcast_cmplx(buf,count,source)
!  **************************************************
!  Call mpi_bcast for complex arrays
!  **************************************************
      integer count,source
      complex buf(*)
      return
      end

      subroutine par_bcast_real(buf,count,source)
!  **************************************************
!  Call mpi_bcast for real arrays
!  **************************************************
      integer count,source
      real buf(*)
      return
      end

      subroutine par_bcast_dc(buf,count,source)
!  **************************************************
!  Call mpi_bcast for double_complex arrays
!  **************************************************
      integer count, source
      complex*16 buf(*)
      return
      end

      subroutine par_bcast_double(buf,count,source)
!  **************************************************
!  Call mpi_bcast for double arrays
!  **************************************************
      integer count, source
      double precision buf(*)
      return
      end

      subroutine MPE_DECOMP1D( n, num_procs, myid, s, e )
!  ******************************************************
!  A routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may 
!  be used in "direct" product decomposition.  The 
!  values returned assume a "global" domain in [1:n]
!  ******************************************************
!  MPE_Decomp1d - Compute a balanced decomposition of
!  a 1-D array
!  ******************************************************
!  Input Parameters:
!  n  - Length of the array
!  num_procs - Number of processors in decomposition
!  myid  - Rank of this processor in the decomposition 
!  (0 <= rank < size)
!  ******************************************************
!  Output Parameters:
!  s,e - Array my_particles are s:e, with the original 
!  array considered as 1:n.  
!  ******************************************************

      integer n, num_procs, myid, s, e
      integer nloc, deficit
 
      nloc  = n / num_procs
      s       = myid * nloc + 1
      deficit = mod(n,num_procs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
        nloc = nloc + 1
      endif
      e = s + nloc - 1
      if (e .gt. n .or. myid .eq. num_procs-1) e = n

      return
      end

      SUBROUTINE SECONDS( W)
!  ***************************************************
!  SECONDS returns the wall clock times for a process
!  in seconds.
!  ***************************************************
 
      REAL*8      W

      W = 0.0

      RETURN
      END
