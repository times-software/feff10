!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_kklist.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2012/01/30 06:01:58 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       THE K-MESH TO SAMPLE THE BRILLOUIN ZONE
!*****************************************************************************
        module kklist
!    Number of k-points for the BZ mesh
      integer nkp
!    BZ mesh size, specified as nkx x nky x nkz
     integer nkx,nky,nkz
!    Use symmetry (1) or not (0) for this mesh
      integer usesym
!    Type of k-mesh
      integer ktype
	  ! ktype=1  :  regular mesh of nkp points for all modules
	  ! ktype=2  :  use nkp points for ldos/fms and nkp/5 points for pot  (significant time savings)
	  ! ktype=3  :  use nkp points for ldos/fms and nkp/5 points for pot (near edge) ; reduce nkp for all modules as we get away from near-edge
!    Rotation matrices for spherical harmonics
      !complex*16 drot(32,32,48,2)  !lx=3
      complex*16 drot(50,50,48,2)   !lx=4
        complex*16, allocatable :: mrot(:,:,:)
!    The k-mesh itself!
      real*8, allocatable :: bk(:,:)
!    Corresponding integration weights
      real*8, allocatable :: weight(:)
!    Sum of the integration weights
      real*8 sumweights
!    Correspondence between wien2k and sprkkr symmetry matrices
      integer symid(2,48)
!    Arrays that code for the relation between full and reduced k-mesh
      integer,allocatable :: intn(:),inti(:,:,:)
!    Which symmetries are actually used for the k-mes
      integer symact(48)
!    More arrays
      real*8,allocatable :: intw(:,:)


      CONTAINS
          subroutine init_kklist(n,nsym) !KJ 6-09
!            use struct,only: nsym !KJ 6-09
            implicit none
            integer,intent(in) :: n,nsym !KJ added nsym 6-09
            allocate(bk(3,n),weight(n))
            allocate(intn(n),inti(n,nsym,2),intw(n,nsym))
            intn=0
            inti=0
            intw=dble(0)
            symact=0
            bk=dble(0)
            weight=dble(0)
            sumweights=dble(0)
            symid=0
          end subroutine init_kklist

          subroutine destroy_kklist
		     if(allocated(bk)) deallocate(bk)
		     if(allocated(weight)) deallocate(weight)
		     if(allocated(inti)) deallocate(inti)
		     if(allocated(intn)) deallocate(intn)
		     if(allocated(intw)) deallocate(intw)
          end subroutine destroy_kklist

        end module kklist
