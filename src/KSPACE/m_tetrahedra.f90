!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_tetrahedra.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module tetrahedra

!    The number of tetrahedra
      integer ntet
!    The tetrahedra
      integer,allocatable :: itet(:,:)
!    Some work array
      integer,allocatable :: iy(:,:)
!    Some work array
      integer,allocatable :: ittfl(:)
!    Needed for dimensioning ittfl
        integer,parameter :: mwrit=101

      CONTAINS
          subroutine init_tetrahedra(n,m)
            integer n,m
        	allocate(itet(4,n*6),iy(4,6*m),ittfl(5*mwrit))
          end subroutine init_tetrahedra

          subroutine destroy_tetrahedra
           deallocate(itet,iy,ittfl)
          end subroutine destroy_tetrahedra

      end module tetrahedra
