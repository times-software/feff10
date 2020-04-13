!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_xstruc.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xstruc
  !--------------------------------------------------------------------
  ! Module for cluster structural information to be used for 
  !   the full multiple scattering calculation
  !
  ! Variables:
  !   xphi:  matrix of angles between the z axis and pairs of atoms
  !   xrat:  xyz coordinates of the atoms in the cluster, the first
  !          npot+1 entries are examples of each unique potential
  !   iphx:  potential indeces of each atom in the cluster, ordered like
  !          xrat
  !
  ! Subroutines:
  !   init_xstruc: Allocates memory for module variables
  !   kill_xstruc: Deallocates memory for module variables
  !--------------------------------------------------------------------
  
  implicit none

  real,    allocatable :: xphi(:,:),xrat(:,:)
  integer, allocatable :: iphx(:)

contains

  subroutine init_xstruc(nclusx)

    implicit none
    integer, intent(in) :: nclusx

    allocate(xphi(nclusx,nclusx))
    allocate(xrat(3,nclusx))
    allocate(iphx(nclusx))

  end subroutine init_xstruc

  subroutine kill_xstruc
    implicit none

    deallocate(xphi,xrat,iphx)
  end subroutine kill_xstruc

end module xstruc
