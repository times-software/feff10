!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_lnlm.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lnlm
  !--------------------------------------------------------------------
  ! Module for legendre polynomial normalization constants
  !
  ! Variables:
  !     xnlm:
  !   sigsqr:
  !
  ! Subroutines:
  !   init_lnlm: Allocates memory for module variables
  !   kill_lnlm: Deallocates memory for module variables
  !--------------------------------------------------------------------

  implicit none

  real, allocatable :: xnlm(:,:) 
  real, allocatable :: sigsqr(:,:)

contains

  subroutine init_lnlm(lx,nclusx)

    implicit none
    integer, intent(in) :: lx, nclusx

    allocate(xnlm(0:lx,0:lx))
    allocate(sigsqr(nclusx,nclusx))
  end subroutine init_lnlm

  subroutine kill_lnlm
    implicit none

    deallocate(xnlm,sigsqr)
  end subroutine kill_lnlm

end module lnlm
