!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_t3j.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module t3j
  !--------------------------------------------------------------------
  ! Module for saving Clebsch-Gordon coefficients: <LS|J>
  !
  ! Variables:
  !   t3jp:
  !   t3jm:
  !
  ! Subroutines:
  !   init_t3j: Allocates memory for module variables
  !   kill_t3j: Deallocates memory for module variables
  !--------------------------------------------------------------------
  
  implicit none

  real, allocatable, dimension(:,:,:) :: t3jp,t3jm

contains

  subroutine init_t3j(lx)
    implicit none
    integer, intent(in) :: lx

    allocate(t3jp(0:lx,-lx:lx,2))
    allocate(t3jm(0:lx,-lx:lx,2))
  end subroutine init_t3j

  subroutine kill_t3j
    implicit none

   deallocate(t3jp,t3jm)
 end subroutine kill_t3j


end module t3j
           
