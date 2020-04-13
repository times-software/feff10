!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_rotx.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rotx  
  !--------------------------------------------------------------------
  ! Module for rotation matricies
  !
  ! Parameters:
  !    jsavx = 150    :
  !   jbmagk = -9999  :
  !   roteps = 1.e-12 :
  !
  ! Variables:
  !     drix: Rotation matrix for the entire cluster
  !   drisav:
  !   betsav:
  !    ldsav:
  !    mdsav:
  !     jsav:
  !
  ! Subroutines:
  !   init_rotx: Allocates memory for module variables
  !   kill_rotx: Deallocates memory for module variables
  !--------------------------------------------------------------------

  implicit none

  ! Rotation matrix for the entire cluster
  complex, allocatable :: drix(:,:,:,:,:,:) 
  ! Should drisav,betsav be complex ???

  ! Matrices for saving rotation matrices between xanes and rotxan
  real,    allocatable :: drisav(:,:,:,:)
  integer, parameter   :: jsavx=150, jbmagk=-9999
  real,    parameter   :: roteps = 1.e-12
  integer jsav

  integer :: ldsav(jsavx),mdsav(jsavx)
  real    :: betsav(jsavx)

contains

  subroutine init_rotx(lx,nclusx)
    implicit none
    integer, intent(in) :: lx, nclusx

    allocate(drix(-lx:lx,-lx:lx,0:lx,0:1,nclusx,nclusx))
    allocate(drisav(-lx:lx,-lx:lx,0:lx,jsavx))
  end subroutine init_rotx

  subroutine kill_rotx
    implicit none

   deallocate(drix,drisav)
 end subroutine kill_rotx


end module rotx
           
