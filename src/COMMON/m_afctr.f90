!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_afctr.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module afctr
  !--------------------------------------------------------------------
  ! Module for factorial calculation w/o overflow issues
  !
  ! Variables:
  !     afac:
  !   flzero:
  !      flg:
  !
  ! Subroutines:
  !   TODO: MOVE XFCTST HERE???
  !--------------------------------------------------------------------

  implicit none

  real :: afac, flzero, flg(0:210)

end module afctr
