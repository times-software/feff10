!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_kinds.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_Kinds

! Definition of the standard kinds used throughout the program

  integer, parameter :: i04 = selected_int_kind(4)
  integer, parameter :: i08 = selected_int_kind(8)
  integer, parameter :: i16 = selected_int_kind(16)

  integer, parameter :: r04 = selected_real_kind(4)
  integer, parameter :: r08 = selected_real_kind(8)
  integer, parameter :: r16 = selected_real_kind(16)

end module m_Kinds
