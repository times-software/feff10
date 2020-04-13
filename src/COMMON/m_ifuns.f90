!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_ifuns.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ifuns
  ! Contains indexing functions for fixvar,fixdsp,fixdsx

  real*8, parameter :: xx00 = 8.8
  
  ! xxx(j) = -xx00 + (j-1)*delta
  ! rrr(j) = exp (-xx00 + (j-1)*delta)
  ! jjj(r) = (log(r) + xx00) / delta + 1

contains

  function xxx(j,delta)

    real*8,  intent(in) :: delta
    integer, intent(in) :: j
    real*8 :: xxx

    xxx = -xx00 + (j-1)*delta

  end function xxx

  function rrr(j,delta)

    real*8,  intent(in) :: delta
    integer, intent(in) :: j
    real*8 :: rrr

    rrr = exp (-xx00 + (j-1)*delta)

  end function rrr

  function jjj(r,delta)

    real*8,  intent(in) :: delta,r
    integer :: jjj

    jjj = (log(r) + xx00) / delta + 1

  end function jjj

end module ifuns
