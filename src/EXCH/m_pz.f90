!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
! Modified for FEFF by Lazaro Calderin (Mar/2015)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

module pz_mod
implicit none

private
public :: pz_vxc

contains

subroutine pz_vxc( rs, vxc )
  implicit none
  real*8, intent(in)  :: rs 
  ! real*8, intent(in),optional :: xmag
  real*8, intent(out) :: vxc

  real*8 :: exc     

  call lda_xc_pz(rs, exc, vxc)

end subroutine pz_vxc


subroutine lda_xc_pz( rs, exc, vxc )
! Perdew-Zunger LDA

  implicit none

  real*8, intent(in)  :: rs
  real*8, intent(out) :: exc, vxc

  real*8 :: ex, ec, vx, vc

  call slater( rs, ex, vx )
  call pz    ( rs, 1, ec, vc )

  exc = ex + ec
  vxc = vx + vc

end subroutine lda_xc_pz

subroutine slater (rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3
  !
!#ifdef __LIBXC
!  use xc_f90_types_m
!  use xc_f90_lib_m
!#endif
  implicit none
  real*8, intent(in) :: rs
  real*8, intent(out):: ex, vx
!#ifdef __LIBXC  
!  real*8:: rho 
!  real*8, parameter :: pi34 = 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
!  integer :: func_id = 1  ! Slater Exchange
!  integer :: size = 1
!  TYPE(xc_f90_pointer_t) :: xc_func
!  TYPE(xc_f90_pointer_t) :: xc_info
!  
!  rho = (pi34/rs)**3
!  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
!  call xc_f90_lda_exc_vxc(xc_func, size, rho ,ex, vx)  
!  call xc_f90_func_end(xc_func)  
!#else
  real*8, parameter  :: f= -0.687247939924714d0, alpha = 2.0d0/3.0d0
  ! f = -9/8*(3/2pi)^(2/3)
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
!#endif
  !
  return
end subroutine slater
!
!
!-----------------------------------------------------------------------
subroutine pz (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     LDA parameterization from Monte Carlo data
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
!#ifdef __LIBXC
!  use xc_f90_types_m
!  use xc_f90_lib_m
!#endif
  implicit none
  real*8,  intent(in) :: rs
  real*8,  intent(out):: ec, vc
  integer, intent(in)  :: iflag
!#ifdef __LIBXC
!  real*8:: rho 
!  real*8, parameter :: pi34 = 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
!  integer :: func_id = 9   ! Perdew & Zunger
!    integer :: size = 1
!  TYPE(xc_f90_pointer_t) :: xc_func
!  TYPE(xc_f90_pointer_t) :: xc_info
!
!  if (iflag.eq.1)  func_id = 9   ! Perdew & Zunger
!  if (iflag.eq.2)  func_id = 11  ! Ortiz & Ballone (PZ)
!
!  rho = (pi34/rs)**3
!  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
!  call xc_f90_lda_exc_vxc(xc_func, size, rho, ec, vc)  
!  call xc_f90_func_end(xc_func)  
!#else
  real*8 :: a (2), b (2), c (2), d (2), gc (2), b1 (2), b2 (2)
  real*8 :: lnrs, rs12, ox, dox
  !
  data a / 0.0311d0, 0.031091d0 /, b / -0.048d0, -0.046644d0 /, &
       c / 0.0020d0, 0.00419d0 /, d / -0.0116d0, -0.00983d0 /
  data gc / -0.1423d0, -0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 / 0.3334d0, 0.27358d0 /
  !
  if (rs.lt.1.0d0) then
     ! high density formula
     lnrs = log (rs)
     ec = a (iflag) * lnrs + b (iflag) + c (iflag) * rs * lnrs + d ( &
          iflag) * rs
     vc = a (iflag) * lnrs + (b (iflag) - a (iflag) / 3.d0) + 2.d0 / &
          3.d0 * c (iflag) * rs * lnrs + (2.d0 * d (iflag) - c (iflag) ) &
          / 3.d0 * rs
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     ox = 1.d0 + b1 (iflag) * rs12 + b2 (iflag) * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 (iflag) * rs12 + 4.d0 / 3.d0 * &
          b2 (iflag) * rs
     ec = gc (iflag) / ox
     vc = ec * dox / ox
  endif
!#endif
  !
  return
end subroutine pz
!


!! SPIN

!-----------------------------------------------------------------------

subroutine lda_xc_pz_spin( rs, zeta, exc, vxc )
! Perdew-Zunger LDA

  implicit none

  real*8, intent(in)  :: rs, zeta
  real*8, intent(out) :: exc, vxc

  real*8 :: ex, ec, vx, vc

  call slater( rs, ex, vx )
  call pz    ( rs, 1, ec, vc )

  exc = ex + ec
  vxc = vx + vc

end subroutine lda_xc_pz_spin

!-----------------------------------------------------------------------
subroutine pz_polarized (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     spin-polarized energy and potential
  !
!  USE kinds, ONLY : DP
  implicit none
  real(8) :: rs, ec, vc
  real(8) :: a, b, c, d, gc, b1, b2
  parameter (a = 0.01555d0, b = - 0.0269d0, c = 0.0007d0, d = &
       - 0.0048d0, gc = - 0.0843d0, b1 = 1.3981d0, b2 = 0.2611d0)
  real(8) :: lnrs, rs12, ox, dox
  REAL(8), PARAMETER :: xcprefact = 0.022575584d0, pi34 = 0.6203504908994d0 
  ! REAL(8) :: betha, etha, csi, prefact
  !
  if (rs.lt.1.0d0) then
     ! high density formula
     lnrs = log (rs)
     ec = a * lnrs + b + c * rs * lnrs + d * rs
     vc = a * lnrs + (b - a / 3.d0) + 2.d0 / 3.d0 * c * rs * lnrs + &
          (2.d0 * d-c) / 3.d0 * rs
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     ox = 1.d0 + b1 * rs12 + b2 * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 * rs12 + 4.d0 / 3.d0 * b2 * rs
     ec = gc / ox
     vc = ec * dox / ox
  endif
  !
!  IF ( lxc_rel ) THEN
!     betha = prefact * pi34 / rs
!     etha = DSQRT( 1 + betha**2 )
!     csi = betha + etha
!     prefact = 1.0D0 - (3.0D0/2.0D0) * ( (betha*etha - log(csi))/betha**2 )**2
!     ec = ec * prefact
!     vc = vc * prefact
!  ENDIF
  return
end subroutine pz_polarized
!
!
!
subroutine pz_spin (rs, zeta, ec, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
!  USE kinds, ONLY : DP
  implicit none
  real(8) :: rs, zeta, ec, vcup, vcdw
  !
  real(8) :: ecu, vcu, ecp, vcp, fz, dfz
  real(8) :: p43, third
  parameter (p43 = 4.0d0 / 3.d0, third = 1.d0 / 3.d0)
  !
  ! unpolarized part (Perdew-Zunger formula)
  call pz (rs, 1, ecu, vcu)
  ! polarization contribution
  call pz_polarized (rs, ecp, vcp)
  !
  fz = ( (1.0d0 + zeta) **p43 + (1.d0 - zeta) **p43 - 2.d0) / &
       (2.d0**p43 - 2.d0)
  dfz = p43 * ( (1.0d0 + zeta) **third- (1.d0 - zeta) **third) &
       / (2.d0**p43 - 2.d0)
  !
  ec = ecu + fz * (ecp - ecu)
  vcup = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (1.d0 - zeta)
  vcdw = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * ( - 1.d0 - &
       zeta)
  !
  return
end subroutine pz_spin
!
!-----------------------------------------------------------------------
subroutine slater_spin (rho, zeta, ex, vxup, vxdw)
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, spin-polarized case
  !
!  USE kinds, ONLY : DP
  implicit none
  real(8) :: rho, zeta, ex, vxup, vxdw
  real(8) :: f, alpha, third, p43
  parameter (f = - 1.10783814957303361d0, alpha = 2.0d0 / 3.0d0)
  ! f = -9/8*(3/pi)^(1/3)
  parameter (third = 1.d0 / 3.d0, p43 = 4.d0 / 3.d0)
  real(8) :: exup, exdw, rho13
  !
  rho13 = ( (1.d0 + zeta) * rho) **third
  exup = f * alpha * rho13
  vxup = p43 * f * alpha * rho13
  rho13 = ( (1.d0 - zeta) * rho) **third
  exdw = f * alpha * rho13
  vxdw = p43 * f * alpha * rho13
  ex = 0.5d0 * ( (1.d0 + zeta) * exup + (1.d0 - zeta) * exdw)
  !
  return
end subroutine slater_spin

end module pz_mod
