module pdw_mod
  use constants
  implicit none
contains

subroutine pdw_vxc(rs, temp, vxc)
  ! Finite-T exchange-correlation potential
  !
  ! rs: average inter-particle spacing (in Bohr)
  ! temp: temperature in Hartree
  !
  ! vxc: (output) exchange-correlation potential
  !
  ! Reference: Perrot, Dharma-Wardana PRA 30, 2619 (1984)
  !   Note: Confident range 0.1 < t < 12
  double precision, intent(in) :: rs, temp
  double precision, intent(out) :: vxc

  double precision t ! unitless temperature scaled by Fermi energy
  double precision efermi

  efermi = 1.0 /  (2 * (4.0/9/pi)**(2/3.0)* rs**2)
  t = temp / efermi

  call pdw_vxc1(rs, t, vxc)
end subroutine

subroutine pdw_vxc1(rs, t, vxc)
  ! Finite-T exchange-correlation potential
  !
  ! rs: average inter-particle spacing (in Bohr)
  ! t: temperature in units of electron gas Fermi energy
  !
  ! vxc: (output) exchange-correlation potential
  !
  ! Reference: Perrot, Dharma-Wardana PRA 30, 2619 (1984)

  double precision, intent(in) :: rs, t
  double precision, intent(out) :: vxc

  double precision efermi

  Integer :: i, ni
  REAL*8 :: ap, t0, rs0

  vxc = pdw_mu_x(rs,t) + pdw_mu_c(rs,t)

  ! rs0 = 2.0d0
  ! t0 = 0.d0
  ! ni = 100
  ! OPEN(UNIT=111, file="pdw_tun.dat", STATUS="REPLACE")
  ! DO i = 1, ni
  !   t0 = DBLE(i-1)*(0.002d0/ni)
  !   ap = pdw_mu_x(rs0, t0) + pdw_mu_c(rs0, t0)
  !   WRITE(111,*), t0, pdw_mu_x(rs, t0), pdw_mu_c(rs, t0)
  ! ENDDO
  ! WRITE(111,*), t, pdw_mu_x(rs,t), pdw_mu_c(rs,t)
  ! CLOSE(111)
end subroutine

function pdw_mu_x(rs,t)
  use constants, only: fa, pi
  ! calculate mu_x(rs,t) / mu_x(rs,0)
  double precision rs, t, pdw_mu_x

  double precision mu_x0

  mu_x0 = -fa / (pi * rs)
  ! print*, "Tanh(1/t)", tanh(1/t), 1-2*EXP(-2*t), t
  if (t < 1d-5) then
    pdw_mu_x = mu_x0
  else
    pdw_mu_x = ((1.0 + 2.83431*t**2 - 0.215120 * t**3 + 5.27586*t**4) / &
                (1 + 3.94309 * t**2 + 7.91379 * t**4)) * tanh(1/t) * mu_x0

  end if

  !print *, "x", pdw_mu_x
end function

function pdw_mu_c(rs,t)
  double precision :: rs, t, pdw_mu_c

  double precision :: c1, c2, c3, c4, mu_c0, mu_ch

  mu_c0 = -0.02545 * log(1.0d0 + 19.0d0 / rs)
  if (t == 0) then
    pdw_mu_c = mu_c0

  else
    c1 = 9.55432/(1.0d0 + 0.06666 * rs)
    c2 = (3.57912 - 5.99065 * rs**0.25 + 1.29722 * rs**0.75) / (1.0d0 + 1.61126*rs**0.25)
    c3 = 4.80217 / (1 + 0.423387 * rs**0.5)
    c4 = 0.29335 + 0.322565 * rs**0.5

    mu_ch = -0.638168 * sqrt(t/rs) * tanh(1/t)

    pdw_mu_c = mu_c0 * (1.0d0 + c1 * t + c2 * t**0.25d0) * exp(-c3 * t) + &
               mu_ch * exp(-c4 / t)

  end if

  !print *, "c", pdw_mu_c
  !print *, c1, c2, c3, c4
end function
end module
