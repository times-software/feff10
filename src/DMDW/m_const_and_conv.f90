!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_const_and_conv.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_Const_and_Conv

  use m_Kinds

! Define the constants and convertion factors

!#######################################################################
! Constants

! amu in ev / mass_unit
    real(kind=r08), parameter :: amu = 9.31478e+08_r08

! Speed of light in Ang/ps
    real(kind=r08), parameter :: caps = 2.99792458e+06_r08

! Boltzmann constant in eV/K
    real(kind=r08), parameter :: K_B = 8.617385e-05_r08

! hbar in eV*ps
    real(kind=r08), parameter :: hbar = 6.582122e-04_r08

! hbar*c in eV*Ang
    real(kind=r08), parameter :: hbarc = 1973.27e+00_r08

! Pi
    real(kind=r08), parameter :: pi = 3.1415926535897932384626433e+00_r08

!#######################################################################
! Convertion factors

! Angstrom to au
  real(kind=r08), parameter :: ang2au = 1.88972666351031921149e+00_r08

! N/m to amu/s^2
  real(kind=r08), parameter :: npm2amups2 = 602.214198280e+00_r08

! au of force to N/m
  real(kind=r08), parameter :: auf2npm = 1556.89279161e+00_r08

! AMU to au of mass
  real(kind=r08), parameter :: amu2au = 1822.88848121e+00_r08

! J/mol to eV
  real(kind=r08), parameter :: Jpmol2eV = 96.4853100e+03_r08

!#######################################################################

end module m_Const_and_Conv
