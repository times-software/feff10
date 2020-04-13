!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_wavelength.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !ROUTINE: wavelength
! !INTERFACE:
      real*8 function mdff_WaveLength(E)
      use constants
! !INPUT/OUTPUT PARAMETERS:
!   e      : energy in eV.
! !DESCRIPTION:
!   Calculates electron wavelength in a.u. from energy in eV using relativistic formula.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP

      REAL*8,intent(in) :: E

!     H / Sqrt (2 Me) in Sqrt(eV) * A0 (because E is in eV and WL in a.u.
!     = 6.626D-34 / Sqrt(2 * 9.1095D-31 * 1.602189D-19) / 0.52917D-10
      
      mdff_WaveLength = HOnSqrtTwoMe/dsqrt(E + E*E/(dble(2)*MeC2))
      RETURN
      END
