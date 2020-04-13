!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: omegaq.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION Omegaq(Wp,q)

  DOUBLE PRECISION,INTENT(IN) :: Wp, q
  DOUBLE PRECISION MM1 ! inverse moment
  DOUBLE PRECISION pi
  PARAMETER (pi = 3.1415926535897932384626433d0)
  COMPLEX*16 xLogx
  EXTERNAL xLogx
  ! Calculate inverse moment. Proportional to zero frequency limit of
  ! 1 minus the Lindhardt inverse dielectric function. Use function xLogx 
  ! to avoid singularity at q = 2.
  MM1 = 4*q*(8*q**2 + 3*wp**2) - 3*(q + 2)*wp**2*(((q - 2)*LOG(2 + q) - DBLE(xLogx(ABS(q - 2)))))
  MM1 = pi/2.d0*(32.d0*q**3/MM1 - 1.d0)

  ! Calculate Omegaq
  Omegaq = SQRT( -(pi*Wp**2 + 2/7*(2.4*q)**2*MM1)/(2*MM1) )

END FUNCTION Omegaq

DOUBLE PRECISION FUNCTION Gamq(gam0,q)
  DOUBLE PRECISION,INTENT(IN) :: gam0, q

  Gamq = SQRT(gam0**2 + 2.4*q**2)

END FUNCTION Gamq
  
  
  
