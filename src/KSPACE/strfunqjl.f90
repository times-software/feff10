!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strfunqjl.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION STRFUNQJL(FACT,J,L)
!
!   ********************************************************************
!   *                                                                  *
!   *        PRE FACTORS FOR THE REAL SPHERICAL HARMONICS              *
!   *                                                                  *
!   ********************************************************************
!

      IMPLICIT NONE

!
! PARAMETER definitions
!
      REAL*8 PI
      PARAMETER ( PI=3.141592653589793238462643D0 )
!
! Dummy arguments
!
      INTEGER J,L
      REAL*8 FACT(0:100)
      REAL*8 STRFUNQJL
!
! Local variables
!
      REAL*8 WERT
!  /D (10)/
      IF ( J.EQ.0 ) THEN
         WERT = 0.5D0
      ELSE
         WERT = FACT(L-J)/FACT(L+J)
      END IF
      STRFUNQJL = SQRT(WERT*(2.0D0*L+1)/(2.0D0*PI))
!
      END 
