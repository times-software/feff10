!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strconfra.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION STRCONFRA(AA,X)
!
!   ********************************************************************
!   *                                                                  *
!   *         CONTINUED FRACTION   /D (18)/                            *
!   *                                                                  *
!   *         perform evaluation of continued fraction                 *
!   *         f(a,x) = [1/x+ (1-a)/1+ 1/x+ (2-a)/1+ ... ]              *
!   *         breaking sequence at (IMAX-a)/1+                         *
!   *                                                                  *
!   ********************************************************************
!

      IMPLICIT NONE

!
! PARAMETER definitions
!
      INTEGER IMAX0
      PARAMETER (IMAX0=100)
!
! Dummy arguments
!
      REAL*8 AA,X
      REAL*8 STRCONFRA
!
! Local variables
!
      INTEGER I,IMAX,J
      REAL*8 WERT,WERT0

      IMAX = IMAX0 - 20
!
 100  CONTINUE
      IMAX = IMAX + 20
!
      WERT = IMAX/X
!
      I = IMAX
      DO J = 2,IMAX
         WERT = 1.0D0 + WERT
         WERT = X + (I-AA)/WERT
         WERT = (I-1)/WERT
         I = I - 1
      END DO
!
      WERT = 1.0D0 + WERT
      WERT = X + (1.0D0-AA)/WERT
      WERT = 1.0D0/WERT
!
      IF ( IMAX.EQ.IMAX0 ) THEN
         WERT0 = WERT
         GOTO 100
      END IF
      IF ( ABS((WERT-WERT0)/WERT).GT.1.0D-10 ) THEN
         WERT0 = WERT
         GOTO 100
      END IF
!
      STRCONFRA = WERT
      END 
