!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cgcrac.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION CGCRAC(FACT,J1,J2,J3,M1,M2,M3)
!   ********************************************************************
!   *                                                                  *
!   *     CLEBSCH GORDAN COEFFICIENTS FOR ARBITRARY                    *
!   *     QUANTUM NUMBERS  J1,J2 ...                                   *
!   *     ACCORDING TO THE FORMULA OF   RACAH                          *
!   *     SEE: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM          *
!   *          EQUATION (3.19)                                         *
!   *          EDMONDS EQ. (3.6.11) PAGE 45                            *
!   *                                                                  *
!   ********************************************************************
!
      IMPLICIT NONE
!
!
! Dummy arguments
!
      REAL*8 J1,J2,J3,M1,M2,M3
      REAL*8 CGCRAC
      REAL*8 FACT(0:100)
!
! Local variables
!
      INTEGER J,N,N1,N2,N3,N4,N5,NBOT,NTOP
      INTEGER NINT
      REAL*8 RFACT
      REAL*8 S,SUM,VF,X,Y
!
! INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
      RFACT(X) = FACT(NINT(X))
!
!
      CGCRAC = 0.0D0
      IF ( ABS(M3-(M1+M2)).GT.1.0D-6 ) RETURN
      IF ( ABS(J1-J2).GT.J3 ) RETURN
      IF ( (J1+J2).LT.J3 ) RETURN
      IF ( ABS(M1).GT.(J1+1.0D-6) ) RETURN
      IF ( ABS(M2).GT.(J2+1.0D-6) ) RETURN
      IF ( ABS(M3).GT.(J3+1.0D-6) ) RETURN
!
      DO J = ABS(NINT(2*(J1-J2))),NINT(2*(J1+J2)),2
         IF ( J.EQ.NINT(2*J3) ) GOTO 100
      END DO
      RETURN
!
!
 100  CONTINUE
      X = (2.0D0*J3+1.0D0)*RFACT(J1+J2-J3)*RFACT(J1-J2+J3)              &
     &    *RFACT(-J1+J2+J3)*RFACT(J1+M1)*RFACT(J1-M1)*RFACT(J2+M2)      &
     &    *RFACT(J2-M2)*RFACT(J3+M3)*RFACT(J3-M3)
!
      Y = RFACT(J1+J2+J3+1)
!
      VF = DSQRT(X/Y)
!
!
      N1 = NINT(J1+J2-J3)
      N2 = NINT(J1-M1)
      N3 = NINT(J2+M2)
      N4 = NINT(J3-J2+M1)
      N5 = NINT(J3-J1-M2)
      NTOP = MIN(N1,N2,N3)
      NBOT = MAX(0,-N4,-N5)
!
      N = NBOT + 1
      IF ( N.EQ.(2*(N/2)) ) THEN
         S = +1.0D0
      ELSE
         S = -1.0D0
      END IF
      SUM = 0.0D0
!
      DO N = NBOT,NTOP
         S = -S
         Y = FACT(N)*FACT(N1-N)*FACT(N2-N)*FACT(N3-N)*FACT(N4+N)        &
     &       *FACT(N5+N)
         SUM = SUM + (S/Y)
      END DO
      CGCRAC = VF*SUM
      END
