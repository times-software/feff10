!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strharpol.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRHARPOL(X,Y,Z)
!
!   ********************************************************************
!   *                                                                  *
!   *   HARMONIC POLYNOMIAL    R**L * YLM(->R)      ->R = (X,Y,Z)      *
!   *                                                                  *
!   *   WORKING FOR   REAL  AND  COMPLEX  SPHERICAL HARMONICS           *
!   *   SEE: WILLIAMS, JANAK AND MORUZZI                               *
!   *        PHYS.REV. B6 P.4509 (1972)                                *
!   *                                                                  *
!   ********************************************************************
!
      use boundaries
        use workstrfacs
        use workstrfacs2
        use controls,only: cplxylm

      IMPLICIT NONE

! Dummy arguments
!
      REAL*8 X,Y,Z
!
! Local variables
!
      REAL*8 CHP(0:LLARR),F,F1,F10,F2,F20,F3,RSQ,SHP(0:LLARR),XY,ZSQ
      INTEGER I,JJ,JJP1,LL,MM0LL
!

!KJ catch 0 argument (will produce NaN in calculation of T(:,:))
      if(dabs(x)+dabs(y)+dabs(z).lt.dble(0.00000001)) then
           xy=hp(1)
           hp=dble(0)
           hp(1)=xy
           return
        endif
!KJ

      SHP(0)=dble(0) ; CHP(0)=dble(1)
!
      XY = X*X + Y*Y
      ZSQ = Z*Z
      RSQ = XY + ZSQ
!
!     /WJM (B16) (B17)/
!
      DO JJ = 0,(LLMAX-1)
         JJP1 = JJ + 1
         CHP(JJP1) = X*CHP(JJ) - Y*SHP(JJ)
         SHP(JJP1) = X*SHP(JJ) + Y*CHP(JJ)
      END DO
!
!
!     T(L,J)
      T(1,0) = Z
!
!     /WJM (B19)/
!
      F1 = Z
      F2 = 0
      F3 = 1.0D0
!
      DO LL = 1,(LLMAX-1)
         F1 = F1 + Z + Z
         F2 = F2 + RSQ
         F3 = F3 + 1.0D0
         T(LL+1,0) = (F1*T(LL,0)-F2*T(LL-1,0))/F3
      END DO
!
!
!     /WJM (B20)/
!
      IF ( XY.GT.ZSQ ) THEN
!
         F20 = Z/XY
         F10 = 1.0D0 + Z*F20
!
         DO JJ = 0,(LLMAX-1)
            JJP1 = JJ + 1
            F1 = F10*(JJ+JJ+1)
            F2 = F20
!
            DO LL = (JJ+2),LLMAX
               F1 = F1 + F10
               F2 = F2 + F20
               T(LL,JJP1) = F1*T(LL-1,JJ) - F2*T(LL,JJ)
            END DO
         END DO
!
      ELSE
!
         F1 = -XY/Z
         F20 = RSQ/Z
!
         DO LL = 2,LLMAX
            JJ = LL
            F2 = F20*(LL+JJ)
            F3 = LL - JJ
!
            DO I = 1,(LL-1)
               JJ = JJ - 1
               F2 = F2 - F20
               F3 = F3 + 1.0D0
               T(LL,JJ) = (F1*T(LL,JJ+1)+F2*T(LL-1,JJ))/F3
            END DO
         END DO
!
      END IF
!
!     /WJM (B10) - B(12)/
!
!     HP(1) = QJLTAB(0,0) * T(0,0)  IN <STRAA>
!
      MM0LL = 1
      DO LL = 1,LLMAX
         MM0LL = MM0LL + LL + LL
         HP(MM0LL) = QJLTAB(0,LL)*T(LL,0)
         DO JJ = 1,LL
              if(.not.cplxylm) then
!                                                        ------
!                                                         REAL
!                                                        ------
               F = QJLTAB(JJ,LL)*T(LL,JJ)
               HP(MM0LL+JJ) = F*CHP(JJ)
               HP(MM0LL-JJ) = F*SHP(JJ)

              else
!                                                      ---------
!                                                       COMPLEX
!                                                      ---------
!               HP(MM0LL+JJ) = 
!     1			 CQMLTAB( JJ,LL)*T(LL,JJ)*DCMPLX(CHP(JJ),-SHP(JJ))
!               HP(MM0LL-JJ) = 
!     1			 CQMLTAB(-JJ,LL)*T(LL,JJ)*DCMPLX(CHP(JJ), SHP(JJ))
               stop 'activate complex code in strharpol'
!               To do that, the declaration of hp must be changed from real*8 to complex*16!  Otherwise, you will still lose
!               the complex part.
            endif
         END DO
      END DO
!

      END
