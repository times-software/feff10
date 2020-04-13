!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: calccgc.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
!   ********************************************************************
!   *                                                                  *
!   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
!   *                                                                  *
!   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
!   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
!   *   IS= 1/2  SPIN DOWN/UP                                          *
!   *                                                                  *
!   ********************************************************************
!
      IMPLICIT NONE
!
!
! Dummy arguments
!
      INTEGER NKMAX,NKMPMAX,NMUEMAX
      REAL*8 CGC(NKMPMAX,2)
      INTEGER KAPTAB(NMUEMAX),LTAB(NMUEMAX),NMUETAB(NMUEMAX)
!
! Local variables
!
      INTEGER IKM,K,KAPPA,M
      REAL*8 J,L,MUE,TWOLP1
!
      IKM = 0
      DO K = 1,(NKMAX+1)
         L = LTAB(K)
         KAPPA = KAPTAB(K)
         J = ABS(KAPPA) - 0.5D0
         MUE = -J - 1.0D0
         TWOLP1 = 2.0D0*L + 1.0D0
!
         IF ( KAPPA.LT.0 ) THEN
!
!     J = L + 1/2
            DO M = 1,NMUETAB(K)
!
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L-MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = DSQRT((L+MUE+0.5D0)/TWOLP1)
            END DO
         ELSE
!     J = L - 1/2
            DO M = 1,NMUETAB(K)
!
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L+MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = -DSQRT((L-MUE+0.5D0)/TWOLP1)
!
            END DO
         END IF
!
!
      END DO
!
      END
