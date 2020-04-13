!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strsmat.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRSMAT(LMAX,CGC,SRREL,NRREL,IRREL,NKMMAX,NKMPMAX)
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
!   *                                                                  *
!   * 25/10/95  HE  proper convention of trans. matrix introduced      *
!   ********************************************************************
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX*16 CI,C1,C0
      PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
!
! Dummy arguments
!
      INTEGER LMAX,NKMMAX,NKMPMAX
      REAL*8 CGC(NKMPMAX,2)
      INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
      COMPLEX*16 SRREL(2,2,NKMMAX)
!
! Local variables
!
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),                 &
     &           RREL(NKMMAX,NKMMAX)
      INTEGER I,IKM,J,JP05,K,L,LAM,LM,LNR,LR,M,MUEM05,MUEP05,NK,NKM,NLM,&
     &        NS1,NS2
      REAL*8 W
!
      NK = 2*(LMAX+1) + 1
      NLM = (LMAX+1)**2
      NKM = 2*NLM
!     ===================================================
!     INDEXING:
!     IKM  = L*2*(J+1/2) + J + MUE + 1
!     LM   = L*(L+1)     +     M   + 1
!     ===================================================
!
! ----------------------------------------------------------------------
! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
! ----------------------------------------------------------------------
      CREL=C0 !KJ CALL ZCOPY(NKMMAX*NKMMAX,C0,0,CREL,1)
!
      LM = 0
      DO LNR = 0,LMAX
         DO M = -LNR,LNR
            LM = LM + 1
!
            IKM = 0
            DO K = 1,NK
               L = K/2
               IF ( 2*L.EQ.K ) THEN
                  JP05 = L
               ELSE
                  JP05 = L + 1
               END IF
!
               DO MUEM05 = -JP05,(JP05-1)
                  MUEP05 = MUEM05 + 1
                  IKM = IKM + 1
!
                  IF ( L.EQ.LNR ) THEN
                     IF ( MUEP05.EQ.M ) CREL(LM,IKM) = CGC(IKM,1)
                     IF ( MUEM05.EQ.M ) CREL(LM+NLM,IKM) = CGC(IKM,2)
                  END IF
!
               END DO
            END DO
!
         END DO
      END DO
!
! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
      RC=C0 !KJ CALL ZCOPY(NKMMAX*NKMMAX,C0,0,RC,1)
!
      W = 1.0D0/SQRT(2.0D0)
!
      DO L = 0,LMAX
         DO M = -L,L
            I = L*(L+1) + M + 1
            J = L*(L+1) - M + 1
!
            IF ( M.LT.0 ) THEN
               RC(I,I) = -CI*W
               RC(J,I) = W
               RC(I+NLM,I+NLM) = -CI*W
               RC(J+NLM,I+NLM) = W
            END IF
            IF ( M.EQ.0 ) THEN
               RC(I,I) = C1
               RC(I+NLM,I+NLM) = C1
            END IF
            IF ( M.GT.0 ) THEN
               RC(I,I) = W*(-1.0D0)**M
               RC(J,I) = CI*W*(-1.0D0)**M
               RC(I+NLM,I+NLM) = W*(-1.0D0)**M
               RC(J+NLM,I+NLM) = CI*W*(-1.0D0)**M
            END IF
         END DO
      END DO
!
! ----------------------------------------------------------------------
! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
! ----------------------------------------------------------------------
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,  &
     &              NKMMAX)
!
!     ---------------------------------------------------
!     store the elements of  RREL
!     ---------------------------------------------------
      DO LAM = 1,NKM
         NS1 = 0
         NS2 = 0
!
         DO LR = 1,2*NLM
            IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
               IF ( LR.LE.NLM ) THEN
                  NS1 = NS1 + 1
                  IF ( NS1.GT.2 ) STOP ' IN <STRSMAT>   NS1 > 2'
                  SRREL(NS1,1,LAM) = RREL(LR,LAM)
                  IRREL(NS1,1,LAM) = LR
               ELSE
                  NS2 = NS2 + 1
                  IF ( NS2.GT.2 ) STOP ' IN <STRSMAT>   NS2 > 2'
                  SRREL(NS2,2,LAM) = RREL(LR,LAM)
                  IRREL(NS2,2,LAM) = LR - NLM
               END IF
            END IF
         END DO
!
         NRREL(1,LAM) = NS1
         NRREL(2,LAM) = NS2
      END DO
!
      END
