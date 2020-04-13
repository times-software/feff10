!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bastrans.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*==bastrmat.f    processed by SPAG 6.05Rc at 11:40 on 21 Jul 2001
      SUBROUTINE BASTRMAT(LMAX,CGC,RC,CREL,RREL,NKMMAX,NKMPMAX)
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    this is a special version of <STRSMAT> passing the            *
!   *    full BASis TRansformation MATrices  RC, CREL and RREL         *
!   *                                                                  *
!   * 13/01/98  HE                                                     *
!   ********************************************************************
!
      IMPLICIT REAL*8(A-H,O-Z)
!
!*** Start of declarations rewritten by SPAG
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
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),                 &
     &           RREL(NKMMAX,NKMMAX)
!
! Local variables
!
      INTEGER I,IKM,J,JP05,K,L,LM,LNR,M,MUEM05,MUEP05,NK,NKM,NLM
      REAL*8 W
!
!*** End of declarations rewritten by SPAG
!

!      write(*,*) 'lmax,nkmmax,nkmpmax',lmax,nkmmax,nkmpmax

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
      crel=c0 !KJ CALL ZCOPY(NKMMAX*NKMMAX,C0,0,CREL,1)

      LM = 0
      DO LNR = 0,LMAX
         DO M = -LNR,LNR
            LM = LM + 1

            IKM = 0
            DO K = 1,NK
               L = K/2
               IF ( 2*L.EQ.K ) THEN
                  JP05 = L
               ELSE
                  JP05 = L + 1
               END IF

               DO MUEM05 = -JP05,(JP05-1)
                  MUEP05 = MUEM05 + 1
                  IKM = IKM + 1

                  IF ( L.EQ.LNR ) THEN
                     IF ( MUEP05.EQ.M ) CREL(LM,IKM) = CGC(IKM,1)
                     IF ( MUEM05.EQ.M ) CREL(LM+NLM,IKM) = CGC(IKM,2)
                  END IF

               END DO
            END DO

         END DO
      END DO
!
! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
      rc=c0 !KJ CALL ZCOPY(NKMMAX*NKMMAX,C0,0,RC,1)
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
!
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,  &
     &           NKMMAX)
!
      END
!*==changerep.f    processed by SPAG 6.05Rc at 11:40 on 21 Jul 2001
      SUBROUTINE CHANGEREP(A,MODE,B,N,M,RC,CREL,RREL,TEXT,LTEXT)
!   ********************************************************************
!   *                                                                  *
!   *   change the representation of matrix A and store in B           *
!   *   according to MODE:                                             *
!   *                                                                  *
!   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
!   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
!   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
!   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
!   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
!   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
!   *                                                                  *
!   *   the non-relat. representations include the  spin index         *
!   *                                                                  *
!   *   for LTEXT > 0 the new matrix  B  is printed                    *
!   *                                                                  *
!   ********************************************************************
!
      IMPLICIT REAL*8(A-H,O-Z)
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
      COMPLEX*16 C1,C0
      PARAMETER (C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
!
! Dummy arguments
!
      INTEGER LTEXT,M,N
      CHARACTER*7 MODE
      CHARACTER*(*) TEXT
      COMPLEX*16 A(M,M),B(M,M),CREL(M,M),RC(M,M),RREL(M,M)
!
! Local variables
!
      INTEGER KEY
      COMPLEX*16 W1(M,M)
!
!*** End of declarations rewritten by SPAG
!
!---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
      IF ( MODE.EQ.'REL>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'REL>CLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'CLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'CLM>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>CLM' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE
         WRITE (*,*) ' MODE = ',MODE
         STOP 'in <ROTATE>  MODE not allowed'
      END IF
!
!KJ      IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-8,6)
!
      END
