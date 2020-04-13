!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strgaunt.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRGAUNT(LMAX,ALAT,FACT,GNT,NGNT,IGNT,CIPWL,NLMAX,     &
     &                      IG123,LGNT12,LGNT123)
!
!   ********************************************************************
!   *                                                                  *
!   *    CALCULATION OF GAUNT COEFFICIENTS                             *
!   *                                                                  *
!   *    G(L1,M1;L2,M2;L3,M3) == INT Y(L1,M1) * Y(L2,M2) * Y(L3,M3)    *
!   *                                                                  *
!   *    CPLXYLM = FALSE : REAL SPHERICAL HARMONICS                    *
!   *              TRUE  : COMPLEX SPHERICAL HARMONICS                 *
!   *                      THE COMPLEX CONJUGATE OF Y(L1,M1) IS        *
!   *                      TAKEN                                       *
!   *    SEE EDMONDS EQS. (4.6.3), (3.7.3) AND (3.6.11)                *
!   *                                                                  *
!   *                                                                  *
!   *    TO SET UP THE GLL' - MATRIX FROM THE DLM'S THE GAUNTS ARE     *
!   *    MULTIPLIED BY THE FACTOR 4*PI                                 *
!   *    THE ADDITIONAL FACTOR 2*PI/A CONVERTS THE GLL' - MATRIX       *
!   *    FROM D.U.'S TO A.U.'S                                         *
!   *                                                                  *
!   *    THE GAUNTS ARE CALCULATED ONLY FOR  FOR LM2 <= LM1            *
!   *    FOR THE LOOP LM1=1,LM1MAX / LM2=1,LM1 / LM3=1,LM3MAX          *
!   *    NGNT     SPECIFIES THE NUMBER OF NON-ZERO GNT'S               *
!   *             FOR THE CURRENT LINEAR INDEX I == (LM1,LM2)          *
!   *    IGNT(I)  GIVES THE  LM3  FOR THE I-TH NON-ZERO GNT            *
!   *             THE ARRAYLENGTH  LGNT123  HAS TO BE CHOSEN           *
!   *             ACCORDINGLY   --  see below                          *
!   *                                                                  *
!   ********************************************************************
!
      use controls,only:cplxylm
      IMPLICIT NONE

!
! PARAMETER definitions
!
      REAL*8 PI
      PARAMETER ( PI=3.141592653589793238462643D0 )
      COMPLEX*16 CI
      PARAMETER (CI=(0.0D0,1.0D0))
!
! Dummy arguments
!
      REAL*8 ALAT
      INTEGER LGNT12,LGNT123,LMAX,NLMAX
      COMPLEX*16 CIPWL((2*NLMAX)**2)
      REAL*8 FACT(0:100),GNT(LGNT123)
      INTEGER IGNT(LGNT123),NGNT(LGNT12)
!
! Local variables
!
      COMPLEX*16 CC1,CC2,CC3
      REAL*8 CGAUNT,PRE,RGAUNT,WZ05,XL1,XL2,XL3,XM1,XM2,XM3
      REAL*8 CGCRAC
      INTEGER I,I1,I2,I3,IG123,J1,L,L1,L2,L3,LM,LM1,LM1LM2,LM2,M,M1,M2, &
     &        M3,MM1
      INTEGER NINT

      EXTERNAL CGCRAC


      gnt=dble(0)


      FACT(0) = dble(1)
      DO I = 1,100
         FACT(I) = FACT(I-1)*I
      END DO
!
      LM = 0
      DO L = 0,2*LMAX
         DO M = -L,L
            LM = LM + 1
            CIPWL(LM) = CI**L
         END DO
      END DO
!
      WZ05 = SQRT(0.5D0)
      PRE = 4*PI*2*PI/ALAT
      LM1LM2 = 0
      IG123 = 0


!
      IF ( CPLXYLM ) THEN
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
         DO L1 = 0,LMAX
            DO MM1 = -L1,L1
               M1 = -MM1
               XL1 = L1
               XM1 = M1
               J1 = ABS(M1)
               LM1 = L1*(L1+1) + M1 + 1
!
               DO L2 = 0,LMAX
                  DO M2 = -L2,L2
                     XL2 = L2
                     XM2 = M2
                     LM2 = L2*(L2+1) + M2 + 1
!
!
!---> STORE ONLY GAUNT COEFFIENTS FOR  LM2 <= LM1
!
!        IF( LM2.GT.LM1) GOTO 20
                     LM1LM2 = LM1LM2 + 1
                     NGNT(LM1LM2) = 0
!
                     DO L3 = ABS(L1-L2),(L1+L2)
                        DO M3 = -L3,L3
                           XL3 = L3
                           XM3 = M3

!
                           IF ( (M1+M2+M3).EQ.0 ) THEN
!
                              CGAUNT = SQRT((2*XL1+1)*(2*XL2+1)         &
     &                                 /(4*PI*(2*XL3+1)))*(-1)**(-M3)   &
     &                                 *CGCRAC(FACT,XL1,XL2,XL3,0.0D0,  &
     &                                 0.0D0,0.0D0)                     &
     &                                 *CGCRAC(FACT,XL1,XL2,XL3,XM1,XM2,&
     &                                 -XM3)
!
                              IF ( ABS(CGAUNT).GT.1.0D-8 ) THEN
                                 IG123 = IG123 + 1
                                 IF ( IG123.GT.LGNT123 ) THEN
                                    WRITE (6,*) '   STOP IN <STRGAUNT> '
                                    WRITE (6,*) '   LGNT123 =',LGNT123, &
     &                                 '  TOO SMALL'
                                    WRITE (6,*)                         &
     &                                  '   INCREASE ARRAY - LENGTH  '
                                    WRITE (6,*)                         &
     &                                  '    100  for  l_max = 2     '
                                    WRITE (6,*)                         &
     &                                  '    400  for  l_max = 3     '
                                    WRITE (6,*)                         &
     &                                  '   1200  for  l_max = 4     '
                                    STOP
                                 END IF
                                 NGNT(LM1LM2) = NGNT(LM1LM2) + 1
                                 GNT(IG123) = PRE*CGAUNT*(-1)**J1
                                 IGNT(IG123) = L3*(L3+1) + M3 + 1
                              END IF
                           END IF
!
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      ELSE
!
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
!
         DO L1 = 0,LMAX
            DO M1 = -L1,L1
               XL1 = L1
               LM1 = L1*(L1+1) + M1 + 1
!
               DO L2 = 0,LMAX
                  DO M2 = -L2,L2
                     XL2 = L2
                     LM2 = L2*(L2+1) + M2 + 1
!
!---> STORE ONLY GAUNT COEFFIENTS FOR  LM2 <= LM1
!
                     IF ( LM2.LE.LM1 ) THEN
                        LM1LM2 = LM1LM2 + 1
                        NGNT(LM1LM2) = 0
!
                        DO L3 = ABS(L1-L2),(L1+L2)
                           DO M3 = -L3,L3
                              XL3 = L3
!
                              RGAUNT = 0.0D0
                              DO I1 = 1,2
                                 IF ( I1.EQ.1 ) THEN
                                    XM1 = -ABS(M1)
                                    IF ( M1.LT.0 ) CC1 = CI*WZ05
                                    IF ( M1.EQ.0 ) CC1 = 1.0D0
                                    IF ( M1.GT.0 ) CC1 = WZ05
                                 ELSE
                                    XM1 = ABS(M1)
                                    IF ( M1.LT.0 ) CC1 = -CI*WZ05*(-1)  &
     &                                 **ABS(M1)
                                    IF ( M1.EQ.0 ) CC1 = 0.0D0
                                    IF ( M1.GT.0 ) CC1 = WZ05*(-1)      &
     &                                 **ABS(M1)
                                 END IF
!
                                 DO I2 = 1,2
                                    IF ( I2.EQ.1 ) THEN
                                       XM2 = -ABS(M2)
                                       IF ( M2.LT.0 ) CC2 = CI*WZ05
                                       IF ( M2.EQ.0 ) CC2 = 1.0D0
                                       IF ( M2.GT.0 ) CC2 = WZ05
                                    ELSE
                                       XM2 = ABS(M2)
                                       IF ( M2.LT.0 )                   &
     &                                    CC2 = -CI*WZ05*(-1)**ABS(M2)
                                       IF ( M2.EQ.0 ) CC2 = 0.0D0
                                       IF ( M2.GT.0 ) CC2 = WZ05*(-1)   &
     &                                    **ABS(M2)
                                    END IF
!
                                    DO I3 = 1,2
                                       IF ( I3.EQ.1 ) THEN
                                         XM3 = -ABS(M3)
                                         IF ( M3.LT.0 ) CC3 = CI*WZ05
                                         IF ( M3.EQ.0 ) CC3 = 1.0D0
                                         IF ( M3.GT.0 ) CC3 = WZ05
                                       ELSE
                                         XM3 = ABS(M3)
                                         IF ( M3.LT.0 )                 &
     &                                      CC3 = -CI*WZ05*(-1)**ABS(M3)
                                         IF ( M3.EQ.0 ) CC3 = 0.0D0
                                         IF ( M3.GT.0 ) CC3 = WZ05*(-1) &
     &                                      **ABS(M3)
                                       END IF
!
                                       IF ( (NINT(XM1)+NINT(XM2)+NINT(  &
     &                                    XM3)).EQ.0 ) THEN
!
                                         CGAUNT = SQRT((2*XL1+1)        &
     &                                      *(2*XL2+1)/(4*PI*(2*XL3+1)))&
     &                                      *(-1)**(-NINT(XM3))         &
     &                                      *CGCRAC(FACT,XL1,XL2,XL3,   &
     &                                      0.0D0,0.0D0,0.0D0)          &
     &                                      *CGCRAC(FACT,XL1,XL2,XL3,   &
     &                                      XM1,XM2,-XM3)
!
                                         RGAUNT = RGAUNT +              &
     &                                      DREAL(CC1*CC2*CC3)*CGAUNT
                                       END IF
                                    END DO
                                 END DO
                              END DO
!
                              IF ( ABS(RGAUNT).GT.1.0D-8 ) THEN
                                 IG123 = IG123 + 1
                                 IF ( IG123.GT.LGNT123 ) THEN
                                    WRITE (6,*) '   STOP IN <STRGAUNT> '
                                    WRITE (6,*) '   LGNT123 =',LGNT123, &
     &                                 '  TOO SMALL'
                                    WRITE (6,*)                         &
     &                                  '   INCREASE ARRAY - SIZE'
                                    STOP
                                 END IF
                                 NGNT(LM1LM2) = NGNT(LM1LM2) + 1
                                 GNT(IG123) = PRE*RGAUNT
                                 IGNT(IG123) = L3*(L3+1) + M3 + 1
!                 WRITE(6,'('' GAUNT'',4I4,F20.10)')
!     &           LM1,LM2,IGNT(IG123),IG123, GNT(IG123)
!
                              END IF
!
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
            END DO
         END DO
!
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
!
      END IF
!
!      WRITE (6,90001) LMAX,IG123
90001 FORMAT(/,10X,'for LMAX =',I2,I6,                                  &
     &' non-zero GAUNTS tabulated',/)
      END
