!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: straa.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRAA(NQ,FACT,IG123)

!   ********************************************************************
!   *                                                                  *
!   *  CALCULATE ALL QUANTITIES INDEPENDENT OF  ENERGY  AND  K         *
!   *                THIS ROUTINE IS CALLED ONLY ONCE                  *
!   *                                                                  *
!   ********************************************************************
!
      use boundaries, NMARR_FINAL=>NMARR  
	  use controls,only: IPRINT
	  use workstrfacs, m_nl=>nl  !The "NL" below is a local variable, I think
	  use workstrfacs2
      implicit none
 
! PARAMETER definitions
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
      COMPLEX*16 CI
      PARAMETER (CI=(0.0D0,1.0D0))
 
 ! Dummy arguments
      INTEGER, intent(in) :: NQ,IG123 
      REAL*8, intent(in) :: FACT(0:100)
! Local variables
      REAL*8 ATVOL,D1TERM1,F,FAKTOR,GX,GY,GZ,Q1,Q2,Q3,RS,RSQUAD,RX,RY,RZ,WZ05
      DOUBLE PRECISION DSQRT
      INTEGER I,I1,IQQP,J,J22,JJ,LL,MM,MMLL,N,NL,S,SMAXMAX,SMAXMIN
      REAL*8 RX0,RY0,RZ0
      REAL*8,EXTERNAL :: STRCONFRA,STRFUNQJL
	  INTEGER NMARR


	INDR(:,:) = 0
	NMARR=0
!KJ June 2013: replaced NMARR approach.  Now count appropriate value here from G/R arrays and set for later use in strbbdd.
! Replaces old approach of precompiled parameter and throwing error here if exceeded.


! ----------------------------------------------------------------------
! pre factors for the  spherical harmonics
! used in   <STRAA> + <STRHARPOL>
      WZ05 = SQRT(1.0D0/2.0D0)
      DO LL = 0,LLMAX
         DO JJ = 0,LL
!                                                       REAL
            QJLTAB(JJ,LL) = STRFUNQJL(FACT,JJ,LL)
!                                                       COMPLEX
            CQMLTAB(-JJ,LL) = STRFUNQJL(FACT,JJ,LL)*WZ05
            CQMLTAB(+JJ,LL) = STRFUNQJL(FACT,JJ,LL)*WZ05*(-1.0D0)**JJ
            IF ( JJ.EQ.0 ) CQMLTAB(0,LL) = STRFUNQJL(FACT,0,LL)
         END DO
      END DO

! ----------------------------------------------------------------------
! calculation of the argument-independent terms of the
! harmonic polynomials     used in   <STRHARPOL>
!
!     T(L,J)
      T(0,0) = 1.0D0
      T(1,1) = 1.0D0
!
      DO LL = 1,LLMAX
         T(LL,LL) = T(LL-1,LL-1)*(2.0D0*LL-1.0D0)
      END DO
!     T(LL,LL) = (2*LL-1)!!
!
      HP(1) = QJLTAB(0,0)*T(0,0)
!
!
!  =====================================================================
!                      ********
!                      * DLM1 *
!                      ********
!
      RX = BRY(2)*BRZ(3) - BRZ(2)*BRY(3)
      RY = BRZ(2)*BRX(3) - BRX(2)*BRZ(3)
      RZ = BRX(2)*BRY(3) - BRY(2)*BRX(3)
!
      ATVOL = ABS(BRX(1)*RX+BRY(1)*RY+BRZ(1)*RZ)*(2*PI)**3
!
      D1TERM1 = -4.0D0*PI/ATVOL
!
      G123MAX = 0
      DO N = 1,NMAX
         G123MAX = MAX(G123MAX,ABS(G1(N)),ABS(G2(N)),ABS(G3(N)))
         IF ( G123MAX.GT.NMARR ) THEN
            !WRITE (6,99001) G123MAX,NMARR
            !STOP
			NMARR=G123MAX
         END IF
!
         GX = G1(N)*BGX(1) + G2(N)*BGX(2) + G3(N)*BGX(3)
         GY = G1(N)*BGY(1) + G2(N)*BGY(2) + G3(N)*BGY(3)
         GZ = G1(N)*BGZ(1) + G2(N)*BGZ(2) + G3(N)*BGZ(3)
!
         F = D1TERM1*EXP(-(GX**2+GY**2+GZ**2)/ETA)
!
         DO IQQP = 1,NQQP
            EXPGNQ(N,IQQP) = F*CDEXP(CI*2*PI*(GX*QQPX(IQQP)+GY*QQPY(IQQP)+GZ*QQPZ(IQQP)))
         END DO
      END DO
!
!  =====================================================================
!                      ********
!                      * DLM2 *
!                      ********
!
      Q1 = -SQRT(ETA/PI)*0.5D0
!
      R123MAX = 0
!
      DO IQQP = 1,NQQP
!
! ---> suppress contribution for ->R=0 and IQQP=1 (IQ=IQP)
!      using the fact that the ->R's have been ordered
         IF ( IQQP.EQ.1 ) THEN
!KJ            SMAX(iqqp) = SMAX(iqqp) - 1   !now doing this in strvecgen so it doesn't get executed every time eta changes 2-2012
            I1 = 2
         ELSE
            I1 = 1
         END IF
         S = 0
!
         DO I = I1,NRTAB

            RX0 = R1(I)*BRX(1) + R2(I)*BRX(2) + R3(I)*BRX(3)  - QQPX(IQQP)
            RY0 = R1(I)*BRY(1) + R2(I)*BRY(2) + R3(I)*BRY(3)  - QQPY(IQQP)
            RZ0 = R1(I)*BRZ(1) + R2(I)*BRZ(2) + R3(I)*BRZ(3)  - QQPZ(IQQP)
            RX = 2*PI*RX0
            RY = 2*PI*RY0
            RZ = 2*PI*RZ0
            RSQUAD = RX**2 + RY**2 + RZ**2
            RS = DSQRT(RX0**2+RY0**2+RZ0**2)

            IF ( RS.LE.RMAX ) THEN
!  ---------------------------------------------------------------------
!                         ->R    ACCEPTED
!  ---------------------------------------------------------------------
               R123MAX = MAX(R123MAX,ABS(R1(I)),ABS(R2(I)),ABS(R3(I)))
               IF ( R123MAX.GT.NMARR ) THEN
                  !WRITE (6,99002) R123MAX,NMARR
                  !STOP
				  NMARR=R123MAX
               END IF
               S = S + 1
               INDR(S,IQQP) = I
               CALL STRHARPOL(RX,RY,RZ)
               ! Q2 = (-ETA/2.0D0)**LL
               Q2 = 1.0D0/(-ETA/2.0D0)
               Q3 = EXP(-ETA*RSQUAD/4.0D0)
               MMLL = 0
               DO LL = 0,LLMAX
                  Q2 = Q2*(-ETA/2.0D0)
                  F = Q1*Q2*Q3
                  DO MM = -LL,LL
                     MMLL = MMLL + 1
                     QQMLRS(MMLL,S,IQQP) = F*HP(MMLL)
                  END DO
               END DO
               !  /D (23)/
               DO LL = 0,LLMAX
                  FAKTOR = 1.0D0
                  DO J22 = 0,J22MAX
                     GGJLRS(J22,LL,S,IQQP) = STRCONFRA((LL-J22+0.5D0),(RSQUAD*ETA/4.0D0)) *FAKTOR
                     FAKTOR = FAKTOR/(ETA*(J22+1.0D0))
                  END DO
               END DO
            END IF
         END DO
         IF ( S.NE.SMAX(IQQP) ) THEN
            WRITE (6,*) ' WARNING FROM <STRAA> '
            WRITE (6,*) ' IQQP = ',IQQP
            WRITE (6,*) ' S    = ',S
            WRITE (6,*) ' SMAX = ',SMAX(IQQP)
         END IF
      END DO
!
!  =====================================================================
!                      ********
!                      * DLM3 *
!                      ********
!  /D (8),(13)/
!
      ALPHA0 = SQRT(ETA)/(2.0D0*PI)
!
!  =====================================================================
      IF ( IPRINT.GE.3 ) THEN
         WRITE (6,'(''  LATTICE VECTORS OF THE RECIPROCAL LATTICE'',/)')
         DO N = 1,NMAX
            WRITE (6,99003) N,G1(N),G2(N),G3(N)
         END DO
!
         WRITE (6,'(//,''  LATTICE VECTORS OF THE DIRECT LATTICE '',/)')
         DO S = 1,NRTAB
            WRITE (6,99004) S,R1(S),R2(S),R3(S)
         END DO
!
         WRITE (6,'(//,''  R-INDEX TABLE FOR OFF-DIAG. ELEMENTS  '',/)')
         DO S = 1,NRTAB
            WRITE (6,99005) S,(INDR(S,IQQP),IQQP=1,NQQP)
         END DO
      END IF
!
      SMAXMIN = 1000000
      SMAXMAX = 0
      DO I = 1,NQQP
         SMAXMIN = MIN(SMAX(I),SMAXMIN)
         SMAXMAX = MAX(SMAX(I),SMAXMAX)
      END DO
!



!  =====================================================================
      NMARR_FINAL=NMARR  !KJ set global value to proper value counted here
      NL = LLMAX/2 + 1
      if(iprint.gt.0) WRITE (6,99006) NMAX,NGRLMAX,NRTAB,NRDLMAX0,SMAXMAX,NRDLMAX,      &
     &                SMAXMIN,NRDLMAX,NL,NLMAX,NQ,NQMAX,NQQP,NQQPMAX,   &
     &                R123MAX,NMARR,G123MAX,NMARR,J22MAX,J22MAX,        &
     &                (NL**2*(NL**2+1)/2),LGNT12,IG123,LGNT123,ETA
 
	 RETURN
!  =====================================================================
99001 FORMAT ('  STOP IN <STRAA> ',/,'G123MAX=',I3,' > NMARR=',I3)
99002 FORMAT ('  STOP IN <STRAA> ',/,'R123MAX=',I3,' > NMARR=',I3)
99003 FORMAT ('  ->GN(',I3,')=  (',I3,',',I3,',',I3,')')
99004 FORMAT ('  ->RS(',I3,') =  (',I3,',',I3,',',I3,')')
99005 FORMAT ('  INDR(',I3,'):',30I4)
99006 FORMAT (/,1X,79('*'),/,12X,                                       &
     &        'array sizes for calculation of structure constants',/,1X,&
     &        79('*'),/,10X,32X,'used',14X,'available',/,10X,           &
     &        'G-vectors table       NMAX   ',I7,8X,'NGRLMAX ',I7,/,10X,&
     &        'R-vectors table       NRTAB  ',I7,8X,'NRDLMAX0',I7,/,10X,&
     &        'R-vectors (max.)      SMAX   ',I7,8X,'NRDLMAX ',I7,/,10X,&
     &        'R-vectors (min.)      SMAX   ',I7,8X,'NRDLMAX ',I7,/,10X,&
     &        'l-max + 1             NL     ',I7,8X,'NLMAX   ',I7,/,10X,&
     &        'sites                 NQ     ',I7,8X,'NQMAX   ',I7,/,10X,&
     &        'q-q''-blocks           NQQP  ',I8,8X,'NQQPMAX ',I7,/,10X,&
     &        'auxiliary array        R123MAX',I7,8X,'NMARR   ',I7,/,10X,&
     &        'auxiliary array        G123MAX',I7,8X,'NMARR   ',I7,/,10X,&
     &        'Davis eq. (22)        J22MAX ',I7,8X,'J22MAX  ',I7,/,10X,&
     &        'Gaunts                IGNT12 ',I7,8X,'LGNT12  ',I7,/,10X,&
     &        'Gaunts                IGNT123',I7,8X,'LGNT123 ',I7,/,10X,&
     &        'Ewald parameter       ETA    ',F7.2,/,1X,79('*'),/)
!    &     10X,'lattice sum  evaluated by continued fractions',/,
      END
