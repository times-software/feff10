!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strvecgen.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRVECGEN(NQ,ALAT)
!   ********************************************************************
!   *                                                                  *
!   *   GENERATE VECTORS OF DIRECT AND RECIPROCAL SPACE FROM           *
!   *   BASIC TRANSLATION VECTORS (BRX,BRY,BRZ)                        *
!   *                                                                  *
!   ********************************************************************

      use energygrid,only : emin,emax
	  use controls,only: IPRINT
	  use boundaries,only : NGRLMAX,NRDLMAX,NRDLMAX0,NQMAX,NQQPMAX,nr,ng,nrdl 
	  use workstrfacs2,only : BRX,BRY,BRZ,BGX,BGY,BGZ,QX,QY,QZ,QQPX, &
            QQPY,QQPZ,G1,G2,G3,R1,R2,R3,init_workstrfacs2_b, &
            init_workstrfacs2_c,exit_workstrfacs2_b,init_workstrfacs2_d
	  use workstrfacs,only : IJQ,NIJQ,NMAX,NQQP,NRTAB,SMAX,RMAX,GMAX,GMAXSQ,init_workstrfacs_b,init_workstrfacs_c
      IMPLICIT REAL*8(A-H,O-Z)
!
! PARAMETER definitions
      REAL*8, parameter :: PI=3.141592653589793238462643D0
	  real*8, parameter :: QQPtolerance=0.00100d0
!
! Dummy arguments
      REAL*8 ALAT 
      INTEGER NQ 
	  
! Local variables
      DOUBLE PRECISION DABS,DBLE,DSQRT
      REAL*8 DD(3),DG(NGRLMAX),DGMIN,DK(3),DQ,DR(NRDLMAX0),DRMIN,DX,    &
     &       EDUMAX,EDUMIN,GA,GMAX1,GX,GY,GZ,KNX,KNY,KNZ,KSQ,P,RA,REDU, &
     &       RMAX1,RX,RY,RZ,SX,SY,SZ,VOL
      REAL*8 DMAX1,DMIN1
      INTEGER I,I1,I2,I3,IE,IETOP,IG,II,INSIDE,IQ,IQQP,IR,ISHIG,ISLOW,J,&
     &        J1,J2,J3,JHIG,JLOW,JQ,K,KK,M,N,NSG(NGRLMAX),NSH,    &
     &        NSR(NRDLMAX),NUMG,NUMGH,NUMR,NUMRH
      INTEGER INT

!
      if(iprint.gt.1) WRITE (6,99001)
!
!     PRIMITIVE VECTORS (BGX,BGY,BGZ) OF RECIPROCAL SPACE
      DO I = 1,3
         I1 = 1 + MOD(I,3)
         I2 = 1 + MOD(I1,3)
         BGX(I) = BRY(I1)*BRZ(I2) - BRZ(I1)*BRY(I2)
         BGY(I) = BRZ(I1)*BRX(I2) - BRX(I1)*BRZ(I2)
         BGZ(I) = BRX(I1)*BRY(I2) - BRY(I1)*BRX(I2)
      END DO
      VOL = DABS(BRX(1)*BGX(1)+BRY(1)*BGY(1)+BRZ(1)*BGZ(1))
!
      DO I = 1,3
         BGX(I) = BGX(I)/VOL
         BGY(I) = BGY(I)/VOL
         BGZ(I) = BGZ(I)/VOL
         if(iprint.gt.1) WRITE (6,99008) BGX(I),BGY(I),BGZ(I)
      END DO

!KJ First calculate the number of unique q_i - q_j :
      call init_workstrfacs2_b(nqqpmax)
      NQQP = 0
      DO IQ = 1,NQ
         DO JQ = 1,NQ
            DO IQQP = 1,NQQP
               IF ( (ABS(QQPX(IQQP)-(QX(IQ)-QX(JQ))).LT.QQPtolerance) .AND.  &
                   (ABS(QQPY(IQQP)-(QY(IQ)-QY(JQ))).LT.QQPtolerance) .AND.  &
                   (ABS(QQPZ(IQQP)-(QZ(IQ)-QZ(JQ))).LT.QQPtolerance) ) THEN

                  GOTO 40
               END IF
            END DO

            NQQP = NQQP + 1
            QQPX(NQQP) = QX(IQ) - QX(JQ)
            QQPY(NQQP) = QY(IQ) - QY(JQ)
            QQPZ(NQQP) = QZ(IQ) - QZ(JQ)
 40      END DO
      END DO
!	  nqqpmax = nqqp  !It's important to get this dimension exactly right because it features in some gigantic arrays
	  call exit_workstrfacs2_b
	  call init_workstrfacs2_b(nqqp)  !reallocate qqpx,qqpy,qqpz with right sizes.
	  call init_workstrfacs_c !allocate ijq,nijq,smax
!KJ end my code to fix nqqpmax	



!
!   calculate radii RA,GA of spheres holding all vectors used in lattice
!   sums. RMAX1 is longest basis vector. GMAX1 is 2* longest vector
!   in Brillouin zone. must be reconsidered in any new applications
!KJ Note: we already know nqqp, I'm leaving the old code here just b/c too lazy to edit it ...
!
      RMAX1 = 0.0D0
      NQQP = 0
!
      DO IQ = 1,NQ
         DO JQ = 1,NQ
            DO IQQP = 1,NQQP
               IF ( (ABS(QQPX(IQQP)-(QX(IQ)-QX(JQ))).LT.QQPtolerance) .AND.  &
     &              (ABS(QQPY(IQQP)-(QY(IQ)-QY(JQ))).LT.QQPtolerance) .AND.  &
     &              (ABS(QQPZ(IQQP)-(QZ(IQ)-QZ(JQ))).LT.QQPtolerance) ) THEN

                  NIJQ(IQQP) = NIJQ(IQQP) + 1
                  IJQ(NIJQ(IQQP),IQQP) = 100*IQ + JQ
                  GOTO 50
               END IF
            END DO
!
            IF ( (NQQP+1).GT.NQQPMAX ) THEN
               WRITE (6,99002) NQQP
               STOP
            END IF
            NQQP = NQQP + 1
            SMAX(NQQP) = 0
            NIJQ(NQQP) = 1
            IJQ(1,NQQP) = 100*IQ + JQ
            QQPX(NQQP) = QX(IQ) - QX(JQ)
            QQPY(NQQP) = QY(IQ) - QY(JQ)
            QQPZ(NQQP) = QZ(IQ) - QZ(JQ)
            DQ = SQRT(QQPX(NQQP)**2+QQPY(NQQP)**2+QQPZ(NQQP)**2)
            IF ( DQ.GE.RMAX1 ) RMAX1 = DQ
!
 50      END DO
      END DO
!
      if(iprint.gt.1) WRITE (6,99009) NQQP


      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,99010)
         DO N = 1,NQQP
! Modified by FDV
! Split line to make compile in Solaris Studio
            WRITE (6,99011) NIJQ(N), QQPX(N), QQPY(N), QQPZ(N), &
              (INT(IJQ(M,N)/100),(IJQ(M,N)-100*INT(IJQ(M,N)/100)),M=1, MIN(5,NIJQ(N)))
            IF ( NIJQ(N).GT.5 ) WRITE (6,99012) (INT(IJQ(M,N)/100),(IJQ(M,N) -100*INT(IJQ(M,N)/100)),M=6, NIJQ(N))
         END DO
      END IF
!
      DO I = 1,3
         DD(I) = DSQRT(BRX(I)**2+BRY(I)**2+BRZ(I)**2)
         DK(I) = DSQRT(BGX(I)**2+BGY(I)**2+BGZ(I)**2)
      END DO
      DRMIN = DMIN1(DD(1),DD(2),DD(3))
      DGMIN = DMIN1(DK(1),DK(2),DK(3))
      GMAX1 = DMAX1(DK(1),DK(2),DK(3))
!
      RA = RMAX + RMAX1*1.001D0
      GA = GMAX + GMAX1
      GA = GMAX
!
      if(iprint.gt.0) WRITE (6,99015) RMAX,GMAX,RMAX1,GMAX1,RA,GA
!
      NUMR = 2*(INT(RA/DRMIN)+1) + 1
      NUMG = 2*(INT(GA/DGMIN)+1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
      if(iprint.gt.0) WRITE (6,99013) NUMR,NUMG
!
! ======================================================================
!                            REAL SPACE
! ======================================================================
!
!  accept a R-point ->R[S] in list if
!
!                  |->R[S] - (->Q[i] - ->Q[i'])| < rmax
!
!  i.e. the center of the convergence sphere of radius rmax  is at
!  (->Q[i] - ->Q[i']) and is therefore different for every [i,i'].
!  the  SMAX(IQQP)  R-points ->R[S] within the sphere are picked
!  out of the list using the pointer INDR(S,IQQP) S=1,..,SMAX(IQQP)
!  INDR(S,IQQP) is set up in <STRAA> and used in <STRBBDD>.
!
      IF ( IPRINT.GT.2 ) WRITE (6,99003)


!KJ fix nr, and the other thing :)
      NR = 0
      DO I1 = -NUMRH, + NUMRH
         DO I2 = -NUMRH, + NUMRH
            DO I3 = -NUMRH, + NUMRH
               RX = I1*BRX(1) + I2*BRX(2) + I3*BRX(3)
               RY = I1*BRY(1) + I2*BRY(2) + I3*BRY(3)
               RZ = I1*BRZ(1) + I2*BRZ(2) + I3*BRZ(3)
               INSIDE = 0
               DO IQQP = 1,NQQP
                  SX = RX - QQPX(IQQP)
                  SY = RY - QQPY(IQQP)
                  SZ = RZ - QQPZ(IQQP)
                  DX = DSQRT(SX*SX+SY*SY+SZ*SZ)
                  IF ( DX.LE.RMAX ) THEN
                     INSIDE = 1
                     SMAX(IQQP) = SMAX(IQQP) + 1
                  END IF
               END DO

               IF ( INSIDE.NE.0 ) THEN
                  NR = NR + 1
               END IF

            END DO
         END DO
      END DO
	  !NRDLMAX0=NR
	  nrdl=0  !NRDLMAX=0
	  do iqqp=1,nqqp
	     if(smax(iqqp).gt.nrdl) nrdl=smax(iqqp)   !nrdlmax) !nrdlmax=smax(iqqp)
	  enddo
	  call init_workstrfacs2_c(nr,nrdl,nqqp)
	  call init_workstrfacs_b  ! iilers
	  SMAX(:)=0  !reset, because we'll count it again right now!
!KJ end my code  


      NR = 0
      DO I1 = -NUMRH, + NUMRH
         DO I2 = -NUMRH, + NUMRH
            DO I3 = -NUMRH, + NUMRH
               RX = I1*BRX(1) + I2*BRX(2) + I3*BRX(3)
               RY = I1*BRY(1) + I2*BRY(2) + I3*BRY(3)
               RZ = I1*BRZ(1) + I2*BRZ(2) + I3*BRZ(3)
               INSIDE = 0
               DO IQQP = 1,NQQP
                  SX = RX - QQPX(IQQP)
                  SY = RY - QQPY(IQQP)
                  SZ = RZ - QQPZ(IQQP)
                  DX = DSQRT(SX*SX+SY*SY+SZ*SZ)
                  IF ( DX.LE.RMAX ) THEN
                     INSIDE = 1
                     SMAX(IQQP) = SMAX(IQQP) + 1
                     IF ( SMAX(IQQP).GT.NRDLMAX ) GOTO 300
                  END IF
               END DO

               IF ( INSIDE.NE.0 ) THEN
                  NR = NR + 1
                  IF ( NR.GT.NRDLMAX0 ) GOTO 100
                  DR(NR) = DSQRT(RX*RX+RY*RY+RZ*RZ)
                  R1(NR) = I1
                  R2(NR) = I2
                  R3(NR) = I3
               END IF

            END DO
         END DO
      END DO
!
! -------------------------------- SORT VECTORS IN ORDER OF INCREASING D
!
      NSH = 1
      NSR(1) = 1
!
      DO II = 2,NR
         I = II - 1
         K = I
         P = DR(I)
!
         DO J = II,NR
            IF ( DR(J).LT.P ) THEN
               K = J
               P = DR(J)
            END IF
         END DO
!
         IF ( K.NE.I ) THEN
            DR(K) = DR(I)
            DR(I) = P
!
            IR = R1(I)
            R1(I) = R1(K)
            R1(K) = IR
            IR = R2(I)
            R2(I) = R2(K)
            R2(K) = IR
            IR = R3(I)
            R3(I) = R3(K)
            R3(K) = IR
         END IF
!
         IF ( I.GT.1 ) THEN
            IF ( ABS(DR(I)-DR(I-1)).GT.1.0D-6 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSR(NSH)
               NSH = NSH + 1
               NSR(NSH) = 1
            ELSE
               NSR(NSH) = NSR(NSH) + 1
            END IF
         END IF
!
         IF ( IPRINT.GT.2 ) WRITE (6,99018) I,R1(I),R2(I),R3(I),DR(I)
!
         IF ( I.EQ.(NR-1) ) THEN
            J = I + 1
            IF ( ABS(DR(J)-DR(J-1)).GT.1.0D-6 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSR(NSH)
               NSH = NSH + 1
               NSR(NSH) = 1
            ELSE
               NSR(NSH) = NSR(NSH) + 1
            END IF
            IF ( IPRINT.GT.2 ) WRITE (6,99018) J,R1(I),R2(I),R3(I),DR(I)
            IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSR(NSH)
         END IF
!
      END DO
!
!
! ======================================================================
!                           RECIPROCAL SPACE
! ======================================================================
!
!  accept a k-point ->G[n] in list if
!
!                     ( (->k+->G[n])**2 - EDU ) < GMAX**2
!
!  i.e. the center of the convergence sphere of radius  GMAX  is at +EDU
!  EDU can vary between  edumin .. edumax  and ->k is within the BZ
!  the set parameter correspond to E=0..3 Ry     !KJ changed them
!

!KJ original SPRKKR statements :
!      EDUMIN = 0.0D0/(2*PI/ALAT)**2
!      EDUMAX = 3.0D0/(2*PI/ALAT)**2
!KJ my new statements  01/2007 :
!KJ I'm choosing here a range that seems to work reasonably for applications I've tested.
!   However, eg. for diamond treated in a primitive cubic box with 8 atoms, one can see
!   the residue shoot up steeply by a factor of 100000 between energies 2.5 and 4 Ry
!   for site-offdiagonal terms of the free structure factor.
!   If a user suspects something fishy may be going on, he can use the multiplicative factors
!   to both the real and reciprocal space sums instead of tampering with the code here, as these
!   also serve to increase the radius of used vectors.
      edumin=abs(2*emin)/(2*PI/ALAT)**2
        edumax=abs(2*emax)/(2*PI/ALAT)**2
        if(iprint.gt.1) write(*,*) 'min',emin,edumin
        if(iprint.gt.1) write(*,*) 'max',emax,edumax

      IETOP = 3
      KK = 1
      GMAXSQ = GMAX**2
      IF ( IPRINT.GT.2 ) WRITE (6,99004)
	  
!KJ  First count the number of vectors needed (NG):
      NG = 0
!
      DO I1 = -NUMGH, + NUMGH
         DO I2 = -NUMGH, + NUMGH
            DO I3 = -NUMGH, + NUMGH
               GX = I1*BGX(1) + I2*BGX(2) + I3*BGX(3)
               GY = I1*BGY(1) + I2*BGY(2) + I3*BGY(3)
               GZ = I1*BGZ(1) + I2*BGZ(2) + I3*BGZ(3)
               INSIDE = 0
!
               DO J1 = -KK, + KK
                  DO J2 = -KK, + KK
                     DO J3 = -KK, + KK
                        KNX = GX + 0.5D0*(J1*BGX(1)+J2*BGX(2)+J3*BGX(3))
                        KNY = GY + 0.5D0*(J1*BGY(1)+J2*BGY(2)+J3*BGY(3))
                        KNZ = GZ + 0.5D0*(J1*BGZ(1)+J2*BGZ(2)+J3*BGZ(3))
                        KSQ = KNX**2 + KNY**2 + KNZ**2
!
                        DO IE = 0,IETOP
                           IF ( IETOP.NE.0 ) THEN
                              REDU = EDUMIN + IE*(EDUMAX-EDUMIN) /DBLE(IETOP)
                           ELSE
                              REDU = 0.0D0
                           END IF
                           IF ( (KSQ-REDU).LE.GMAXSQ ) THEN
                              INSIDE = 1
                              GOTO 61
                           END IF
                        END DO
                     END DO
                  END DO
               END DO

 61            CONTINUE
               IF ( INSIDE.NE.0 ) NG = NG + 1
            END DO
         END DO
      END DO
	  !NGRLMAX=NG
	  call init_workstrfacs2_d(NG,NQQP)
!KJ end	  
	  
      NG = 0
!
      DO I1 = -NUMGH, + NUMGH
         DO I2 = -NUMGH, + NUMGH
            DO I3 = -NUMGH, + NUMGH
               GX = I1*BGX(1) + I2*BGX(2) + I3*BGX(3)
               GY = I1*BGY(1) + I2*BGY(2) + I3*BGY(3)
               GZ = I1*BGZ(1) + I2*BGZ(2) + I3*BGZ(3)
               INSIDE = 0
!
               DO J1 = -KK, + KK
                  DO J2 = -KK, + KK
                     DO J3 = -KK, + KK
                        KNX = GX + 0.5D0*(J1*BGX(1)+J2*BGX(2)+J3*BGX(3))
                        KNY = GY + 0.5D0*(J1*BGY(1)+J2*BGY(2)+J3*BGY(3))
                        KNZ = GZ + 0.5D0*(J1*BGZ(1)+J2*BGZ(2)+J3*BGZ(3))
                        KSQ = KNX**2 + KNY**2 + KNZ**2
!
                        DO IE = 0,IETOP
                           IF ( IETOP.NE.0 ) THEN
                              REDU = EDUMIN + IE*(EDUMAX-EDUMIN) /DBLE(IETOP)
                           ELSE
                              REDU = 0.0D0
                           END IF
                           IF ( (KSQ-REDU).LE.GMAXSQ ) THEN
                              INSIDE = 1
                              GOTO 60
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
!
 60            CONTINUE
               IF ( INSIDE.NE.0 ) THEN
                  NG = NG + 1
                  IF ( NG.GT.NGRLMAX ) GOTO 200
                  DG(NG) = DSQRT(GX*GX+GY*GY+GZ*GZ)
                  G1(NG) = I1
                  G2(NG) = I2
                  G3(NG) = I3
               END IF
            END DO
         END DO
      END DO
!
! -------------------------------- SORT VECTORS IN ORDER OF INCREASING D
!
      NSH = 1
      NSG(1) = 1
!
      DO II = 2,NG
         I = II - 1
         K = I
         P = DG(I)
!
         DO J = II,NG
            IF ( DG(J).LT.P ) THEN
               K = J
               P = DG(J)
            END IF
         END DO
!
         IF ( K.NE.I ) THEN
            DG(K) = DG(I)
            DG(I) = P
!
            IG = G1(I)
            G1(I) = G1(K)
            G1(K) = IG
            IG = G2(I)
            G2(I) = G2(K)
            G2(K) = IG
            IG = G3(I)
            G3(I) = G3(K)
            G3(K) = IG
         END IF
!
         IF ( I.GT.1 ) THEN
            IF ( ABS(DG(I)-DG(I-1)).GT.1.0D-7 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSG(NSH)
               NSH = NSH + 1
               NSG(NSH) = 1
            ELSE
               NSG(NSH) = NSG(NSH) + 1
            END IF
         END IF
!
         IF ( IPRINT.GT.2 ) WRITE (6,99018) I,G1(I),G2(I),G3(I),DG(I)
!
         IF ( I.EQ.(NG-1) ) THEN
            J = I + 1
            IF ( ABS(DG(J)-DG(J-1)).GT.1.0D-7 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSG(NSH)
               NSH = NSH + 1
               NSG(NSH) = 1
            ELSE
               NSG(NSH) = NSG(NSH) + 1
            END IF
            IF ( IPRINT.GT.2 ) WRITE (6,99018) J,G1(I),G2(I),G3(I),DG(I)
            IF ( IPRINT.GT.2 ) WRITE (6,99017) NSH,NSG(NSH)
         END IF
!
      END DO
!
      NMAX = NG
      NRTAB = NR
!
      if(iprint.gt.0) WRITE (6,99014) NRTAB,NMAX
      IF ( NQQP.GT.1 ) THEN
         IF ( IPRINT.GT.0 ) THEN
            WRITE (6,99016) (I,SMAX(I),I=1,NQQP)
         ELSE
            JHIG = 0
            JLOW = 1000000
            DO I = 1,NQQP
               IF ( SMAX(I).GT.JHIG ) THEN
                  ISHIG = I
                  JHIG = SMAX(I)
               END IF
               IF ( SMAX(I).LT.JLOW ) THEN
                  ISLOW = I
                  JLOW = SMAX(I)
               END IF
            END DO
            if(iprint.gt.0)WRITE (6,99016) ISLOW,JLOW
            if(iprint.gt.0)WRITE (6,99016) ISHIG,JHIG
         END IF
      END IF
!

      !KJ This used to be in straa; moved it here so it wouldn't get executed every time ETA changes.  2-2012
	  SMAX(1) = SMAX(1) - 1



      RETURN
!
 100  CONTINUE
      WRITE (6,*) NR,DX,RA,RMAX
      WRITE (6,99005) NRDLMAX0,DX,RA,RMAX
      STOP
 200  CONTINUE
      WRITE (6,99006) NGRLMAX,DX,GA,GMAX
      STOP
 300  CONTINUE
      WRITE (6,*) NR,DX,RA,RMAX
      WRITE (6,99007) NRDLMAX0,SMAX
      STOP
99001 FORMAT (/,10X,'primitive vectors of reciprocal space',/,10X,      &
     &        '       in units of 2*pi/a',/)
99002 FORMAT (/,10X,'STOP in <STRVECGEN> ',/,'NQQP =',I3,' > NQQPMAX')
99003 FORMAT (10X,'result from <STRVECGEN> for real space vectors',//,  &
     &        14X,'NO',5X,'SX',8X,'SY',8X,'SZ',8X,'D',/)
99004 FORMAT (10X,'result from <STRVECGEN> for reciprocal vectors',//,  &
     &        14X,'NO',5X,'KX',8X,'KY',8X,'KZ',8X,'D',/)
99005 FORMAT (' *** NR > NRDLMAX0 = ',I5,/,'   DX,RA,RMAX: ',3F10.5)
99006 FORMAT (' *** NG > NGRLMAX  = ',I5,/,'   DX,GA,GMAX: ',3F10.5)
99007 FORMAT (' *** SMAX > NRDLMAX0 =',I5,/,(10X,14I5))
99008 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
99009 FORMAT (/,10X,'number of inequivalent QQP-block-matrices',        &
     &        '  NQQP =',I3)
99010 FORMAT (10X,'NIJQ (   ->Q[IQ] -  ->Q[JQ]   ) [IQ,JQ] ...')
99011 FORMAT (10X,I3,2X,'(',F7.3,',',F7.3,',',F7.3,' )',1X,             &
     &        5(:,'[',I2,',',I2,']'))
99012 FORMAT ((43X,5(:,'[',I2,',',I2,']')))
99013 FORMAT (10X,'NUMR  =',I10,10X,'NUMG  =',I10)
99014 FORMAT (10X,'NRTAB =',I10,10X,'NMAX  =',I10)
99015 FORMAT (//,10X,'RMAX  =',F10.4,10X,'GMAX  =',F10.4,/,10X,         &
     &        'RMAX1 =',F10.4,10X,'GMAX1 =',F10.4,/,10X,'RA    =',F10.4,&
     &        10X,'GA    =',F10.4,/)
99016 FORMAT (10X,'IQQP=',I4,': SMAX =',I4)
99017 FORMAT (' ',//,15X,'shell number',I5,' with',I5,' points',//)
99018 FORMAT (' ',10X,I5,3I10,F10.6)
      END
