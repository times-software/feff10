      SUBROUTINE define_kpath(BRAVAIS,KPATH,NKDIR,LBLKDIR,KA,KE,BGX,BGY,BGZ)
! **********************************************************************
! *                                                                    *
! *   generate a path in the reciprocal lattice according to           *
! *   BRAVAIS    and    KPATH                                          *
! *   the path may consists of several segments which do not have      *
! *   to join continously                                              *
! *   NKDIR   number of k-path segments created                        *
! *   LBLKDIR name of each segment - contains name of first, last      *
! *           and intermediate k-point (symmetry point in BZ)          *
! *                                                                    *
! *   in multiples of  2 * PI/A                                        *
! *                                                                    *
! *--------------------------------------------------------------------*
! *  BRAVAIS                             KPATH    NKDIR   LBLKDIR      *
! *                                                                    *
! *      1 triclinic   primitive                                       *
! *      2 monoclinic  primitive                                       *
! *      3 monoclinic  base centered                                   *
! *      4 orthorombic primitive                                       *
! *      5 orthorombic base-centered                                   *
! *      6 orthorombic body-centered                                   *
! *      7 orthorombic face-centered                                   *
! *      8 tetragonal  primitive                                       *
! *      9 tetragonal  body-centered                                   *
! *     10 trigonal    primitive                                       *
! *     11 hexagonal   primitive                                       *
! *     12 cubic       primitive                                       *
! *     13 cubic       face-centered                                   *
! *     14 cubic       body-centered                                   *
! *                                                                    *
! **********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
! PARAMETER definitions
      REAL*8 R0,R1,R2,R3,R4,R5,R8
      PARAMETER (R0=0.0D0,R1=1.0D0,R2=2.0D0,R3=3.0D0,R4=4.0D0,R5=5.0D0,R8=8.0D0)
	  integer,parameter :: nsegmax=20
! Dummy arguments
      INTEGER, intent(in)  :: BRAVAIS
	  integer,intent(inout):: KPATH
	  integer, intent(out) :: NKDIR
      REAL*8, intent(in)   :: BGX(3),BGY(3),BGZ(3)
	  real*8, intent(inout) :: KA(3,nsegmax),KE(3,nsegmax)
      CHARACTER*8 LBLKDIR(nsegmax)
! Local variables
      REAL*8 AVEC(3),BVEC(3),CVEC(3),DVEC(3),EVEC(3),G1VEC(3),G2VEC(3), &
            G3VEC(3),GAMV(3),HVEC(3),KVEC(3),LVEC(3),MVEC(3),NVEC(3), &
            PVEC(3),RVEC(3),SVEC(3),TVEC(3),UVEC(3),WVEC(3),XVEC(3), &
            YVEC(3),ZVEC(3)
      INTEGER I,ID

!      G1VEC(1) = BGX(1)
!      G1VEC(2) = BGY(1)
!      G1VEC(3) = BGZ(1)
!      G2VEC(1) = BGX(2)
!      G2VEC(2) = BGY(2)
!      G2VEC(3) = BGZ(2)
!      G3VEC(1) = BGX(3)
!      G3VEC(2) = BGY(3)
!      G3VEC(3) = BGZ(3)
!KJ changed the meaning of the dummy argument:
!   now given the reciprocal basis vectors directly.
	  g1vec = bgx
	  g2vec = bgy
	  g3vec = bgz
!
      CALL LC3RVEC(GAMV,R0,R0,R0,G1VEC,G2VEC,G3VEC,3)
!
      NKDIR = 0
      ID = 0
!
!----------------------------------------------- 1 triclinic   primitive
      IF ( BRAVAIS.EQ.1 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!----------------------------------------------- 2 monoclinic  primitive
      ELSE IF ( BRAVAIS.EQ.2 ) THEN
         CALL LC3RVEC( BVEC,-R1/R2,  R0  ,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( YVEC,  R0  , R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( ZVEC,  R0  ,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( CVEC,  R0  , R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( DVEC,-R1/R2,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( AVEC,-R1/R2, R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( EVEC,-R1/R2, R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
!
         STOP '<KDIRTAB>: Bravais lattice not treated '
!------------------------------------------- 3 monoclinic  base centered
      ELSE IF ( BRAVAIS.EQ.3 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!----------------------------------------------- 4 orthorombic primitive
      ELSE IF ( BRAVAIS.EQ.4 ) THEN
         CALL LC3RVEC( YVEC,  R0  , R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( XVEC, R1/R2,  R0  ,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( ZVEC,  R0  ,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( UVEC, R1/R2,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( TVEC,  R0  , R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( SVEC, R1/R2, R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( RVEC, R1/R2, R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
!
         IF ( KPATH.EQ.1 ) NKDIR = 12
         IF ( KPATH.EQ.2 ) NKDIR = 7
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) THEN
            NKDIR = 3
            GOTO 50
         END IF
         IF ( KPATH.EQ.5 ) NKDIR = 1
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 1
            GOTO 50
         END IF
         IF ( KPATH.EQ.7 ) THEN
            NKDIR = 2
            GOTO 50
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-X '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -G -U '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -A -Z '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,ZVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'Z -GL-GG'
         CALL DCOPY(3,ZVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
 50      CONTINUE
         IF ( KPATH.EQ.7 ) THEN
            ID = ID + 1
            LBLKDIR(ID) = 'X -GS-GG'
            CALL DCOPY(3,XVEC,1,KA(1,ID),1)
            CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-Y '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,YVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'Y -H -T '
         CALL DCOPY(3,YVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'T -B -Z '
         CALL DCOPY(3,TVEC,1,KA(1,ID),1)
         CALL DCOPY(3,ZVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -D -S '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,SVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'S -C -Y '
         CALL DCOPY(3,SVEC,1,KA(1,ID),1)
         CALL DCOPY(3,YVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -P -R '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,RVEC,1,KE(1,ID),1)
         ID = ID + 1  
         LBLKDIR(ID) = 'R -E -T '
         CALL DCOPY(3,RVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'S -Q -T '
         CALL DCOPY(3,SVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
!------------------------------------------- 5 orthorombic base-centered
      ELSE IF ( BRAVAIS.EQ.5 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!------------------------------------------- 6 orthorombic body-centered
      ELSE IF ( BRAVAIS.EQ.6 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!------------------------------------------- 7 orthorombic face-centered
      ELSE IF ( BRAVAIS.EQ.7 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!----------------------------------------------- 8 tetragonal  primitive
      ELSE IF ( BRAVAIS.EQ.8 ) THEN
         CALL LC3RVEC( MVEC, R1/R2, R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( ZVEC,  R0  ,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( AVEC, R1/R2, R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( RVEC,  R0  , R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( XVEC,  R0  , R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
!
         STOP '<KDIRTAB>: Bravais lattice not treated '
!------------------------------------------- 9 tetragonal  body-centered
      ELSE IF ( BRAVAIS.EQ.9 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!---------------------------------------------- 10 trigonal    primitive
      ELSE IF ( BRAVAIS.EQ.10 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
!---------------------------------------------- 11 hexagonal   primitive
      ELSE IF ( BRAVAIS.EQ.11 ) THEN
         CALL LC3RVEC( MVEC,  R0  , R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( AVEC,  R0  ,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( LVEC,  R0  , R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( KVEC,-R1/R3, R2/R3,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( HVEC,-R1/R3, R2/R3, R1/R2, G1VEC,G2VEC,G3VEC,3)
!
         IF ( KPATH.EQ.1 ) NKDIR = 9
         IF ( KPATH.EQ.2 ) NKDIR = 7
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) NKDIR = 1
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 1
            GOTO 100
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-M '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -T''-K '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,KVEC,1,KE(1,ID),1)
 100     CONTINUE
         ID = ID + 1
         LBLKDIR(ID) = 'K -T -GG'
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-A '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,AVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'A -R -L '
         CALL DCOPY(3,AVEC,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -S''-H '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'H -S -A '
         CALL DCOPY(3,HVEC,1,KA(1,ID),1)
         CALL DCOPY(3,AVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -U -L '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'K -P -H '
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
!---------------------------------------------- 12 cubic       primitive
      ELSE IF ( BRAVAIS.EQ.12 ) THEN
         CALL LC3RVEC( XVEC,  R0  , R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( MVEC, R1/R2, R1/R2,  R0  , G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( RVEC, R1/R2, R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
!
         NKDIR = 1
         IF ( KPATH.EQ.1 ) NKDIR = 5
         IF ( KPATH.EQ.2 ) NKDIR = 4
         IF ( KPATH.EQ.3 ) NKDIR = 3
         IF ( KPATH.EQ.4 ) NKDIR = 2
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 3
            DO ID=1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL DCOPY(3,R0,0,KE(1,ID),1)
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 0.5D0            
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 0.5D0            
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 0.5D0            
            ID = NKDIR
         END IF
!
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-X '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -Y -M '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -V -R '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,RVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'R -GL-GG'
         CALL DCOPY(3,RVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-M '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
!------------------------------------------ 13 cubic       face-centered
      ELSE IF ( BRAVAIS.EQ.13 ) THEN
         CALL LC3RVEC( XVEC, R1/R2,  R0  , R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( LVEC, R1/R2, R1/R2, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( WVEC, R3/R4, R1/R4, R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( KVEC, R3/R4, R3/R8, R3/R8, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( UVEC, R5/R8, R1/R4, R5/R8, G1VEC,G2VEC,G3VEC,3)
!
         IF ( KPATH.EQ.0 ) KPATH = 4
         IF ( KPATH.EQ.1 ) NKDIR = 9
         IF ( KPATH.EQ.2 ) NKDIR = 5
         IF ( KPATH.EQ.3 ) NKDIR = 2
         IF ( KPATH.EQ.4 ) THEN
            NKDIR = 1
            ID = ID + 1
            LBLKDIR(ID) = 'GG-GD-X '
            CALL DCOPY(3,GAMV,1,KA(1,ID),1)
            CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         END IF
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 1
            GOTO 150
         END IF
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 3
            DO ID=1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL DCOPY(3,R0,0,KE(1,ID),1)
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 1D0            
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 1D0            
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 1D0            
            ID = NKDIR
         END IF
!
         ID = ID + 1
         LBLKDIR(ID) = 'X -GD-GG'
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
 150     CONTINUE
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GL-L '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -Q -W '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,WVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'W -N -K '
         CALL DCOPY(3,WVEC,1,KA(1,ID),1)
         CALL DCOPY(3,KVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'K -GS-GG'
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -M -U '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -S -X '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -Z -W '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,WVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'W -B -U '
         CALL DCOPY(3,WVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
!------------------------------------------ 14 cubic       body-centered
      ELSE IF ( BRAVAIS.EQ.14 ) THEN
         CALL LC3RVEC( HVEC, R1/R2, R1/R2,-R1/R2, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( PVEC, R1/R4, R1/R4, R1/R4, G1VEC,G2VEC,G3VEC,3)
         CALL LC3RVEC( NVEC, R1/R2,  R0  ,  R0  , G1VEC,G2VEC,G3VEC,3)
!
         IF ( KPATH.EQ.0 ) KPATH = 5
         IF ( KPATH.EQ.1 ) NKDIR = 6
         IF ( KPATH.EQ.2 ) NKDIR = 5
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) NKDIR = 3
         IF ( KPATH.EQ.5 ) NKDIR = 1
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 3
            DO ID=1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL DCOPY(3,R0,0,KE(1,ID),1)
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 1D0            
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 1D0            
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 1D0            
            ID = NKDIR
         END IF
!
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-H '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'H -G -N '
         CALL DCOPY(3,HVEC,1,KA(1,ID),1)
         CALL DCOPY(3,NVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'N -GS-GG'
         CALL DCOPY(3,NVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GL-P '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,PVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'P -F -H '
         CALL DCOPY(3,PVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'N -D -P '
         CALL DCOPY(3,NVEC,1,KA(1,ID),1)
         CALL DCOPY(3,PVEC,1,KE(1,ID),1)
      ELSE
         WRITE (*,*) '<KDIRTAB> called for BRAVAIS=',BRAVAIS
         STOP
      END IF
! ......................................................................
      WRITE (6,99001) KPATH,NKDIR
      WRITE (6,99002) '->G(1): ',G1VEC
      WRITE (6,99002) '->G(2): ',G2VEC
      WRITE (6,99002) '->G(3): ',G3VEC
      WRITE (6,*) ' '
      DO ID = 1,NKDIR
         WRITE (6,99002) LBLKDIR(ID),(KA(I,ID),I=1,3),(KE(I,ID),I=1,3)
      END DO
      WRITE (6,*) ' '
99001 FORMAT (/,1X,79('*'),/,35X,'<KDIRTAB>',/,1X,79('*'),//,5X,'for KPATH =',I2,I10,' k-directions created',/)
99002 FORMAT (5X,A,2X,'(',3F8.4,' )',:,'   ...   (',3F8.4,' )')
      END



      SUBROUTINE LC3RVEC(V,C1,C2,C3,V1,V2,V3,L)
! *     ->V = C1 * ->V1 + C2 * ->V2 + C3 * ->V3                        *

      IMPLICIT REAL*8(A-H,O-Z)
! Dummy arguments
      REAL*8,  intent(in) :: C1,C2,C3
      INTEGER, intent(in) :: L
      REAL*8, intent(out) :: V(L)
	  real*8, intent(in)  :: V1(L),V2(L),V3(L)
! Local variables
      INTEGER I

      DO I = 1,L
         V(I) = C1*V1(I) + C2*V2(I) + C3*V3(I)
      END DO
      END






      subroutine  dcopy(n,dx,incx,dy,incy)
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!        code for both increments equal to 1
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
