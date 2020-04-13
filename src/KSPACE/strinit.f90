!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strinit.f90,v $:
! $Revision: 1.11 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRINIT(ETOP,NQ,ALAT,FACT,CGC,gmulti,rmulti)

!   ********************************************************************
!   *                                                                  *
!   *               KKR - structure constant routines                  *
!   *                                                                  *
!   *                       ===============                            *
!   *                          <STRINIT>                               *
!   *                       ===============                            *
!   *                                                                  *
!   *  call once to initialize all                                     *
!   *  <STRINIT> calls <VECGEN>, <GAUNT>, <INITS>, <STRAA>             *
!   *  supply:                                                         *
!   *  ARGUMENTS       ETA, RMAX, GMAX, QX,QY,QZ                       *
!   *                                                                  *
!   *                  NT,NQ, NL,NLM,NLMQ, NK,NKM,NLIN, NKKR           *
!   *                     ^   ^                                        *
!   *                  ALAT                                            *
!   *                  ^                                               *
!   *  <STRAA>  calculates all energy- and k-indep. terms of D[L,M]    *
!   *  <INITS>  creates coefficients to transform real non-relat.      *
!   *           G[L,L'] to relat. G[Lam,Lam']                          *
!   *                                                                  *
!   *                                                                  *
!   *                       ===============                            *
!   *                           <STRCC>                                *
!   *                       ===============                            *
!   *                                                                  *
!   *  call once per energy-value                                      *
!   *  <STRC> calls NO other routines                                  *
!   *  supply:                                                         *
!   *                  ERYD,P                                          *
!   *                  ^                                               *
!   *  ERYD is converted to  EDU  used in <STRCC> and <STRBBDD>        *
!   *                                                                  *
!   *                       ===============                            *
!   *                          <STRSET>                                *
!   *                       ===============                            *
!   *  supply:                                                         *
!   *  arguments: DLLMMKE, KX,KY,KZ                                    *
!   *  KX,KY,KZ   in multiples of 2*pi/a  i.e.  in  d.u.'s             *
!   *                                                                  *
!   *  call once per k-point                                           *
!   *  <STRBBDD> is called to get the  D[L,M]'s                        *
!   *  the G[L,L']'s are then set up in a.u.                           *
!   *  conversion is managed by multiplying the gaunts with  2*pi/a    *
!   *  the factor  4*pi  is also included in the gaunts                *
!   *                                                                  *
!   *                       ===============                            *
!   *                          <STRBBDD>                               *
!   *                       ===============                            *
!   *                                                                  *
!   *  call once per k-point                                           *
!   *  supply:                                                         *
!   *  arguments: DLLMMKE, KX,KY,KZ                                    *
!   *  KX,KY,KZ   in multiples of 2*pi/a  i.e.  in  d.u.'s             *
!   *  D[L,M]     is created in  d.u.'s                                *
!   *                                                                  *
!   *  NOTE:   all str-routines work internally in d.u.'s              *
!   *                                                                  *
!   ********************************************************************
!   *                                                                  *
!   *  INITIALIZE CALCULATION OF STRUCTURE CONSTANTS                   *
!   *                                                                  *
!   ********************************************************************
!   *                                                                  *
!   *   ETA    :EWALD PARAMETER                                        *
!   *   RMAX   :RADIUS OF CONVERGENCE SPHERE IN REAL SPACE             *
!   *   GMAX   :RADIUS OF CONVERGENCE SPHERE IN RECIPROCAL SPACE       *
!   *   BRX                                                            *
!   *   BRY    :PRIMITIVE VECTORS IN REAL SPACE  (UNITS OF A)          *
!   *   BRZ                                                            *
!   *   NL     :NUMBER OF L-QUANTUM NUMBERS,LMAX=NL-1                  *
!   *   NQ     :NUMBER OF ATOMS IN THE PRIMITIVE CELL                  *
!   *   QX                                                             *
!   *   QY     :BASIS VECTORS IN REAL SPACE      (UNITS OF A)          *
!   *   QZ                                                             *
!   *                                                                  *
!   ********************************************************************
!   *                                                                  *
!   *   NPT    :NUMBER OF K-POINTS IN IRREDUCIBLE ZONE                 *
!   *   KX                                                             *
!   *   KY     :K-VECTOR     (UNITS OF 2*PI/A)   RECTANGULAR COORD.    *
!   *   KZ                                                             *
!   *                                                                  *
!   ********************************************************************
!   *                                                                  *
!   * NOTE: in contrast to the previous version of <STRINIT> the       *
!   *       subroutine is supplied with ALL necessary structural       *
!   *       information; i.e.          BRX,BRY,BRZ, NQ, QX,QY,QZ       *
!   *                                                                  *
!   *       to avoid conflicts the variable names in the argument list *
!   *       have 'P' added at the beginning                            *
!   *                                                                  *
!   ********************************************************************
!
      use workstrfacs
      use workstrfacs2
      use boundaries
      use controls,only:iprint
      implicit none

! PARAMETER definitions
      REAL*8 PI,SMALL,ALIM
      PARAMETER (PI=3.141592653589793238462643D0,SMALL=1D-15,ALIM=225D0)
!
! Dummy arguments
      REAL*8, intent(in) :: ALAT,ETOP
      REAL*8, intent(in) :: FACT(0:100),CGC(NKMPMAX,2)
      real*8, intent(in) :: rmulti,gmulti
! Local variables
      REAL*8 B,PRETA,VUC,PRETA2
      INTEGER NQ
      INTEGER I,IG123,IQ,LMAX


      if(iprint.gt.0) WRITE (6,99001)
! -------------------------------- set ETA according to AKAI's algorithm
         VUC = BRX(1)*(BRY(2)*BRZ(3)-BRY(3)*BRZ(2)) + BRX(2)      & 
               *(BRY(3)*BRZ(1)-BRY(1)*BRZ(3)) + BRX(3)             &
               *(BRY(1)*BRZ(2)-BRY(2)*BRZ(1))
         VUC = ABS(VUC)
         B = -LOG(SMALL)
         PRETA = 1D0/VUC**(2D0/3D0)/PI
         PRETA = PRETA*0.75D0
         PRETA2 = MAX(PRETA,ABS(ETOP)/ALIM)
      IF (gmax.lt.1.0)  GMAX = SQRT(B*PRETA2)
      IF (rmax.lt.1.0)  RMAX = SQRT(B/PRETA2)/PI
      IF (rmulti.gt.1d-3) rmax=rmax*rmulti  !KJ allows to rescale easily
      IF (gmulti.gt.1d-3) gmax=gmax*gmulti  !KJ
      IF ( ETA.LT.1D-3 ) THEN
         PRETA = 5D-1*SQRT(PRETA)
         ETA = 5D-1*SQRT(PRETA2)
      END IF

      if(eta0.lt.0.01) eta0=eta

!    Remark that, according to the prescription above,
!    The number of r-vectors is directly proportional to the volume of the unit cell
!    The number of g-vectors is inversely proportional to the volume of the unit cell !KJ

! ---------------------- generate the  SMAX  and  NMAX  shortest vectors of real and reciprocal space
      if(.not.allocated(EXPGNQ)) &   !Meaning: only on the first run of STRINIT.  I should really find a cleaner way of doing this.
            CALL STRVECGEN(NQ,ALAT)

! ----------------------------------------- calculate Gaunt coefficients
      LMAX = NL - 1
      LLMAX = 2*LMAX
      CALL STRGAUNT(LMAX,ALAT,FACT,GNT,NGNT,IGNT,CIPWL,NLMAX,IG123,LGNT12,LGNT123)

! ------------------------------------ calculate transformation matrices
      CALL STRSMAT(LMAX,CGC,SRREL,NRREL,IRREL,NKMMAX,NKMPMAX)



      CALL STRAA(NQ,FACT,IG123)





      RETURN
!
99001 FORMAT (/,1X,79('*'),/,35X,'<STRINIT>',/,1X,79('*'),//,10X,'parameters for calculation of structure constants:')
99002 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
99003 FORMAT (/,10X,'primitive vectors for Bravais lattice',/)
      END
