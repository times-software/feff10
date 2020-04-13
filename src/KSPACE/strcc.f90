!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strcc.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/10/23 18:41:32 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRCC(ERYD,ALAT,RESTART)
!   ********************************************************************
!   *                                                                  *
!   *    CALCULATE ALL QUANTITIES WHICH DEPEND ON THE ENERGY           *
!   *                                                                  *
!   *    the  terms  QQMLRS( MMLL, S, IQQP )  set up in  <STRAA>       *
!   *    are multiplied by  IILERS(LL,S,IQQP)  to save storage         *
!   *    this term is stored and removed at the next call of           *
!   *    <STRCC>                                                       *
!   *                                                                  *
!   ********************************************************************
!
      use boundaries
        use workstrfacs
        use workstrfacs2
        use workstrfacssimple,only: copy_workstrfacssimple
        use controls,only: singleprec

      IMPLICIT NONE


!
! PARAMETER definitions
!
      REAL*8 PI
      PARAMETER ( PI=3.141592653589793238462643D0 )
      COMPLEX*16 CI,C0
      PARAMETER (CI=(0.0D0,1.0D0), C0=(0.0D0,0.0D0))
	  real*8,parameter :: Ewald_terms_threshold=1.d14
	  ! For SnO2 Sn M4 edge FSR, 1.d8 gives much improved but somewhat broadened spectrum at the tail end;
	  ! 1.d30 gives unimproved spectrum where tail is completely wrong due to divergent Ewald summation.
	  ! 1.d20 gives only marginal improvement.
	  ! 1.d15-16 gets most of the correction but still lacking in a few energy points.
	  ! 1.d05-14 gets the complete improvement.
	  
	  

! Dummy arguments
!
      REAL*8 ALAT
      COMPLEX*16 ERYD
      LOGICAL RESTART
!
! Local variables
!
      REAL*8 ALPHA,test
      COMPLEX*16 EHOCHJ,EPWMLLH(0:LLARR),PDU,ealpha,d300b 
      INTEGER ICALL,IQQP,J13,J22,LL,MM,MMLL,S
        logical,parameter :: debug=.false.
      SAVE ICALL
!
      DATA ICALL/0/
!
      IF (RESTART) ICALL = 0
      
      IF (.NOT. RESTART) THEN



700   continue !KJ 1-2012 come here after changing Ewald parameter eta at end of this routine
      EDU = ERYD/(2*PI/ALAT)**2
      PDU = CDSQRT(EDU)
!
!  ===============================================================
!                      ********
!                      * DLM1 *
!                      ********
!  /D (6)/
!
      EPWMLLH(0) = 1.0D0
      D1TERM3(0) = CDEXP(EDU/ETA)
!
      DO LL = 1,LLMAX
         EPWMLLH(LL) = EPWMLLH(LL-1)*CI/PDU
         D1TERM3(LL) = D1TERM3(LL-1)/PDU
      END DO
!
!  ===============================================================
!                      ********
!                      * DLM2 *
!                      ********
!
!  ---------------------------------------------------------------
      IF ( ICALL.EQ.0 ) THEN
         IILERS=C0
      ELSE
!     remove the energy-dependent factor IILERS from the last run
!
         DO IQQP = 1,NQQP
            DO S = 1,SMAX(IQQP)
!
               MMLL = 0
               DO LL = 0,LLMAX
!
                  DO MM = -LL, + LL
                     MMLL = MMLL + 1
                     QQMLRS(MMLL,S,IQQP) = QQMLRS(MMLL,S,IQQP) /IILERS(LL,S,IQQP)
                  END DO
               END DO
            END DO
         END DO
      END IF
!  ---------------------------------------------------------------
!
      ICALL = ICALL + 1
!        open(127,file='testtest.txt',position='append')
!
      DO IQQP = 1,NQQP
         DO S = 1,SMAX(IQQP)
            DO LL = 0,LLMAX
               EHOCHJ = 1.0D0
               IILERS(LL,S,IQQP) = 0.0D0
               DO J22 = 0,J22MAX
                  IILERS(LL,S,IQQP) = IILERS(LL,S,IQQP)   + EHOCHJ*GGJLRS(J22,LL,S,IQQP)
                  EHOCHJ = EHOCHJ*EDU
               END DO
!               test=abs(ehochj/edu *ggjlrs(j22max,ll,s,iqqp))  /abs(iilers(ll,s,iqqp))
!               if(test.gt.0.00000001) then
!                   write(127,*) iqqp,s,ll,j22max
!                   write(127,'(4e12.5)') iilers(ll,s,iqqp),test
!                 endif
!
!     /D (20) AND (22)/
               IILERS(LL,S,IQQP) = IILERS(LL,S,IQQP)*EPWMLLH(LL)
            END DO
         END DO
      END DO

!        close(127)
!
!  ---------------------------------------------------------------
!     multiply the energy-dependent factor IILERS for the
!     current energy  ERYD
!
      DO IQQP = 1,NQQP
         DO S = 1,SMAX(IQQP)
!
            MMLL = 0
            DO LL = 0,LLMAX
!
               DO MM = -LL, + LL
                  MMLL = MMLL + 1
                  QQMLRS(MMLL,S,IQQP) = QQMLRS(MMLL,S,IQQP)  *IILERS(LL,S,IQQP)
               END DO
            END DO
         END DO
      END DO
!
!  ===============================================================
!                      ********
!                      * DLM3 *
!                      ********
!    /D (13)/
!

!      IF (DEBUG) THEN
!         if(dble(eryd).gt.dble(8.0)) then
!	     open(71,file='d300terms.txt')
!         D300 = 0.0D0
!         EHOCHJ = 1.0D000
!         ALPHA = ALPHA0
!	     ealpha = ehochj*alpha
!	     d300b = 0.0d0
!         J13 = -1
! 110     CONTINUE
!         J13 = J13 + 1
!         D300 = ALPHA*EHOCHJ + D300
!	     d300b = ealpha + d300b 
!	     write(71,'(i5,7(e25.16,x))') j13,dble(d300),alpha,dble(ehochj),dble(alpha*ehochj),dble(ealpha),dble(d300b) 
!         ALPHA = ALPHA*(2.0D0*J13-1.0D0)  /(ETA*(J13+1.0D0)*(2.0D0*J13+1.0D0))
!         EHOCHJ = EHOCHJ*EDU
!	     ealpha = ealpha *( edu *(2.0D0*J13-1.0D0)  /(ETA*(J13+1.0D0)*(2.0D0*J13+1.0D0)))
!        prevent floating point underflow
!         IF ( ABS(EHOCHJ).LT.1D-50 ) EHOCHJ = C0
!
!         IF ( CDABS(ALPHA*EHOCHJ/D300).GT.1.0D-10 ) GOTO 110
!         IF ( J13.LT.J13MIN ) GOTO 110
!	     close(71)
!	     write(*,*) 'ERYD IS ',eryd
!	     write(*,*) 'ETA IS ',eta
!	    ! stop
!	     else
!	     write(*,*) 'ERYD IS ',eryd
!	     endif
!	  ENDIF !DEBUG
	  
      D300 = 0.0D0
      EHOCHJ = 1.0D000
      ALPHA = ALPHA0
      J13 = -1
 100  CONTINUE
      J13 = J13 + 1
      D300 = ALPHA*EHOCHJ + D300
      ALPHA = ALPHA*(2.0D0*J13-1.0D0)  /(ETA*(J13+1.0D0)*(2.0D0*J13+1.0D0))
      EHOCHJ = EHOCHJ*EDU
!     prevent floating point underflow
      IF ( ABS(EHOCHJ).LT.1D-50 ) EHOCHJ = C0

      IF ( CDABS(ALPHA*EHOCHJ/D300).GT.1.0D-10 ) GOTO 100
      IF ( J13.LT.J13MIN ) GOTO 100


      if(singleprec) call copy_workstrfacssimple



!KJ 1-2012 If terms are exploding, adjust Ewald parameter to knock them down.
        if ( dble(max(abs(d300),abs(qqmlrs(1,1,1)),abs(d1term3(1)))) .gt. Ewald_terms_threshold ) then
		   call change_eta(eryd)
		   ICALL=0
		   ! recalculate everything
		   goto 700
		endif
		



      END IF
! .NOT.RESTART

      if(.not.restart.and.debug) then
         write(101,*) 'eryd in strcc',eryd
         write(101,*) 'edu',edu
         write(101,*) 'alat',alat
         write(101,*) 'pi',pi
         write(101,*) 'd300',d300
         write(101,*) 'd1term3',d1term3
         write(102,*) 'qqmlrs',qqmlrs
  !      stop
      endif

      END
