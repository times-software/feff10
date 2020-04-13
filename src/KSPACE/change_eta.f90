
! eta (the Ewald parameter) is declared in m_workstrfacs.
! In kprep / kpreppot, an optional user input streta is copied to eta
! Then, at the beginning of the module, strinit sets a value for eta (if no user input) and calls straa, which uses eta to calculate E and k independent parts of the Ewald sum.
! Also strcc, which contains energy-dependent data and is called once for each energy point by fmskspace, depends on eta because it uses the results of straa (eta does not
! occur explicitly in strcc).
! Finally, for each energy point and k-point in the FMS algorithm, structurefactor calls strset which calls strbbdd, which again uses eta to calculate the whole Ewald sum for a
! given E and k.

! So, in order to "change eta" for a calculation, one only needs to call strinit again.  Below is essentially a stripped-down version of kprep:

subroutine change_eta(e)

        use struct,lat=>alat
        use workstrfacs2,only: eta
        use boundaries,only: nmuemax,nkmpmax,nkmax		
		implicit none
		
real*8 :: eta_new
  complex*16,intent(in) :: e

  ! LOCALS
  ! For structure factors :
  integer i
  integer,allocatable :: LTAB(:),KAPTAB(:),NMUETAB(:)
  real*8 ALAT,FACT(0:100),etop,gmulti,rmulti  !gmulti,rmulti are dangerous: if anyone ever changes them in kprep/kpreppot, there'll be trouble ...
  real*8,allocatable :: CGC(:,:)
character*512 slog

real*8,parameter :: eta_increase_factor = 1.4d0  ! Changing this from 1.3 to 1.5 didn't matter for SnO2 Sn M4 edge FSR
real*8,parameter :: eta_max = 3.0d0


eta_new = eta * eta_increase_factor
write(slog,'(a,f7.3,a,f7.3,a,f8.4,a)') 'Changing Ewald parameter eta from ',eta,' to ',eta_new,' at energy ',dble(e),'.'
eta=eta_new
call wlog(slog)
if (eta .gt. eta_max) then
! Modified by FDV
! Split line to make compile in Solaris Studio
   write(slog,*) 'This is larger than the hardwired maximum of ', eta_max, &
                 ' .  Please make sure there is nothing wrong in your calculation.  Quitting now.'
   call wlog(slog)
   stop
endif


   !#######################################################################
  !       FIRST:    INITIALIZE DIMENSIONS 
  !#######################################################################


  !* STRUCTURE
  ALAT = lat(1)
  !  Next variables allow to expand the clusters used for the calculation of the structure constants
  !  in real (rmulti) and reciprocal (gmulti) space.
  rmulti = dble(1)
  gmulti = dble(1)

  ! allocate some locals
  allocate(LTAB(NMUEMAX),KAPTAB(NMUEMAX),NMUETAB(NMUEMAX),CGC(NKMPMAX,2) )

  FACT(0) = 1.0D0
  DO I=1,100
     FACT(I) = FACT(I-1) * DBLE(I)
  END DO

  !   ********************************************************************
  !   *                                                                  *
  !   *   rel. quantum numbers    up to                                  *
  !   *                                                                  *
  !   ********************************************************************
  !                     s   p   p   d   d   f   f   g   g   h
  !     DATA LTAB    /  0,  1,  1,  2,  2,  3,  3,  4,  4,  5 /
  !     DATA LBTAB   /  1,  0,  2,  1,  3,  2,  4,  3,  5,  4 /
  !     DATA KAPTAB  / -1,  1, -2,  2, -3,  3, -4,  4, -5,  5 /
  !     DATA NMUETAB /  2,  2,  4,  4,  6,  6,  8,  8, 10, 10 /

  DO I = 1,NMUEMAX
     LTAB(I) = I/2
     IF( 2*LTAB(I) .EQ. I ) THEN
        KAPTAB(I) = LTAB(I)
     ELSE
        KAPTAB(I) = - LTAB(I) - 1
     END IF
     NMUETAB(I) = 2*ABS( KAPTAB(I) )
  END DO

  ! calculate some gaunt symbols before ltab and nmuetab get overwritten for non-sprel calcul.
  CALL CALCCGC( LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX )

           !sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
           !                         structure constants
           !sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! ETOP influences eta,gmax,rmax in strinit if those values are very small.
! However, since they have already been set before (from kprep/kpreppot), those instructions will not be executed.
! Therefore, at this point ETOP doesn't matter.
ETOP=dble(0)

 CALL STRINIT(ETOP,nats,ALAT,FACT,CGC,gmulti,rmulti ) !KJ eliminating arguments that I think are unnecessary
!CALL STRINIT( ETA,RMAX,GMAX,ETOP,BRX,BRY,BRZ,R1,R2,R3,NRDLTAB,QX,QY,QZ,BGX,BGY,BGZ,G1,G2,G3,NGRLTAB,nats,ALAT,FACT,CGC,gmulti,rmulti )


! kill some locals
deallocate(LTAB,KAPTAB,NMUETAB,CGC)


return
end
	  
