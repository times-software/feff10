!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strbbdd.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/11/23 22:57:43 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRBBDD(DLLMMKE,KX,KY,KZ)
!   ********************************************************************
!   *                                                                  *
!   *  calculate the missing   ->k - dependent quantities              *
!   *  then set up   DLM(K,E)   by performing the lattice sum          *
!   *  the DLM's for the various (IQ,IQ')-blocks of  G  are            *
!   *  evaluated in parallel - as far as possible                      *
!   *                                                                  *
!   ********************************************************************
!
      use boundaries
        use workstrfacs
        use workstrfacs2

      IMPLICIT NONE

! PARAMETER definitions
!
      REAL*8 PI
      PARAMETER ( PI=3.141592653589793238462643D0 )
      COMPLEX*16 CI,C0,CI2PI
      PARAMETER (CI=(0.0D0,1.0D0),C0=(0.0D0,0.0D0),CI2PI=CI*2*PI)
!
! Dummy arguments
!
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQPMAX)
      complex*16,allocatable :: dllmmkeb(:,:,:)  !KJ for debugging  (nllmmmax,nqqpmax,3)      
!
! Local variables
!
      COMPLEX*16 DENOM,EXIK1,EXIK2,EXIK3,EXIKQQP,EXIKRS,F1,F2,F3
      complex*16, allocatable :: PWEXIK1(:),PWEXIK2(:), PWEXIK3(:)
      REAL*8 EX2K1,EX2K2,EX2K3
      REAL*8 EX2KGN,F0,KNX,KNY,KNZ
	  real*8, allocatable :: PWEX2K1(:),PWEX2K2(:),PWEX2K3(:)	  
      INTEGER I,IPW,IQQP,LL,MM,MMLL,N,S
	  
	  logical,parameter :: debug=.false.

      allocate(PWEXIK1(-NMARR:NMARR),PWEXIK2(-NMARR:NMARR),PWEXIK3(-NMARR:NMARR), &
	           PWEX2K1(-NMARR:NMARR),PWEX2K2(-NMARR:NMARR),PWEX2K3(-NMARR:NMARR))
      PWEX2K1(0)=1.0D0
	  PWEX2K2(0)=1.0D0
	  PWEX2K3(0)=1.0D0
      PWEXIK1(0)=1.0D0
	  PWEXIK2(0)=1.0D0
	  PWEXIK3(0)=1.0D0

      if (debug) then
         open(22,file='hohoho.txt')
         write(22,*) 'r123max=',r123max
         write(22,*) 'g123max=',g123max
         write(22,*) 'smax=',smax(1:2)
         write(22,*) 'k= ',kx,ky,kz
         write(22,*) 'd300= ',d300
         write(22,*) 'd1term3= ',d1term3
         do i=1,nl**2
         do ipw=1,min(smax(1),smax(2))
         write(22,'(2i6,4(e14.5,2x))') i,ipw,qqmlrs(i,ipw,1:2)
         enddo
         enddo
         close(22)
         open(22,file='huhaho.txt')
      endif

      dllmmke=c0
      if(debug) then
	     allocate(dllmmkeb(nllmmmax,nqqpmax,3))
		 dllmmkeb=c0
	  endif

!  ===============================================================
!                      ********
!                      * DLM1 *
!                      ********
!
!     TABLE TO CALCULATE         EXP( -2 * ->KK[N] * ->K )
!
      EX2K1 = EXP(-2.0D0*(BGX(1)*KX+BGY(1)*KY+BGZ(1)*KZ)/ETA)
      EX2K2 = EXP(-2.0D0*(BGX(2)*KX+BGY(2)*KY+BGZ(2)*KZ)/ETA)
      EX2K3 = EXP(-2.0D0*(BGX(3)*KX+BGY(3)*KY+BGZ(3)*KZ)/ETA)
!
      DO IPW = 1,G123MAX
         PWEX2K1(IPW) = PWEX2K1(IPW-1)*EX2K1
         PWEX2K1(-IPW) = 1.0D0/PWEX2K1(IPW)
!
         PWEX2K2(IPW) = PWEX2K2(IPW-1)*EX2K2
         PWEX2K2(-IPW) = 1.0D0/PWEX2K2(IPW)
!
         PWEX2K3(IPW) = PWEX2K3(IPW-1)*EX2K3
         PWEX2K3(-IPW) = 1.0D0/PWEX2K3(IPW)
      END DO

     if(debug) then
        write(22,*) 'pwex2k1= ',pwex2k1(0:g123max)
	    write(22,*) 'pwex2k2= ',pwex2k2(0:g123max)
        write(22,*) 'pwex2k3= ',pwex2k3(0:g123max)
        write(22,*) 'RECIPROCAL LATTICE SUM :'
    endif
!
!  ---------------------------------------------------------------
!                                           RECIPROCAL LATTICE SUM
      F0 = EXP(-(KX*KX+KY*KY+KZ*KZ)/ETA)
      if (debug) then
	      write(22,*) 'f0= ',f0
          open(12,file='d1term3.dat',position='append')
	      write(12,1098) edu,d1term3(1)
	      close(12)
		  open(12,file='expgnq.dat',position='append')
	      write(12,1098) edu,expgnq(10,1)
	      close(12)
      endif
 
      DO N = 1,NMAX
!
!     ->K(N) = ->K + ->KK(N)
!
         KNX = KX + G1(N)*BGX(1) + G2(N)*BGX(2) + G3(N)*BGX(3)
         KNY = KY + G1(N)*BGY(1) + G2(N)*BGY(2) + G3(N)*BGY(3)
         KNZ = KZ + G1(N)*BGZ(1) + G2(N)*BGZ(2) + G3(N)*BGZ(3)
!
         DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
!
!
         IF ( DREAL(DENOM).LE.GMAXSQ ) THEN
!  ...............................................................

            EX2KGN = PWEX2K1(G1(N))*PWEX2K2(G2(N))*PWEX2K3(G3(N))
            F1 = F0*EX2KGN/DENOM
            CALL STRHARPOL(KNX,KNY,KNZ)

            DO IQQP = 1,NQQP

               F2 = F1*EXPGNQ(N,IQQP)
               MMLL = 0
               DO LL = 0,LLMAX
                  F3 = F2*D1TERM3(LL)

                  DO MM = -LL, + LL
                     MMLL = MMLL + 1

      if(debug) write(22,*) 'n,knx,kny,knz,denom,ex2kgn,f1= ',n,knx,kny,knz,denom, ex2kgn,f1
	  if(debug) write(22,*) 'iqqp,f2,expgnq,ll,mm,f3,hp= ',iqqp,f2,expgnq(n,iqqp), ll,mm,f3,hp(mmll)

                     DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)  + F3*HP(MMLL)
                  END DO
               END DO
            END DO
         END IF

!  ...............................................................
      END DO


!
!  MULTIPLY THE MISSING PHASE FACTOR    EXP( I ->K * (->Q[I] - ->Q[J]) )
!
      DO IQQP = 2,NQQP
! kvec 2pi / alat   brav* qqpx  ; brav == alat/sqrt2  alat/sqrt2  0 ; etc.
		EXIKQQP = CDEXP( &    !dsqrt(dble(2))* !dble(8.884)*
     &              CI2PI*(KX*QQPX(IQQP)+KY*QQPY(IQQP)+KZ*QQPZ(IQQP)))
        if(debug.and.(iqqp.eq.nqqp))  write(97,1097)kx,ky,kz,qqpx(iqqp),qqpy(iqqp),qqpz(iqqp),exikqqp
1097    format(100f12.4)
		DLLMMKE(:,IQQP) = EXIKQQP*DLLMMKE(:,IQQP)
      END DO
      if(debug) then
         open(12,file='dlm1.txt',position='append')
         write(12,1098) edu,dllmmke(1,1)
         close(12)
         denom=dllmmke(1,1)
         dllmmkeb(:,:,1)=dllmmke
  	  endif
!      
!
!  ---------------------------------------------------------------
!
!
!  ===============================================================
!                      ********
!                      * DLM2 *
!                      ********
!
!     TABLE TO CALCULATE         EXP( I ->K * ->R[S] )
!
      EXIK1 = CDEXP(CI2PI*(BRX(1)*KX+BRY(1)*KY+BRZ(1)*KZ))
      EXIK2 = CDEXP(CI2PI*(BRX(2)*KX+BRY(2)*KY+BRZ(2)*KZ))
      EXIK3 = CDEXP(CI2PI*(BRX(3)*KX+BRY(3)*KY+BRZ(3)*KZ))
!
      DO IPW = 1,R123MAX
         PWEXIK1(IPW) = PWEXIK1(IPW-1)*EXIK1
         PWEXIK1(-IPW) = 1.0D0/PWEXIK1(IPW)
!
         PWEXIK2(IPW) = PWEXIK2(IPW-1)*EXIK2
         PWEXIK2(-IPW) = 1.0D0/PWEXIK2(IPW)
!
         PWEXIK3(IPW) = PWEXIK3(IPW-1)*EXIK3
         PWEXIK3(-IPW) = 1.0D0/PWEXIK3(IPW)
      END DO
      if(debug) then
         write(22,*) 'pwexik1= ',pwex2k1(0:r123max)
         write(22,*) 'pwexik2= ',pwex2k2(0:r123max)
         write(22,*) 'pwexik3= ',pwex2k3(0:r123max)
      endif
!
!  ---------------------------------------------------------------
!                                               DIRECT LATTICE SUM
!
       if(debug) then
         open(12,file='qqmlrs.dat',position='append')
         write(12,1098) edu,qqmlrs(1,10,1)
         close(12)
       endif
      DO IQQP = 1,NQQP
         DO S = 1,SMAX(IQQP)

            I = INDR(S,IQQP)
            EXIKRS = PWEXIK1(R1(I))*PWEXIK2(R2(I))*PWEXIK3(R3(I))

            dllmmke(:,iqqp)=dllmmke(:,iqqp)+exikrs*qqmlrs(:,s,iqqp)

	
         END DO
      END DO
      if(debug) then
         open(12,file='dlm2.txt',position='append')
         write(12,1098) edu,dllmmke(1,1)-denom
         close(12)
         dllmmkeb(:,:,2)=dllmmke-dllmmkeb(:,:,1)
      endif
!  ===============================================================
!                      ********
!                      * DLM3 *
!                      ********
!
!             ADD THE MISSING TERM   D300

      DLLMMKE(1,1) = DLLMMKE(1,1) + D300

      if(debug) then
          open(12,file='d300.txt',position='append')
          write(12,1098) edu,d300
          close(12)
1098      format(1000e12.4)
      endif
      
	  
      deallocate(PWEXIK1,PWEXIK2,PWEXIK3,PWEX2K1,PWEX2K2,PWEX2K3)	   	  
	  
 	  if(debug) then
          open(12,file='strset.txt',position='append')
          write(12,1098) edu,dllmmke(:,1)
          close(12)
          dllmmkeb(:,:,3)=dllmmke-dllmmkeb(:,:,1)-dllmmkeb(:,:,2)
          open(12,file='dllmmkeb.txt',position='append')
          write(12,1099) edu,dllmmkeb(1,1,:)
		  1099 format(2f10.4,1000(f30.5,2x))
          close(12)
      endif       

      END
