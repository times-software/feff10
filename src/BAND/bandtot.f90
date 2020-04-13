      subroutine bandtot    

  use DimsMod, only: nex, nspx=>nspu, ltot, nphx=>nphu, lx
  use boundaries,only : maxl,msize
  use constants
  use par
  use fms_inp,only: lmaxph,rclust=>rfms2,sigma2=>sig2g
  use global_inp, only: ispin
  use fitting
  use band_inp
  use struct 

  implicit none
  
  character*512 :: slog
  complex*16 ene,refene
  character*6, allocatable  :: potlbl(:)
  integer :: ie !,nsp
  complex*16, allocatable  :: em(:), eref(:,:)
  complex*16 :: dck
  integer :: kinit,linit

  real*8 yr,yi
  real*8,allocatable :: efit(:)
  
! the band structure  
  complex,allocatable :: eival(:),eigen(:,:,:)
  integer nbandsmin,nbandsmax
  integer,allocatable :: n_pos(:,:)
  integer n_pos_eigenval
  real*8,allocatable :: bands(:,:)
  integer,allocatable :: nbands(:)

! reading xsph output
  integer :: i,ik,ik0,lmaxp1,ib,n,i0,isp,ill,ipp
  real*8  :: rnrmav,xmu,edge
  complex*16, allocatable :: ph(:,:,:,:)
  integer, allocatable    :: lmax(:,:)
  complex*16, allocatable :: rkk(:,:,:)
  complex, allocatable, dimension(:,:,:)   :: xphase
  complex ck(nspx)
  real rpart,aipart
  integer ne, ne1, ne3, ihole
  integer :: iz(0:nphx)

! creating the k-path
  integer,parameter :: nsegmax = 20
  character*8 label_k_segments(nsegmax)
  real*8 ka(3,nsegmax),ke(3,nsegmax)
  integer n_k_segments,nkp_segment(nsegmax),indkdir(nsegmax)
  integer bravais
  real*8, allocatable :: bk(:,:),kp(:)  ! The k-mesh
  
  integer,external :: ibravais
  real*8 ebroad,sumd,del(3),sum,a2pi,refere,refeim
  integer iq,it,id,nkp0,ikd,ieverb
  real*8,external :: DNRM2
  logical,parameter :: debug=.true.

  ! Allocate local variables
  allocate(em(nex), eref(nex, nspx))
  allocate(potlbl(0:nphx))
  allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(lmax(nex, 0:nphx))
  allocate(rkk(nex,8,nspx))
  allocate(xphase(nspx, -lx:lx, 0:nphx))


!CCCCCCCCCC READ SCATTERING PHASE SHIFTS CALCULATED BY XSPH : CCCCCCCCCCCCC

      xphase = cmplx(0,0)
	  em=dcmplx(0,0)
	  eref=dcmplx(0,0)
      call rdxsph (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      call setkap (ihole, kinit, linit)
	  
	  write(*,*) 'ne,ne1,ne3',ne,ne1,ne3
!	  write(*,*) 'em',dble(em(:))
	  write(*,*) 'min ',em(1)
	  write(*,*) 'max ',em(ne1)
!	  stop

      call kprep(em,ne,nex,.true.)
      nsp = 1
      if (abs(ispin).eq.1 ) nsp = nspx


!CCCCCCCCCC SET UP ENERGY GRID TO SEARCH FOR BANDS : CCCCCCCCCCCCC
    if(debug) then
        open(99,file='emesh.dat')
        do ie=1,nex
           write(99,'(i4,x,4(e12.4,x,e12.4,5x))') ie,em(ie),eref(ie,:)
        enddo
        close(99)
	endif

    call wlog('Solving band structure.')
!   Cut off grid to range for which phase shifts are available :
	allocate(efit(ne1))
	efit=dble(em(1:ne1))
	emin=emin/hart !eV->Ha
	emax=emax/hart
	estep=estep/hart
	emax=min(emax,dble(em(ne1))-xmu)
	emin=max(emin,dble(em(1))-xmu)
!   Use broadening of em-grid
	ebroad=dimag(em(1))*dble(0)
!	ebroad=dble(0.01)
	write(*,'(a,x,f8.4)') 'Broadening ebroad= ',ebroad
	write(*,'(a,x,f8.4)')'Fermi energy= ',xmu
!   Construct remaining variables :
	nep=nint((emax-emin)/estep)+1
	estep=(emax-emin)/dble(nep-1)
	write(*,'(2(a,x,f9.3,x),a,x,i5,x,a,x,f9.3)') 'Energy range ',emin*hart,' to ', emax*hart,' eV using ',nep,'steps of size',estep*hart


!CCCCCCCCCC SET UP K-SPACE GRID TO SEARCH FOR BANDS : CCCCCCCCCCCCC

    a2pi= alat(1)/(2*pi)
    n_k_segments=0 ; ka=0.d0 ; ke=0.d0 ; label_k_segments=' '
	bravais=ibravais(sgroup,lattice)
    call define_kpath(bravais,ikpath,n_k_segments,label_k_segments,ka,ke,b1*a2pi,b2*a2pi,b3*a2pi)
	write(*,*) 'ibravais ',bravais
	write(*,*) 'kpaths: ',n_k_segments,label_k_segments(1:n_k_segments)


      IF ( n_k_segments.EQ.0 ) THEN
         n_k_segments = 1
		 ka(:,:)=0.d0
		 ke(:,:)=0.d0
         KE(1,1) = 1.0D0
         nkp_segment(1) = nkp
         label_k_segments(1) = 'GG-x -1 '
      ELSE
         SUMD = 0.0D0
         DO ID = 1,n_k_segments
            DO I = 1,3
               DEL(I) = KE(I,ID) - KA(I,ID)
            END DO
            SUM = DNRM2(3,DEL,1)
            nkp_segment(ID) = INT(1000000*SUM)
            SUMD = SUMD + SUM
         END DO
         nkp0 = 0
         DO ID = 1,n_k_segments
            nkp_segment(ID) = INT(nkp*nkp_segment(ID)/nint(SUMD*1000000))
            nkp_segment(ID) = MAX(2,nkp_segment(ID))
            nkp0 = nkp0 + nkp_segment(ID)
         END DO

         nkp = nkp0
      END IF

      INDKDIR(1) = nkp_segment(1)
      DO ID = 2,n_k_segments
         INDKDIR(ID) = INDKDIR(ID-1) + nkp_segment(ID)
      END DO

    allocate(bk(3,nkp),kp(nkp))
    bk=0.d0

!----------------------------------- set up k-points and linear array KP
      IK = 0
      KP(1) = 0D0
      DO ID = 1,n_k_segments

         DO I = 1,3
            DEL(I) = (KE(I,ID)-KA(I,ID))/DBLE(nkp_segment(ID)-1)
         END DO

         DO IKD = 1,nkp_segment(ID)
            IK = IK + 1
            DO I = 1,3
               bk(I,IK) = KA(I,ID) + DEL(I)*DBLE(IKD-1)
            END DO
            IF ( IKD.NE.1 ) THEN
               KP(IK) = KP(IK-1) + DNRM2(3,DEL,1)
            ELSE IF ( IK.EQ.1 ) THEN
               KP(IK) = 0.0D0
            ELSE
               KP(IK) = KP(IK-1)
            END IF
         END DO
      END DO

      IF ( debug ) THEN
	     open(39,file='kpath.dat',form='formatted')
         IK = 0
         DO ID = 1,n_k_segments
            DO IKD = 1,nkp_segment(ID)
               IK = IK + 1
               WRITE (39,'(A,4I5,3F8.3,F10.3)') ' ->K ',ID,nkp_segment(ID),IKD,IK,bk(:,IK),KP(IK)
            END DO
         END DO
		 close(39)
      END IF
	  
	  write(*,'(a,i5,a)') 'Using ',nkp,' k-points.'

      allocate(eival(msize))
	  allocate(eigen(nep,nkp,msize))
	  allocate(n_pos(nep,nkp))
	  
goto 1413	  
	  nkp=1
	  bk(:,1)=0.d0
	  bk(2,1)=0.5d0
      open(44,file='mssq.txt',form='formatted')
      open(45,file='zeig.txt',form='formatted')
1413 continue

!CCCCCCCCCCCCCCCCCCCCCC GRIDS READY CCCCCCCCCCCCCCCCCCCCCCCCCCC

      ieverb=nint(dble(nep)*0.05d0)
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCC START LOOP OVER K-POINTS CCCCCCCCCCCCCCCCCCCCCCCCCC
      do ie=1,nep
         if(mod(ie,ieverb).eq.1) write(*,'(a,i5)') 'energy point',ie

         !CCCCCCCCCCCCCCCCCCCCCCCCCCCC LOOP OVER ENERGY POINTS CCCCCCCCCCCCCCCCCCCCCCCCCCC
         do ik=1,nkp

	       ene=dcmplx(emin+(ie-1)*estep,ebroad) !+xmu
           do  isp = 1, nsp
			   call terp(efit,dble(eref(:,isp)),ne1,3,dble(ene),refere)
			   call terp(efit,dimag(eref(:,isp)),ne1,3,dble(ene),refeim)
			   refene=dcmplx(refere,refeim)
               dck=sqrt(2*(ene-dcmplx(refere,refeim)))
!               dck=sqrt(2*(ene-eref(1,isp)))  !Assuming reference energy to be constant w.r.t. ie
               rpart  = real( dble(dck))
               aipart = real(dimag(dck))
               ck(isp) = cmplx(rpart, aipart)
           enddo

!          Get the phase shifts for this energy by interpolation
           do ipp = 0,nph
           do isp = 1, nsp
           do ill = -lmaxph(ipp), lmaxph(ipp)
               call terp(efit,dble(ph(:,ill,isp,ipp)),ne1,3,dble(ene),yr)
               call terp(efit,dimag(ph(:,ill,isp,ipp)),ne1,3,dble(ene),yi)
               xphase(isp,ill,ipp)=cmplx(yr,yi)
           enddo
	       enddo
	       enddo

!          Get the eigenvalues of [t^-1 - G]   for bk(:,ik) and ene(ie), returned in eival :
           call fmsband(ispin,xphase,ene-refene,eival(:),bk(:,ik)/a2pi,freeprop)

           eigen(ie,ik,:)=eival(:)


!          So far, so good; but what metric do we use to determine an eigenstate?
!          copied from SPRKKR : count the number of positive eigenvalues; an increase means we're at an E(k).
!          Is it done this way because broadening spreads the band??

           n_pos_eigenval=0
		   do i=1,msize
		      if(dble(eival(i)).gt.0) n_pos_eigenval=n_pos_eigenval+1
		   enddo
		   n_pos(ie,ik)=n_pos_eigenval
		   
	     enddo  ! ik
         !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC END OF LOOP OVER K-POINTS CCCCCCCCCCCCCCCCCCCC
          
	enddo  ! ie
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC END OF LOOP OVER ENERGY POINTS  CCCCCCCCCCCCCCCCCCCCCCCC




  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  IDENTIFY BANDS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  allocate(bands(nep,nkp),nbands(nkp))  !first dimension of bands will index bands ; "nep" is merely an absolute upper limit for this.
  bands=0.d0
  nbands=0
  
  do ik=1,nkp
     do ie=2,nep
	    if(n_pos(ie,ik).gt.n_pos(ie-1,ik)) then
		   ! We found another E(k) !!
		   do i=1,(n_pos(ie,ik)-n_pos(ie-1,ik))
		      nbands(ik)=nbands(ik)+1
		      bands(nbands(ik),ik)=emin+(ie-1)*estep
		   enddo
		endif
	 enddo
  enddo
  
  ! Count minimum/maximum number of bands
  nbandsmin=100000
  nbandsmax=0
  do ik=1,nkp
     if(nbands(ik).gt.nbandsmax) nbandsmax=nbands(ik)
	 if(nbands(ik).lt.nbandsmin) nbandsmin=nbands(ik)
  enddo

  WRITE (*,99009) nbandsmax

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  WRITE OUTPUT FILES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
		   
  ! Output to bandstructure.dat
  open(42,file='bandstructure.dat',form='formatted')
  write(42,*) '# grid of ',nkp,' k-points.'
  write(42,*) '# grid of ',ne,' energy points  emin= ',emin,' , emax= ',emax,' , estep= ',estep
  write(42,*) '# Found between ',nbandsmin,' and ',nbandsmax,' number of bands.'
  do ik=1,nkp
     write(42,'(i5,x,3(f8.4,x),i4,x,300(f8.4,x))') ik,bk(:,ik),nbands(ik),bands(1:nbands(ik),ik)
  enddo
  close(42)

! Output to SPRKKR format output file
! This file needs to be in Ry units -- multiply energies by 2
! Header
      OPEN (43,FILE="bandstructure_sprkkr.dat",form="formatted")

      WRITE (43,9006) 'KEYWORD   ','DISPERSION'
      WRITE (43,9001) 'TITLE     ','FEFF9 calculation for Po'
      WRITE (43,9001) 'SYSTEM    ','Po'
      WRITE (43,9003) 'NQ        ',nats
      WRITE (43,9003) 'NT        ',nph
      WRITE (43,9003) '          '
      WRITE (43,9003) 'NE        ',0
      WRITE (43,9003) 'IREL      ',3

      WRITE (43,9004) 'EFERMI    ',xmu*2.d0
      WRITE (43,9002) 'INFO      ','FEFF9      ',300
      WRITE (43,9002) '          '

      WRITE (43,'(''   IQ  NLQ '')')
      DO IQ = 1,nats
         WRITE (43,FMT='(2I5)') IQ,lpot(iq)
      END DO

      WRITE (43,'(''   IT       TXTT      CONC  NAT IQAT'')')
      DO IT = 1,nph
!         WRITE (43,9005) IT,TXTT(IT),CONC(IT),NAT(IT), (IQAT(IA,IT),IA=1,NAT(IT))
         WRITE (43,9005) it,label(it),1.d0,natom(it), 1
      END DO
9001 FORMAT (A10,A)
9002 FORMAT (A10,A,I5)
9003 FORMAT (A10,I10)
9004 FORMAT (A10,F10.5)
9005 FORMAT (I5,1X,A4,6X,F10.5,I5,10I3,:,(41X,10I3))
9006 FORMAT (A10,A10)

! Data
      WRITE (43,99008)
      WRITE (43,99005) 'NKDIR     ',n_k_segments
      DO I = 1,n_k_segments
         WRITE (43,99006) label_k_segments(I)
      END DO
      WRITE (43,99008)
      DO I = 1,n_k_segments
         WRITE (43,99005) 'INDKDIR   ',INDKDIR(I)
      END DO
      WRITE (43,99008)
      WRITE (43,99005) 'NKTAB      ',nkp
      DO IK = 1,nkp
         WRITE (43,99007) KP(IK)
      END DO
      WRITE (43,99008)
!      WRITE (43,99005) 'NBAND     ',NKKR
      WRITE (43,99005) 'NBAND     ',nbandsmax

      if(debug) open(44,file='ejk.dat')

      DO IK = 1,nkp

         WRITE (43,'(2I5)') IK,nbands(IK)
         WRITE (43,'(10F8.4)') (bands(I,IK)*2.d0,I=1,nbands(IK))

         IF ( debug ) THEN
            WRITE (44,99010) IK,(bk(I,IK),I=1,3),KP(IK)
            WRITE (44,99011) nbands(IK),EMIN*2.d0,EMAX*2.d0,(bands(I,IK)*2.d0,I=1,nbands(IK))
         END IF

      END DO
	  close(43)
	  if(debug) close(44)
	  
99005 FORMAT (A10,I10)
99006 FORMAT (10X,A)
99007 FORMAT (8F10.4)
99008 FORMAT (80('#'))
99009 FORMAT (/,5X,'maximum number of eigenvalues ',I3,/)
99010 FORMAT (' -> K ',I3,3F7.3,F10.3)
99011 FORMAT (I3,' eigenvalues between ',F9.4,'  and ',F9.4,/,(8F10.4))


      return
      end
	  
	  
	  
	  
	  
	  
	  
		DOUBLE PRECISION FUNCTION DNRM2 ( N, X)
      INTEGER                           N
      DOUBLE PRECISION                  X( * )
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!     DNRM2 := sqrt( x'*x )

      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER               IX
      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
      INTRINSIC             ABS, SQRT

      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE

         DO  IX = 1, N
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
         enddo
         NORM  = SCALE * SQRT( SSQ )
      END IF

      DNRM2 = NORM
      RETURN

      END

