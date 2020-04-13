!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: eels.f90,v $:
! $Revision: 1.13 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sub-program eelsmod
      program eelsmod
!     subroutine eelsmod
!     !KJ : This routine finishes calculating the EELS spectrum.
      
      use program_control
      use qvectors
      use eels_inp !KJ 7-09 replaces old module input throughout program
      use work
      use spectrum
      use par
      use constants,only: pi,hbarc_eV,hbarc_atomic
	  use errorfile
      implicit none

      real*8,external :: wavelength
!     
      integer i,j,ios,iq,p,j1,j2,iheader
      complex*16,allocatable :: qve(:,:,:)
      complex*16,parameter :: ie = (0.0d0,1.0d0)
      real*8 factor,qfac
      character*512 slog
	  character*100 header(5)

      call par_begin
      if (worker) go to 400
      call OpenErrorfileAtLaunch('eels')

!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='logeels.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logeels.dat', 'feff')

      call init_control
      call eels_init
      call eels_read
      if(relat.eq.0) relatQ=.false. ! relatQ was initialized to true by init_control
      if(eels.eq.0) goto 400 
      call wlog ('Calculating EELS spectra ...')

!     Read spectra from file - output of ffmod6.
      !call wlog('Reading spectra from file.')
      call readsp(header,iheader)


!     write task description to logeels.dat
      if(eels.eq.9) then
          call wlog('Printing GOS tables for Jo Verbeeck.')
      else
          !call wlog('Calculating ELNES spectrum.')
      endif

      if(relatQ) then
         call wlog('Using relativistic theory.')
      else
         call wlog('Using nonrelativistic calculation - may be wrong for anisotropic measurements.')
      endif
      write(slog,'(a,f8.2,a)') 'Beam energy=',ebeam/1000.d0,' keV'
      call wlog(slog)
      write(slog,'(a,3(f6.3,x),a)') 'Beam direction=',xivec(1:3),' in coordinate frame of feff.inp'
      call wlog(slog)
      write(slog,'(a,f6.2,a,2x,a,f6.2,a)') 'Collection semiangle=', &
        acoll*1000.d0,' mrad','convergence semiangle=',aconv*dble(1000),' mrad'
      call wlog(slog)
      write(slog,'(a,i5,a,a,i4,a,i4,a,i3,a)') 'Integration mesh for q-vectors=', nqr*nqf*nqr,' points','(',nqr,'x',nqr,'x',nqf,')'
      call wlog(slog)
      write(slog,'(a,2(f6.2,a,x))') 'Detector position=', thetax*dble(1000),'mrad',thetay*dble(1000),'mrad'
      if((1000.d0*(dabs(thetax)+dabs(thetay))).gt.0.05) call wlog(slog)
      write(11,'(/,/)') 


!     set some variables to default values :
      th0=dble(0.05)/dble(1000) ! value in rad!!
      if(eels.eq.9) then
         qmodus='1'
      else
         qmodus='L'
      endif
!     initialize working variables
      call init_work(acoll,aconv,nqr,nqf,qmodus)
      call make_Qvecs2(NPos,1)  !KJ  mult(natom))
      allocate (qve(3,1,npos))


      !call wlog ('Converting XAS to EELS.')
!     to convert cross-section from XAS to EELS, we need to include a prefactor
!     factor = 4 gamma^2 / a0^2 * kf/ki * 4*c*epsilon_0 / e^2 / Omega
!     we put everything in a.u.
!     The factor 4 clearly has to go; I'm not sure why. 12/2010
      do i=1,ne
         factor=wavelength(ebeam)/wavelength(ebeam-s(i,1)) &! wavevectors are inversely proportional to wavelengths
            * (dble(1) + ebeam/dble(511004))**2 &! gamma, provided ebeam is in eV
            / pi                   &! in atomic units, 4 pi epsilon_0 = 1
            * hbarc_atomic / s(i,1)               ! s is energy loss in eV ; s/hbar is omega !KJ 12-2010 atomic: recalibration
!           / dble(1)^2  ! e=1 in atomic units
         s(i,2:2+nip)=s(i,2:2+nip)*factor
      enddo


      call wlog ('Setting up k_i and k_f vectors and calculating integration weights.')
!     now, the k'-mesh is constructed and the q-integration weights are prepared
      call calculateweights

     
	  if (headers) then
		  !call wlog ('Creating headers.')
		  open(7,file='eels.dat',form='formatted',status='unknown')
	!     Write header to eels.dat file :
		  if(aver.eq.1) then
			write(7,'(a,a,f6.0,a)') '# Orientation averaged EELS ', 'calculation - beam energy = ',ebeam/1000,'keV'
		  else
			write(7,'(a,f6.0,a)') '# Orientation sensitive EELS calculation - beam energy = ',ebeam/1000,'keV'
			write(7,'(a,3(f8.3,x))') '# Sample to beam orientation : ',xivec
		  endif
		  write(7,'(a,2(f10.3,x),a,i5,a,i2)') &
                    '# Collection and convergence semiangle: ', &
                    acoll*1000,aconv*1000,'  ; # points: ',nqr,' x',nqf
		  write(7,'(a,2(f10.4,x))') '# Detector position: ',thetax,thetay
		  write(7,'(a)') '# Units are a_0^2 / eV.  Multiply by 28.00 10^-18  to get cm^-2 / eV.  Or by 28 to get Mbarn / eV.'
		  if(relat.eq.1.and.cross.eq.1) then
			write(7,'(a)') '# Relativistic and cross-terms.'
		  elseif(relat.eq.1.and.cross.eq.0) then
			write(7,'(a)') '# Relativistic, no cross-terms.'
		  elseif(relat.eq.0.and.cross.eq.1) then
			write(7,'(a)') '# Nonrelativistic and cross-terms.'
		  else
			write(7,'(a)') '# Nonrelativistic, no cross-terms.'
		  endif
		  do i=1,iheader
		     write(7,'(a)') header(i)
		  enddo
		  if(aver.eq.0) then
! Modified by FDV
! Split line to compile on Solaris Studio
		  write(7,'(a,a,a)') &
                    '#  Energy       total         atomic-bg     ' // &
                    'fine-struct   xx            xy            xz' // &
                    '            yx            yy            yz' // &
                    '            zx            zy            zz'
          elseif(aver.eq.1) then
		  write(7,'(a,a,a)') '#  Energy       total         atomic-bg     fine-struct'
		  endif
      endif ! if headers

      
      x(:)=dcmplx(0)
      xpart(:,:)=dcmplx(0)
      

      write(slog,'(a,i5,a,i5,a)') 'Calculating absorption for',ne,' energy points and',npos,' impulse transfer vectors.'
      call wlog (slog)
!     now, we have a loop over energy loss :
      do i=1,ne

!     for this energy loss, the q-vectors are calculated
         call qmesh(ebeam-s(i,1),npos)
!     for feff compatibility, we transform them from a cartesian representation to l,m representation
!     qve(1,:,:)=(qv(1,:,:)-ie*qv(2,:,:))/dsqrt(dble(2))
!     qve(2,:,:)=dcmplx(qv(3,:,:))
!     qve(3,:,:)=-(qv(1,:,:)+ie*qv(2,:,:))/dsqrt(dble(2))
         qve=qv                 !to stay in carthesian coordinates
         

         do iq=1,npos
!           in dipole approximation, now a q-dependent factor (q'_i q'_j) / (q^2-(E/hbar c)^2 has to be added
!           we combine this with the integration over alfa and beta :

!           Spectrum =   SUM (iq=1,nq)  integration_weight(iq) * SUM (i=1,3 ; j=1,3)  Spectrum(iq;i,j) * q-dependent_factor(iq;i,j2)

            if(relatQ) then
               qfac=(QLenVClas(1,iq)**2-(s(i,1)/hbarc_eV)**2)**2
            else
               qfac=QLenVClas(1,iq)**4
            endif

            do j1=1,3
               do j2=1,3
                  p=3*(j1-1)+j2+1 ! 2-10  
                  x(i)=x(i)+WeightV(iq)/qfac*qve(j1,1,iq)*qve(j2,1,iq)*s(i,p)
				  if(j1.eq.j2) bg(i)=bg(i)+WeightV(iq)/qfac*qve(j1,1,iq)*qve(j2,1,iq)*s(i,11)
                  xpart(i,p-1)=xpart(i,p-1)+WeightV(iq)/qfac*qve(j1,1,iq)*qve(j2,1,iq)*s(i,p)
               enddo
            enddo       

         enddo
		 
!         write(12,'(22g14.6)') s(i,1),dble(x(i)),factor,qfac,qlenv(1,1),qlenv(1,1)**2 / qfac
      enddo

      call wlog('Writing EELS to eels.dat .')
	  do i=1,ne
	     if(aver.eq.1) then  
            write(7,'(22g14.6)') s(i,1),dble(x(i)),dble(bg(i)),dble(x(i)-bg(i))  ! partial terms are meaningless for averaged calculation		 
		 elseif(aver.eq.0) then
            write(7,'(22g14.6)') s(i,1),dble(x(i)),dble(bg(i)),dble(x(i)-bg(i)),(dble(xpart(i,j)),j=1,9)		 
		 endif
	  enddo
      close(7)


      if(eels.eq.9) then
         call writeangulardependence3
      endif


      if(magic.eq.1) then
         call wlog('Producing angular output for the magic angle in magic.dat .')
         open(59,file='magic.dat',form='formatted',status='unknown')
         call writeangulardependence2
         close(59)
      endif


 111  continue
      call wlog('Done with module: EELS.'//char(13)//char(10))

 400  call par_barrier
      call par_end
      if(master)call WipeErrorfileAtFinish
!     sub-program eelsmod
      stop
!     return
      end

