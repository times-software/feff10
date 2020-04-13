!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_eels.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/05/15 22:57:34 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sub-program mdffmod
      program mdffmod
!     subroutine mdffmod
!     !KJ : This routine finishes calculating the EELS spectrum.
!     It is currently set up to use two q-vectors (q and qq below).  For these vectors, it will calculate sigma_q, sigma_qq, and sigma_q,qq -
!     the former two being given by a DFF(q,E) and the latter being given by MDFF(q,qq,E) and MDFF(qq,q,E).
!     This approach can easily be generalized to an arbitrary number of q-vectors.
      
      use mdff_program_control
      use mdff_qvectors
      use eels_inp
      use mdff_work
      use mdff_spectrum
      use par
      use constants,only: pi,hbarc_eV,MeC2,hbarc_atomic
      implicit none

      real*8,external :: mdff_wavelength
!     
      integer i,j,ios,p,j1,j2,iq,iqq,iqqq
	  integer imdff !usually in global_inp but I don't want to use that whole routine ...  It contains too many 'confusing' variables for another MDFF calculation paradigm
      real*8,allocatable :: qve(:,:)
      complex*16,parameter :: ie = (0.0d0,1.0d0)
      real*8 factor,qfac,qqfac,qqqfac,beta
      character*512 slog
	  
	  !The following are dedicated MDFF input and should come from an input file:
	  integer nq !NOT to be confused with the variable in global_inp ...
	  integer,parameter:: qinput=2,itask=1
	  real*8 q(3),qq(3)
	  complex*16,allocatable :: aq(:)
	  complex*16 term
	  logical :: normalize_waves = .true.
	  

      call par_begin
      if (worker) go to 400
	  open(5,file='global.inp',form='formatted',status='old') ; do i=1,15 ;read(5,*) ; enddo ;read(5,*) i,imdff ; close(5) ; if(imdff.ne.3) goto 400 ! 3 is the dedicated value for MDFF calculation.
	  
	  
!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='logmdff.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logmdff.dat', 'feff')

      call init_control
      call eels_init
!     to avoid duplication, we read several parameters from the eels input file, and then forget about some of them.	  
      call eels_read
      if(relat.eq.0) relatQ=.false. ! relatQ was initialized to true by init_control

!     set some variables to default values :
      th0=dble(0.05)/dble(1000) ! value in rad!!
      if(eels.eq.9) then
         qmodus='1'
      else
         qmodus='L'
      endif

	  open(15,file='mdff.inp',form='formatted',status='old')
	  read(15,*); read(15,*) itask,qinput
	  if(qinput.eq.1) then
	     call wlog("Calculating MDFF for user-specified q,q' - e.g. for plotting")
!        add some initialization here - much of this stuff will move to an input file later.
         nq = 2
	     allocate(aq(nq))
	     aq(1)= dcmplx(1,0)
	     aq(2)= dcmplx(0.8,-0.2)
	     q(:)=0 ; qq(:) = 0
	     q(1)=0.
	     q(2)=0.
	     q(3)=-0.03755
	     qq(1)=0.
	     qq(2)=-0.23240
	     qq(3)=-0.03755
		 if(nq.ne.2) call wlog('Warning - these calculations designed for nq=2.  Stuff Might Look Funky.')
	  elseif(qinput.eq.2) then
	     call wlog("Calculating MDFF for given experimental parameters - e.g. for simulating an EELS experiment")
!        initialize working variables
         call init_work(acoll,aconv,nqr,nqf,qmodus)
         call make_Qvecs2(NPos,1)  !KJ  mult(natom))
         allocate (qve(3,1,npos))
	  else
	     call par_stop('strange value for qinput - exiting.')
	  endif

	  ! The beam is given by psi_i(r,t)=SUM(iq=1,nq) aq(iq) exp(i k_iq . r) exp(i w_i t)
	  ! Where w_i is the beam energy, and k_iq = k_f + q_iq
	  ! Where k_f is defined by the energy loss and the detector position.

		 if(normalize_waves) then
		    do iq=1,nq
			   aq(iq)=aq(iq)/dsqrt(dble(aq(iq))**2+dimag(aq(iq))**2)
			enddo
		 endif

	  	    
      call wlog ('Starting MDFF module.')

!     Read spectra from file - output of module ffmod6/ff2x.
      call wlog('Reading Sigma tensor from file.')
      call mdff_readsp


!     write task description to logmdff.dat
      if (itask.eq.1) then
          call wlog('Calculating EELS cross-section.')
      elseif (itask.eq.2) then
	      call wlog('Calculating MDFF.')
      elseif (itask.eq.3) then
	      call wlog('Calculating CAMDFF.')
      else
	      call wlog('Unknown task request.  Quitting now.') ; stop
      endif

      if(relatQ) then
         call wlog('Relativistic corrections are switched on.')
      else
         call wlog('Nonrelativistic calculation - results may be wrong.')
      endif
	  
      write(slog,'(a,f11.1)') 'Beam energy in eV : ',ebeam
      call wlog(slog)
      write(slog,'(a,a,3(f8.3,x))') 'Beam direction w.r.t. coordinate', ' frame of feff.inp : ',xivec
      call wlog(slog)
      write(slog,'(a,2(f6.2,x))') 'Detector position in mrad : ', thetax*dble(1000),thetay*dble(1000)
	  call wlog(slog)
	  if(qinput.eq.1) then
	     write(slog,'(a,a)') 'There is no collection/convergence angle integration angle in this program, and all corresponding variables will be ignored.', &
	       '  This includes acoll, aconv, nqr, nqf, qmodus, th0, imagic.'
         write(11,'(/,/)')
	  elseif(qinput.eq.2) then
         write(slog,'(a,f6.2)') 'collection semiangle in mrad : ', acoll*dble(1000)
         call wlog(slog)
         write(slog,'(a,f6.2)') 'convergence semiangle in mrad : ', aconv*dble(1000)
         call wlog(slog)
         write(slog,'(a,2(i5,x))') 'Integration mesh for q-vectors : ', nqr,nqf
         call wlog(slog)	  
	  endif
	  write(slog,'(a,i5,a)') 'Beam composition : ',nq,' beams'
	  call wlog(slog)
	  write(slog,'(a)') 'Beam amplitudes and phases - i_beam, Re{A_beam},Im{A_beam}'
	  call wlog(slog)
	  do iq=1,nq
		 write(slog,'(i5,2x,f8.3,2x,f8.3)') iq,dble(aq(iq)),dimag(aq(iq))
		 call wlog(slog)
      enddo
	  

!     process q-mesh further
      call make_qvecs1(nq)
	  call make_qvecs2(nq,1)
      allocate(qve(3,nq))
      call mdff_AngularMesh(thetax,thetay) ! sets final wave vectors == q-vectors
      qve=0.
	  npos=nq !This is tentative and should be checked more carefully
	  if (qinput.eq.1) then
	     write(slog,'(a)') 'You have selected manual input of q-vectors.'
		 call wlog(slog)
		 write(slog,'(a)') 'FEFF will not adjust q for energy, nor transform from lab to crystal basis.'
		 call wlog(slog)
		 if(relatq) then
		    write(slog,'(a)') 'However, FEFF will apply the relativistic correction since you asked for it.  Since this will produce garbage if the beam and crystal z are not parallel - use caution.'
		    call wlog(slog)
		 endif
	     qve(:,1)=q
	     qve(:,2)=qq
		 beta=dsqrt((2+ebeam/MeC2)/(2+ebeam/MeC2+MeC2/ebeam))
         do iq=1,nq
	        QLenVClas(1,iq)=dsqrt(qve(1,iq)**2+qve(2,iq)**2+qve(3,iq)**2)
			if(relatq) qve(3,iq)=qve(3,iq)*(dble(1)-beta*beta)
	     enddo
	  else
	     write(slog,'(a)') 'Automatic determination of q-vectors based on user input: beam direction and distance between q-vectors in diffraction plane.'
	  endif
	  	 



      call wlog ('Converting XAS to EELS.')
!     to convert cross-section from XAS to EELS, we need to include a prefactor
!     factor = 4 gamma^2 / a0^2 * kf/ki * 4*c*epsilon_0 / e^2 / Omega
!     we put everything in a.u.
!     The factor 4 clearly has to go; I'm not sure why. 12/2010
      do i=1,ne
         factor=mdff_wavelength(ebeam)/mdff_wavelength(ebeam-s(i,1)) &! wavevectors are inversely proportional to wavelengths
            * (dble(1) + ebeam/dble(511004))**2 &! gamma, provided ebeam is in eV
            / pi                   &! in atomic units, 4 pi epsilon_0 = 1
            * hbarc_atomic / s(i,1)               ! s is energy loss in eV ; s/hbar is omega !KJ 12-2010 atomic: recalibration
!           / dble(1)^2  ! e=1 in atomic units
         s(i,2:10)=s(i,2:10)*factor
      enddo

     
	  if (headers) then
		  call wlog ('Creating headers.')
		  open(7,file='mdff.dat',form='formatted',status='unknown')
	!     Write header to eels.dat file :
		  if(aver.eq.1) then
			write(7,'(a,a,f6.0,a)') '# Orientation averaged EELS ', 'calculation - beam energy = ',ebeam/1000,'keV'
		  else
			write(7,'(a,f6.0,a)') '# Orientation sensitive EELS calculation - beam energy = ',ebeam/1000,'keV'
			write(7,'(a,3(f8.3,x))') '# Sample to beam orientation : ',xivec
		  endif
		  write(7,'(a,2(f10.3,x),a,i5,a,i2)') '# Collection and convergence semiangle: ',acoll,aconv,'  ; # points: ',nqr,' x',nqf
		  write(7,'(a,2(f10.4,x))') '# Detector position: ',thetax,thetay
		  if(relat.eq.1.and.cross.eq.1) then
			write(7,'(a)') '# Relativistic and cross-terms.'
		  elseif(relat.eq.1.and.cross.eq.0) then
			write(7,'(a)') '# Relativistic, no cross-terms.'
		  elseif(relat.eq.0.and.cross.eq.1) then
			write(7,'(a)') '# Nonrelativistic and cross-terms.'
		  else
			write(7,'(a)') '# Nonrelativistic, no cross-terms.'
		  endif
		  write(7,'(a,a,a)') '#  Energy       total         xx            ', &
			'xy            xz            yx            yy            yz',   &
			'            zx            zy            zz'
      endif ! if headers

      
      call allocate_spectrum_2(ne,9,nq)
      
      call wlog ('Entering big loop over energy.')            
!     now, we have a loop over energy loss :
      do i=1,ne

!        If qinput=1, I don't recalculate q for each energy (this implies that the detector aperture and kf shift with each energy ...)
!        If qinput=2, the q-vectors are calculated for this energy loss
         if(qinput.eq.2) then
		    call mdff_qmesh(ebeam-s(i,1),npos)   
			qve(:,:)=qv(:,1,:)
		 endif

         do iq=1,nq
		 do iqq=1,nq
!           in dipole approximation, now a q-dependent factor (q'_i q'_j) / (q^2-(E/hbar c)^2 has to be added
!           we combine this with the integration over alfa and beta :

!           Spectrum =   SUM (iq=1,nq)  integration_weight(iq) * SUM (i=1,3 ; j=1,3)  Spectrum(iq;i,j) * q-dependent_factor(iq;i,j2)

            if(relatQ) then
               qfac=(QLenVClas(1,iq)**2-(s(i,1)/hbarc_eV)**2)
			   qqfac=(QLenVClas(1,iqq)**2-(s(i,1)/hbarc_eV)**2)
            else
               qfac=QLenVClas(1,iq)**2
			   qqfac=QLenVClas(1,iqq)**2
            endif
			iqqq=iqq+(iq-1)*nq+1
			qqqfac=dble(1)/(qfac*qqfac)

            do j1=1,3
               do j2=1,3
                  p=3*(j1-1)+j2+1 ! 2-10 
				  term= qqqfac*qve(j1,iq)*qve(j2,iqq)*s(i,p)*aq(iq)*dconjg(aq(iqq))
                  x(i,1)=x(i,1)+term
				  x(i,iqqq)=x(i,iqqq)+term
                  xpart(i,p-1,1)=xpart(i,p-1,1)+term
				  xpart(i,p-1,iqqq)=xpart(i,p-1,iqqq)+term
               enddo
            enddo       

         enddo
		 enddo
         write(7,'(22g14.6)') s(i,1),((x(i,j)),j=1,1+nq*nq)
!         write(7,'(22g14.6)') s(i,1),(dble(x(i,j)),j=1,1+nq*nq)
      enddo
      close(7)


      call wlog('Module mdff is finished.  Exiting.')

 400  call par_barrier
      call par_end

!     sub-program mdffmod
      stop
	  call WipeErrorfileAtFinish
!     return
      end

