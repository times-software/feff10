!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2x.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sub-program exchange
      program ffmod6
!     subroutine ffmod6 (iabs)

!     final calculations for various spectroscopies
!     (EXAFS, XANES, FPRIME, DANES, XES)
!     written by a.ankudinov 2000
!     modified by a.ankudinov 2001 for new I/O structure

!     INPUT: mod6.inp global.dat xsect.bin fms.bin list.dat and feff.bin
!     OUTPUT: xmu.dat (chi.dat for EXAFS)

	  use par
	  use eels_inp,only: elnes=>eels,ipmin,ipmax,ipstep
      use ff2x_inp
	  use global_inp
	  use nrixs_inp
	  use errorfile
      use dimsmod, only: init_dimensions
      implicit none
	  integer ios,iabs


      call par_begin
      if (worker) go to 400
      call OpenErrorfileAtLaunch('ff2x')
      call init_dimensions

      iabs = 1
!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log6.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log6.dat', 'feff')

!     read  input files
      call reff2x

      if (mchi .eq. 1)  then
         call wlog('Calculating XAS spectra ...')

		 if (do_nrixs .eq. 1) then					!NRIXS
			if (ispec.gt.0 .and. ispec.lt.3) then
!				using FMS+Paths method
				call ff2xmujas (ispec, ipr6, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, kfinmax, xivec, ldecmx,abs(ldecmx))
			elseif (ispec.eq.3 .or. ispec.eq.4) then 
!				using FMS+Paths method
				call ff2afsjas ( ipr6, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, kfinmax, xivec, ldecmx,abs(ldecmx))
			else
!				using MS Paths expansion
				call ff2chijas (ispec, ipr6, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, kfinmax, xivec, ldecmx,abs(ldecmx))
			endif
         else										!REGULAR FEFF
			if (ispec.gt.0 .and. ispec.lt.3) then 
!				using FMS+Paths method
            call ff2xmu (ispec, ipr6, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, absolu,  &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, elnes,ipmin,ipmax,ipstep,iGammaCH)    
			elseif (ispec.eq.3 .or. ispec.eq.4) then 
!				using FMS+Paths method
				call ff2afs ( ipr6, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, absolu,  &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, elnes,ipmin,ipmax,ipstep)        
			else
!				using MS Paths expansion
				call ff2chi (ispec, ipr6, idwopt, critcw, s02, sig2g,  tk, thetad, mbconv, absolu,  &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, elnes,ipmin,ipmax,ipstep)      
			endif
		 endif
         call wlog('Done with module: XAS spectra (FF2X: DW + final sum over paths).'//char(13)//char(10))
      endif
      close (unit=11)

  400 call par_barrier
      call par_end
	  if(master)call WipeErrorfileAtFinish
!     sub-program exchange
      stop  
!     return

      end
