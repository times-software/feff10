!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: reff2x.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reff2x  !KJ stuffed everything into modules 7-09

	  use constants, only: hart
	  use global_inp
	  use eels_inp
	  use ff2x_inp
      use nrixs_inp
      implicit none
	  integer ios


      call global_read
	  call eels_read
      call ff2x_read

	  if (do_nrixs .eq. 1) then  ! it may be that this and the next line mean the same, but I'm not sure !KJ 7-09
	      if (ldecmx.gt.0) then 
             open (unit=1,file='xsecl.bin',status='unknown',iostat=ios)
             call chopen(ios,'xsecl.bin','ff2xmuq')
             read(1,*) kfinmax
             close(1)
          else
             if (ispec.gt.0 .and. ispec.lt.3) then 
                ldecmx=-1
                kfinmax=-1
             elseif (ispec.eq.3 .or. ispec.eq.4) then 
				kfinmax=0.00d0
             else
				kfinmax=-1
             endif
          end if
		  call nrixs_init
	  endif


!     transform energies to atomic units
      vrcorr = vrcorr / hart
      vicorr = vicorr / hart

      return
      end



