!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: regenf.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/24 09:17:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine regenf 
 
	  use constants
	  use global_inp
	  use genfmt_inp
	  use eels_inp
	  use nrixs_inp
      use pdata, only: init_pdata
      use str, only: init_str
      implicit none

      call global_read  ! read global.inp
	  call genfmt_read  ! read genfmt.inp/mod5.inp
      call eels_read    ! read eels.inp

	  if (do_nrixs .eq. 1) then
	     IF(mfeff.ne.0) call nrixs_init ! Josh - only call nrixs_init if genfmt is going to run.
         if (elpty.lt.0) then
			call wlog('Spherically averaged NRIXS in module genft - setting jinit=jmax.')
			jinit=jmax
         endif
      end if

      if(mfeff.eq.1) then
          call init_pdata
          call init_str
      endif

      return
      end
