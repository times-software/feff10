!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wrtall.f90,v $:
! $Revision: 1.14 $
! $Author: bmattern $
! $Date: 2012/02/09 18:04:57 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrtall

!     all necessary input files for other modules.
!     version 1.0 written by Alexei Ankudinov, March 2001
!     Torn to shreds Kevin Jorissen 7-09
!     Now everything repackaged into modules COMMON/m_inpmodules.f90

!     Adding new variables now only involves declaring them in one module
!     in said file, and making sure they are both in the corresponding
!     input and output routine of the same module.  Done!

      use par, only: master
      use geometry_inp
	  use global_inp
	  use reciprocal_inp
	  use potential_inp
	  use ldos_inp
      use opcons_inp
	  use screen_inp 
	  use xsph_inp
	  use fms_inp
	  use paths_inp
	  use genfmt_inp
	  use ff2x_inp
	  use sfconv_inp
	  use eels_inp
      use compton_inp
	  use band_inp
      use rixs_inp
	  use hubbard_inp
	  use crpa_inp
      use fullspectrum_inp
	  implicit none

      if (.not. master) return

      call geometry_write_atoms
	  call global_write(.true.)
	  call reciprocal_write
	  call potential_write
      call ldos_write
      call opcons_write
	  call screen_write
	  call xsph_write
	  call fms_write
	  call paths_write
	  call genfmt_write
	  call ff2x_write
	  call sfconv_write
	  call eels_write
      call compton_write
	  call band_write
      call rixs_write
	  call hubbard_write
	  call crpa_write
      call fullspectrum_write
! screen.inp is optional ; it is written by screen_inp_parse_and_write in the presence of a SCREEN card in feff.inp

     
      return
      end
