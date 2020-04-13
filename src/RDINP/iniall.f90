!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: iniall.f90,v $:
! $Revision: 1.13 $
! $Author: bmattern $
! $Date: 2012/02/09 18:04:57 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniall

!     initializes all input variables
!     written by Alexei Ankudinov , march 2001.
!     complete overhaul Kevin Jorissen, july 2009

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

!  called in order of appearance in case one initialization affects the next
!  (modules get variables in first come first served order!)
!  although I don't think such dependencies occur at initialization level.      
	  call geometry_init
	  call global_init
	  call reciprocal_init
	  call potential_init
	  call ldos_init
      call opcons_init
	  call screen_init
	  call xsph_init
	  call fms_init
	  call paths_init
	  call genfmt_init
	  call ff2x_init
      call sfconv_init
	  call eels_init
      call compton_inp_init
      call band_init
      call rixs_init
      call hubbard_init
      call crpa_init
      call fullspectrum_init

      return
      end
