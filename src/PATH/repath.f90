!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: repath.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/24 09:17:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine repath !KJ 7-09 I put everything in modules

	  use constants
      use atoms_inp
	  use paths_inp
	  use global_inp
	  use eels_inp
	  use nrixs_inp
      implicit none


      call atoms_read  ! read geom.inp
      call global_read ! read global.inp
      call paths_read  ! read paths.inp
	  call eels_read   ! read eels.inp
          
	  IF(mpath.ne.0) call nrixs_init  ! getsizes only if path module is going to run. Josh
            

      return
      end
