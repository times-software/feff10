!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_m_program_control.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module mdff_program_control
!     Use relativistic corrections to the cross section formula?
      logical          RelatQ
!     Controls the amount of output given by EELS (0-3)
      integer          verbosity
!     Messy fix - remove later
      logical        writeqmesh
!     Write headers or not
      logical        headers

      CONTAINS
            subroutine init_control
!            Initialize the variables for program flow
             RelatQ=.true.! use shortened q-vector
             verbosity=2 ! give a lot of output
             writeqmesh=.true.
             headers=.true.
            end subroutine init_control
      end module mdff_program_control
