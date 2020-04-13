!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_program_control.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/07/02 05:32:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module program_control
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
             writeqmesh=.false.
             headers=.true.
            end subroutine init_control
      end module program_control
