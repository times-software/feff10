!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_wigner3j.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module wigner3j
! Wigner 3j symbols needed to set up FEFF's t-matrix in routine fms2.

        real,allocatable :: t3jp(:,:,:),t3jm(:,:,:)

        contains
           subroutine init_wigner3j(l)
             implicit none
             integer,intent(in) :: l

             allocate(t3jp(0:l,-(l+1):(l+1),2),t3jm(0:l,-(l+1):(l+1),2))
! m-field is allowed to go up to l+1, in order to zero out m2=m1+1 contributions for m1=l
! in the setup of off-diagonal elements of the t-matrix in subr fms(2).
             t3jp=real(0)
             t3jm=real(0)
           end subroutine init_wigner3j
        end module wigner3j
