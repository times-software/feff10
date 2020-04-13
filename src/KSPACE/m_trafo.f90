!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_trafo.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module trafo
      complex*16, allocatable :: crel(:,:),rrel(:,:),rc(:,:)
        real*8 mrotr(3,3,48)

      contains
           subroutine init_trafo(n)
             implicit none
             integer n
             allocate(crel(n,n),rrel(n,n),rc(n,n))
             crel=dcmplx(0,0)
             rrel=dcmplx(0,0)
             rc=dcmplx(0,0)
           end subroutine init_trafo
        end module trafo
