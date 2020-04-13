!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_energygrid.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*************************************************************************
      module energygrid
        complex*16,allocatable :: egrid(:)
        complex*16 emin,emax
        integer nene

        contains
           subroutine init_energygrid(n)
              implicit none
              integer n
              nene=n
              allocate(egrid(n))
              egrid=dcmplx(0,0)
           end subroutine init_energygrid
        end module energygrid
