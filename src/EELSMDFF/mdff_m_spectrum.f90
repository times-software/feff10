!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_m_spectrum.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2010/12/14 00:22:37 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module mdff_spectrum
      
! The sigma tensor of xmu??.dat :
      real*8, allocatable :: s(:,:)
! The number of energy loss values considered :
      integer ne
! The EELS spectrum :
      complex*16, allocatable :: x(:,:)
! The partial EELS spectra :        
      complex*16, allocatable :: xpart(:,:,:)

      CONTAINS
       subroutine allocate_spectrum_1(n,nip)
         implicit none
         integer n,nip
         allocate(s(n,1+nip))
         s=dble(0)
       end subroutine allocate_spectrum_1
	   subroutine allocate_spectrum_2(n,nip,nq)
	     implicit none
		 integer n,nip,nq
		 allocate(x(n,1+nq*nq),xpart(n,nip,1+nq*nq))
		 x=dcmplx(0)
		 xpart=dcmplx(0)
	   end subroutine allocate_spectrum_2
      end module mdff_spectrum      
