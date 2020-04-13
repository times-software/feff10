!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_spectrum.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/06/23 04:00:14 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module spectrum
      
! The sigma tensor of xmu??.dat :
      real*8, allocatable :: s(:,:)
! The number of energy loss values considered :
      integer ne
! The EELS spectrum :
      complex*16, allocatable :: x(:)
! The "atomic background" :
      complex*16, allocatable :: bg(:)	  
! The partial EELS spectra :        
      complex*16, allocatable :: xpart(:,:)

      CONTAINS
       subroutine allocate_spectrum(n,nip)
         implicit none
         integer n,nip
         allocate(s(n,1+nip+1),x(n),xpart(n,nip),bg(n))
         s=dble(0)
         x=dcmplx(0)
         xpart=dcmplx(0)
		 bg=dcmplx(0)
       end subroutine allocate_spectrum
         
      end module spectrum      
