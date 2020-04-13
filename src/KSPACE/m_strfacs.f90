!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_strfacs.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/11/23 22:57:43 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       THE KKR STRUCTURE FACTORS
!*****************************************************************************
        module strfacs
!    Ewald parameter
      real*8 streta
!    R-space cutoff
      real*8 strrmax
!    K-space cutoff
      real*8 strgmax
!    Energy broadening
      real*8 eimag
!    Structure factor
        complex,allocatable :: gk(:,:,:)

      CONTAINS
          subroutine init_strfacs
           !streta=dble(0)
           !strrmax=dble(0)
           !strgmax=dble(0)
           eimag=dble(0)
!	   open(98,file='eimag.txt',err=1010)
!	   read(98,*,err=1010,end=1010) eimag
!	   close(98)
           return
1010       eimag=dble(0);return
          end subroutine init_strfacs
        subroutine init_gk(n,j)
           implicit none
           integer n,j
           allocate(gk(n,n,j))  !nkkrmax,nkkrmax,nemax,nktabmax)
           gk=cmplx(0,0)
          end subroutine init_gk
          subroutine exit_gk
           deallocate(gk)
          end subroutine exit_gk

        end module strfacs
