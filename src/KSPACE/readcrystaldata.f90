!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: readcrystaldata.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine cross(a,b,c)
!    Vector product :  c = a x b
      implicit none
      real*8 a(3),b(3),c(3)
      integer i,j,k
      do i=1,3
        j=i+1
        if (j.gt.3) j=j-3
        k=j+1
        if (k.gt.3) k=k-3
        c(i)=a(j)*b(k)-a(k)*b(j)
      end do
      return
      end
!--------------------------------------------
      function dotthree(b,a)
!    Scalar product : dotthree = a . b
      implicit none
      real*8 dotthree,a(3),b(3),temp
      integer i
      temp=0.d0
      do i=1,3
        temp=temp+b(i)*a(i)
      end do
      dotthree=temp
      return
      end






