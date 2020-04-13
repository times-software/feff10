!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: change_car.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine change_car(bvs,avs,mtrx,imtrx,r)

!     calculates   r = bvs' * mtrx(imtr) * avs
      implicit none
      double precision bvs(1:3,1:3),avs(1:3,1:3),r(1:3,1:3)
      integer mtrx(48,1:3,1:3),imtrx
      integer i,j,k,l

        r=dble(0)
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  r(i,j)=r(i,j)+bvs(l,i)*                               &
     &                 dble(mtrx(imtrx,l,k))*avs(k,j)
               end do
            end do
         end do
      end do


      
      return 
      end 

