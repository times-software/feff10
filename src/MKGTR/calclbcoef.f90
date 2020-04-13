!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: calclbcoef.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calclbcoef(lx,jlmax,mjlmax,clbcoef)
!KJ 7-09 moved a lot of original content to COMMON/m_nrixs.f90  -  not much left now.

      implicit none
      integer,intent(in) :: lx,jlmax,mjlmax
      real*8,intent(out) :: clbcoef(1:mjlmax,1:jlmax,0:1,0:lx)
      integer jnow,mj,ms,ll,lnow,ii,im,is
      real*8 res3j
      real*8,external :: cwig3j
      
      clbcoef(:,:,:,:)=0.0d0
      do ll=0,lx
         lnow=2*ll
         do is=0,1
            ms=2*is-1
            do ii=1,jlmax
               jnow=2*ii-1
               if (jnow.le.2*ll+1) then 
                  do im=1,2*ii
                     mj=-jnow+2*(im-1)
!     Here we calculate Clebsch-Gordan coefficient
!     <j_f m_f | l_f m_f s ms_f > 
                     res3j=cwig3j(1,jnow,lnow,ms,-mj,2)
                     if (mod((lnow+mj-1)/2,2).ne.0) res3j=-res3j
                     clbcoef(im,ii,is,ll) = res3j
                  end do
               end if
            end do
         end do
      end do

      return 
      end
