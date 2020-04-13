!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: iniptz.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniptz(ptz,iptz,modus)
        !KJ This routine rewrites the ptz-matrix.

      implicit none
!     which polarization tensor to create      
      integer iptz
!     two ways of working
      integer modus
!     the polarization tensor
      complex*16 ptz(-1:1,-1:1)

      complex*16 zero,one,coni
      parameter (zero=(0,0),one=(1,0),coni=(0,1))
      integer i,j
      real*8 unity(3,3)



      if (iptz.lt.1.or.iptz.gt.10) then
          call wlog('Inieln sees weird iptz - returns without           &
     &      changing ptz - danger of calculating nonsense !!')
      endif


      do i=1,3
      do j=1,3
        unity(i,j)=dble(0)
      enddo
      unity(i,i)=dble(1)/dble(3)
      enddo
      do i=-1,1
      do j=-1,1
        ptz(i,j)=zero
      enddo
      enddo
      

      if (modus.eq.1) then
! work in spherical coordinates

         if(iptz.eq.10) then
        do i=-1,1
        do j=-1,1
           ptz(i,j)=unity(i+2,j+2)
        enddo
        enddo
         else
            i=(iptz-1)/3+1  ! row index
            j=iptz-3*(i-1)  ! column index
            i=i-2 !shift from 1..3 to -1..1
            j=j-2
            ptz(i,j)=one
         endif  


      elseif(modus.eq.2) then
! work in carthesian coordinates      

      if (iptz.eq.10) then ! orientation averaged spectrum
        do i=-1,1
        do j=-1,1
           ptz(i,j)=unity(i+2,j+2)
        enddo
        enddo
        
        elseif (iptz.eq.1) then   ! x x*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=-one/dble(2)
        ptz(1,-1)=-one/dble(2)
        elseif (iptz.eq.5) then ! y y*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=one/dble(2)
        ptz(1,-1)=one/dble(2)
         elseif (iptz.eq.9) then ! z z*
        ptz(0,0)=one
        elseif (iptz.eq.2) then ! x y*
          ptz(1,1)=one*coni/dble(2)
        ptz(-1,-1)=-one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.4) then ! x* y
          ptz(1,1)=-one*coni/dble(2)
        ptz(-1,-1)=one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.3) then ! x z*
          ptz(-1,0)=one/dsqrt(dble(2))
        ptz(1,0)=-one/dsqrt(dble(2))
        elseif (iptz.eq.7) then ! x* z
          ptz(0,-1)=one/dsqrt(dble(2))
        ptz(0,1)=-one/dsqrt(dble(2))
        elseif (iptz.eq.6) then ! y z*
          ptz(-1,0)=-one*coni/dsqrt(dble(2))
        ptz(1,0)=-one*coni/dsqrt(dble(2))
        elseif (iptz.eq.8) then ! y* z
          ptz(0,-1)=one*coni/dsqrt(dble(2))
        ptz(0,1)=one*coni/dsqrt(dble(2))
      endif
      
      
        else
          stop 'alien modus in inieln'
        endif


      return
      end


