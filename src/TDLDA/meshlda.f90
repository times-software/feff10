!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: meshlda.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine meshlda (xkstep, ne, ne1, ne3, em, ik0)
      !KJ Feb 14 added ne3, currently set=0
     
      use dimsmod, only: nex
	  use constants
      implicit double precision (a-h, o-z)
     
      complex*16 em(nex)

!     constant step near Ef
!  manual input
!c for 3d elements
!     eleft = -15.0 / hart
!     eright = 15.0 / hart
!     eext = 60.0 / hart
!c for 4d, 5d elements
!     eleft = -15.0 / hart
!     eright = 75.0 / hart
!     eext = 150.0 / hart
!c for Xe and diamond and MgO
!     eleft = -15.0 / hart
!     eright = 45.0 / hart
!     eext = 135.0 / hart
!c for W and Ta
      eleft = -20.0 / hart
      eright = 200.0 / hart
      eext = 450.0 / hart

      ne = 100
      next = 20

      step = (eright - eleft)/(ne-1)
      step1 = (eext - eright)/(next-1)

      nk = int(sqrt((eright-eleft)/(2*xkstep)))
      next = 20
        
      ne1 = ne + next


      em(1) = eleft
      do i = 2, ne
        em(i) = em(i-1) + step 
      enddo
   
      do i = ne+1, ne1
          em(i) = em(i-1) + step1 
      enddo
!      do 20 i = ne+1, ne1
!         em(i) = em(i-1) + 4*((xkstep*(i - ne + 2))**2)
!   20 continue

!      do 25 i = 2, ne1
!         em(i) = em(i-1) + ((xkstep*i)**2)/2
!   25 continue
 
     
!     don't need ik0
      ik0 = 0

      !KJ Feb 14  No horizontal auxiliary grid?
      ne3=0

      return
      end

