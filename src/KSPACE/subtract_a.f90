!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: subtract_a.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! given a vector r (in cartesian coord) and bvs (reciprocal vectors)
! subroutine calculates reduced vector r_red (in units of lattice vectors)
! if there was a subtraction l is not zero on exit


      subroutine subtract_a(bvs,r,r_red,l)
        use constants,only: pi2
      implicit none
      real*8 bvs(1:3,1:3),r(1:3),r_red(1:3)
      integer l
      
      integer i,j,chng
      l=0

        r_red=dble(0)
      do i=1,3
         do j=1,3
            r_red(i)=r_red(i)+bvs(i,j)*r(j)
         end do
         r_red(i)=r_red(i)/pi2
      end do
      do i=1,3
         chng=nint(r_red(i))
         l=l+abs(chng)
         r_red(i)=r_red(i)-dble(chng)
      end do
      return 
      end




      subroutine reduce(avs,bvs,r,r_red,l)
        use constants,only: pi2
      implicit none
        real*8,intent(in) :: avs(3,3),bvs(3,3),r(3)
        real*8,intent(out) :: r_red(3)
      real*8 rr(3)
      integer,intent(inout) :: l
!	logical debug
        real*8,parameter :: eps=-0.00000001d0
      
      integer i,j,chng
      l=0

!      if(debug) write(11,*) 'r',r
!	if(debug) write(11,*) 'bvs',bvs

        r_red=dble(0)
      do i=1,3
         do j=1,3
            r_red(i)=r_red(i)+bvs(i,j)*r(j)
         end do
         r_red(i)=r_red(i)/pi2
      end do

!      if(debug) write(11,*) 'r_red',r_red
      do i=1,3
         chng=nint(r_red(i))
!	if(debug) write(11,*) i,r_red(i)-dble(1)
         l=l+abs(chng)
!	   if(debug) write(11,*) i,r_red(i)
         r_red(i)=r_red(i)-dble(chng)
!	   if(debug) write(11,*) i,r_red(i)
!    Here we must take care to deal with numerical errors ( ~ 10^-16)
!    If r is somewhere in the middle of the cell, this noise is not important.
!    If r is just above 0 or just above 1 (or any integer), it will be reduced to 0+noise - that's fine.
!    BUT if r is just below 0 or 1 (or any integer), then it will be reduced to 0-noise, and then by the
!    following statement increased to 1-noise   -   while obviously the correct result (no numerical noise) is 0!!
!    Therefore, allow for a small margin of error where negative numbers are brought to 0 instead of 1.
         if(r_red(i).lt.dble(0).and.r_red(i).gt.eps) then
            r_red(i)=dble(0)
           elseif(r_red(i).lt.dble(0)) then
              r_red(i)=r_red(i)+dble(1)
         endif
!	   if(debug) write(11,*) 'i,chng,l',i,chng,l
      end do
!	if(debug) write(11,*) 'r_red',r_red

! back to normal coordinates
      rr=dble(0)
      do i=1,3
         do j=1,3
            rr(i)=rr(i)+avs(i,j)*r_red(j)
         end do
      end do
      r_red=rr
!      if(debug) write(11,*) 'r_red',r_red

      return 
      end
