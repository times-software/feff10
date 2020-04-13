!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fixlinenow.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2012/09/04 23:15:20 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fixlinenow(w,nw)
      use dimsmod,only:nwordx
      implicit none 
!      integer nwordx     
!      parameter (nwordx = 20)
      character*120 w(nwordx),k  !KJ 9-2012 raised 20->120 in accordance with rdinp.f90 .  Solved some bugs.
      integer nw,nwold,i
      logical iscomm
      external iscomm
 
      nwold=nw
      do i=1,nw
        k=w(i)
        call untab(k)
        call triml(k)
        w(i)=k
!	write(*,'(a1,a20,a1,i1)') '-',w(i),'-',iscomm(w(i))
        if (iscomm(w(i))) then
          nw=i-1
          exit
        endif
      enddo
      
      if (nw.ne.nwold) then
        do i=nw+1,nwold
        
          w(i)='                    '
        enddo
      endif
      
      
      return
      end
      
