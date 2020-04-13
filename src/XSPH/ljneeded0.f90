!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ljneeded0.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  ljneeded0(ljmax,ljneeded,kfinmax,indmax,ljind,icalc,indmap,indcalc)
      implicit none
      integer ljmax,kfinmax
      integer ljneeded(0:ljmax)
      integer indmax,icalc
      integer indmap(kfinmax)
      integer indcalc(kfinmax,3)
      integer ljind(kfinmax)
      
      integer ii,ll
      
      ljneeded(:)=0

      do ii=1,indmax
         if (abs(indmap(ii)).eq.icalc) then
            ll=ljind(ii)
            if (ll.gt.ljmax) stop "Something wrong in ljneeded0, ll.gt.ljmax"
            ljneeded(ll)=1
         end if
      end do
      

      return
      end
