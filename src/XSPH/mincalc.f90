!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mincalc.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     This subroutine attempts to minimize the number of
!     matrix element calculations needed.

      subroutine mincalc(kfinmax,indmax,kind,lind,ljind,indcalc,indmap,ljj,ncalc)
      implicit none
      integer kfinmax,indmax
      integer kind(kfinmax),lind(kfinmax),ljind(kfinmax)
      integer indcalc(kfinmax,3)
      integer indmap(kfinmax)
      integer ljj,ncalc

      integer k0,know,lnow,ljnow
      integer ikk,ikmap,ifound

!     initialize
      ncalc=1
      indcalc(1,1)=kind(1)  ! given kfin
      indcalc(1,2)=ljind(1) ! maximum lj for given kfin
      indcalc(1,3)=lind(1)  ! l corresponding to kfin
      indmap(1)=1
      ljj=ljind(1)

      do ikk=2,indmax
         know=kind(ikk)
         lnow=lind(ikk)
         ljnow=ljind(ikk)
         
         if (ljj.lt.ljnow) ljj=ljnow
         ifound=-1 
         ikmap=0
         do while(ikmap.lt.ncalc .and. ifound.lt. 0)
            ikmap=ikmap+1
            if (know.eq.indcalc(ikmap,1)) then
               ifound=1
               indmap(ikk)=-ikmap
               if (ljnow.gt.indcalc(ikmap,2)) indcalc(ikmap,2)=ljnow  
            end if
         end do
!
!     if ifound<0 add new calculation
!
         if (ifound.lt.0) then
            ncalc=ncalc+1
            indmap(ikk)=ncalc
            indcalc(ncalc,1)=know
            indcalc(ncalc,2)=ljnow
            indcalc(ncalc,3)=lnow
         end if
         
      end do
      return
      end
