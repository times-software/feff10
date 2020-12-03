!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: muatcc.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine muatcc(xnval) 
!               * angular coefficients *
!        sous programmes utilises  cwig3j
!

      use dimsmod, only: ltot
      implicit double precision (a-h,o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      dimension xnval(41)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabc/afgkc
      dimension afgkc(-ltot-1:ltot,41,0:4)
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
!#mn
       external cwig3j

      do 511 i=-ltot-1,ltot
      do 511 j=1,41
      do 511 k=0,3
 511  afgkc(i,j,k)=0.0d00
 601  do 701 ikap=-ltot-1,ltot
         if (ikap .eq. 0) go to 701
         li= abs(ikap)*2-1
         do 700 j=1,norb-1
            lj= abs(kap(j))*2-1
            kmax=(li+lj)/2
            kmin= abs(li-lj)/2
            if ((ikap*kap(j)).lt.0) kmin=kmin+1
            if (xnval(j) .gt. 0.0d0) goto 700
! calculate b_k(i,j)
            do 675 k = kmin, kmax,2
               index=(k-kmin)/2
               afgkc(ikap,j,index)=xnel(j)*(cwig3j(li,k*2,lj,1,0,2)**2)
 675        continue
 700     continue
 701  continue
      return
      end
