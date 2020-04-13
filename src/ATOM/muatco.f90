!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: muatco.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine muatco(xnval) 
!               * angular coefficients *
!        sous programmes utilises  cwig3j
!
      implicit double precision (a-h,o-z)
      dimension xnval(30)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabk/afgk
      dimension afgk(30,30,0:3)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
!#mn
       external cwig3j

      do 511 i=1,30
      do 511 j=1,30
      do 511 k=0,3
 511  afgk(i,j,k)=0.0d00
 601  do 701 i=1,norb
         li= abs(kap(i))*2-1
         do 701 j=1,i
            lj= abs(kap(j))*2-1
            kmax=(li+lj)/2
            kmin= abs(li-lj)/2
            if ((kap(i)*kap(j)).lt.0) kmin=kmin+1
! calculate a_k(i,j)
            m=0
            if (j.eq.i .and. xnval(i).le.0.0d0) m=1
!           use to test SIC
!           if (j.eq.i) m=1

            afgk(j,i,0)=afgk(j,i,0)+xnel(i)*(xnel(j)-m)
            if (xnval(i).gt.0.0d0 .and. xnval(j).gt.0.0d0) goto 700
! calculate b_k(i,j)
            b=afgk(j,i,0)
            if (j.eq.i .and. xnval(i).le.0.0d0) then
               a=li
               b=-b*(a+1.0d00)/a
               kmin = kmin+2
            endif
            do 675 k = kmin, kmax,2
               afgk(i,j,k/2)=b*(cwig3j(li,k*2,lj,1,0,2)**2)
 675        continue

 700        continue
 701  continue
      return
      end
