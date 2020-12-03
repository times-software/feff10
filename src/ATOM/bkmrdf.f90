!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bkmrdf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine bkmrdf (i,j,k)
!     angular coefficients for the breit term. i and j are the numbers
!     of orbitals and  k is the value of k in uk(1,2)
!        this programm uses cwig3j
!     coefficients for magnetic interaction  are in cmag
!     and those for retarded term are in cret
!     the order correspond to -1 0 and +1
 
! Josh Kas - Changed array dimensions from 30 to 41 for high Z elements
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      common/tabre/cmag(3),cret(3)
!#mn
       external cwig3j
 
      do 12 l=1,3
        cmag(l)=0.0d00
 12     cret(l)=0.0d00
      ji=2* abs(kap(i))-1
      jj=2* abs(kap(j))-1
      kam=kap(j)-kap(i)
      l=k-1
      do 51 m=1,3
         if (l.lt.0) go to 51
         a=cwig3j(ji,jj,l+l,-1,1,2)**2
         if (a.eq.0.0d00) go to 51
         c=l+l+1
         if (m-2) 14,16,17
 14      cm=(kam+k)**2
         cz=kam*kam-k*k
         cp=(k-kam)**2
         n=k
 15      l1=l+1
         am=(kam-l)*(kam+l1)/c
         az=(kam*kam+l*l1)/c
         ap=(l+kam)*(kam-l1)/c
         d=n*(k+k+1)
         go to 31

 16      d=k*(k+1)
         cm=(kap(i)+kap(j))**2
         cz=cm
         cp=cm
         go to 41

 17      cm=(kam-l)**2
         cz=kam*kam-l*l
         cp=(kam+l)**2
         n=l
         c=-c
         go to 15

 31      c= abs(c)*d
         if (c.ne.0.0d00) c=n/c
         cret(1)=cret(1)+a*(am-c*cm)
         cret(2)=cret(2)+(a+a)*(az-c*cz)
         cret(3)=cret(3)+a*(ap-c*cp)
 41      if (d.eq.0.0d00) go to 51
         a=a/d
         cmag(1)=cmag(1)+cm*a
         cmag(2)=cmag(2)+cz*(a+a)
         cmag(3)=cmag(3)+cp*a
 51      l=l+1
      return
      end
