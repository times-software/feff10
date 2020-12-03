!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: potrdf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine potrdf (ia)
!        this programm uses akeato(bkeato),aprdev,multrk,yzkrdf
      implicit double precision (a-h,o-z)
      common cg(251,41), cp(251,41), bg(10,41), bp(10,41),              &
     &        fl(41), fix(41), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),dv(251),av(10), &
     &              eg(251),ceg(10),ep(251),cep(10)
!     dg,dp to get data from yzkrdf, dv,eg,ep -output for soldir
      dimension at(251),bt(251)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      ! JK - 435 may need to be changed to 820 for high Z elements.
      common/scrhf1/eps(820),nre(41),ipl
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)
!#mn
       external akeato, bkeato, aprdev
 
      do 9 i=1,ndor
         cep(i)=0.0d00
         ceg(i)=0.0d00
 9       av(i)=anoy(i)
      do 11 i=1,idim
         at(i)=0.0d00
         bt(i)=0.0d00
         ep(i)=0.0d00
         eg(i)=0.0d00
 11      dv(i)=0.0d00
 
!     coulomb terms
      jia=2* abs(kap(ia))-1
      k=0
 21   do 25 i=1,idim
 25   dg(i)=0.0d00
      do 31 i=1,ndor
 31   ag(i)=0.0d00
      max0=0
      do 51 j=1,norb
         do 33 i = 1,10
            bgj(i) = bg(i,j)
 33         bpj(i) = bp(i,j)
         m=2* abs(kap(j))-1
         if (k.gt.m) go to 51
         a=akeato(ia,j,k)/xnel(ia)
         if (a.eq.0.0d00) go to 51
         m=nmax(j)
         do 35 i=1,m
 35         dg(i)=dg(i)+a*(cg(i,j)*cg(i,j)+cp(i,j)*cp(i,j))
         n=2* abs(kap(j))-k
         l=ndor+2-n
         if (l.le.0) go to 51
!        quick fix of development coefficients
         a = a * fix(j)**2
         do 41 i=1,l
            m=n-2+i
 41         ag(m)=ag(m)+a*(aprdev(bgj,bgj,i)+                           &
     &            aprdev(bpj,bpj,i))
 51      max0= max(max0,nmax(j))
      call yzkrdf (0,max0,k)
      do 61 i=1,ndor
         l=k+i+3
         if (l.gt.ndor) go to 61
         av(l)=av(l)-ag(i)
 61   continue
      do 81 i=1,idim
 81   dv(i)=dv(i)+dg(i)
      k=k+2
      if (k.le.ndor) av(k)=av(k)+ap(1)
      if (k.lt.jia) go to 21
 
!     exchange terms
      if (method.eq.0) go to 411
      do 201 j=1,norb
         if (j-ia) 105,201,105
 105     max0=nmax(j)
         jj=2* abs(kap(j))-1
         kma=(jj+jia)/2
         k= abs(jj-kma)
         if ((kap(j)*kap(ia)).lt.0) k=k+1

 111     a=bkeato(j,ia,k)/xnel(ia)
         if (a.eq.0.0d00) go to 151
         call yzkrdf (j,ia,k)
         do 121 i=1,max0
            eg(i)=eg(i)+a*dg(i)*cg(i,j)
 121        ep(i)=ep(i)+a*dg(i)*cp(i,j)
         n=k+1+ abs(kap(j))- abs(kap(ia))
         if (n.gt.ndor) go to 141
         do 135 i=n,ndor
            ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1) *fix(j)/fix(ia)
 135        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1) *fix(j)/fix(ia)
 141     i=2* abs(kap(j))+1
         if (i.gt.ndor) go to 151
         do 143 ix = 1,10
            bgj(ix) = bg(ix,j)
 143        bpj(ix) = bp(ix,j)
         do 145 n=i,ndor
            ceg(n)=ceg(n)-a*aprdev(ag,bgj,n+1-i) *fix(j)**2
 145        cep(n)=cep(n)-a*aprdev(ag,bpj,n+1-i) *fix(j)**2
 151     k=k+2
         if (k.le.kma) go to 111
 201  continue
 411  if (ipl.eq.0) go to 511
      do 481 j=1,norbsc
         if (kap(j).ne.kap(ia).or.j.eq.ia) go to 481
         if (nre(j).lt.0.and.nre(ia).lt.0) go to 481
         m= max(j,ia)
         i= min(j,ia)+((m-1)*(m-2))/2
         a=eps(i)*xnel(j)
         max0=nmax(j)
         do 461 i=1,max0
            at(i)=at(i)+a*cg(i,j)
 461        bt(i)=bt(i)+a*cp(i,j)
         do 471 i=1,ndor
            ceg(i)=ceg(i)+bg(i,j)*a
 471        cep(i)=cep(i)+bp(i,j)*a
 481  continue
 
! addition of nuclear potential and division of potentials and
!       their development limits by speed of light
 511  do 527 i=1,ndor
         av(i)=av(i)/cl
         cep(i)=cep(i)/cl
 527     ceg(i)=ceg(i)/cl
      do 531 i=1,idim
         dv(i)=(dv(i)/dr(i)+dvn(i))/cl
         ep(i)=(ep(i)+bt(i)*dr(i))/cl
 531     eg(i)=(eg(i)+at(i)*dr(i))/cl
      return
      end
