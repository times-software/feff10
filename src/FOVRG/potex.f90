!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: potex.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine potex( ps, qs, aps, aqs, jri, p2)
!        this programm uses bkeato,aprdec,multrk,yzkrdc
      use dimsmod, only: nrptx,ltot  !KJ 2014 added "only" to protext 'nx'
      implicit double precision (a-h,o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      complex*16 aprdec, p2
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,41),cp(nrptx,41),bg(10,41),bp(10,41),         &
     &             fl(41), fix(41), ibgp
      complex*16 dg,ag,dp,ap,dv,av,eg,ceg,ep,cep
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),dv(nrptx),  &
     &          av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
!     dg,dp to get data from yzkrdc, dv,eg,ep -output for soldir
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      common/mulabc/afgkc
      dimension afgkc(-ltot-1:ltot,41,0:4)
      dimension bgj(10),bpj(10)
!#mn
       external aprdec
 
!     ia=norb
      jia=2* abs(kap(norb))-1
      do 9 i=1,10
         cep(i)=0.0d00
 9       ceg(i)=0.0d00
      do 11 i=1,idim
         ep(i)=0.0d00
 11      eg(i)=0.0d00
 
!     exchange terms
      do 201 j=1,norb-1
 105     jj=2* abs(kap(j))-1
         kma=(jj+jia)/2
         k= abs(jj-kma)
         if ((kap(j)*kap(norb)).lt.0) k=k+1
         kmin = k
!        kma=min(kma,15)
!        if (k.lt.kma) goto 201

!111     a=bkeato(j,ia,k)/xnel(ia)
 111     a=afgkc(kap(norb),j,(k-kmin)/2)
         if (a.eq.0.0d00) go to 151
         call yzkrdc (j,k,fl(norb),ps,qs,aps,aqs, p2, norb)
         do 121 i=1,idim
            eg(i)=eg(i)+a*dg(i)*cg(i,j)
 121        ep(i)=ep(i)+a*dg(i)*cp(i,j)
         n=k+1+ abs(kap(j))- abs(kap(norb))
!         differrent for irregular solution
         if (fl(norb) .lt.0.0) n=k+1+ abs(kap(j)) + abs(kap(norb))
         if (n.gt.ndor) go to 141
         do 135 i=n,ndor
            ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1)*fix(j)/fix(norb)
 135        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1)*fix(j)/fix(norb)
 141     i=2* abs(kap(j))+1
         if (i.gt.ndor) go to 151
         do 143 ix = 1,10
            bgj(ix) = bg(ix,j)
 143        bpj(ix) = bp(ix,j)
         do 145 n=i,ndor
            nx = n + 1 - i
            ceg(n) = ceg(n) - a * aprdec(ag,bgj,nx)*fix(j)**2
 145        cep(n) = cep(n) - a * aprdec(ag,bpj,nx)*fix(j)**2
 151     k=k+2
         if (k.le.kma) go to 111
 201  continue
 
!    division of potentials and
!    their development limits by speed of light
      do 527 i=1,ndor
         cep(i)=cep(i)/cl
 527     ceg(i)=ceg(i)/cl
      do 531 i=1,jri
         ep(i)=ep(i)/cl
 531     eg(i)=eg(i)/cl
      do 532 i=jri+1,nrptx
         ep(i)=0.0d0
 532     eg(i)=0.0d0

      return
      end
