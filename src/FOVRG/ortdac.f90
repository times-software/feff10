!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ortdac.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ortdac(ikap,ps,qs,aps,aqs)
!        * orthogonalization by the schmidt procedure*
! the ia orbital is orthogonalized toa all orbitals of the same
! symmetry if ia is positive, otherwise all orbitals of the same
! symmetry are orthogonalized
!        this program uses dsordc
 
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)
      complex*16 dsordc
      complex*16 ps(nrptx), qs(nrptx), aps(10),aqs(10)
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),     &
     &             fl(30), fix(30), ibgp
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      complex*16 a
 
      do 51 j=1,norb-1
         if (kap(j).ne.ikap .or. xnel(j).le.0) go to 51
         a = dsordc(j,fl(norb),ps,qs,aps,aqs)
         do 41 i=1,idim
            ps(i)=ps(i)-a*cg(i,j)
 41         qs(i)=qs(i)-a*cp(i,j)
         do 42 i=1,ndor
            aps(i)=aps(i)-a*bg(i,j)
 42         aqs(i)=aqs(i)-a*bp(i,j)
 51   continue
      return
      end
