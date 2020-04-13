!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fdrirk.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function fdrirk (i,j,l,m,k)
!                       * calculate radial integrales rk *
!        rk = integral of f(r) * uk(r,s) * g(s)
! uk(r,s) = rinf**k / rsup**(k+1)    rinf=min(r,s)   rsup=max(r,s)
!        if nem=0  f(.)=cg(.,i)*cg(.,j)+cp(.,i)*cp(.,j)
!                  g(.)=cg(.,l)*cg(.,m)+cp(.,l)*cp(.,m)
!        if nem non zero f(.)=cg(.,i)*cp(.,j)
!                        g(.)=cg(.,l)*cp(.,m)
!                  cg (cp) large (small) componenents of the orbitales
! moreover if nem > or =0 the integration is made from 0 to infinity,
! and otherwise from 0 to r.
!        this programm uses yzkrdf and dsordf
 
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
! comdir is used just to exchange variables between dsordf,yzkrdf,fdrirk
      dimension hg(251)
      common/inelma/nem
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      save
 
      fdrirk=0.0d00
      if (i.le.0.or.j.le.0) go to 201
      call yzkrdf (i,j,k)
      nn= abs(kap(i))+ abs(kap(j))
      nn=max(nn-k,1)
      a=k+1
      do 21 n=1,ndor
 21   hg(n)=0.0d00
      do 31 n=1,ndor
         if (nn.gt.ndor) go to 31
         hg(nn)=-ag(n)
 31      nn=nn+1
      do 41 n=1,ndor
 41      ag(n)=hg(n)
      ag(1)=ag(1)+ap(1)

 201  if (l.le.0.or.m.le.0) return
      n=-1
      if (nem.ne.0) n=-2
      fdrirk=dsordf(l,m,-1,n,a)
      return
      end
