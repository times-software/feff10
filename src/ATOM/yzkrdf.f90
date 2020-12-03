!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yzkrdf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine yzkrdf (i,j,k)
!       * calculate  function yk *
! yk = r * integral of f(s)*uk(r,s)
! uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
! f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)      if nem=0
! f(s)=cg(s,i)*cp(s,j)                      if nem is non zero
! f(s) is constructed by the calling programm  if i < or =0
! in the last case a function f (lies in the block dg) is supposedly
! tabulated untill point dr(j), and its' devlopment coefficients
! at the origin are in ag and the power in r of the first term is k+2

! the output functions yk and zk are in the blocks dp and dg.
! at the origin  yk = cte * r**(k+1) - developement limit,
! cte lies in ap(1) and development coefficients in ag.
!        this programm uses aprdev and yzkteg
 
      implicit double precision (a-h,o-z)
! Josh Kas - Changed array dimensions from 30 to 41 (and others) for high Z elements
! according to Pavlo Baranov's changes.
      common cg(251,41), cp(251,41), bg(10,41), bp(10,41),              &
     &         fl(41), fix(41), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
      dimension chg(10)
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      common/inelma/nem
      dimension bgi(10),bgj(10),bpi(10),bpj(10)
!#mn
       external aprdev
 
      if (i.le.0) go to 51
!     construction of the function f
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
        bgj(l) = bg(l,j)
        bpi(l) = bp(l,i)
  5     bpj(l) = bp(l,j)
      id= min(nmax(i),nmax(j))
      ap(1)=fl(i)+fl(j)
      if (nem.ne.0) go to 31
      do 11 l=1,id
 11   dg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
      do 21 l=1,ndor
 21   ag(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
      go to 55

 31   do 35 l=1,id
 35   dg(l)=cg(l,i)*cp(l,j)
      do 41 l=1,ndor
 41   ag(l)=aprdev(bgi,bpj,l)
      go to 55

 51   ap(1)=k+2
      id=j
 55   call yzkteg (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim)
      return
      end
