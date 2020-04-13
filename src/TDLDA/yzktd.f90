!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yzktd.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine yzktd(i,k,flps,ps,qs,aps,aqs,ykgr,j)
! nesvi - don't want to use common blocks to pass ykgr, otherwise it's
!         the same subroutine as yzkrdc

!       * calculate  function yk *
! yk = r * integral of f(s)*uk(r,s)
! uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
! j=norb for photoelectron
! f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)
! f(s) is constructed by the calling programm  if i < or =0
! in the last case a function f (lies in the block dg) is supposedly
! tabulated untill point dr(j), and its' devlopment coefficients
! at the origin are in ag and the power in r of the first term is k+2

! the output function ykgr.
! at the origin  yk = cte * r**(k+1) - developement limit,
! cte lies in ap(1) and development coefficients in ag.
!        this programm uses aprdec and yzktec
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)

      complex*16 aprdec, dyzk
!     complex*16 a1,a2,b1,b2,coni
!     complex*16 xck, temp, ck, phx
      parameter (coni=(0.d0,1.d0))
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),      &
     &             fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,bidcom, chg(10)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),            &
     &   bidcom(3*nrptx+30)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bpi(10)
      complex*16 bgj(10),bpj(10)

      complex*16 ykgr(nrptx)
!#mn
       external aprdec

!     construction of the function f
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
  5     bpi(l) = bp(l,i)

      if (j.eq.0) then
        id=min(nmax(i),np)
        ap(1)=fl(i)+flps
        do 11 l=1,id
 11     dg(l)=cg(l,i)*ps(l)+cp(l,i)*qs(l)
        do 12 l = id+1,idim
 12      dg(l) = 0.0d0
        do 13 l=1,ndor
 13     ag(l) = aprdec(aps,bgi,l) + aprdec(aqs,bpi,l)
      else
        do 15 l= 1,ibgp
          bgj(l) = bg(l,i)
 15       bpj(l) = bp(l,i)
        id=min(nmax(i),nmax(j))
        ap(1)=fl(i)+fl(j)
        do 21 l=1,id
 21     dg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
        do 22 l = id+1,idim
 22      dg(l) = 0.0d0
        do 23 l=1,ndor
 23     ag(l) = aprdec(bgj,bgi,l) + aprdec(bpj,bpi,l)
      endif

      dyzk = 0

      call yzktec (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim, dyzk)

      do 777 l=1,nrptx
!       yk is in dg
        ykgr(l) = dg(l)
  777 continue
 
      return
      end
