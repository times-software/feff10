!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cdos.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cdos (iorb, p2, ikap, rmt, jri,jlast, dx, ri, vxc,     &
     &   iz, dgcn, dpcn, adgc, adpc, dos, iph) !KJ iph
!     calculate central atom DOS with wronskian normalization
!     which can be used to find very deep core levels.
!     coded by a.ankudinov 2004

!     input:
!        p2  - energy
!        ikap    quantum number kappa for the orbital
!        rmt     muffin-tin radius
!        jlast   last point for integration of Dirac eq.
!        jri     first interstitial grid point (imt + 1)
!        dx      dx in loucks' grid (usually .05)
!        ri(nr)  loucks' position grid, r = exp ((i-1)*dx - 8.8)
!        vxc(nr) coulomb+xc potential for total density
!     output:
!        dos  - dos at energy p2

      use dimsmod, only: nrptx, ltot
	  use constants
      implicit double precision (a-h, o-z)

      complex*16 vxc(nrptx), p2
      complex*16 vxcp(nrptx)
      dimension ri(nrptx)
      complex*16 ph0, pu, qu, vm(nrptx)
      complex*16 ps(nrptx), qs(nrptx), aps(10),aqs(10), xnorm1
      complex*16  api(10),aqi(10)
!     all atoms' dirac components and their development coefficients
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      dimension adgc(10,30), adpc(10,30)

!     new arrays for TDLDA
!     solution of homogeneous equations (needed for normalization)
      complex*16 ph(nrptx), qh(nrptx)
      complex*16 pir(nrptx), qir(nrptx), xirf

      complex*16 dum1, factor, pun, qun
      complex*16 wronsk, hpp, jpp
!     need to save pun and qun to start irregular solutions

!     bessel functions
      complex*16 ck, xkmt, tl
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 j0(ltot+2), n0(ltot+2)
      complex*16 jlp(0:ltot+1), hl(0:ltot+1)
      complex*16 jl0(0:ltot+1), hl0(0:ltot+1)

!     iph atom's dirac components and their development coefficients
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),      &
     &             fl(30), fix(30), ibgp
!     fl power of the first term of development limits.
!     ibgp first dimension of the arrays bg and bp (=10)

      complex*16 gg,gp,ag,ap,dv,av,bid
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),            &
     &              dv(nrptx),av(10),bid(2*nrptx+20)
!      gg,gp are the input and output for solout
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/messag/dlabpr,numerr
      character*8 dlabpr
!      xnel here - number of core electrons only
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idm

      external besjn

! Debug: Fer
!     print *, 'Entering cdos'
!     print *, 'p2: ', p2

!     get orbitals(dgc,kappa), their energy(eorb) and occupation(xnval)
! Debug: Fer
! These variables don't seem to be used
!     ipr1 = 0
!     iph = 0
!     nph =1
!     ispinr = 0

      xion = 0
      iunf = 1

      nmax = 0
      ndor = 3
      cl = alpinv
      numerr = 0
      nz = iz
      hx = dx
      ibgp=10
      x0 = 8.8d0
      dr = ri

      idm= 1 + nint(250*0.05/dx)
      if (idm .gt. nrptx) idm = nrptx
      if (mod(idm,2) .eq. 0) idm=idm-1

!     copy information into common's of atomic code
      bg = adgc
      bp = adpc
      cg = dgcn
      cp = dpcn

! initialize
      iholep = 0
      kkap = 0
      call inmuac (iholep, xion, iunf, kkap, iph) !KJ iph
!     norb = norbp
      nmax(norb)=jlast
! Debug
!     print *, 'nmax in cdos: ', nmax
!     note that here norb correspond to photoelectron

!     numerical integration of Dirac eq. works if you have 6 grid points
!     for one period of oscillations, switch to analytical expression
!     for a steplike potential  at large distances

      ck   = sqrt (2*p2+ (p2*alphfs)**2)
      rwkb = 0.5d0 / dx / abs(ck) 
      iwkb = (log(rwkb) + x0) / dx  +  2
      if ( iwkb .gt. idm )     iwkb = idm
      if ( iwkb .lt. 10 )      iwkb = 10
      if ( iwkb .ge. jlast-1 ) iwkb = idm

!     note that here norb correspond to photoelectron

!     calculate initial photoelectron orbital using lda
      call diff (vxc,ri,ikap,cl,hx,jri,vm)
      vm(jri+1:nrptx) = 0.0d0
      ic3 = 0

!     set complex momentum and call bessel functions
      ck =  sqrt (2*p2+ (p2*alphfs)**2)
      rmtx =  10.d0 / abs(dimag(ck))
      if (rmtx.gt.rmt) rmtx = rmt
      jrip = (log(rmtx) + x0) / dx + 2
      rmtp = ri(jrip) - 1.d-20

      vxcp(1:jrip) = vxc(1:jrip)
      vxcp(jrip+1:nrptx) = vxc(jri+1)

      xkmt = rmtp * ck
      call besjn (xkmt, jl, nl)
      ix = ltot
      xkmt = ri(jrip) * ck
      call besjh (xkmt, ix, jl0, hl0)

!     set iwkb
      rwkb = 0.5d0 / dx / abs(ck) 
      iwkb = (log(rwkb) + x0) / dx  +  2
      if ( iwkb .gt. idm )     iwkb = idm
      if ( iwkb .lt. 10 )      iwkb = 10
      if ( iwkb .ge. jlast-1 ) iwkb = idm

      kfin = ikap
      kap(norb) = kfin
      irr = -1
      call wfirdc (p2, kap, nmax, vxcp, ps, qs,aps,aqs,irr,ic3,        &
     &           vm, rmtp, jrip, iwkb)
      ph(1:idm) = ps(1:idm)
      qh(1:idm) = qs(1:idm)

!     wronskian normalization of regular solution done below
!     find irregular solution first that matches H_l outside
      il = abs(kfin) 
      if ( kfin .lt. 0) il  = il - 1
      ilp = il - 1
      if ( kfin .lt. 0) ilp = il + 1
      sign = -1.0
      if ( kfin .gt. 0 ) sign = 1.0
      factor = ck*alphfs
      factor = sign * factor/(1+sqrt(1+factor**2))
      dum1 = 1/ sqrt(1+factor**2)
      api(1) = hl0(il)  * rmtp * dum1
      aqi(1) = hl0(ilp) * rmtp * dum1 * factor
      irr = 1
      call wfirdc (p2,kap,nmax,vxcp,ps,qs,api,aqi,irr,ic3,vm,          &
     &       rmtp, jrip, iwkb)

!     set irregular solution outside jrip
      il = abs(kfin)
      if ( kfin .lt. 0) il = il - 1
      ilp = il - 1
      if ( kfin .lt. 0) ilp = il + 1
      ilx = il + 1
      if ( ilp .gt. il) ilx = ilp + 1
      ilx = ltot
      do i = jrip+1, idm
        xkmt = ck * ri(i)
        call besjh( xkmt, ilx, jlp, hl)
        dum2 = ri(i) / ri(jrip)
        ps(i) = ps(jrip) * dum2 * hl(il) /hl0(il)
        qs(i) = qs(jrip) * dum2 * hl(ilp)/hl0(ilp)
      end do

!     calculate wronskian at jrip
!     hpp = (ps(jrip+1)/ri(jrip+1) - ps(jrip-1)/ri(jrip-1) ) /          &
!    &               (ri(jrip+1) - ri(jrip-1) )
!     jpp = (ph(jrip+1)/ri(jrip+1) - ph(jrip-1)/ri(jrip-1) ) /          &
!    &               (ri(jrip+1) - ri(jrip-1) )
!     hpp = qs(jrip)/ri(jrip)
!     jpp = qh(jrip)/ri(jrip)
!     wronsk = ri(jrip)* ( ps(jrip) * jpp - ph(jrip)*hpp)
      wronsk = 2*alpinv*(ps(jrip)*qh(jrip) - ph(jrip)*qs(jrip))

!c    project green's function on the atomic orbital

!     put wronskian normalization into regular solution
      wronsk = 2 /wronsk
      do 90 i=1, idm 
        ph(i) = ph(i) * wronsk
        qh(i) = qh(i) * wronsk
! Debug: Fer
!       print *, ' i, cg, cp: ', i, cg(i,iorb), cp(i,iorb) 
        qir(i) = cg(i,iorb) * ph(i) + cp(i,iorb) * qh(i)
        if (i.gt.1) then
          pir(i)= pir(i-1) + (qir(i) + qir(i-1))*(ri(i)-ri(i-1))
        else
          pir(i)= qir(i)*ri(i)
        endif
  90  continue
      do 100 i=1,idm
        ph(i) = cg(i,iorb) * pir(i)
        qh(i) = cp(i,iorb) * pir(i)
  100 continue
      do 110 i=1, idm 
        qir(i) = ps(i) * ph(i) + qs(i) * qh(i)
  110 continue

      xirf = 0
      call csomm2 (ri, qir, dx, xirf, ri(jrip), jrip+1)
      gamb = dimag(p2)
      dos = dimag(xirf) * gamb

! Debug: Fer
!     print *, 'cdos: ', dos
!     print *, 'Leaving  cdos'

      return
      end
