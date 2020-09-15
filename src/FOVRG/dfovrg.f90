!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dfovrg.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfovrg (ncycle, ikap, rmt, jlast, jri, p2, dx,         &
                       ri, vxc, vxcval, dgcn, dpcn, adgc, adpc,        &
                       xnval, pu, qu, ps, qs,                          &
                       iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
!     Dirac equation solver for complex energy
!     coded by a.ankudinov 1996
!     modified by a.ankudinov 1997 to get irregular solution 

!     fully relativistic version of subroutine fovrg.f
!     input:
!        ncycle  times to calculate photoelectron wave function
!                with nonlocal exchange
!        ikap    quantum number kappa for photoelectron
!        rmt     muffin-tin radius
!        jri     first interstitial grid point (imt + 1)
!        jlast   last point for integration of Dirac eq.
!        p2      current complex energy
!        dx      dx in loucks' grid (usually .05)
!        ri(nr)  loucks' position grid, r = exp ((i-1)*dx - 8.8)
!        vxc(nr) coulomb+xc potential for total density
!        vxcval  coulomb+xc potential for valence density
!        both vxc and vxcval include coulomb and nuclear potential
!        dgcn(dpcn) large(small) dirac components for 'iph' atom
!        adgc(adpc) their development coefficients
!     work space:
!        must be dimensioned in calling program.  coded like this
!        to make using different r-grids with different nrmax easy.
!
!     output:
!        pu, qu  upper and lower components at muffin tin
!        ps and qs are  upper and lower components for photoelectron
      use dimsmod, only: nrptx, ltot
	  use constants
      implicit double precision (a-h, o-z)
!Changed the dimensions to 41 to account for superheavy elements. Pavlo Baranov 07/2016

      integer, intent(in) :: ncycle, ikap, jri, iz, ihole, iunf, irr, ic3, iph
      integer, intent(out) :: jlast
      real*8, intent(in) :: xnval(41)
      complex*16, intent(out) :: vxc(nrptx), vxcval(nrptx)
      real*8, intent(in) :: ri(nrptx)
      complex*16, intent(out) :: ps(nrptx), qs(nrptx), pu, qu
!     all atoms' dirac components and their development coefficients
      real*8, intent(in) :: dgcn(nrptx,41), dpcn(nrptx,41)
      real*8, intent(in) :: adgc(10,41), adpc(10,41)
      real*8, intent(in) :: dx, rmt
      complex*16, intent(in) :: p2

      complex*16 aps(10),aqs(10)
      complex*16 ph0, amp, vu, vm(nrptx)

!     iph atom's dirac components and their development coefficients
      common/dff/cg(nrptx,41), cp(nrptx,41), bg(10,41), bp(10,41),      &
     &             fl(41), fix(41), ibgp
!     fl power of the first term of development limits.
!     ibgp first dimension of the arrays bg and bp (=10)

      complex*16 gg,gp,ag,ap,dv,av,bid
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),            &
     &              dv(nrptx),av(10),bid(2*nrptx+20)
!      gg,gp are the input and output for solout
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabc/afgkc
      real*8 afgkc(-ltot-1:ltot,41,0:4)
      common/messag/dlabpr,numerr
      character*8 dlabpr
!      xnel here - number of core electrons only
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      ! JK - 435 below may need to be 820 for superheavies.
      common/scrhf1/eps(820),nre(41),ipl
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim


!     initialize the data and test parameters
      ndor = 3
      cl = alpinv
      if (irr.gt.0) then
!        for irregular solution
         ndor=2
         aps(1) =  pu
         aqs(1) =  qu
         do 5 i=1, jri
           gg(i) = ps(i)
           gp(i) = qs(i)
 5       continue
      endif
      do 9 i = jri+1,nrptx
         vxc(i)=vxc(jri+1)
 9       vxcval(i)=vxc(jri+1)
      ibgp=10
      numerr = 0
      nz = iz
      hx = dx
      idim= 1 + nint(250*0.05/dx)
      if (idim .gt. nrptx) idim = nrptx
      if (mod(idim,2) .eq. 0) idim=idim-1
      
!     numerical integration of Dirac eq. works if you have 6 grid points
!     for one period of oscillations, switch to analytical expression
!     for a steplike potential  at large distances
      aa = 0.5
!     if (irr.gt.0) aa = 0.05
      rwkb = aa / dx / sqrt(abs(2*p2+(p2/cl)**2))
      x0 = 8.8
      iwkb= (log(rwkb) + x0) / dx  +  2
      if (iwkb.gt.idim) iwkb = idim
      if (iwkb.lt. 10) iwkb = 10
      
!     copy information into common's of atomic code
      do 13 j=1,41
      do 13 i=1,10
         bg(i,j)=adgc(i,j) 
 13      bp(i,j)=adpc(i,j) 
      do 15 j=1,41
      do 15 i=1,idim
         cg(i,j)=dgcn(i,j) 
 15      cp(i,j)=dpcn(i,j) 

      call inmuac (ihole, xion, iunf, ikap, iph) !KJ iph
      nmax(norb)=jlast
      if (iwkb.ge. jlast-1) iwkb = idim
!     note that here norb correspond to photoelectron

!     calculate initial photoelectron orbital using lda
      call diff (vxc,ri,ikap,cl,hx,jri,vm)
      do 18 i = jri, nrptx
  18  vm(i)=0.0d0
      call wfirdc (p2,kap,nmax,vxc,ps,qs,aps,aqs,irr,ic3,vm,            &
     &             rmt,jri, iwkb)

      if (numerr .ne. 0) call par_stop('error in wfirdc')
      if (ncycle .eq. 0) go to 999

!     to get orthogonalized photo e w.f., use alternative exit below
!     in general it should not be orthogonolized. Use for testing only 
!     ala

!     further need only core electrons for exchange term
      do 40 i=1, norb-1
  40  xnel(i) = xnel(i) - xnval(i)
!     take vxcval at the origin as vxcval=vcoul +const1 + i*const2
      av(2)=av(2)+(vxcval(1)-vxc(1))/cl
      do 50 i=1,iwkb
  50  dv(i)=vxcval(i)/cl
!     keep dv=vxc/cl above iwkb

      nter=0
 
!     angular coefficients 
      call muatcc(xnval)

!     no orthogonalization needed. Looking for g.f., not w.f.
!     if (ipl.ne.0) call ortdac (ikap,ps,qs,aps,aqs)
!     ortdac orthogonalizes photoelectron orbital to core orbitals
!     have to use exchange 5 card to exit here; also want vxc=vxcval
!     if (ncycle .eq. 0) go to 999

!     iteration over the number of cycles
 101  continue
         nter=nter+1
!        calculate exchange potential
         jriwkb = min (jri, iwkb)
         call potex( ps, qs, aps, aqs, jriwkb, p2)

!        resolution of the dirac equation
         if (irr.lt.0) then
            call solout (p2, fl(norb), aps(1), aqs(1), ikap, rmt,       &
     &        jri, nmax(norb), ic3, vm, iwkb)
         else
            call solin (p2, fl(norb), pu, qu, ikap, rmt,                &
     &        jri, nmax(norb), ic3, vm, iwkb)
         endif

!     no orthogonalization needed. Looking for g.f., not w.f.
!        if (ipl.ne.0) call ortdac (ikap,gg,gp,ag,ap)

!        acceleration of the convergence 
         scc(norb)=1.0d0
         do 151 i=1,idim
            ps(i)=gg(i)
 151        qs(i)=gp(i)
         do 155 i=1,ndor
            aps(i) =ag(i) 
 155        aqs(i) =ap(i) 

      if (nter.le.ncycle) go to 101

 999  if (numerr .eq. 0) then
        if (irr.lt.0 ) then
!c        need pu, qu for regular solution
!c        want to have vxc(jri)-smooth and vxc(jri+1)=v_mt
!c        assume no exchange beyond jri 
           vu=vxc(jri+1)
           call flatv                                                   &
     &     (ri(jri), rmt, ps(jri), qs(jri), p2, vu, ikap, pu, qu)
           jlast = nmax(norb)
!          jlast might change on very rare occasion
        endif
      else
        call par_stop('error in dfovrg.f')
      endif

      return
      end

      subroutine flatv (r1, r2, p1, q1, en, vav, ikap, p2, q2)
      use constants
	  use dimsmod, only: ltot
      implicit double precision (a-h, o-z)
!     solution of Dirac equation for flat potential for ikap is known
!     exactly (see e.g. in Loucks T.L. eq. 4-19)
!     given p1 and q1 at point r1 this subrotuine finds p2, q2 at r2
!     for given energy(en) and average potential (vav)
!     en and vav in hartrees
      external besjn, atancc

      complex*16 p1, q1, en, vav, p2, q2
      complex*16 ck, xkr, jl(ltot+2), nl(ltot+2), a,b, factor 

!     initialize staff
      ck = sqrt(2*(en-vav) + (alphfs*(en-vav))**2)
      xkr = ck*r1
      if (ikap.lt.0) then
        isign = -1
        lp = -ikap - 1
        lq = lp + 1
      else
        isign = 1
        lp = ikap
        lq = lp - 1
      endif
      a = ck * alphfs
      factor = isign*a/(1+sqrt(1+a**2))

!     find a and b that p1 = r1*(a*jl+b*nl), q1=factor*r1*(a*jl'+b*nl')
      call besjn (xkr, jl, nl)
      a = isign*ck*xkr* (p1*nl(lq+1) - q1*nl(lp+1)/factor)
      b = isign*ck*xkr* (q1*jl(lp+1)/factor - p1*jl(lq+1))

!     get values at r2
      xkr = ck * r2
      call besjn (xkr, jl, nl)
      p2 =  r2 * (jl(lp+1)*a + nl(lp+1)*b)
      q2 =  r2* factor * (jl(lq+1)*a + nl(lq+1)*b)

      return
      end

