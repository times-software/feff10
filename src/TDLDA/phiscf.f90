!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phiscf.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine phiscf (ifxc, rmt, jlast, jri, p2, edge, emu, dx,      &
     &      ri, vxc, edens, dgcn, dpcn, adgc, adpc, iz, ihole, neg, eng,&
     &      rhoj, kappa, norbp, fscf, yvec, maxsize, matsize, sfun, iph) !KJ iph
!     Zangwill-Soven effective field calculation
!     coded by a.ankudinov 2003

!     input:
!        rmt     muffin-tin radius
!        jlast   last point for integration of Dirac eq.
!        jri     first interstitial grid point (imt + 1)
!        p2      current complex energy
!        edge    shifted Fermi level
!        dx      dx in loucks' grid (usually .05)
!        ri(nr)  loucks' position grid, r = exp ((i-1)*dx - 8.8)
!        vxc(nr) coulomb+xc potential for total density
!        dgcn(dpcn) large(small) dirac components for 'iph' atom
!        adgc(adpc) their development coefficients
!     output:
!        fscf  - self-consistent radiation field
      use dimsmod, only: nrptx,ltot,nex
      use errorfile
	  use constants
      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

      complex*16 vxc(nrptx), p2, p2p
      complex*16 vxcp(nrptx)
      dimension ri(nrptx), edens(nrptx), fxc(nrptx)
      complex*16 pu, qu, vm(nrptx)
      complex*16 ps(nrptx), qs(nrptx), aps(10),aqs(10), xnorm1
      complex*16  api(10),aqi(10)
!     DOS information
      dimension neg(41), eng(nex,41), rhoj(nex,41), kappa(41)

!     new arrays for TDLDA
!     solution of homogeneous equations (needed for normalization)
      complex*16 ph(nrptx), qh(nrptx)
      complex*16 pir(nrptx), qir(nrptx)

      complex*16 dum1, factor, fscf(nrptx), yvec(nrptx, maxsize)
      complex*16 wronsk, hpp, jpp
      complex*16 chik(251,251), ww
      complex cchik(251,251)

!     bessel functions
      complex*16 ck, xkmt, ff, tl
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 j0(ltot+2), n0(ltot+2)
      complex*16 jlp(0:ltot+1), hl(0:ltot+1)
      complex*16 jl0(0:ltot+1), hl0(0:ltot+1)


!     all atoms' dirac components and their development coefficients
      dimension dgcn(nrptx,41), dpcn(nrptx,41)
      dimension adgc(10,41), adpc(10,41)
 
!     iph atom's dirac components and their development coefficients
      common/dff/cg(nrptx,41), cp(nrptx,41), bg(10,41), bp(10,41), fl(41), fix(41), ibgp
!     fl power of the first term of development limits.
!     ibgp first dimension of the arrays bg and bp (=10)

      complex*16 gg,gp,ag,ap,dv,av,bid
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),av(10),bid(2*nrptx+20)
!      gg,gp are the input and output for solout
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/messag/dlabpr,numerr
      character*8 dlabpr
!      xnel here - number of core electrons only
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),nq(41),kap(41),nmax(41)
      ! JK - 435 below may need to be 820 for superheavies.
      common/scrhf1/eps(820),nre(41),ipl
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idm

      external besjn
     ! call wlog('in phiscf')
      iholep = 0
      xion = 0
      iunf = 1
      nmax(:) = 0
      cl = alpinv
      do 9 i = jri+2,nrptx
 9      vxc(i)=vxc(jri+1)
      ibgp=10
      numerr = 0
      nz = iz
      hx = dx
      idm= 1 + nint(250*0.05/dx)
      if (idm .gt. nrptx) idm = nrptx
      if (mod(idm,2) .eq. 0) idm=idm-1
! initialize
      fscf(:) = 1
      cchik(:,:) =  0

!      itest = 1
!      if (itest .eq.1) goto 999

!     xc model: ifxc=0 --> RPA; ifxc=1 -->Znagwil soven adiabatic fxc
      do i= 1, nrptx
!       also calculate the local exchange term fxc = vxc*r_s**3/r**2
        if (edens(i).le.0) then
          rs = 100
        else
          rs = (4*pi*edens(i)/3)**(-third)
        endif
!       vvbh from Von Barth Hedin paper, 1971
!       see eq.60-61 of Gross&Kohn, Adv. Quant. Chem. 21, p.255(1990).
!       eq.60 fxc = d V_xc / d rho * (2\ell + 1) /(4*pi*r**2)
!       where second factor comes from \delta(r-r')
!       for dipole field \ell=1 
!       fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.22177412/rs -1.512/(30+rs))
!c      below are coefficients in Zangwill/Soven paper
        fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.222/rs -0.75924/(11.4+rs))
        if (ifxc.eq.0) fxc(i) = 0
      enddo


!     numerical integration of Dirac eq. works if you have 6 grid points
!     for one period of oscillations, switch to analytical expression
!     for a steplike potential  at large distances
      aa = 0.5d0
      ck =  sqrt (2*p2+ (p2*alphfs)**2)
      rwkb = aa / dx / abs(ck) 
      x0 = 8.8d0
      do 12 i =1, nrptx
12    dr(i) = exp(-x0 + (i-1)*dx)
      iwkb= (log(rwkb) + x0) / dx  +  2
      if (iwkb.gt.idm) iwkb = idm
      if (iwkb.lt. 10) iwkb = 10
      
!     copy information into common's of atomic code
      do 13 j = 1, 41
      do 13 i = 1, 10
         bg(i,j) = adgc(i,j) 
 13      bp(i,j) = adpc(i,j) 
      do 15 j = 1, 41
      do 15 i = 1, nrptx
         cg(i,j) = dgcn(i,j) 
 15      cp(i,j) = dpcn(i,j) 

      ikap = 0
      call inmuac (iholep, xion, iunf, ikap, iph) !KJ iph
      nmax(norb)=jlast
      if (iwkb.ge. jlast-1) iwkb = idm
!     note that here norb correspond to photoelectron

!     calculate initial photoelectron orbital using lda
      vm(:)=0.0d0
      ic3 = 0

      if (ihole.ne.0) then
!       true x-ray photon energy, calcualtions are done for
        wp = dble(p2 + emu - edge)
!       wp = p2 + emu - edge -0.5
!       single electron energy for which chi0 is calculated
        ww = p2 - eng(1,ihole)

!       wp = 1
!       first power is needed to satisfy f-sumrule
        wp = (dble(ww)/wp)
!       wp = (dble(ww)/wp)**2

!       cycle over atomic orbitals and energy points
        do 155 iorb = 1, norbp
        do 155 ieg = 1, neg(iorb)
          kinit = kappa(iorb) 
          xx = rhoj(ieg, iorb)

!         cycle over 2 poles
          do 150 ind = 1, 2
!           print*, iorb, dble(ww)*hart, dble(eng(ieg,iorb)) * hart
            p2p = eng(ieg,iorb) + ww
            if (ind.eq.2)  p2p = eng(ieg,iorb) -dble(ww) +coni*dimag(ww)
            if (dble(p2p).lt.edge) then
!              no contibution to Im chi0. Only Re chi0 is important
!              use large broadening to remove  possible resonant behaviour
!              when E_i - E_j = omega for occupied i and j
               yy = eng(ieg, iorb) + dble(ww) 
               if (ind.eq.2) yy = eng(ieg, iorb) 
!              gamb = max( dimag(ww) , (edge - yy)/100)
               gamb = max( dimag(ww) , (edge - yy)/10)
!              gamb = dimag(ww) 
               p2p = dble(p2p) + coni*gamb
            endif
             
!           set complex momentum and call bessel functions
            ck =  sqrt (2*p2p+ (p2p*alphfs)**2)
            rmtx =  10.d0 / abs(dimag(ck))
            if (rmtx.gt.rmt) rmtx = rmt
            jrip = (log(rmtx) + x0) / dx + 2
            rmtp = ri(jrip) - 1.d-20
            do 133 i = 1, jrip
133         vxcp(i) = vxc(i)
            do 134 i = jrip+1, nrptx
134         vxcp(i) = vxc(jri+1)
            xkmt = rmtp * ck
            call besjn (xkmt, jl, nl)
            ix = ltot
            xkmt = ri(jrip) * ck
            call besjh (xkmt, ix, jl0, hl0)

!           set iwkb
            aa = 0.5d0
            rwkb = aa / dx / abs(ck) 
            iwkb= (log(rwkb) + x0) / dx  +  2
            if (iwkb.gt.idm) iwkb = idm
            if (iwkb.lt. 10) iwkb = 10
            if (iwkb.ge. jlast-1) iwkb = idm

!           cycle over dipole selection rules
            do 149 ik = -1,1
              kfin = kinit + ik
              if (ik.eq.0) kfin = -kfin
              if (kfin.eq.0) goto 149
!  test : final f states only
!             if (kfin.ne.-4 .and. kfin.ne.3 .and.kinit.ne.-4 .and.
!    1            kinit.ne.3) goto 149
!             if (kfin.ne.-3 .and. kfin.ne.2) goto 149
!  end test

              kap(norb) = kfin
              irr = -1
              call wfirdc (p2p, kap, nmax, vxcp, ps, qs,aps,aqs,irr,ic3,vm, rmtp, jrip, iwkb)
              do 130 i=1, idm
                ph(i) = ps(i)
                qh(i) = qs(i)
  130         continue
!             wronskian normalization of regular solution is below
!             find irregular solution first
              il = abs(kfin) 
              if (kfin.lt. 0) il = il - 1
              ilp = il - 1
              if (kfin .lt. 0) ilp = il + 1
              sign = -1.0
              if (kfin.gt.0) sign = 1.0
              factor = ck*alphfs
              factor = sign * factor/(1+sqrt(1+factor**2))
              dum1 = 1/ sqrt(1+factor**2)
              api(1) = hl0(il)  * rmtp * dum1
              aqi(1) = hl0(ilp) * rmtp * dum1 * factor
              irr = 1
              call wfirdc (p2p,kap,nmax,vxcp,pir,qir,api,aqi,irr,ic3,vm,rmtp, jrip, iwkb)

!             set irregular solution outside jrip
              il = abs(kfin)
              if (kfin.lt. 0) il = il - 1
              ilp = il - 1
              if (kfin .lt. 0) ilp = il + 1
              ilx = il+1
              if (ilp.gt.il) ilx=ilp+1
              ilx = ltot

!             get irregular solutions at jrip+1
              i = jrip+1
              xkmt = ck * ri(i)
              call besjh( xkmt, ilx, jlp, hl)
              dum2 = ri(i) / ri(jrip)
              pir(i) = pir(jrip) * dum2 * hl(il) /hl0(il)
              qir(i) = qir(jrip) * dum2 * hl(ilp)/hl0(ilp)

!             calculate wronskian at jrip
              wronsk= 2*alpinv* (pir(jrip)*qh(jrip) -ph(jrip)*qir(jrip))
              wronsk = 2 /wronsk

!             put wronskian normalization into irregular solution
              do i=1, jrip
                ph(i) = ph(i) * wronsk
                qh(i) = qh(i) * wronsk
              enddo
              xkmt = ck * ri(jrip)
              tl = (ph(jrip)/(2*xkmt) - jl0(il)) /hl0(il)
              do i = jrip+1, idm
                xkmt = ck * ri(i)
                call besjh( xkmt, ilx, jlp, hl)
                dum2 = ri(i) / ri(jrip)
                pir(i) = pir(jrip) * dum2 * hl(il) /hl0(il)
                qir(i) = qir(jrip) * dum2 * hl(ilp)/hl0(ilp)
                ph(i) = 2*xkmt* (jlp(il)+tl*hl(il))
                qh(i) = 2*xkmt* (jlp(ilp)+tl*hl(ilp))
              enddo

!             make K(r,r') * chi0(r'r'') = chik(r,r")
              call lipman (iorb, ph,qh, pir,qir, fxc, jrip, nmax(norb), chik)
              jfin2 = 2*abs(kfin) - 1
              jin2 = 2*abs(kinit) - 1
              aa= -cwig3j( jfin2,2,jin2, 1,0,2)**2 *(jfin2+1)*(jin2+1)/3
!             xx - fraction of the shell occupation
              aa =  aa * xx
!             wp - correction due to inequality of photon energy
!             and single electron energy differences
              aa = aa * wp
!             separation function (between ZS and PM)
              aa = aa * sfun
         
              do 161 i = 1,251
              do 161 j = 1,251
               itest = 0
               if (itest.eq.0) then
                 if (ind.eq.1 .and. dble(p2p).gt. edge) then
                  cchik(i,j) = cchik(i,j) + real(aa) *  cmplx( real(dble(chik(i,j))), real(dimag(chik(i,j))))
                 else
                  cchik(i,j) = cchik(i,j) + real(aa) *  cmplx( real(dble(chik(i,j))), 0 )
                 endif
               else
                 if (ind.eq.1) then
                  cchik(i,j) = cchik(i,j) + real(aa) *  cmplx( real(dble(chik(i,j))), real(dimag(chik(i,j))))
                 else
                  cchik(i,j) = cchik(i,j) + real(aa) *  cmplx( real(dble(chik(i,j))), -real(dimag(chik(i,j))))
                 endif
               endif
  161         continue
  149       continue
  150     continue
  155   continue

        itest = 0
        if (itest.eq.0) then
!         invert matrix and multiply by r
          nx = idm/5 + 1
          call chiklu(nx, dr, cchik, fscf, yvec, maxsize, matsize)
        elseif (itest.eq.1) then
!         do one iteration test
          do 191 i = 1,251
            i1 = 1+5*(i-1)
!           qir(i) =  fscf(i1) * dr(i1)
            qir(i) =  dr(i1)
            do 192 j=1,251
              j1 = 1+5*(j-1)
              qir(i) = qir(i) + cchik(i,j)*dr(j1)
!             qir(i) = qir(i) - cchik(i,j)*fscf(j1)*dr(j1)
  192       continue
            fscf(i1) = qir(i) / dr(i1)
  191     continue
        else
!         get f' and f" = r*chi0*r'
!         replace lipman.f with lipman.new also
          ff = 0
          do 194 i = 1,251
            i1 = 1+5*(i-1)
            qir(i) =  0
            do 193 j=1,251
              j1 = 1+5*(j-1)
!             qir(i) = qir(i) + chik(i,j)*dr(j1)
              qir(i) = qir(i) + cchik(i,j)*dr(j1)
  193       continue
            if (i.gt.1) then
              ff = ff + (qir(i) + qir(i-1)) * (dr(i1) - dr(i1-5))/2
            endif
  194     continue
!          do 196 i = 1,251
!           i1 = 1+5*(i-1)
!           write(56, *) dr(i1) , real(cchik(i,i)), aimag(cchik(i,i))
! 196      continue
          write(56, *) dble(ww)*hart, dble(ff), dimag(ff)
!         print*, dble(ww)*hart, dble(ff), dimag(ff)
        endif

!       write out scf field
        p2ev = dble(p2)*hart
!       if (p2ev.gt.14.d0 .and. p2ev.lt.15.d0) then
!       if (p2ev.gt.27.d0 .and. p2ev.lt.28.d0) then
!       if (p2ev.gt.41.d0 .and. p2ev.lt.42.d0) then
!       if (p2ev.gt.54.d0 .and. p2ev.lt.56.d0) then
        if (p2ev.gt.99.d0 .and. p2ev.lt.101.d0) then
!          print*, 305, idm
!c        write out scf field Re and Im parts
          do 305 i = 1,251
            i1 = 1 +5*(i-1)
            write(43,306) dr(i1), dr(i1)*dble(fscf(i1)),                &
     &        dr(i1)*dimag(fscf(i1))
!           write(43,306) dr(i1), dble (dgcn(i1, ihole))
  306       format(3f11.5)
! 306       format(3e14.5)
  305     continue
          close (unit=43)
          itest = 0
          if (itest.gt.1) then
		     stop 'test finished'
			 call WipeErrorfileAtFinish
		  endif
        endif
      endif
  !    call wlog('done with phiscf')

 999  continue
      return
      end
