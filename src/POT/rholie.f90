!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rholie.f90,v $:
! $Revision: 1.6 $
! $Author: bmattern $
! $Date: 2012/07/01 16:29:35 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rholie ( ri05, nr05, dx, x0, ri, em,                   &
     &                  ixc, rmt, rnrm,                                 &
     &                  vtot, vvalgs, xnval, dgcn, dpcn, eref,          &
     &                  adgc, adpc, xrhole, xrhoce, yrhole, yrhoce, ph, &
     &                  iz, xion, iunf, ihole, lmaxsc, iph) !KJ iph

      use DimsMod, only: nrptx, lx
      use constants
      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

!     INPUT
!     dx, x0, ri(nr)
!                  Loucks r-grid, ri=exp((i-1)*dx-x0)
!     ne, em(ne)   number of energy points,  complex energy grid
!     ixc          0  Hedin-Lunqist + const real & imag part
!                  1  Dirac-Hara + const real & imag part
!                  2  ground state + const real & imag part
!                  3  Dirac-Hara + HL imag part + const real & imag part
!                  5  Dirac-Fock exchange with core electrons +
!                     ixc=0 for valence electron density
!     rmt          r muffin tin
!     rnrm         r norman
!     vtot(nr)     total potential, including gsxc, final state
!     dgcn(dpcn)   large (small) dirac components for central atom
!     adgc(adpc)   their development coefficients
!
!     OUTPUT
!     xrhole(0:lx)  integral over r of density function
!     xrhoce(0:lx)  the same integral for embedded atom only
!     yrhole(251,0:lx)   density function
!     yrhoce(251)        density function for embedded atom

!     max number allowed in xsect r-grid
      parameter (nrx = nrptx)

!     output
      complex*16, intent(inout) :: xrhole(0:lx)
      complex*16, intent(inout) :: xrhoce(0:lx)
      complex*16, intent(inout) :: yrhole(251,0:lx), yrhoce(251)
      complex*16, intent(inout) :: ph(lx+1)

      dimension ri(nrptx), ri05(251)
      dimension  vtot(nrptx), vvalgs(nrptx)
      complex*16 vtotc(nrptx), vvalc(nrptx)
      dimension xnval(41), dgcn(nrptx,41), dpcn(nrptx,41)
      dimension adgc(10,41), adpc(10,41)

!     energy grid in complex e-plane
      complex*16 em, eref

!     work space for dfovrg: regular and irregular solutions
      complex*16 pr(nrx), qr(nrx), pn(nrx), qn(nrx)

      complex*16  p2, xkmt, ck, xck
      complex*16  pu, qu
      complex*16  xfnorm, xirf
      complex*16  temp,  phx, tempc

      complex*16 jl,jlp1,nl,nlp1
      complex*16  xpc(nrx)

!     initialize
      lmax=lmaxsc
      if (lmax.gt.lx) lmax = lx
      do 20 i = 1, nrptx
         vtotc(i)=vtot(i)
         vvalc(i)= vvalgs(i)
  20  continue
!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt  = (log(rmt) + x0) / dx  +  1
      jri  = imt+1
      if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
      inrm = (log(rnrm) + x0) / dx  +  1
      jnrm = inrm+1

!     set limits for tabulations
      nr05= (log(rnrm) + x0) / 0.05d0 + 5
      if (nr05.gt.251) nr05 = 251
!     ilast is the last integration point
!     it is larger than jnrm for better interpolations
      ilast = nint( (nr05-1) *0.05d0 / dx ) + 1
      if (ilast.gt.nrptx) ilast=nrptx

      do 10 lll = 0, lx
      do 10 j = 1, 251
         yrhole(j,lll) = 0
  10  continue
      do 30 j = 1, 251
  30  yrhoce(j) = 0

!     p2 is 0.5*(complex momentum)**2 referenced to energy dep xc
!     need hartree units for dfovrg
      p2 = em - eref
      if (mod(ixc,10) .lt. 5) then
        ncycle = 0
      else
        ncycle = 3
      endif
      ck = sqrt(2*p2 + (p2*alphfs)**2)
      xkmt = rmt * ck

      do 200 lll=0,lx
        if (lll.gt.lmax) then
           ph(lll+1) = 0
           xrhoce(lll) = 0
           xrhole(lll) = 0
           do 110 i = 1,251
  110      yrhole(i,lll) = 0
           goto 200
        endif

!       may want to use ihole=0 for new screening. 
!       don't want ro use it now
!       ihole = 0
        ikap = -1-lll
        irr = -1
        ic3 = 1
        if (lll.eq.0) ic3 = 0
        call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,            &
     &                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,         &
     &                xnval, pu, qu, pn, qn,                            &
     &                iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        call exjlnl (xkmt, lll, jl, nl)
        call exjlnl (xkmt, lll+1, jlp1, nlp1)
        call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,         &
     &                  phx, temp)
        ph(lll+1)=phx
!     Normalize final state  at rmt to
!     rmt*(jl*cos(delta) - nl*sin(delta))
        xfnorm = 1 / temp
!     normalize regular solution
        do 133  i = 1,ilast
          pr(i)=pn(i)*xfnorm
          qr(i)=qn(i)*xfnorm
  133   continue

!      find irregular solution
        irr = 1
        pu = ck*alphfs
        pu = - pu/(1+sqrt(1+pu**2))
!       set pu, qu - initial condition for irregular solution at ilast
!       qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt
!       pu = (nl*cos(phx)+jl*sin(phx)) *rmt
        qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt 
        pu = (nl*cos(phx)+jl*sin(phx)) *rmt 

        call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,             &
     &              ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,            &
     &              xnval, pu, qu, pn, qn,                              &
     &              iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
!c      set N- irregular solution , which is outside
!c      N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
!c      N = i*R - H*exp(i*ph0)
        temp = exp(coni*phx)
!       calculate wronskian
        qu = 2 * alpinv * temp * ( pn(jri)*qr(jri) - pr(jri)*qn(jri) )
        qu = 1 /qu / ck
!       qu should be close to 1
        do i = 1, ilast
          pn(i) = coni * pr(i) - temp * pn(i)*qu
          qn(i) = coni * qr(i) - temp * qn(i)*qu
        enddo

!     ATOM,  dgc0 is large component, ground state hole orbital
!     .      dpc0 is small component, ground state hole orbital
!     FOVRG, p    is large component, final state photo electron
!     .      q    is small component, final state photo electron

            
!    combine all constant factors to temp
!    add relativistic correction to normalization and factor 2*lll+1
        pu = ck*alphfs
        pu = - pu/(1+sqrt(1+pu**2))
        temp = (2*lll+1.0d0)/(1+pu**2) /pi *ck * 2
!    also scale by appropriate step in complex energy
        do 190  i = 1, ilast
          xpc(i) = pr(i) * pr(i) + qr(i) * qr(i) 
 190    continue
          
        do 191 ir=1,nr05
           call terpc(ri, xpc, ilast, 3, ri05(ir), tempc)
           tempc = tempc * temp
           yrhole(ir,lll)= tempc
 191    continue

        xirf = lll*2 + 2
!       i0 should be less or equal to  ilast
        i0=jnrm+1
        call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
!       print out xirf for Bruce
        xrhole(lll) = xirf*temp

!     only central atom contribution needs irregular solution
        do 195  i = 1, ilast
          xpc(i) = pn(i)*pr(i)-coni*pr(i)*pr(i)                         &
     &           + qn(i)*qr(i)-coni*qr(i)*qr(i)
!         yrhoce(i)=yrhoce(i) - temp*xpc(i)
 195    continue
        do 196 ir=1,nr05
           call terpc(ri, xpc, ilast, 3, ri05(ir), tempc)
           yrhoce(ir)=yrhoce(ir) - temp*tempc
 196    continue

        xirf =  1
        call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
        xrhoce(lll) =  - xirf* temp
 200  continue 

      return
      end
