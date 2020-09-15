!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xsect.f90,v $:
! $Revision: 1.17 $
! $Author: jorissen $
! $Date: 2012/03/31 04:30:33 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Josh - argument iPl has been added to arguments of xsect
!     Josh - added NPoles and Eps0 3/8/2010
      subroutine xsect (ipr2, dx, x0, ri, ne, ne1, ik0, em, edge,       &
                       ihole, emu, corr, dgc0, dpc0, jnew,             &
                       ixc, lreal, rmt, rnrm, xmu,                     &
                       vi0, iPl, NPoles, Eps0, EGap, gamach,           &
                       vtot, vvalgs, edens, dmag, edenvl,              &
                       dgcn, dpcn, adgc, adpc, xsec, xsnorm, rkk,      &
                       iz, xion, iunf, xnval,                          &
                       izstd, ifxc, eorb, kappa, iorb, l2lp,           &
                       ipol, ispin, le2, angks, ptz, iph) !KJ

      USE IOMOD
      use constants
      use DimsMod, only: nex, nrptx, lx, MxPole, nspx=>nspu
      USE SelfEnergyMod
!     right now the same self-energy is used for calculation of the central atom part (xsec) and dipole m.e. for
!     scattering (rkk). You may want to run xsect separately for xsec and for rkk, if you want to use different self-energy
!     for central and scattering parts.  ala. fix later

      !implicit double precision (a-h, o-z)
      implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      integer, intent(in) :: ihole, ipr2, ne, ne1, ik0, l2lp, le2, ispin, ipol, iPl, &
        ixc, lreal, iunf, izstd, iz, ifxc
      integer, intent(inout) :: NPoles, iph  !KJ I hope that is as intended??
      real*8, intent(in) :: dx, x0, rmt, rnrm, xmu, vi0, gamach, corr, xion, emu
!     INPUT
!     dx, x0, ri(nr)
!                  Loucks r-grid, ri=exp((i-1)*dx-x0)
!     ne, em(ne)   number of energy points, real energy grid
!     edge         chemical potential (energy for k=0)
!     ihole        hole code
!     emu          position of chemical potential in absorption specrum
!     dgc0(nr)     dirac upper component, ground state hole orbital
!     dpc0(nr)     dirac lower component, ground state hole orbital
!     ixc          0  Hedin-Lunqist + const real & imag part
!                  1  Dirac-Hara + const real & imag part
!                  2  ground state + const real & imag part
!                  3  Dirac-Hara + HL imag part + const real & imag part
!                  5  Dirac-Fock exchange with core electrons +
!                     ixc=0 for valence electron density
!     lreal        logical, true for real phase shifts only
!     rmt          r muffin tin
!     xmu          fermi level
!     vi0          const imag part to add to complex potential
!     gamach       core hole lifetime
!     vtot(nr)     total potential, including gsxc, final state
!     edens(nr)    density, hole orbital, final state
!     dmag(251)     density magnetization
!     edenvl      valence charge density
!     dgcn(dpcn)   large (small) dirac components for central atom
!     adgc(adpc)   their development coefficients
!
!     OUTPUT
!     xsec(ne)    atomic absorption cross section to multiply \chi
!                 (atomic background for XMCD)
!     xsnorm(ne)  atomic  absorption cross section (norm for XMCD)
!     rkk(ne, 8)  normalized reduced matrix elements for construction
!                 of termination matrix in genfmt.

      complex*16, intent(in) :: ptz(-1:1, -1:1)
      complex*16, intent(in) :: em(nex)
      real*8, intent(in) :: ri(nrptx), vtot(nrptx), edens(nrptx),dmag(nrptx)
      real*8, intent(in) :: dgc0(nrptx), dpc0(nrptx), vvalgs(nrptx), edenvl(nrptx)
      real*8, intent(in) :: dgcn(nrptx,41), dpcn(nrptx,41), eorb(41)
      real*8, intent(in) :: adgc(10,41), adpc(10,41), xnval(41)
      integer, intent(in) :: iorb(-5:4),kappa(41)
      complex*16, intent(out) :: rkk(nex, 8), xsec(nex)
      real*8, intent(out) :: xsnorm(nex)

      integer kiind(8), lind(8) !KJ renamed kind to kiind to avoid conflict with reserved name and its use in modules
      real*8 xp(nrptx), xq(nrptx)

!     work space for xcpot
      real*8 vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      real*8 vvxcrm(nrptx), vvxcim(nrptx)

!     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx), fscf(nrptx)
      complex*16 pp(nrptx), qp(nrptx), pnp(nrptx), qnp(nrptx)
!     storage for calculation of cross term (SPIN 1 only)
      complex*16 xrcold(nrptx) , xncold(nrptx), yvec(nrptx,1)

      complex*16  p2, ck, xkmt, xkmtp
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf, xirf1
      complex*16  temp, aa, bb, cc, rkk1, rkk0, phold
      complex*16  phx(8), ph0
      complex*16  eref, xm1, xm2, xm3, xm4

      complex*16 jl,jlp1,nl,nlp1
      complex*16  v(nrptx), vval(nrptx)
      complex*16  xrc(nrptx), xnc(nrptx)
      character*512 slog
      logical ltrace
      complex*16 xrhoce(nex), xrhopr(nex), chia(nex), cchi(nex)
      real*8 omega1(nex), bf(0:2, nrptx)

      real*8 pat(nrptx),qat(nrptx)
      complex*16 intr(nrptx),var(nrptx) 
!     to pass energy levels and projected DOS
      integer neg(41)
      real*8 eng(nex, 41), rhoj(nex,41)
!     Josh - Added iPl switch for PLASMON card
!          - and WpCorr = Wi/Wp, Gamma, AmpFac
!          - to describe Im[eps^-1]
      integer ipole, NData
      real*8, allocatable :: Energy(:), Loss(:)
      real*8 WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole), Eps0, EGap
      character(LEN=10) ColumnLabels(20)
!     Josh END
      integer imt, jri, jri1, inrm, jnrm, i, index, ifirst, ie, ncycle, ilast, ll, &
         kinit, linit, isp, mult, kx, ks, id, ifl, j, ilp, lfin, ic3p, ic3, iold, &
         ikap, ind, kdif, ios, i0, k1, kdif1, jproj, isel, irr, matsize, maxsize, &
         norbp, jnew
      real*8 xinorm, del, angks, p2f, edge, omega, xk0, x, sinx, cosx, ww, wse, &
         dum, sign, ratio, ratiop, xsedge, edg50, prefac, vrcorr, vicorr, xnorm, sfun

      complex*16, allocatable :: bmat(:,:,:,:,:,:)
      logical,parameter :: correction_josh = .false.
      logical CorrectOrbitalEnergies ! Set to true until this is done, then set to false.
!     Debug: FDV
!     dimension sh_eng(nex, 41)
      
      CorrectOrbitalEnergies = .TRUE.
      allocate(bmat(-lx:lx,0:1,8, -lx:lx,0:1,8))

      call setkap(ihole, kinit, linit)

!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

! Debug: FDV
!     call correorb(iz, ihole, rmt, jri, dx,ri, p2f,edge, v, dgcn, dpcn, adgc, adpc, eorb, neg, sh_eng, rhoj, kappa, norbp)
!     e_chsh = -sh_eng(1,1)
!     write(6,fmt='(a,f12.8)') 'e_chsh: ', e_chsh

!     nesvi: define jnrm
      inrm = (log(rnrm) + x0) / dx + 1
      jnrm = inrm + 1
!write(*,*) 'x0,dx,rnrm,inrm,jnrm',x0,dx,rnrm,inrm,jnrm
!     We'll need <i|i> later to normalize dipole matrix elements
!     <i|r|f>.  NB, dgc and dpc are r*wave_fn, so use '0' in somm to get integral  psi**2 r**2 dr.
!     Square the dgc0 and dpc0 arrays before integrating.
!     <i|i> == xinorm.
!     dgc and dpc should be normalized <i|i>=1, check this here
      do i = 1, nrptx
         xp(i) = dpc0(i)**2
         xq(i) = dgc0(i)**2
      enddo
!     nb, xinorm is used for exponent on input to somm
      xinorm = 2*linit + 2
      call somm (ri, xp, xq, dx, xinorm, 0, jnrm)
      del = abs (abs(xinorm) - 1)
      if (del .gt. 1.e-2) then
         write(slog,'(a,i8,1p2e13.5)') ' ihole, xinorm ', ihole , xinorm
         call wlog(slog)
!        if using real phase shifts, don't expect great results
         if (lreal.lt.2)  then
           call wlog(' There may be convergence problems.')
           call wlog(' Xinorm should be 1. If you set the RGRID, minor ' // &
             'interpolation errors that will not affect final results may occur.')
         endif
      endif

!     use ixc for testing
      index = ixc
!     Always use ground state self energy for xsection, quick fix
!     JJR, Jan 93
!     change for testing broadened plasmon pole 6/93
!     index = 2
!     ALA found that it is better to use index=ixc and real part of
!     self-energy for atomic xsection. 12/96
      ltrace = .true.
      call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, kiind, lind, bmat)
!     set spin index to use bmat
      isp = 0
      if (ispin.eq.1) isp = nspx - 1

!     zero rkk and phx
      phx(:)=dcmplx(0)
	  rkk(:,:)=dcmplx(0)

      ifirst = 0
!     Josh - if PLASMON card is set, and using HL exc,
!          - read loss function from loss.dat and find poles etc.
      ! Josh - Changing to read directly from loss.dat : 3/9/2010
      IF ( (iPl.gt.0).and.(ixc.eq.0) ) THEN
         WpCorr(:) = -1.d30
         CALL OpenFl('loss.dat', FileStatus = 'OLD')
         NData = NumberOfLines('loss.dat')

         ALLOCATE(Energy(NData), Loss(NData))

         CALL ReadArrayData('loss.dat', Double1 = Energy, Double2 = Loss)
  
         !NPoles = 100 ! 100 poles should be enough for anything. Can add a card for this later.
         !Eps0 = -2.d0 ! This will not do anything. Will add a card to set Eps0 later.
         CALL MkExc(Energy, Loss, Eps0, WpCorr, AmpFac, NPoles)

         Gamma(:)  = Gamma(:)/hart
         WpCorr(:) = (WpCorr(:)/hart) / SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)
         CALL CloseFl('loss.dat')
         DEALLOCATE(Energy,Loss)
      END IF
      ! Write wp as calculated from density to sigma.dat
      ColumnLabels(:) = ' '
      ColumnLabels(1) = 'Rs'
      ColumnLabels(2) = 'wp (eV)'
      CALL WriteData('mpse.dat',                                                       &
              Double1 = (3 / (4*pi*edens(jri+1))) ** third,                          &
              Double2 = SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart,     &
              ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE.,               &
              Headers = (/ 'This file contains information about the self-energy.' /))

!     Josh END


!     ********************* START OF BIG LOOP OVER ENERGY MESH **************
      do ie = 1, ne
         iph = 0
!        Josh - xcpot now has new arguments:
!             - iPl, WpCorr, Gamma, AmpFac, EGap         
         call xcpot (iph, ie, index, lreal, ifirst, jri, em(ie), xmu, vtot, vvalgs, edens, dmag, edenvl,  &
                     eref, v, vval, iPl, WpCorr, Gamma, AmpFac, EGap, vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)

!       set the method to calculate atomic cross section
!       p2 is (complex momentum)**2 referenced to energy dep xc
        p2 = em(ie) - eref
        p2f = edge - dble(eref)
        ck = sqrt (2*p2 + (p2*alphfs)**2)
        xkmt = rmt * ck
        if (mod(index,10) .lt. 5) then
           ncycle = 0
        else
!          fix later. maybe ncycle can be less
           ncycle = 3
        endif
        omega = (dble(em(ie)) - edge) + emu
        omega = max (omega, 0.001d0 / hart)
!       nesvi: add omega1(ie)- need it later
        omega1(ie) = omega

!       remember the bessel functions for multipole matrix elements
        xk0 = omega * alphfs
        ilast = jnrm+6
        if (ilast.lt.jnew) ilast = jnew
        if (ilast.gt.nrptx) ilast = nrptx
        do i = 1, ilast
          temp = xk0 * ri(i)
          if (abs(temp).lt.1.d0) then
!           use series expansion
            do ll = 0,2
              call bjnser(temp,ll, xirf, dum1,1)
              bf(ll,i) = dble(xirf)
            enddo
          else
!           use formula
            x = dble(temp)
            sinx = sin(x)
            cosx = cos(x)
            bf(0,i) = sinx/x
            bf(1,i) = sinx/x**2 - cosx/x
            bf(2,i) = sinx*(3/x**3-1/x) - 3*cosx/x**2
          endif
        enddo

!       notice for spin-dep case xsnorm and xsec are spin-dep
!       and kept separately (see call to xsect in subroutine xsph)
        xsnorm(ie) = 0 
        xsec(ie) = 0
        if (dble(em(ie)).lt.-10.d0) cycle !goto 400
        if (dimag(p2).le.0.d0 .and. dble(p2).le.0.d0) cycle !goto 400

!       matrix elements for multipole (E1,E2,M1) transitions
!       The terms up to (k/c)^2 for absorption are kept.
!       L3 edge: kdif=1 (3d5/2)      kdif=2 (3d3/2), kdif=3(4s)
!       L2 edge: kdif=1 (no transition), 2 (4s),      3 (3d3/2)
        do 350 mult = 0, 2
          if (mult.eq.0) then
!           E1 transitions
            kx = 1
            ks = 2
          else
!           M1 transitions
            kx = 1
            ks = 6
!           E2 transitions
            if (mult.eq.2) kx = 2
          endif 
!         skip unnecessary calculations
          if (mult.gt.0 .and. (mult.ne.le2)) goto 350
 
!         set ilast larger than jri for better interpolation for pu
!         also need 5 points after jri for irregular solution
          ilast = jnrm + 6
          if (ilast.lt.jnew) ilast = jnew

!c        calculate screened dipole field
          ww = dble(emu+p2-edge)

          if (mult.eq.0 .and. izstd.gt.0) then
            ! In f' calculations this gets skipped by line 262
            ! Change to logical first.
            !if (ie.eq.1) then
            if(CorrectOrbitalEnergies) then
               CorrectOrbitalEnergies = .FALSE.
               call correorb(iz, ihole, rmt, jri, dx,ri, p2f,edge, v, dgcn, &
                 dpcn, adgc, adpc, eorb, neg, eng, rhoj, kappa, norbp, iph) !KJ
               if(correction_josh) eng(1,:) = eng(1,:) + p2f - edge ! Josh added this line so that
               ! correorb calculates orbital energies relative to absolute zero.
            end if
            maxsize = 1
            matsize = 0
            sfun = 1.d0
            call phiscf (ifxc, rmt, ilast, jri, p2, p2f, emu, dx, ri, v, &
              edens, dgcn, dpcn, adgc, adpc, iz, ihole, neg, eng, rhoj, &
              kappa, norbp, fscf, yvec, maxsize, matsize, sfun, iph) !KJ
            wse = dble(p2-eng(1,ihole))
          else
            fscf(:) = 1.d0
            wse = ww
          endif
if(.false.) then
write(*,*) 'ie',ie
write(*,*) 'p2',p2
write(*,*) 'p2f',p2f
write(*,*) 'em(ie)',em(ie)
write(*,*) 'eref',eref
write(*,*) 'alphfs',alphfs
write(*,*) 'ck',ck
write(*,*) 'omega',omega
write(*,*) 'eng(1,ihole)',eng(1,ihole)
write(*,*) 'edge',edge
write(*,*) 'emu',emu
write(*,*) 'ww',ww
write(*,*) 'wse',wse
endif


!         ww = 1
!         ww = wse / ww
          ww = sqrt(wse/ww)

          do 300 kdif = -kx, kx
            if (omega.le.0.0) goto 300
            ind = kdif + ks
            ikap = kiind (ind)
            if (ikap .eq. 0) goto 300
!           use l2lp =0 to include both transitions l-->l+/-1
!           if (l2lp.ne.0) only dipole transitions are calculated.
!            l-->l+1 transitions
            if (l2lp.eq. 1 .and. ((kinit.lt.0 .and. ind.ge.3) .or. (kinit.gt.0 .and. ind.ne.3)) ) goto 300
!            l-->l-1 transitions
            if (l2lp.eq.-1 .and. ((kinit.lt.0 .and. ind.ne.3) .or. (kinit.gt.0 .and. ind.ge.3)) ) goto 300

            iold = 0
            ic3=0
!           start cycle  do ic3=0,1
!           return for ic3=1 calculations only for |ispin|=1
!           where the central atom cross terms are needed
  100       continue
            irr = -1
!           ic3p=1 to test K2Cr2O7  L3 XES 
            ic3p = ic3
            call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx, ri, v, vval, dgcn, dpcn, adgc, adpc,                        &
                     xnval, pu, qu, p, q, iz, ihole, xion, iunf, irr, ic3p, iph) !KJ iph
            lfin = lind (ind)
            ilp = lfin - 1
            if (ikap .lt. 0) ilp = lfin + 1
            call exjlnl (xkmt, lfin, jl, nl)
            call exjlnl (xkmt, ilp, jlp1, nlp1)
            call phamp(rmt,pu,qu, ck, jl,nl,jlp1,nlp1, ikap, ph0,temp)

            sign = -1.0
            if (ikap.gt.0) sign = 1.0
            factor = ck*alphfs 
            factor = sign * factor/(1+sqrt(1+factor**2))
            dum1 = 1/ sqrt(1+factor**2)
            xfnorm = 1 / temp *dum1
!           normalization factor
!           xfnorm = dum1*rmt*(jl*cos(delta) - nl*sin(delta))/ Rl(rmt)
!           dum1 is relativistic correction to normalization
!           normalize regular solution
            do i = 1,ilast
              p(i)=p(i)*xfnorm
              q(i)=q(i)*xfnorm
            enddo

!           calculate xirf including fscf - TDLDA result
            do id = 1, 2
               if((izstd.le.0).AND.(id.eq.2)) EXIT ! Josh added to fix glitch in xmcd.
               if (id.eq.1) then
                 do j = 1,ilast
                   pp(j)  = p(j)*dble(fscf(j))
                   qp(j)  = q(j)*dble(fscf(j))
                 enddo
               else ! id.eq.2
                 do j = 1,ilast
                   pp(j)  = p(j)*dimag(fscf(j))
                   qp(j)  = q(j)*dimag(fscf(j))
                 enddo
               endif
               ifl = 1
               if (izstd.gt.0) ifl = -1
               xirf1 = 0
               call radint(ifl, mult, bf, kinit, dgc0,dpc0, ikap, pp,qp, pn,qn,ri,dx, ilast,iold, xrc,xnc, xrcold,xncold, xirf1)
!              if (ifl.lt.0) xirf1 = xirf1 * xk0 * ww
               if (ifl.lt.0) xirf1 = xirf1 * xk0
               if (id.eq.1) then
                 xirf = xirf1
               else
                 if (abs(xirf) .eq. 0.d0) then
                   xirf = xirf1
                 elseif (abs(xirf1) .eq. 0.d0) then
                   xirf = xirf
                 elseif (abs(xirf1) .lt. abs(xirf)) then
                   dum = abs(xirf1) / abs(xirf)
                   xirf = xirf * sqrt(1.d0 + dum**2)
                 else
                   dum = abs(xirf) / abs(xirf1)
                   xirf = xirf1 * sqrt(1.d0 + dum**2)
                 endif
               endif
            enddo

!           note that for real potential  xirf is real or reduced matrix
!           element for dipole transition is pure imaginary.
            if (ic3.eq.0) then
!              can remember only E2 or M1 matrix elements
               if (mult.eq.0 .or. le2.eq.mult) then
                 rkk(ie,ind)=xirf 
                 phx(ind) = ph0
               endif
!              for f' want to include both E2 and M1 for xsnorm and xsec
!              but now only one of them is included (fix later)
               xsnorm(ie)=xsnorm(ie) +  ( dble(xirf)**2 + dimag(xirf)**2 )/ (2*kx+1)
               aa =  - coni*xirf**2
               xsec(ie) = xsec(ie) -  aa * bmat(0,isp,ind, 0,isp,ind)
            elseif (iold.eq.1) then
                rkk1=xirf
                phold = ph0
            elseif (iold.eq.2) then
                rkk0=xirf
            endif

!           get irregular solution and atomic cross-section xsec
!           find irregular solution
            if(dimag(em(ie)).gt.0.d0) then
              irr = 1
!             set pu, qu - initial condition for irregular solution 
              pu = (nl*cos(ph0)+jl*sin(ph0)) *rmt * dum1
              qu=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
             
!             test on bessel functions
!             if (ikap.gt.0) print*,'test1',xkmt**2*(jl*nlp1-nl*jlp1)

              call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx, ri, v,vval, dgcn, dpcn, adgc, adpc,  &
                    xnval, pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3p, iph) !KJ iph
!             set N- irregular solution , which is outside
!             N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
!             N = i*R - H*exp(i*ph0)
              temp = exp(coni*ph0)
              do i = 1, ilast
                pn(i) = coni * p(i) - temp * pn(i)
                qn(i) = coni * q(i) - temp * qn(i)
              enddo
            else
              pn(1:ilast) = 0
              qn(1:ilast) = 0
            endif

!           combine regular and irregular solution into the central atom absorption coefficient xsec (mu = dimag(xsec))
!           thus for real energy dimag(xsec)=xsnorm
!           also include TDLDA effects
            do id = 1, 2
               if((izstd.le.0).AND.(id.eq.2)) EXIT ! Josh - added to avoid glitch in XMCD.
               if (id.eq.1) then
                 do j = 1,ilast
                   pp(j)  = p(j)*dble(fscf(j))
                   qp(j)  = q(j)*dble(fscf(j))
                   pnp(j)  = pn(j)*dble(fscf(j))
                   qnp(j)  = qn(j)*dble(fscf(j))
                 enddo
               else
                 do j = 1,ilast
                   pp(j)  = p(j)*dimag(fscf(j))
                   qp(j)  = q(j)*dimag(fscf(j))
                   pnp(j)  = pn(j)*dimag(fscf(j))
                   qnp(j)  = qn(j)*dimag(fscf(j))
                 enddo
               endif

!              TDLDA theory is written for the r-form of matrix elements
!              so one might want to use ifl=-1,-2 for these calculations
!              on the other hand want ifl=1,2 for DANES calculations
!              since it is more reliable at high energies and gives better results for Cu test.
               ifl = 2
               if (izstd.gt.0) ifl = -2

               call radint(ifl,mult, bf, kinit, dgc0, dpc0, ikap, pp, qp, &
                 pnp, qnp, ri,dx, ilast,iold, xrc, xnc, xrcold, xncold, xirf1)
               if (ifl.lt.0) xirf1 = xirf1 * xk0**2 * ww**2
               if (id.eq.1) then
                 xirf = xirf1
                 isel=0
               else
                 if (abs(xirf) .eq. 0.d0) then
                   xirf = xirf1
                   isel=1
                 elseif (abs(xirf1) .eq. 0.d0) then
                   xirf = xirf
                   isel=2
                 elseif (abs(xirf1) .lt. abs(xirf)) then
                   dum = abs(xirf1) / abs(xirf)
                   xirf = xirf * sqrt(1.d0 + dum**2)
                   isel=3
                 else
                   dum = abs(xirf) / abs(xirf1)
                   xirf = xirf1 * sqrt(1.d0 + dum**2)
                   isel=4
                 endif
               endif
               !write(*,'(3i4,2e12.4,i4,6e12.4)') ie,id,ic3,xirf1,isel,xirf,ww,xk0
            enddo
            if (ic3.eq.0) xsec(ie) = xsec(ie) - xirf * bmat(0,isp,ind, 0,isp,ind)


!           ------start of density of states part------------- 
!           added by nesvi 07/12/00
!           Calculate rhoc00 and rho_projected on the same grid as xsect. Need this to calculate the smooth atomic ratio rho_0/mu_0 or rho_proj/mu_0.
!           The atomic functions are normalized to 1 inside Norman radius. This procedure can be called 'Renormalized atomic sphere method".
!           It gives very reasonable numbers for hole counts. In order to get Mulliken counts one can extend integration limits till very
!           large R, but it's currently not recommended because of the problems with the wave function's tails above Rnm.
 

            jproj =  iorb(ikap)
            if (jproj.eq.0 .and. ikap.lt.0) jproj = iorb(-ikap-1)
            kdif1 = -1
            if(kinit.gt.0) kdif1 =  1
                
            if (kdif .eq. kdif1 .and. ic3.eq.0 .and. jproj.gt.0) then
!              calculate rhoc00 (rho_0)

               temp = (2*lfin+1.0d0) / (1+factor**2) /pi *4*ck /hart
               do i = 1, ilast
                 xrc(i) = pn(i)*p(i) - coni*p(i)*p(i) + qn(i)*q(i) - coni*q(i)*q(i)
               enddo
               xirf = 1
!              integration is till Norman radius, not Rmt as in xsect
               i0 = jnrm + 1
               call csomm2 (ri, xrc, dx, xirf, rnrm, i0)
               xrhoce(ie) = - xirf * temp
            
!              calculate rho_projected:              

!              pat, qat - atomic functions that we make projection on.
               do i=1,nrptx
                 pat(i) = dgcn(i,jproj)
                 qat(i) = dpcn(i,jproj)
               enddo

!             normalize pat and qat in the Norman radius sphere: <n|n>=1, (renormalized atomic sphere method)
               do i = 1, ilast
                  xp(i) = pat(i)**2 + qat(i)**2
                  xq(i) = 0
               enddo
              !nb, xinorm is used for exponent on input to somm
               xinorm = 2*lfin + 2
               call somm2 (ri, xp, dx, xinorm, rnrm, 0, i0)
!              call somm (ri, xp, xq, dx, xinorm, 0, jnrm)
      
               xinorm = sqrt(xinorm)
               do i=1,nrptx
                  pat(i) = pat(i) / xinorm
                  qat(i) = qat(i) / xinorm
               enddo
  
              !calculate overlap integral between f and atomic function (integral Rl(r)*Psi_at(r)dr from 0 till r').
              !intr(i) is that overlap integral. Later it will be multiplied by pr(i)*Psi_at(r') and integrated till Norman radius.
               do i=1,ilast
                  var(i)=pat(i)*p(i)+qat(i)*q(i)
!                 factor of 2 -integration r< r>  -->2 r r'
               enddo

!              integration by trapezoid method
               intr(1)=var(1)*ri(1)
               do i=2,ilast
                  intr(i)=intr(i-1)+ (var(i)+var(i-1))*(ri(i)-ri(i-1))
               enddo

!              now calculate rho_projected - xrhopr
               temp = (2*lfin+1.0d0) / (1+factor**2) /pi *4*ck /hart
!              temp = abs(ikap) / (1+factor**2) /pi *4*ck /hart
               do i = 1, ilast
                 xrc(i) = pn(i)*pat(i)*intr(i)+  qn(i)*qat(i)*intr(i)
                 xrc(i) = xrc(i) - coni*(p(i)*pat(i)*intr(i) +  q(i)*qat(i)*intr(i))
               enddo

               xirf =  1
               call csomm2 (ri, xrc, dx, xirf, rnrm, i0)
               xrhopr(ie) = - xirf * temp
    
            endif
!           ----------end of density of states part---


            if (iold.gt.0) then
!             calculate cross term contribution to XMCD
!             in both cases coupling between neighbors 
!             need to remove SO interaction (ic3=1) in order to avoid unphysical peak in Gd XMCD. a.l. ankudinov
              k1 = ind - 1
              if (k1.ge.1 .and.k1.le.8) then
              if (lind(k1).eq.lind(ind) .and. lind(k1).gt.0) then
                aa = exp( coni*(ph0 - phold))
                bb = 1/aa
                cc = - ( bmat(0,isp,k1, 0,isp,ind) + bmat(0,isp,ind, 0,isp,k1) ) / 2.d0
                xsec(ie) = xsec(ie) - coni * rkk1 * rkk0 * (bb+aa) * cc
!               combine regular and irregular solution into the central atom absorption coefficient (mu=dimag(xsec))
!               thus for real energy dimag(xsec)=xsnorm
                call radint (3, mult, bf, kinit, dgc0, dpc0, ikap, p, q, pn, &
                  qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold, xirf)
                xsec(ie) = xsec(ie) + xirf * cc * bb
  
                call radint (4, mult, bf, kinit, dgc0, dpc0, ikap, p, q, pn, &
                  qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold, xirf)
                xsec(ie) = xsec(ie) + xirf * cc * aa
              endif
              endif
            endif
!           end of |ispin=1| cross term calculations

!           prepare for ic3=1 cross term calculations if needed
            if (ic3.eq.0 .and. abs(ispin).eq.1) then
              iold = 0
              ! Josh - Write equivalent in if elseif (easier to see what is happening).
              ! if orbital angular momentum (lind) is not zero
              if (lind(ind).gt.0) then
                 if(ind.eq.1) then
                    ! If this is the first transition of this angular momentum, save xrc and xnc from radint. 
                    k1 = ind + 1
                    if (lind(k1).eq.lind(ind)) iold = 1
                 else
                    ! Otherwise, use the old result since xrc and xnc should be the same.
                    k1 = ind - 1
                    if(lind(k1).eq.lind(ind)) iold = 2
                 end if
              end if
              !if (ind.lt.8 .and. lind(ind).gt.0) then
              !  k1 = ind + 1
              !  if (lind(k1).eq.lind(ind)) iold = 1
              !endif
              !if (ind.gt.1 .and. lind(ind).gt.0) then
              !  k1 = ind - 1
              !  if (lind(k1).eq.lind(ind)) iold = 2
              !endif
!             need to remove SO interaction to calculate cross term
!             big effect for Gd XMCD calculations
              if (iold.gt.0) then
!               repeat calculation for current kdif with SO turned off
                ic3 = 1
                goto 100
              endif
            endif

  300     continue
  350   continue
        if (omega.gt.0.0) then
!         prefac = (8 * pi / 3)  * alphfs * omega  -- nonrelativistic
!         relativistic is (for alpha form)
          prefac = 4 * pi * alpinv / omega * bohr**2
          xsnorm(ie) =  xsnorm(ie) * prefac * 2*abs(ck) 
          xnorm= sqrt( xsnorm(ie) )
          xsec(ie) = xsec(ie) * prefac* 2*ck

!         put complex sqrt(prefactor) into reduced matrix elements rkk
          ck = sqrt ( prefac * (2*ck))
!         guarantee that we have the right root
          if (dimag(ck) .lt. 0) ck = -ck
!         add central atom phase shift here. 
          do kdif = 1 , 8
             rkk(ie,kdif)= rkk(ie,kdif) * ck/xnorm * exp(coni*phx(kdif))
          enddo
        endif
      enddo
!     ********************* END OF BIG LOOP OVER ENERGY MESH **************


      CALL CloseFl('mpse.dat')

!     *** optionally write ratio.dat and ratiop.dat *****
      if (ipr2.ge.3) then
!       calculate mu_0/rho_0 for XMCD normalization.
        chia(1:ne) = 0
        vrcorr = 0
        vicorr = 0
        call xscorr(1, em, ne1, ne, ik0, xrhoce,xsnorm,chia,vrcorr, vicorr, cchi)
        do ie = 1, ne1
            xrhoce(ie)  = coni* dimag(xrhoce(ie)+cchi(ie))
        enddo
        call xscorr(1, em, ne1, ne, ik0, xrhopr,xsnorm,chia,vrcorr, vicorr, cchi)
        do ie = 1, ne1
            xrhopr(ie)  = coni* dimag(xrhopr(ie)+cchi(ie))
        enddo
        call xscorr(1, em, ne1, ne, ik0, xsec,xsnorm,chia,vrcorr, vicorr, cchi)
        do ie = 1, ne1
            cchi(ie)  = coni* dimag(xsec(ie)+cchi(ie))
        enddo
        open(unit=3,file='ratio.dat',status='unknown', iostat=ios)
        open(unit=4,file='ratiop.dat',status='unknown', iostat=ios)
!       normalize to xsec at 50 ev above edge
        edg50 = emu +50.0 / hart
        call terp (omega1, xsnorm, ne1, 1, edg50, xsedge)
        write(3,440) xsedge, emu * hart 
  440   format ('# Normalization factor:', e12.4,' Angstrom**2. Fermi level at ', f7.1, ' eV.')
        write(3,450)
  450   format ('#   Energy      rho_0        mu_0       rho_0/mu_0 ')
        write(4,440) xsedge, emu * hart 
        write(4,455)
  455   format ('#   Energy      rho_proj      mu_0      rho_proj/mu_0','    mu_deloc ')

        do ie=1,ne1
           if (dimag(cchi(ie)).eq.0.d0 .and. ie.lt.ik0) then
              cchi(ie)=cchi(ik0)
              xrhoce(ie)=xrhoce(ik0)
              xrhopr(ie)=xrhopr(ik0)
           endif
           ratio = dimag(xrhoce(ie)) / dimag(cchi(ie)) * xsedge
           ratiop = dimag(xrhopr(ie)) / dimag(cchi(ie)) * xsedge

           write(3,460)  dble(em(ie))*hart, dimag(xrhoce(ie)), dimag(cchi(ie))/xsedge, ratio*corr
!          corr is the ratio N_av/N_j, responsible for difference in
!          counts due to variation of wave function due to spin-orbit
  460      format(f12.6, 2x, e12.6,2x,e12.6,2x,e12.6,1x,e12.6)      
           write(4,465)  dble(em(ie))*hart, dimag(xrhopr(ie)), dimag(cchi(ie))/xsedge, ratiop, dimag(xrhoce(ie)-xrhopr(ie))/ratio
!          also write contribution to mu_0 from delocalized states defined as
!          (rho-rho_proj)/ratio
  465      format(f12.6, 2x, e12.6,2x,e12.6,2x,e12.6,1x,e12.6,2x,e12.6)    
      
        enddo
        close(unit=3)
        close(unit=4)
      endif 

      deallocate(bmat)

      return
      end
