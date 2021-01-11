!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xsectd.f90,v $:
! $Revision: 1.11 $
! $Author: jorissen $
! $Date: 2011/12/10 23:17:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Josh - added argument iPl to control many pole self energy
      subroutine xsectd (ipr2, dx, x0, ri, ne, ne1, ik0, em, edge,      &
     &                  ihole, emu, corr, dgc0, dpc0, jnew,             &
     &                  ixc, lreal, rmt, rnrm, xmu,                     &
     &                  vi0, iPl, NPoles, Eps0, EGap, gamach,           &
     &                  vtot, vvalgs, vch, edens, dmag, edenvl,         &
     &                  dgcn, dpcn, adgc, adpc, xsec, xsnorm, rkk,      &
     &                  iz, xion, iunf, xnval,                          &
     &                  ipmbse, ifxc, ibasis, eorb, kappa, iorb, l2lp,  &
     &                  ipol, ispin, le2, angks, ptz, itdlda, iph) !KJ

      USE IOMOD
      use dimsmod, only: nrptx, lx_xsph, MxPole, nspx=>nspu
      use constants
      USE SelfEnergyMod
!     right now the same self-energy is used for calculation
!     of the central atom part (xsec) and dipole m.e. for
!     scattering (rkk). You may want to run xsect separately
!     for xsec and for rkk, if you want to use different self-energy
!     for central and scattering parts.  ala. fix later

      !implicit double precision (a-h, o-z)
      implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

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
!     iPl - Josh: added to control many pole self energy
!
!     OUTPUT
!     xsec(ne)    atomic absorption cross section to multiply \chi
!                 (atomic background for XMCD)
!     xsnorm(ne)  atomic  absorption cross section (norm for XMCD)
!     rkk(ne, 8)  normalized reduced matrix elements for construction
!                 of termination matrix in genfmt.
      complex*16 ptz(-1:1, -1:1)
      complex*16,parameter :: conr = (1,0)
      integer, parameter :: nex = 500 ! JJK - set a separate nex here
                                      ! Should allocate in future.
                             ! since nex = 2000 is way too large and causes problems
                             ! with memory in this routine.

      real*8 ri(nrptx), vtot(nrptx), edens(nrptx),dmag(nrptx)
      real*8 dgc0(nrptx), dpc0(nrptx), vvalgs(nrptx), edenvl(nrptx)
      real*8 dgcn(nrptx,41), dpcn(nrptx,41)
      real*8 dgcnp(nrptx,41), dpcnp(nrptx,41)
      real*8 adgc(10,41), adpc(10,41), xnval(41)
      integer iorb(-5:4)
      real*8 eorb(41)
      integer kappa(41)
      complex*16 rkk(nex, 8), xsec(nex)
      complex*16 bmat(-lx_xsph:lx_xsph,0:1,8, -lx_xsph:lx_xsph,0:1,8)
      dimension kiind(8), lind(8) !KJ renamed kind to kiind to avoid conflicts with reserved name
     
      real*8 xp(nrptx), xq(nrptx), vch(nrptx)

!     work space for xcpot
      real*8 vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      real*8 vvxcrm(nrptx), vvxcim(nrptx)

!     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx)
!     storage for calculation of cross term (SPIN 1 only)
      complex*16 xrcold(nrptx) , xncold(nrptx)

      complex*16  p2, ck
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf
      complex*16  aa, bb, cc, rkk1, rkk0, phold
      complex*16  phx(8), ph0
      complex*16  eref, xm1, xm2, xm3, xm4

      complex*16  v(nrptx), vval(nrptx)
      complex*16  xrc(nrptx), xnc(nrptx)
      character*512 slog
      logical ltrace
!     nesvi:  
      complex*16 xrhoce(nex), xrhopr(nex), chia(nex), cchi(nex)
      real*8 omega(nex), bf(0:2, nrptx)

      complex*16 xrhoce1(nex), xsec1(nex)
      real*8 emr(nex), sfun(nex)
      real*8 chil2(nex), chil3(nex), chil4(nex), chil5(nex)
      complex*16 em(nex), emx
      real*8 pat(nrptx),qat(nrptx)

      integer, parameter :: maxsize = 78
      integer jinit(maxsize), minit(maxsize), jfin(maxsize), mfin(maxsize)
      integer kinitm(maxsize), kfinm(maxsize), nph(maxsize)
      real*8 chi0r(nex,maxsize, maxsize), chi0im(nex,maxsize, maxsize)
      integer ncore(maxsize)
      real*8 dipmatl(nex,maxsize), dipmat(nex,maxsize)
      complex*16 chi0br(nex), chi0(nex,maxsize, maxsize)
      real*8 refsh(maxsize)
      real*8 dml(maxsize), dmdl(maxsize), chi(maxsize, maxsize)
      real*8 phf(maxsize)
      complex*16 wmat(nex, maxsize,maxsize), wm(maxsize,maxsize)
      complex*16 xkmat(nex, maxsize,maxsize), xkm(maxsize,maxsize)
      complex*16 xkmatp( maxsize,maxsize, nex)
      complex*16 dipscf(nex,maxsize), ctemp
      
      real*8 xsnorml3(nex), xsnorml2(nex), xsnorml5(nex), xsnorml4(nex)
      real*8 xsscfl3(nex), xsscfl2(nex), xsscfl5(nex), xsscfl4(nex)
      real*8 gammab(maxsize)
      real*8 ckl3(nex), ckl2(nex)
      complex*16 dmlbr(nex), dmdlbr(nex)
      complex*16 xl3br(nex), xl2br(nex), xl5br(nex), xl4br(nex)

!     Josh - Added iPl switch for PLASMON card
!          - and WpCorr = Wi/Wp, Gamma, AmpFac
!          - to describe Im[eps^-1]
      integer iPl, ipole, NPoles
      real*8, allocatable :: Energy(:), Loss(:)
      real*8 WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole), Eps0, EGap
      character(LEN=10) ColumnLabels(20)
!     Josh END      

!     for implicit none:
      integer ihole, kinit, linit, iz, iholep, nlp, nlm, nlpoc, nlmoc, ibasis
      real*8 gamach, gaml2, gaml3, rmt, x0, dx, rint, rnrm
      integer nch
      integer imt, matsize, jri, jri1, inrm, jnrm, iint, jint, ie, ne1, im
      real*8 deltaso,temp
      real*8 xmu
      real*8 xion,xinorm,del
      real*8 angks,e1,e2,e12,de,xx,edge,emu,xsnorm
      integer ixc,ndata,ixcp,jj,iph,lreal,ifirst,iunf,i,ipol,le2,ispin,kiind
      integer lind,isp,ipmbse,nelast,ifxc,imp,imi,imj,ne,imx,ios
      real*8 prefac
      real*8 prefacl2,prefacl3,vi0,corr
      integer ipr2,ik0,jnew,l2lp,itdlda
      !implicit none

      call setkap(ihole, kinit, linit)
      if (kinit.ge.0) then
        call wlog('  Initial state kappa should be negative: ')
        call wlog('  E.g.: HOLE K,L3,M5 and  not HOLE L2, M4, etc.')
        stop
      endif
!  begin   manual input
!KJ Feb 2014:  gaml2,gaml3,deltaso are now set automatically below.  Manual input not needed anymore.
!     set splitting between 2 edges and their lifetimes
! V
!     deltaso = 8.620 /hart
!     gaml2 = 0.458/hart/2.0
!     gaml3 = 0.275/hart/2.0
      call setgam (iz, ihole, gamach)
      gaml3 = gamach / hart / 2.d0
      if (kinit.lt.-1) then
        iholep = ihole - 1
        call setgam (iz, iholep, gamach)
        gaml2 = gamach / hart / 2.d0
      else
        gaml2 = gaml3
      endif
       
!  Co 
!     deltaso = 15.24 /hart
!     gaml2 = 1.120/hart/2.0
!     gaml3 = 0.480/hart/2.0
!  Ni
!     deltaso = 17.5 /hart
!     gaml2 = 1.40/hart/2.0
!     gaml3 = 0.550/hart/2.0
! Diamond
!     deltaso = 0.00/hart
!     gaml2 = 0.087/hart/2.0
!     gaml3 = 0.087/hart/2.0
! Mg
!     gaml2 = 0.34/hart/2.0
!     gaml3 = 0.34/hart/2.0
!  Xe
!     deltaso = 1.96/hart
!     gaml2 = 0.109/hart/2.0
!     gaml3 = 0.109/hart/2.0
!  W  
!     deltaso = 62.2/hart
!     gaml2 = 3.430/hart/2.0
!     gaml3 = 1.985/hart/2.0
! Ta
!     deltaso = 58.1/hart
!     gaml2 = 3.175/hart/2.0
!     gaml3 = 1.885/hart/2.0

! set basis set size for L+1 and L-1 orbitals
      nlp = 3
      nlm = 0
! set number of completely occupied L+1 and L-1 orbitals
!     want to know how many nodes basis w.f. should have
      nlpoc = 1
      nlmoc = 0
!KJ ibasis is now read from input (PMBSE card):
! choose the basis
!       choose ibasis = 0 to use occupied orbitals for the basis
!       ibasis = 1 - read basis orbitals from file
!       ibasis = 2 - calculate them by requireing 0 at r_int
!     ibasis = 2
! end of  manual input for xsectd
      if (ibasis.eq.0) then
!       cannot have more than one orbital for projections
        if (nlp .gt. 0) nlp = 1
        if (nlm .gt. 0) nlm = 1
      endif

!     set number of channels (Number of edges involved times number
!     of final l-channels. E.g. for K-edge nch=1, for L2,3 edges with
!     final d only nch=2, with final d,s in basis set nch=4)
      nch = 2
      if (kinit .eq.-1) then
        nch = 1
        nlm = 0
      endif
      if (nlm.gt.0) nch = 4


!     set matrix indices and size
      call getmat(ihole, linit, nlp, nlm, jinit, minit, kinitm,         &
     &     jfin, mfin, kfinm, ncore, nph, matsize, kappa, xnval, ibasis)
!      nph(im) - index of orbital to project on within basis set
!          if nph>0 index refer to array dgcn(dpcn), i.e.
!          this is one of the partly filled atomic orbitals
!          if nph<0 index refer to array dgcnp(dpcnp) which
!          are calculated by getwf subroutines and represent
!          completely unoccupied orbitals

!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (dlog(rmt) + x0) / dx  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
!     nesvi: define jnrm
!     test - increase rnrm
!      rnrm = rnrm*5
      inrm = (log(rnrm) + x0) / dx + 1
      jnrm = inrm + 1

!     set the cutoff radius for integrations
       rint = rnrm
!      rint = 2.2d0 / bohr
       iint = (log(rint) + x0) / dx + 1
       jint = iint + 1

!     read in the fine structure
      do ie = 1, ne1
        emr(ie) = dble(em(ie))
      enddo
      call ridxmu(kinit, ne1, emr, chil2, chil3, chil4, chil5, deltaso)
      do 5 im = 1, matsize
         if (kinitm(im) .gt. 0) then  
             refsh(im) = - deltaso   
             gammab(im)= gaml2
         else
             refsh(im) = 0.0
             gammab(im)= gaml3
         endif
    5 continue


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
      ColumnLabels(:) = ' '
      ColumnLabels(1) = 'Rs'
      ColumnLabels(1) = 'wp (eV)'
      !write(*,*) 'writing mpse.dat in XSECTD'
      CALL WriteData('mpse.dat',Double1 = (3 / (4*pi*edens(jri+1))) ** third, &
           &    Double2 = SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart,     &
           &    ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE., &
           &    Headers = (/ 'This file contains information about the' //&
           &    ' self-energy.' /))
      
      !     Josh END
        
      ixcp  = ixc
!      set atomic orbitals basis set
!      0)  by default the calculated partially occupied orbitals are used
!      1) basis set constructed by having fixed number of nodes
!         ar R_int. Note the user should  specify the energies and
!         check that node is there by plotting fort.78
!      2) natural atomic orbitals can be read from file; see inside
!          subroutine getwf for details

       do 10 jj = 0, nlm+nlp-1
!       manual input for case 1
!c     Xe, nlp=4, nlm=2, nlpoc=0, nlmoc=4
!       temp = 14.8/hart
!       if (jj.eq.1) temp = 31.5/hart
!       if (jj.eq.2) temp = 67.0/hart
!       if (jj.eq.3) temp = 118.0/hart
!       if (jj.eq.4) temp = 25.0/hart
!       if (jj.eq.5) temp = 75.0/hart

!c     tungsten, nlp=5, nlm=0
        temp =  36.5 /hart
        if (jj.eq.1) temp = 110.0/hart
        if (jj.eq.2) temp = 220.0/hart
        if (jj.eq.3) temp = 370.0/hart
        if (jj.eq.4) temp = 550.0/hart

!c     Mg, nlp = 2, nlm=0, nlpoc = 1, nlmoc = 0
!       temp = -5.0/hart
!       if (jj.eq.2) temp = 32.0/hart
!       if (jj.eq.3) temp = 92.0/hart
!       if (jj.eq.4) temp = 160.0/hart
!       temp = 17.0/hart
!       if (jj.eq.1) temp = 103.0/hart

        iph = 0
        ie = 1
        
!       Josh - added arguments: iPl, WpCorr, Gamma, AmpFac
!              for calculation of many pole self energy.        
        ctemp = DCMPLX(temp)
        call xcpot (iph, ie, ixcp , lreal, ifirst, jri, ctemp, xmu,     &
     &               vtot, vvalgs, edens, dmag, edenvl,                 &
     &               eref, v, vval, iPl, WpCorr, Gamma, AmpFac,0.d0,         & !KJ 12-2011 argument "EGap" missing: I added "0.0"
     &               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)
        call getwf (ibasis, jj, nlp, nlm, nlpoc, nlmoc, rmt, rint, jri, &
     &           jint, temp, eref, dx, x0, ri, v, vval, pat, qat,       &
     &           dgcn, dpcn, adgc, adpc, dgcnp, dpcnp, xnval,           &
     &           iz, ihole, xion, iunf, kinitm, kfinm, nph, matsize, iph) !KJ iph
!c      write out orbitals in fort.78 for visual check of nodes
        do i = 1, jint
           write(78,777) ri(i)*bohr, pat(i), qat(i)
        enddo
  10   continue
!     itest=2
!     if (itest.eq.2) stop

!     We'll need <i|i> later to normalize dipole matrix elements
!     <i|r|f>.  NB, dgc and dpc are r*wave_fn, so use '0' in somm to
!     get integral  psi**2 r**2 dr.
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
!       Always use ground state self energy for xsection, quick fix
!       JJR, Jan 93
!       change for testing broadened plasmon pole 6/93
!       ixcp  = 2
!   ALA found that it is better to use ixcp =ixc and real part of 
!   self-energy for atomic xsection. 12/96
      ltrace = .true.
      call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, kiind, lind, bmat)
!     set spin index to use bmat
      isp = 0
      if (ispin.eq.1) isp = nspx - 1

!     zero rkk and phx
      rkk(:,:) = 0
      phx(:) = 0

      ifirst = 0

!     define s-function that goes from 0 at E1 to 1 at E2 smoothly
      e1 = 100.d0 / hart
      e2 = 150.d0 / hart
      e12 = (e1+e2) / 2
      de = abs (e2 - e12)
      if (de.lt.0.05/hart) de = 0.05/hart
      do 35 ie = 1, ne1
        xx = (emr(ie) - e12)/de
        if (xx.lt.-1) then
          sfun(ie) = 0
        elseif (xx.ge.1.d0) then
          sfun(ie) = 1.d0
        else
          sfun(ie) = 0.25d0 * ( 2+3*xx-xx**3)
        endif
!       reset sfun to  0 for PMBSE only; 
!              or 1 for combined (LF via TDDFT, corehole via PMBSE)
        if (ipmbse.eq.1 .or. ipmbse.eq.3) then
          sfun(ie) = 1.d0
        elseif (ipmbse.eq.2) then
          sfun(ie) = 0.d0
        endif
   35 continue

!     if ( itdlda.eq.1) then
!       do 343 ie = 1, ne1
!         chil2(ie) = 1.d0
!         chil3(ie) = 1.d0
!         chil4(ie) = 1.d0
!         chil5(ie) = 1.d0
!343    continue
!     else
!       do 344 ie = 1, ne1
!         read(15,*) dum, chil3(ie), chil2(ie)
!         if (nch.gt.2) read(16,*) dum, chil5(ie), chil4(ie)
!344    continue
!     endif

!     calculate Im chi0, and matrix elements
      call wlog('energy loop for chi0')
      nelast = ne1
      do 400 ie =1, nelast
        iph = 0
        emx = emr(ie)
        !write(*,*) 'energy ie= ',ie
        
!       Josh - added arguments: iPl, WpCorr, Gamma, AmpFac
!              for calculation of many pole self energy.
        call xcpot (iph, ie, ixcp , lreal, ifirst, jri, emx, xmu, vtot, vvalgs, edens, dmag, edenvl, &
                  eref, v, vval, iPl, WpCorr, Gamma, AmpFac,0.d0,vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)
                  !KJ 12-2011 argument "EGap" missing: I added "0.0"


!c       set the method to calculate atomic cross section
!c       p2 is (complex momentum)**2 referenced to energy dep xc
         p2 = emr(ie) - eref
         ckl3(ie) = sqrt (2*p2 + (p2*alphfs)**2)
         p2 = emr(ie) - (eref + deltaso)
         ckl2(ie) = sqrt (2*p2 + (p2*alphfs)**2)

        omega(ie) = (emr(ie) - edge) + emu
        omega(ie) = max (omega(ie), 0.1d0 / hart)

        xsnorml3(ie) = 0
        xsnorml2(ie) = 0   
        xsscfl3(ie) = 0
        xsscfl2(ie) = 0
        xsnorml5(ie) = 0
        xsnorml4(ie) = 0   
        xsscfl5(ie) = 0
        xsscfl4(ie) = 0
   
        if (emr(ie).lt.-10.d0) goto 400

     !   call wlog('entering getchi0')
        call getchi0 (ie, emr(ie), eref, edge, emu, refsh, omega(ie),   &
     &           ipmbse, ifxc, sfun(ie),                                &
     &           rmt, rint, jri, jint, dx, x0, ri, edens, v, vval, vch, &
     &           eorb, kappa, dgcn, dpcn, adgc, adpc,                   &
     &           dgcnp, dpcnp, xnval, iz, ihole, xion, iunf, nlp, nlm,  &
     &           jinit,minit,kinitm,jfin, mfin, kfinm, ncore, nph, matsize,&
     &           dml, dmdl, chi,wm, xkm,xkmatp(1,1,ie),phf, iph) !KJ 

        do 112 im = 1, matsize
            dipmatl(ie,im) = dml(im)
            dipmat(ie,im) = dmdl(im)
            do 113 imp = 1, matsize
              if (kinitm(im).gt.0) then
               if (abs(kfinm(im)+1).gt.abs(kinitm(im)+1)) then
!               jinit = linit - 1/2; lfin = linit + 1
                chi0im(ie,im,imp) = chi(im,imp)*chil2(ie)
               else
!               jinit = linit - 1/2; lfin = linit - 1
                chi0im(ie,im,imp) = chi(im,imp)*chil4(ie)
               endif
              else
               if (abs(kfinm(im)+1).gt.abs(kinitm(im)+1)) then
!               jinit = linit + 1/2; lfin = linit + 1
                chi0im(ie,im,imp) = chi(im,imp)*chil3(ie)
               else
!               jinit = linit + 1/2; lfin = linit - 1
                chi0im(ie,im,imp) = chi(im,imp)*chil5(ie)
               endif
              endif
  113       continue
  112   continue


        do imi = 1, matsize
        do imj = 1, matsize
           xkmat(ie,imi,imj) = xkm(imi,imj)
           wmat(ie,imi,imj) = wm(imi,imj)
        enddo
        enddo

!       test - write matrix K into file
!       if (ie .eq. 60) then
!        open(unit=3,file='xkmat.dat',status='unknown', iostat=ios)
!        do 471 i=1,matsize
!        do 471 j =1, matsize 
!              write(3,461)  i, j, xkmat(ie,i,j)
!  461         format(i6, 2x, i6,2x,e12.6, 2x, e12.6)      
!  471   continue  
!        close(unit=3)
!      endif
        

  400 continue
!     Josh - Close sigma.dat
      close(45)
!     Josh END

!     multiply by factor 2 in nonrelativistic case
!     do 410 ie =1, ne1
!     do 410 im = 1, matsize
!     do 410 imp = 1, matsize
!           chi0im(ie,im,imp) = 2* chi0im(ie,im,imp)
! 410 continue
  !    call wlog('adding broadening')
!     add broadening to chi0
!     gamma=0.2/hart/2.0
      do 430 im = 1, matsize
      do 430 imp = 1, matsize
        do 420 ie = 1, ne1
          chi0br(ie) = chi0im(ie,im,imp)*conr
  420   continue
                    
        call conv(emr,chi0br,ne1,gammab(im))

        do 425 ie =1, ne1
           chi0im(ie,im,imp) = dble(chi0br(ie))
  425   continue

  430 continue   
 
                
!     calculate real part of chi0 
!        print*,'kkchi'
      call wlog('calculating kkchi')
      call kkchi (emu, edge, refsh, kinitm, kfinm, matsize, emr, ne, ne1, chi0im, chi0r)
       
      nelast = ne1
      do 510 ie =1, ne1
        do 34 im = 1, matsize
        do 34 imp = 1, matsize
          chi0(ie,im,imp) = (chi0r(ie,im,imp)+ chi0im(ie,im,imp)*coni)
   34   continue
  510 continue
!     change chi0 to Zangwill-Soven response matrix
!     chi^zs = chi^0 * wmat * chi^0
!      use 
      do 540 ie = 1, ne1
        do im = 1, matsize
         do imp = 1, matsize
          xkm(im,imp) = 0
          do imx = 1, matsize
            xkm(im,imp) = xkm(im,imp) + wmat(ie,im,imx)*chi0(ie,imx,imp)
          enddo
         enddo
        enddo

        do im = 1, matsize
         do imp = 1, matsize
          wm(im,imp) = 0
          do imx = 1, matsize
            wm(im,imp) = wm(im,imp) + chi0(ie,im,imx)*xkm(imx,imp)
          enddo
         enddo
        enddo

        do im = 1, matsize
         do imp = 1, matsize
!           chi0(ie,im,imp) = chi0(ie,im,imp) + wm(im,imp)
         enddo
        enddo
  540 continue
               
           
!     calculate screened matrix element
      call dmscf(emr, nelast,nex, matsize, chi0, dipmatl, xkmat, dipscf)    

      im = 5
      do ie =1, nelast
        write(77, 777)  emr(ie)*hart, dble(chi0(ie, 1,1)), dimag(chi0(ie, 1,1)), dble(chi0(ie,im,im)), dimag(chi0(ie,im,im))
  777   format( 6e12.4)
      enddo
               

      do 500 ie =1, nelast
!c      test case - test KK transform for a Lorentzian
!       gamma = 2.0
!       ener = emr(ie) - 500.0/hart
!       do 34 im = 1, matsize
!       do 34 imp = 1, matsize
!         chi0im(im,imp) = - ener / (ener**2 + gamma**2)
!         chi0(ie,im,imp) = chi0r(im,imp)*conr + chi0im(im,imp)*coni
!   34  continue

!        -- nonrelativistic
!       prefac = - 8 * pi / 3 * alphfs * omega(ie) *  bohr**2
!        -- relativistic is (for alpha form)
        prefac = - 4 * pi * alpinv / omega(ie) * bohr**2 * 100
!       last factor 100 transforms to Mbarn from A**2

!       prefactor with  - 2*ck   
        prefacl3 = - 2 * dble(ckl3(ie))* prefac
        prefacl2 = - 2 * dble(ckl2(ie))* prefac

!       sum over m.m' 
        do 305 im = 1, matsize
          dipmat(ie,im) = dipmat(ie,im)
          if (kinitm(im) .lt. 0) then
!           single-electron approximation
            if (im.le.15) then
              xsnorml3(ie) =  xsnorml3(ie) + dipmat(ie,im)**2
            else
              xsnorml5(ie) =  xsnorml5(ie) + dipmat(ie,im)**2
            endif
!           including TDLDA effect
            temp = dipmat(ie,im)
            do 555 imp = 1, matsize
            do 555 imx = 1, matsize
  555       temp=temp+xkmatp(im,imp,ie)*chi0(ie,imp,imx)*dipscf(ie,imx)
            if (im.le.15) then
              xsscfl3(ie) =  xsscfl3(ie) + (abs(temp))**2
            else
              xsscfl5(ie) =  xsscfl5(ie) + (abs(temp))**2
            endif
          else
!           single electron approximation
            if (im.le.15) then
              xsnorml2(ie) =  xsnorml2(ie) + dipmat(ie,im)**2
            else
              xsnorml4(ie) =  xsnorml4(ie) + dipmat(ie,im)**2
            endif
!           including TDLDA effect
            temp = dipmat(ie,im)
            do 556 imp = 1, matsize
            do 556 imx = 1, matsize
  556       temp=temp+xkmatp(im,imp,ie)*chi0(ie,imp,imx)*dipscf(ie,imx)
            if (im.le.15) then
              xsscfl2(ie) =  xsscfl2(ie) + (abs(temp))**2
            else
              xsscfl4(ie) =  xsscfl4(ie) + (abs(temp))**2
            endif
          endif
            
  305   continue 

        xsnorml3(ie) = xsnorml3(ie) * prefacl3
        xsscfl3(ie) = xsscfl3(ie) * prefacl3
        if (nch.gt.1) then
          xsnorml2(ie) = xsnorml2(ie) * prefacl2
          xsscfl2(ie) = xsscfl2(ie) * prefacl2
          if (nch.gt.2) then
            xsnorml4(ie) = xsnorml4(ie) * prefacl2
            xsnorml5(ie) = xsnorml5(ie) * prefacl3
            xsscfl5(ie) = xsscfl5(ie) * prefacl3
            xsscfl4(ie) = xsscfl4(ie) * prefacl2
          endif
        endif

 500  continue
!c    end of energy cycle


!     add broadening 
      do 820 ie = 1, ne1
        if (emr(ie) .lt. edge ) then
          xl3br(ie) = 0.0             
          xl5br(ie) = 0.0             
        else
          xl3br(ie) = xsscfl3(ie)*conr
          if (nch.gt.2) xl5br(ie) = xsscfl5(ie)*conr
        endif

        if (nch.gt.1) then
        if (emr(ie) .lt. edge + deltaso ) then
          xl2br(ie) = 0.0             
          if (nch.gt.2) xl4br(ie) = 0.0             
        else
          xl2br(ie) = xsscfl2(ie)*conr
          if (nch.gt.2) xl4br(ie) = xsscfl4(ie)*conr
        endif        
        endif        
  820 continue
                    
      call conv(emr,xl3br,ne1,gaml3)
      if (nch.gt.1) call conv(emr,xl2br,ne1,gaml2)
      if (nch.gt.2) call conv(emr,xl5br,ne1,gaml3)
      if (nch.gt.2) call conv(emr,xl4br,ne1,gaml2)

      do 825 ie =1, ne1
        xsscfl3(ie) = dble(xl3br(ie))
        if (nch.gt.1) xsscfl2(ie) = dble(xl2br(ie))
        if (nch.gt.2) xsscfl5(ie) = dble(xl5br(ie))
        if (nch.gt.2) xsscfl4(ie) = dble(xl4br(ie))
  825 continue

!     add broadening 
      do 821 ie = 1, ne1
        if (emr(ie) .lt. edge) then
          xl3br(ie) = 0.0             
          if (nch.gt.2) xl5br(ie) = 0.0             
        else
          xl3br(ie) = xsnorml3(ie)*conr
          if (nch.gt.2) xl5br(ie) = xsnorml5(ie)*conr
        endif

        if (nch.gt.1) then
        if (emr(ie) .lt. (edge + deltaso)) then
          xl2br(ie) = 0.0             
          if (nch.gt.2) xl4br(ie) = 0.0             
        else
          xl2br(ie) = xsnorml2(ie)*conr
          if (nch.gt.2) xl4br(ie) = xsnorml4(ie)*conr
        endif        
        endif        
  821 continue
                    
      call conv(emr,xl3br,ne1,gaml3)
      if (nch.gt.1) call conv(emr,xl2br,ne1,gaml2)
      if (nch.gt.2) call conv(emr,xl5br,ne1,gaml3)
      if (nch.gt.2) call conv(emr,xl4br,ne1,gaml2)

      do 826 ie =1, ne1
        xsnorml3(ie) = dble(xl3br(ie))
        if (nch.gt.1) xsnorml2(ie) = dble(xl2br(ie))
        if (nch.gt.2) xsnorml5(ie) = dble(xl5br(ie))
        if (nch.gt.2) xsnorml4(ie) = dble(xl4br(ie))
  826 continue

!--------- make final output ---------------------------------------
      open(unit=4,file='xsedge.dat',status='unknown', iostat=ios)
      do 460 ie=1,nelast 
        xsnorml3(ie) = xsnorml3(ie)*chil3(ie)
        if (nch.gt.1) xsnorml2(ie) = xsnorml2(ie)*chil2(ie)
        if (nch.gt.2) xsnorml4(ie) = xsnorml4(ie)*chil4(ie)
        if (nch.gt.2) xsnorml5(ie) = xsnorml5(ie)*chil5(ie)
        xsscfl3(ie) = xsscfl3(ie)*chil3(ie)
        if (nch.gt.1) xsscfl2(ie) = xsscfl2(ie)*chil2(ie)
        if (nch.gt.2) xsscfl4(ie) = xsscfl4(ie)*chil4(ie)
        if (nch.gt.2) xsscfl5(ie) = xsscfl5(ie)*chil5(ie)
  460 continue
      do 470 ie=1,nelast 
        if (nch.gt.2) then
          write(4,465)  dble(emr(ie)+emu) * hart,                       &
     &    dble(xsnorml3(ie)+xsnorml2(ie)+xsnorml5(ie)+xsnorml4(ie)),    &
     &    dble(xsscfl3(ie) +xsscfl2(ie) +xsscfl5(ie) +xsscfl4(ie) ),    &
     &    dble(xsnorml3(ie)+xsnorml5(ie)),                              &
     &    dble(xsnorml2(ie)+xsnorml4(ie)),                              &
     &    dble(xsscfl3(ie)+ xsscfl5(ie)),                               &
     &    dble(xsscfl2(ie)+ xsscfl4(ie))
        elseif (nch.gt.1) then
          write(4,465)  dble(emr(ie)+emu) * hart,                       &
     &    dble(xsnorml3(ie)+xsnorml2(ie)),                              &
     &    dble(xsscfl3(ie) +xsscfl2(ie) ),                              &
     &    dble(xsnorml3(ie)),                                           &
     &    dble(xsnorml2(ie)),                                           &
     &    dble(xsscfl3(ie)),                                            &
     &    dble(xsscfl2(ie))
        else
          write(4,465)  dble(emr(ie)+emu) * hart,                       &
     &    dble(xsnorml3(ie)),                                           &
     &    dble(xsscfl3(ie)) 
        endif

  465   format(f10.5, 2x, e10.5,2x,e10.5,2x,e10.5,1x,e10.5,1x,e10.5,1x,e10.5)
  470 continue  
      close(unit=4)

      return
      end
