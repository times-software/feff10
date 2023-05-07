!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phase.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Josh - argument iPl has been added to arguments of xsect
!     Josh - added NPoles and Eps0
      subroutine phase (iph, dx, x0, ri, ne, ne1, ne3, em,              &
     &                  ixc, nsp, lreal, rmt, xmu,                      &
     &                  vi0, iPl, NPoles, Eps0, EGap, gamach,           &
     &                  vtot, vvalgs, edens, dmag, edenvl,              &
     &                  dgcn, dpcn, adgc, adpc, eref, ph, lmax,         &
     &                  iz, ihole, xion, iunf, xnval, ispin)

      USE IOMOD
      use DimsMod, only:  nex, nrptx, ltot, MxPole
      use constants
      USE SelfEnergyMod
      USE xsph_inp, only : PrintRl

      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

!     INPUT
!     iph          unique pot index (used for messages only)
!     dx, x0, ri(nr)
!                  Loucks r-grid, ri=exp((i-1)*dx-x0)
!     ne, em(ne)   number of energy points, real energy grid
!     ixc          0  Hedin-Lunqist + const real & imag part
!                  1  Dirac-Hara + const real & imag part
!                  2  ground state + const real & imag part
!                  3  Dirac-Hara + HL imag part + const real & imag part
!                  4, 5, 6, see rdinp or xcpot
!     lreal        1 for real self energy and 2 for real phase shifts 
!     rmt          r muffin tin
!     xmu          fermi level
!     vi0          const imag part to add to complex potential
!     gamach       core hole lifetime
!     vtot(nr)     total potential, including gsxc
!     vvalgs(nr)   overlap Coulomb+gsxc potential for valence electrons
!     edens(nr)    density
!     dmag(nr)     density magnetization
!     edenvl(nr)  valence charge density
!     dgcn(dpcn)   large (small) dirac components for 'iph' atom
!     adgc(adpc)   their development coefficients
!
!     OUTPUT
!     eref(ne)     complex energy reference including energy dep xc
!     ph(nex,ltot+1) complex scattering phase shifts
!     lmax         max l (lmax = kmax*rmt)

      complex*16 em(nex)
      dimension  ri(nrptx), vtot(nrptx), edens(nrptx)
      dimension  dmag(nrptx), vvalgs(nrptx), edenvl(nrptx)
      dimension  adgc(10,41), adpc(10,41), xnval(41)
      dimension  dgcn(nrptx,41), dpcn(nrptx,41)
      complex*16  eref(nex)
      complex*16  ph(nex,-ltot:ltot)
      integer ispin
      character(len=50) :: fname
!     work space for xcpot
      dimension   vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      dimension   vvxcrm(nrptx), vvxcim(nrptx)
!     p and q were needed in xsect to calc. matrix elements.
      complex*16 p(nrptx), q(nrptx)

      complex*16  p2, p2EC, ck, ckEC, xkmtp, xkmt, xkmtEC, temp, pu, qu, quAmp
      complex*16 jl(ltot+2), nl(ltot+2), jlp(ltot+2), nlp(ltot+2)
      complex*16 jlEC(ltot+2), nlEC(ltot+2), jlpEC(ltot+2), nlpEC(ltot+2)

      complex*16 v(nrptx), vval(nrptx)
      character*512 slog
!     Josh - Added iPl switch for PLASMON card
!          - and WpCorr = Wi/Wp, Gamma, AmpFac
!          - to describe Im[eps^-1]
      integer iPl, ipole, NPoles
      double precision, allocatable :: Energy(:), Loss(:)
      double precision WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole), Eps0, EGap
      character(LEN=10) ColumnLabels(20)
!     Josh END

!{#mn: g77 (and other compilers) have an intrinsic function besjn, 
!      so besjn should be declared  external 
      external besjn
!#mn}

!     zero phase shifts (some may not be set below)
      xkmax = 0
      ne12 = ne - ne3  !KJ this appears to give rubbish numbers sometimes but for now I'm not touching it.
      !write(*,*) 'In phase ne,ne3,ne12 ',ne,ne3,ne12
	  ph(:,:)=dcmplx(0)
      do ie = 1, ne
         if (ie.le.ne12 .and. xkmax.lt.dble(em(ie))) xkmax= dble(em(ie))
      enddo
      xkmax = sqrt(xkmax * 2)

!     Use kmax to find accurate l-points
!     limit l, lmax = prefac* kmax * rmt
!     prefac is set not to have warning message for Cu metal for kmax=20
      prefac = 0.7d0
      lmax = prefac * rmt * xkmax
      lmax = max(lmax, 5)
      if (lmax.gt.ltot) then
        ik = nint( ltot / rmt / bohr / prefac )
        write (slog, 110) ik
  110   format('      Phase shift calculation is accurate to k=', i2)
        call wlog(slog)
        write (slog, 120)
  120   format('      See FEFF document to increase the range.')
        call wlog(slog)
      endif
      lmax = min (lmax, ltot)
!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

      ! Josh - Print out rl.dat for RIXS
      IF(iph.EQ.0.AND.PrintRl) THEN
         CALL WriteData('rl.dat', Double1 = rmt, Int2 = lmax, Int3 = jri)
         CALL WriteData('rl.dat', Double1 = dx, Double2 = x0)
      END IF
      ifirst = 0
      index = ixc
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

!     Write wp as calculated from density to sigma.dat
!         open(file='mpse.dat', unit=45, status='replace', iostat=ios)
!         call chopen(ios, 'sigma.dat', 'ffmod2(phase)')
!         write(45,*) '# ', 'rs      wp(eV)'
!         write(45,*) '# ', (3 / (4*pi*edens(jri+1))) ** third,          &
!     &        SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart
!         write(45,*) '# mu (eV)'
!         write(45,*) '# ', xmu
!         write(45,'(a)')                                                &
!     &         '# E-EFermi (eV)   Re[Sigma(E)] (eV)   Im[Sigma(E)] (eV)'&
!     &       // '   Re[Z]   Im[Z]   Mag[Z]   Phase[Z]   Lambda(E) (/A)'
      ColumnLabels(:) = ' '
      ColumnLabels(1) = 'Rs'
      ColumnLabels(1) = 'wp (eV)'
!KJ         write(*,*) 'writing mpse.dat in PHASE'      
      fname='mpse'   
      CALL WriteData(TRIM(fname)//'.dat',Double1 = (3 / (4*pi*edens(jri+1))) ** third, &
           &    Double2 = SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart,     &
           &    ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE., &
           &    Headers = (/ 'This file contains information about the' //&
           &    ' self-energy.' /))
!     Josh END
!KJ     write(*,*) 'iPl in PHASE',iPl 
!     calculate phase shifts
      do 220 ie = 1, ne12

!        Josh - xcpot now has new arguments:
!             - iPl, WpCorr, Gamma, AmpFac         
         call xcpot (iph, ie, index, lreal, ifirst, jri,                &  !here
     &               em(ie), xmu,                                       &
     &               vtot, vvalgs, edens, dmag, edenvl,                 &
     &               eref(ie), v, vval, iPl, WpCorr, Gamma, AmpFac, EGap, &
     &               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, fname)

         if (dble(em(ie)).lt.-10.d0 .or. dble(em(ie)) .gt.3.d2) goto 220
!        p2 is (complex momentum)**2 referenced to energy dep xc
!        notice that constant Im part (gamach/2+vi0) is cancelled,
!        since it is also present in v and vval.
         p2 = em(ie) - eref(ie) 
         p2EC = p2 - vtot(1)
         if (lreal.gt.1 .and. ie.le.ne1) p2 = dble(p2)
         ck =  sqrt (2*p2+ (p2*alphfs)**2)
         ckEC = sqrt (2*p2EC+ (p2EC*alphfs)**2)
         xkmt = rmt * ck
         xkmtEC = rmt*ckEC
         if (dble(p2).le.0.d0 .and. dimag(p2) .le.0.d0) goto 220

         call besjn (xkmt, jl, nl)
         call besjn (xkmtEC, jlEC, nlEC)

      !    if (mod(ixc,10) .lt. 5) then
         IF ((ixc.EQ.0).OR.(ixc.EQ.1).OR.(ixc.EQ.2).OR.(ixc.EQ.3).OR.(ixc.EQ.6).OR.(ixc.EQ.7)) THEN ! Replace mod(ixc,10) .lt. 5
             ncycle = 0
         else
             ncycle = 3
         endif

         do 210  ll = -lmax, lmax
            il = abs(ll) +  1
!           nonlocal exchange is unstable for high il.
!           need to do integrals instead of diff. eq. fix later
!           use local xc for high il
            if (il*dx.gt.0.50) then
               ncycle=0
            endif

!  v should be V_N+V_COUL+V_XCtotal-V_mt, vval= V_N+V_COUL+V_XCVAL-V_mt
            ikap = ll - 1
            if ( ll.gt.0 ) ikap=ll
            ilp = il + 1
            if (ikap.gt.0) ilp = il - 1
            ic3 = 0

            if(nsp.eq.1 .and. ispin.eq.0) then
!              remove spin-orbit interaction
!              otherwise, get wrong results e.g. for Pt metal
               if (ll.ne.0) ic3 = 1
               ikap = -1 - abs(ll)
               ilp = il + 1
            endif

!_lz  add term (C L_z) (p.32 of Ankoudinov's thesis) 
!     currently just add constant potential only within mt radius
!     keep intersitial level the same
! OPC for U for jj coupling
!           if (ll.eq.3 .and. iph.eq.1) then
!              clz = -0.5d0 / hart
!              if (ispin.lt.0) clz = -clz
!              do 180 i = 1, jri
!                 v(i) = v(i) + clz
!                 vval(i) = vval(i) + clz
!180           continue
!           endif
! OPC for U for LS coupling
            if (abs(ll).eq.3 .and. iph.eq.1 .and. ispin.eq.1) then
               clz = -0.0d0 / hart
               if (ikap.lt.0) clz = -clz
               do 180 i = 1, jri
                  v(i) = v(i) + clz
                  vval(i) = vval(i) + clz
 180           continue
            endif

!           never use irr=0, only positive or negative
            irr = -1
            if(iz.gt.0) then
               ! This cell contains a nucleus. Use numerical solution.
               call dfovrg (ncycle, ikap, rmt, jri, jri, p2, dx,           &
                    &       ri, v,vval, dgcn, dpcn, adgc, adpc,            &
                    &       xnval, pu, qu, p, q,                           &
                    &       iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
               ! JK try using perturbation in vx for non-local corehole vx
               !IF(.TRUE.) THEN
    
            else
               ! This is an empty cell. Use bessel function solution.
               isign = 1
               IF(ikap.lt.0) isign = -1
               quAmp = ck*alphfs
               quAmp = isign*quAmp/(1.d0+SQRT(1.d0+quAmp**2))
               pu = jlEC(il)*rmt
               qu = jlEC(ilp)*rmt*quAmp
            end if
!        restore potential for clz=0
! OPC for U for jj coupling
!           if (ll.eq.3 .and. iph.eq.1) then
! OPC for U for LS coupling
            if (abs(ll).eq.3 .and. iph.eq.1 .and. ispin.eq.1) then
               do 190 i = 1, jri
                  v(i) = v(i) - clz
                  vval(i) = vval(i) - clz
 190           continue
            endif
            call phamp (rmt, pu, qu, ck, jl(il), nl(il), jl(ilp), nl(ilp), ikap, ph(ie,ll), temp)
            ! Josh - debugging output.
            if(iz.eq.0.and.ll.eq.1) CALL WriteData('phaseEC.debug', DComplex1 = em(ie), DComplex2 = ph(ie,1))

            IF(PrintRl.AND.((iph.EQ.0).AND.(ll.LE.0)).AND.(ABS(ll).LE.lmax)) THEN
               ! Josh Kas - form rkap for RIXS.
               !sign = -1.0
               !if (ikap.gt.0) sign = 1.0
               !factor = ck*alphfs 
               !factor = sign * factor/(1+sqrt(1+factor**2))
               !dum1 = 1/ sqrt(1+factor**2)
               !dum1 = 1.d0
               DO i = 1,jri
                  p(i)=p(i) / temp
                  q(i)=q(i) / temp
               END DO
               ! Josh Kas - Write out rl.dat for RIXS stuff.
               CALL WriteData('rl.dat', DComplex1 = em(ie), Int2 = -ll, DComplex3 = ph(ie,ll))
               CALL WriteArrayData('rl.dat', DComplex1 = p(1:jri), DComplex2 = q(1:jri))
            END IF
!           cut phaseshift calculation if they become too small
            if (abs(ph(ie,ll)) .lt. 1.0e-6 .and. ll.ge.4)  goto 220
!           new cut function introduced by Rivas
            if(abs(exp((0,2)*ph(ie,ll))-1.).lt.1.0e-5) ph(ie,ll)=0
            if (abs(ph(ie,ll)) .lt. 1.0e-5 .and. ll.ge.4)  goto 220

  210    continue

  220 continue
!     Josh - Close mpse.dat, rl.dat
      IF((iph.EQ.0).AND.PrintRl) CALL CloseFl('rl.dat')
      CALL CloseFl(TRIM(fname)//'.dat')
!     Josh END

      do ie = max(1,ne12+1), ne  !KJ Feb 14 - added "max(1" to respect array bounds
         eref(ie) = eref(ne1)
      enddo

      return
      end
