!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phase.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine phase_h (iph, dx, x0, ri, ne, ne1, ne3, em,              &
                       ixc, nsp, lreal, rmt, xmu,                      &
                       vi0, iPl, NPoles, Eps0, EGap, gamach,           &
                       vtot, vvalgs, edens, dmag, edenvl,              &
                       dgcn, dpcn, adgc, adpc, eref, ph, lmax,         &
                       iz, ihole, xion, iunf, xnval, ispin, is_p, Vnlm, &
                       aph)

      USE IOMOD
      use DimsMod, only: nex, nrptx, ltot, MxPole, lx, nphx=>nphu
      use constants
      USE SelfEnergyMod
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

      real*8, intent(in) ::   Vnlm(0:lx,(lx+1)**2)
      complex*16, intent(out) :: aph(nex, lx+1, (lx+1)**2)
      complex*16 em(nex)
      dimension  ri(nrptx), vtot(nrptx), edens(nrptx)
      dimension  dmag(nrptx), vvalgs(nrptx), edenvl(nrptx)
      dimension  adgc(10,41), adpc(10,41), xnval(41)
      dimension  dgcn(nrptx,41), dpcn(nrptx,41)
      complex*16  eref(nex)
      complex*16  ph(nex,-ltot:ltot)
      integer ispin

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
!     Hubbard Definitions
      real*8 fermi_shift
!      dimension  Vnlm(0:lx,(lx+1)**2,2)
      complex*16 vtotc_tmp(nrptx), vvalc_tmp(nrptx)
      complex*16 ph_m(nex, -ltot:ltot)
!      g77 (and other compilers) have an intrinsic function besjn,
!      so besjn should be declared  external 
      external besjn
      character*30 fname, f1name
!  Hubbard Input

      
!       write(fname,"('hubbard', i2.2, '.dat')") iph
!       open (unit=24, file=fname, status='old', iostat=ios)
!      ! write(*,*) 'For potential', iph, 'iostat= ', ios
!       call chopen (ios,'hubbard.dat','phase')
!       read(24,*)
!       do iss=1,2
!       do lll=0,lx
!       do mmm=(lll**2+1),(lll+1)**2
!         read(24,*) aa,aa,aa,aa,Vnlm(lll,mmm,iss)
!       enddo
!       enddo
!       enddo
!       close(24)

! End of Hubbard Input

!     zero phase shifts (some may not be set below)
      xkmax = 0
      ne12 = ne - ne3
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

      ColumnLabels(:) = ' '
      ColumnLabels(1) = 'Rs'
      ColumnLabels(1) = 'wp (eV)'
      CALL WriteData('mpse.dat',Double1 = (3 / (4*pi*edens(jri+1))) ** third, &
           &    Double2 = SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart,     &
           &    ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE., &
           &    Headers = (/ 'This file contains information about the' //&
           &    ' self-energy.' /))

!     calculate phase shifts
      do 220 ie = 1, ne12

!        Josh - xcpot now has new arguments:
!             - iPl, WpCorr, Gamma, AmpFac         
         call xcpot (iph, ie, index, lreal, ifirst, jri,                &  !here
                    em(ie), xmu,                                       &
                    vtot, vvalgs, edens, dmag, edenvl,                 &
                    eref(ie), v, vval, iPl, WpCorr, Gamma, AmpFac, EGap, &
                    vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)

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

         if (mod(ixc,10) .lt. 5) then
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
            if (abs(ll).eq.3 .and. iph.eq.1 .and. ispin.eq.1) then
               do 190 i = 1, jri
                  v(i) = v(i) - clz
                  vval(i) = vval(i) - clz
 190           continue
            endif
            call phamp (rmt, pu, qu, ck, jl(il), nl(il), jl(ilp), nl(ilp), ikap, ph(ie,ll), temp)
!  Hubbard Modifications
           if(ll.ge.0.and.ll.le.lx) then
             do 195 imm=(ll**2+1),(ll+1)**2
             if(is_p.eq.1) then 
               do 191 ir00=1,nrptx
!                vtotc_tmp(ir00)=v(ir00)-1.0*abs(Vnlm(ll,imm,is_p))
!                vvalc_tmp(ir00)=vval(ir00)-1.0*abs(Vnlm(ll,imm,is_p))
                vtotc_tmp(ir00)=v(ir00)-1.0*abs(Vnlm(ll,imm))
                vvalc_tmp(ir00)=vval(ir00)-1.0*abs(Vnlm(ll,imm))
 191           continue
             elseif(is_p.eq.2) then  
               do 192 ir00=1,nrptx
                vtotc_tmp(ir00)=v(ir00)+1.0*abs(Vnlm(ll,imm))
                vvalc_tmp(ir00)=vval(ir00)+1.0*abs(Vnlm(ll,imm))
!                vtotc_tmp(ir00)=v(ir00)+1.0*abs(Vnlm(ll,imm,is_p))
!                vvalc_tmp(ir00)=vval(ir00)+1.0*abs(Vnlm(ll,imm,is_p))
 192           continue
              endif

                call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,   &
     &               ri, vtotc_tmp,vvalc_tmp, dgcn, dpcn, adgc, adpc, & 
     &               xnval, pu, qu, p, q,                             &
     &               iz, ihole, xion, iunf, irr, ic3, iph)
                call phamp (rmt, pu, qu, ck, jl(il), nl(il),          &
     &                  jl(ilp), nl(ilp), ikap, ph_m(ie,ll), temp)

                aph(ie,ll+1,imm)=ph_m(ie,ll) 
 195         continue
            endif
! End of Hubbard Modification

!           cut phaseshift calculation if they become too small
            if (abs(ph(ie,ll)) .lt. 1.0e-6 .and. ll.ge.4)  goto 220
!           new cut function introduced by Rivas
            if(abs(exp((0,2)*ph(ie,ll))-1.).lt.1.0e-5) ph(ie,ll)=0
            if (abs(ph(ie,ll)) .lt. 1.0e-5 .and. ll.ge.4)  goto 220

  210    continue
! End of Hubbard Output
  220 continue
! Hubbard Output to aphase_up.dat and aphase_down.dat
!       if(is_p.eq.1) then 
! 2011  format('aphase_up', i2.2, '.dat')
!       write(f1name,2011) iph
!       open (unit=27, file=f1name, status='unknown', iostat=ios)   
!       call chopen (ios,'aphase.dat','phase')
!       do 225 ie=1,ne
!       do 225 lll=0,lx
!       do 225 mmm=(lll**2+1),(lll+1)**2
!         write(27,*) dble(em(ie))*hart,lll,mmm,is_p, aph(ie,lll+1,mmm,is_p)
! 225   continue
!       close(27)       
!       elseif(is_p.eq.2) then 
! 2012  format('aphase_down', i2.2, '.dat')
!       write(f1name,2012) iph
!       open (unit=27, file=f1name, status='unknown', iostat=ios)   
!       call chopen (ios,'aphase.dat','phase')
!       do 226 ie=1,ne
!       do 226 lll=0,lx
!       do 226 mmm=(lll**2+1),(lll+1)**2
!         write(27,*) dble(em(ie))*hart,lll,mmm,is_p, aph(ie,lll+1,mmm,is_p)
! 226   continue
!       close(27)       
!       endif

      CALL CloseFl('mpse.dat')

      do 230 ie = ne12+1, ne
  230 eref(ie) = eref(ne1)

      return
      end
