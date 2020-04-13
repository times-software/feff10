!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xcpot.f90,v $:
! $Revision: 1.20 $
! $Author: hebhop $
! $Date: 2012/03/29 22:52:37 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xcpot (iph, ie, index, lreal, ifirst, jri, em, xmu,                                        &
     &                 vtot, vvalgs, densty, dmag, denval,               &
     &                  eref, v, vval, iPl, WpCorr, Gamma, AmpFac, EGap, &
     &                  vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)

      USE IOMOD
      use dimsmod, only: nrptx, MxPole
	  use constants
      implicit double precision (a-h, o-z)
!     calculate self-energy correction
!     first coded j. mustre de leon
!     last modified a.ankudinov 1996 for non-local self-energies
!     Ankudinov, Rehr, J. Physique IV, vol. 7, C2-121, (1997).

!     INPUT
!     iph, ie used only for debug and labels.
!     index       0  Hedin-Lunqvist + const real & imag part
!                 1  Dirac-Hara + const real & imag part
!                 2  ground state + const real & imag part
!                 3  Dirac-Hara + HL imag part + const real & imag part
!                 4  See rdinp for comment
!     lreal       not equal zero for real self energy
!     ifirst      first entry flag, set to zero before first call for
!                 each unique potential, see vxcrmu and vxcimu below
!     jri         index of first interstitial point in current
!                 Loucks r grid
!     em          current energy grid point
!     xmu         fermi level
!     vi0         const imag part to subtract from potential
!     gamach      core hole lifetime
!     vtot(nr)    total potential (coulomb and gs exchange corr)
!     vvalgs(nr)  total coulomb + gs xc potential from valence electrons
!     densty(nr)  electron density
!     dmag(nr)    density magnetization
!     denval(nr)  valence electron density
!     iPl         Control for many pole self energy (Josh)
!
!     OUTPUT
!     eref        complex energy reference for current energy
!     v(nr)       complex potential including energy dep xc
!     vval(nr)    as above,but xc from valence electrons only
!     em          current energy
!
!     WORKSPACE
!     vxcrmu and vxcimu are calculated only on first entry for a
!     particular unique potential, re-used on subsequent entries.
!     vxcrmu(nr)  real part of xc at fermi level
!     vxcimu(nr)  imag part of xc at fermi level
!     gsrel(nr)   ratio of gs xc potentials with and without magnetization
!     vvxcrm(nr)  real part of xc at fermi level from valence electrons
!     vvxcim(nr)  imag part of xc at fermi level from valence electrons



      dimension   vtot(nrptx), vvalgs(nrptx), densty(nrptx)
      dimension   dmag(nrptx), denval(nrptx)
      complex*16  em, eref, v(nrptx), vval(nrptx)
      dimension   vxcrmu(nrptx), vxcimu(nrptx)
      dimension   vvxcrm(nrptx), vvxcim(nrptx), gsrel(nrptx)

!     Josh added variables:
!     ZRnrm      - renormalization constant
!     csig       - control for using many pole self energy
!     NRPts      - Number of points to use inside atom.
!                  Other points are linearly interpolated.
!     WpCorr     - Array of frequencies for many pole self energy.
!     rsTmp      - Temp var for rs
!     WpTmp      - Temp var for Wp
!     AmpFac     - g_i (pole strengths)
!     Gamma      - pole broadening
!     RsInt      - Rs in the intersitial
!     DRs        - RsCore - RsInt
!     delrHL     - Re[delta sigma] for many pole self energy
!     deliHL     - Im[delta sigma] for many pole self energy
!     Rs1(NRPts) - Array of Rs points for interpolation

      character*512 slog      
      logical csig
      integer NRPts, NPoles
      parameter (tol=0.0004)
      parameter (NRPts=10)
      complex*16  delta, deltav, ZRnrm, cmu, SigF(NRPts), deltaHL(NRPts)
      double precision WpCorr(MxPole), rsTmp, WpTmp,                    &
     &     AmpFac(MxPole), Gamma(MxPole)
      double precision RsInt, DRs, delrHL(NRPts), deliHL(NRPts),        &
     &     Rs1(NRPts), EGap, RsMin, RsMax
      character(LEN=10) ColumnLabels(20)
      save SigF
!    Josh END

!     First calculate vxc to correct the local momentum dispersion
!     relation, delta = vxc(e,k) - vxc(mu,k), and
!               p^2 = k^2 -mu + kf^2 - delta.
!     In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
!                          = vcoul(r) + vxcgs(r) + delta(e,r).


!    at jri potential is smooth continuation of potential to r(jri)
!    at this point potential jumps to interstitial value at jri+1
      csig=.false.
      jri1 = jri + 1
      nmax=1
      nul=0
      ibp = index / 10
      ixc = mod(index,10)
      ixcTmp=ixc
      ZRNrm = 1.d0
      if((iPl.gt.0).and.(ixc.eq.0)) then
         csig=.true.
      end if
      if (ixc .eq. 2 .or. dble(em).le.xmu)  then
         do 10  i = 1, jri1
            v(i) = vtot(i)
            vval(i) = vvalgs(i)
   10    continue
!        Ground state exchange, no self energy calculation
         goto 888
      endif
!     Josh - Added CSigma to calculate HL Sigma with broadening and such.
!     Calculate Rs at the core and interstitial densities.
      if (densty(jri1).le.0) then
         RsInt =10
      else
         RsInt = (3 / (4*pi*densty(jri1))) ** third
      endif
      if (densty(1).le.0) then
         rscore =101.d0
      else
         rscore = (3 / (4*pi*densty(1))) ** third
      endif
      if (MAXVAL(densty(1:jri1)).le.0) then
         RsMin =RsInt*1.d-3
      else
         RsMin = (3 / (4*pi*MAXVAL(densty(1:jri1)))) ** third
      endif
      if (MINVAL(densty(1:jri1)).le.0) then
         RsMax=RsInt*2.d0
      else
         RsMax = (3 / (4*pi*MINVAL(densty(1:jri1)))) ** third
      endif
      !RsMax = RsInt

      DRs = (RsMax-RsMin)/(NRPts-2)
      IF(iPl.EQ.2) THEN
         Rs1(1) = RsMin
         DO i= 2, NRPts
            IF(RsMax.GT.RsInt) THEN
               IF(Rs1(i-1).LT.RsInt) THEN
                  Rs1(i) = RsMin + DBLE(i-1)*DRs
                  IF(Rs1(i).GT.RsInt) Rs1(i) = RsInt
               ELSE  
                  IF(i.GT.2) THEN
                     Rs1(i)=RsMin+DBLE(i-2)*DRs
                  END IF
               END IF
            ELSE
               Rs1(i) = RsMin + DBLE(i-1)*DRs
            END IF
         END DO
      ELSE
         Rs1(:) = RsInt
      END IF
      !WRITE(66,*) Rs1
!     Now calculate delta sigma as a function of Rs and energy
      if (csig) then
         ! Count the number of poles
         do i = 1, MxPole
            if(WpCorr(i).lt.-1.d0) then
               NPoles = i - 1
               EXIT
            end if
         end do
         ! Self energy at the fermi level is calculated once only since it is independent of
         ! energy.
         if(ifirst .eq. 0)  then            
            cmu = xmu*1.00001d0
            do i= 1, NRPts
               SigF(i) = 0.d0
               !Rs1(i)=RsMin+DBLE(i-1)*DRs
               
               ! If iPl = 1, use r independent sigma with parameters calculated
               ! from the interstitial density, i.e. bulk self-energy
               ! If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RInt)]
               if((iPl.eq.2).or.(i.eq.NRPts)) then
                  ! The next line will not work with ipl = 2. Get rid of ipl = 2 in future.
                  WpCorr(:) = WpCorr(:)*SQRT(3.d0/Rs1(i)**3)
                  call CSigZ(cmu,xmu,Rs1(i),SigF(i),ZRnrm,WpCorr,Gamma,AmpFac,EGap,NPoles,.TRUE.,.FALSE.)
                  WpCorr(:) = WpCorr(:)/SQRT(3.d0/Rs1(i)**3)
               end if
            end do
            if (iPl.eq.1) then
               SigF(1:NRPts-1) = SigF(NRPts)
            end if            
         end if

         do i= 1, NRPts
            deltaHL(i) = 0.d0
            !Rs1(i)=RsMin+DBLE(i-1)*DRs
            ! If iPl = 1, use r independent sigma with parameters calculated
            ! from the interstitial density, i.e. bulk self-energy
            ! If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RInt)]
            if((iPl.eq.2).or.(i.eq.NRPts)) then
               WpCorr(:) = WpCorr(:)*SQRT(3.d0/Rs1(i)**3)
               call CSigZ(em,xmu,Rs1(i),deltaHL(i),ZRnrm,WpCorr,Gamma,AmpFac,EGap,NPoles,.TRUE.,.FALSE.)
               WpCorr(:) = WpCorr(:)/SQRT(3.d0/Rs1(i)**3)
               deltaHL(i) = ZRnrm*(deltaHL(i) - SigF(i))
            end if
         end do
         if (iPl.eq.1) then
            deltaHL(1:NRPts-1) = deltaHL(NRPts)
         end if
         delrHL(:) = DBLE(deltaHL(:))
         deliHL(:) = DIMAG(deltaHL(:))

      end if

!     END Josh
      
!     Add the self energy correction
      do 20  i =  jri1,1,-1
         niter = 0
         if (densty(i).le.0) then
            rs =10
         else
            rs = (3 / (4*pi*densty(i))) ** third
         endif
!         write(22,*) 1.d0*exp(dble(i)*0.01), densty(i)         
!        Josh - If csigma is turned on, interpolate onto rs.
!        Then skip to 15 (skip other calculations and self
!        consistency)
         ! Test with constant SE.
         if(.FALSE.) THEN
            delrHL(:) = 0.d0/hart
            deliHL(:) = -5.d0/hart
         end if
         
         if(csig) then
            IF(iPl.EQ.2) THEN
               IF((rs.LT.RsMin).OR.(rs.GT.RsMax)) THEN
                  delrHL = 0.d0
               ELSE
                  call terp (Rs1, delrHL, NRPts, 1, rs, delr)
                  call terp (Rs1, deliHL, NRPts, 1, rs, deli)
               END IF
            ELSE
               delr = delrHL(1)
               deli = deliHL(1)
            END IF
            goto 15
         end if
!        END Josh
         
!        xf = 1.9191.../rs
         xf = fa / rs
         rsm = rs / (1+dmag(i))**third
         xfm = fa / rsm

         if (ixc.eq.5) then
            if ( denval(i) .gt. 0.00001) then
               rsval = (3 / (4*pi*denval(i))) ** third
               if (rsval.gt.10.0) rsval=10.0
            else
               rsval = 10.0
            endif
            xfval = fa / rsval
         elseif (ixc.ge.6) then
            if (densty(i) .le. denval(i) ) then
               rscore = 101.0
            else
               rscore = (3 / (4*pi*(densty(i)-denval(i)))) ** third
            endif
         endif

         if (ifirst .eq. 0)  then
!           vxc_mu indep of energy, calc only once
!           Calculate vxc at fermi level e = mu, j.m. 1/12/89
            xk = xf * 1.00001
            gsrel(i) = 1.0d0
            if (ixc .lt. 5) then
              call sigma(ixc, ibp,rs,rscore,xk,vxcrmu(i),vxcimu(i))
              if (index .eq. 0) then
!  do not need 4 following lines for gs difference in potential
!                xmag = 1.0d0+ dmag(i)
!                call vbh(rs,xmag,v1)
!                call vbh(rs, 1.0d0,v0)
!                if (v0 .ne. 0) gsrel(i) = v1/v0
              endif
            else
              call sigma(nul,ibp, rs, rscore,xk,vxcrmu(i),vxcimu(i))
            endif
            if (ixc.eq.5 ) then
               xkpp = xfval * 1.00001
               call sigma                                               &
     &         (ixc, ibp, rsval, rscore, xkpp, vvxcrm(i),vvxcim(i))
               if (ixc.eq.5 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            elseif (ixc .ge. 6) then
               call sigma                                               &
     &         (ixc, ibp, rs, rscore, xk, vvxcrm(i), vvxcim(i))
               if (ixc.eq.6 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            else
               vvxcrm(i) = 0.0d0
               vvxcim(i) = 0.0d0
            endif
         endif

!        xk2 is the local momentum squared, p^2 = k^2 - 2*mu + kf^2,
!        k^2 represents energy measured from vacuum.
!        See formula 2.15 in Lee and Beni's paper with the last 2
!        terms neglected.  (complete reference?)
         xk2 = 2 * (dble(em) - xmu) + xf**2
         xk = sqrt(xk2)
         xkm2 = 2 * (dble(em) - xmu) + xfm**2
!        quick fix
         if (xkm2.lt.0) xkm2=xk2
         xkm = sqrt(xkm2)

!        find \delta_1
         if (ixc .lt. 5) then
            call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
         else
            call sigma (nul, ibp, rs, rscore, xk, vxcr, vxci)
         endif
         del1r = gsrel(i) * (vxcr - vxcrmu(i))

!        Correct local momentum according to the formula
!        p^2 = k^2 - 2*mu + kf^2 - 2*delta.  Note that imag part
!        of delta is ignored, since xk2 is a real quantity.

!        find xk(em) by iterative solution of dyson equation
  50     continue
         xk2 = 2*(dble(em) - xmu - del1r) + xf**2
         if (xk2 .lt. 0)  then
            write(slog,'(1pe13.5, 3i8, a)') xk2, i, ie, iph, ' xk2, i, ie, iph'
            call wlog(slog)
            call wlog(' em, xf**2, xmu, delta')
            write(slog,'(1p, 5e13.5)') dble(em), xf**2, xmu, del1r
            call wlog(slog)
            call par_stop('XCPOT-2')
         endif
         xk = sqrt (xk2)

!        calculate \delta_2 and \delta_v,2 with the corrected
!        local momentum
         call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
!        delta corrected calculated with new local momentum
         delr = gsrel(i) * (vxcr - vxcrmu(i))
         deli = vxci-vxcimu(i)

         if (ixc.ge.5 .and. i.eq.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5 .or. ixc.eq.6) then
               delvr = delr
               delvi = deli
            endif
         endif

         if (niter.lt.nmax) then
            del1r=delr
            niter=niter+1
            go to 50
         endif

         if (ixc .ge. 5 .and. i.lt.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5) then
               xkpp=sqrt(xk**2-xf**2+xfval**2)
               call sigma (ixc, ibp, rsval,rscore,xkpp,vxcvr,vxcvi)
            else
               call sigma (ixc, ibp, rs, rscore, xk, vxcvr, vxcvi)
            endif
            delvr = vxcvr-vvxcrm(i)
            delvi = vxcvi-vvxcim(i)
         endif

!        Josh - Skip SC loop if CSigma is called. CSigma calculates self consistently.
 15      continue
         
         delta = dcmplx(delr,deli)

!	 Josh - write out delta sigma at interstitial level to sigma.dat.
         if(i.eq.jri1) then
!            write(45,'(X,20e14.6)') (DBLE(em) - xmu)*hart, delr*hart,   &
!     &                        deli*hart, DBLE(ZRnrm), DIMAG(ZRnrm),     &
!     &                        SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2),     &
!     &                        ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),          &
!     &                        SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr
         ColumnLabels(:) = ' '
         ColumnLabels(1) = 'E-Mu'
         ColumnLabels(2) = 'Re[Sigma]'
         ColumnLabels(3) = 'Im[Sigma]'
         ColumnLabels(4) = 'Re[Z]'
         ColumnLabels(5) = 'Im[Z]'
         ColumnLabels(6) = '|Z|'
         ColumnLabels(7) = 'Phase[Z]'
         ColumnLabels(8) = 'IMFP'
            CALL WriteData('mpse.dat',                                  & ! Specify file name.
                 & Double1 = (DBLE(em) - xmu)*hart,                     & ! Specify 1st col.
                 & DComplex2 = delta*hart,                              & ! Specify 2nd col.
                 & DComplex3 = ZRnrm,                                   & ! Specify 3nd col.
                 & Double4 = SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2),      & ! Specify 4rd col.
                 & Double5 = ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),           & ! 5th col.
                 & Double6 = SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr,    & ! 6th col
                 & ColumnLabels = ColumnLabels)
         end if
!	 Josh END

         if (ixc .eq. 5) delta = dcmplx(delr,delvi)
         v(i) = vtot(i) + delta
         if (ixc .ge. 5) then
            deltav = dcmplx(delvr,delvi)
            vval(i) = vvalgs(i) + deltav
         endif
 20   continue
 25   continue

      ifirst = 1

!     Reference the potential with respect to mt potential, ie,
!     first interstitial point.  v(jri1) = 0

!     Note that the reference does not contain the core hole lifetime
!     since the total atomic potential should have it. However in the
!     perturbation  deltav = v - vmt it cancels out.
!     ( deltav = vat - igamma - (vatmt-igamma) ).

 888  eref = v(jri1)
      do 910 i = 1, jri1
  910 v(i) = v(i) - eref
      if (ixc.ge.5) then
         do 920 i = 1, jri1
  920    vval(i) = vval(i) - eref
      else
         do 930 i = 1, jri1
  930    vval(i) = v(i)
      endif

!     Real self energy, zero imag part
      if (lreal.gt.0)  then
         do 950  i = 1, jri1
            v(i) = dble(v(i))
            if (ixc.gt.4)  vval(i) = dble(vval(i))
  950    continue
         eref = dble(eref)
      endif


      return
      end

      subroutine sigma (ixc, ibp, rs, rscore, xk, vr, vi)
      implicit double precision (a-h, o-z)

      if ((ixc.eq.0 .or. ixc.ge.5) .and. ibp .eq. 0) then
         call rhl (rs, xk, vr, vi)
      elseif ((ixc.eq.0.or. ixc.ge.5) .and. ibp .eq. 1) then
         call rhlbp (rs, xk, vr, vi)
      elseif (ixc .eq. 1) then
         vi = 0
         call edp(rs,xk,vr)
      elseif (ixc .eq. 3) then
         call edp(rs,xk,vr)
         call imhl (rs,xk,vi,icusp)
      endif

      if (ixc .ge. 6) then
         call edp(rscore,xk,vrp)
         vr = vr - vrp
      endif

      return
      end

