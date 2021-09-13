!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xcpot.f90,v $:
! $Revision: 1.20 $
! $Author: hebhop $
! $Date: 2012/03/29 22:52:37 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xcpot(iph, ie, index, lreal, ifirst, jri, em_in, xmu,      &
  &                 vtot, vvalgs, densty, dmag, denval,               &
  &                 eref, v, vval, iPl, WpCorr, Gamma, AmpFac, EGap,  &
  &                 vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, fname)

  USE IOMOD
  USE dimsmod, only: nrptx, MxPole
  USE constants
  USE potential_inp, ONLY: scf_temperature
  IMPLICIT NONE

  !  implicit double precision (a-h, o-z)
  !  calculate self-energy correction
  !  first coded j. mustre de leon
  !  last modified a.ankudinov 1996 for non-local self-energies
  !  Ankudinov, Rehr, J. Physique IV, vol. 7, C2-121, (1997).

  !  INPUT
  !  iph, ie used only for debug and labels.
  !  index       0  Hedin-Lunqvist + const real & imag part
  !              1  Dirac-Hara + const real & imag part
  !              2  ground state + const real & imag part
  !              3  Dirac-Hara + HL imag part + const real & imag part
  !              4  See rdinp for comment
  !  lreal       not equal zero for real self energy
  !  ifirst      first entry flag, set to zero before first call for
  !              each unique potential, see vxcrmu and vxcimu below
  !  jri         index of first interstitial point in current
  !              Loucks r grid
  !  em          current energy grid point
  !  xmu         fermi level
  !  vi0         const imag part to subtract from potential
  !  gamach      core hole lifetime
  !  vtot(nr)    total potential (coulomb and gs exchange corr)
  !  vvalgs(nr)  total coulomb + gs xc potential from valence electrons
  !  densty(nr)  electron density
  !  dmag(nr)    density magnetization
  !  denval(nr)  valence electron density
  !  iPl         Control for many pole self energy (Josh)
  INTEGER, INTENT(INOUT) :: iph, ie, index, lreal, jri, iPl
  INTEGER, INTENT(INOUT) :: ifirst
  COMPLEX*16, INTENT(INOUT) :: em_in
  REAL*8, INTENT(INOUT) :: xmu
  REAL*8, DIMENSION(nrptx), INTENT(INOUT) :: vtot, vvalgs, dmag, denval, densty
  !  OUTPUT
  !  eref        complex energy reference for current energy
  !  v(nr)       complex potential including energy dep xc
  !  vval(nr)    as above,but xc from valence electrons only
  !  em          current energy
  COMPLEX*16, INTENT(INOUT) :: eref
  COMPLEX*16, DIMENSION(nrptx), INTENT(INOUT) :: v, vval

  !  WORKSPACE
  !  vxcrmu and vxcimu are calculated only on first entry for a
  !  particular unique potential, re-used on subsequent entries.
  !  vxcrmu(nr)  real part of xc at fermi level
  !  vxcimu(nr)  imag part of xc at fermi level
  !  gsrel(nr)   ratio of gs xc potentials with and without magnetization
  !  vvxcrm(nr)  real part of xc at fermi level from valence electrons
  !  vvxcim(nr)  imag part of xc at fermi level from valence electrons
  REAL*8, DIMENSION(nrptx), INTENT(INOUT) :: gsrel
  REAL*8, DIMENSION(nrptx), INTENT(INOUT) :: vxcrmu, vxcimu
  REAL*8, DIMENSION(nrptx), INTENT(INOUT) :: vvxcrm, vvxcim
  CHARACTER(LEN=50) :: fname

  REAL*8 :: xkpp, xkm2, xkm, xk2, xk, xfval, xf, xfm
  REAL*8 :: vxcvi, vxcvr, vxci, vxcr, rsm, rsval, rs, rscore
  REAL*8 :: delr, delvi, delvr, deli, del1r, tol
  COMPLEX*16 :: em
  INTEGER :: nul, nmax, jri1, niter, ixc, ixcTmp, i, ibp
  !  Josh added variables:
  !  ZRnrm      - renormalization constant
  !  csig       - control for using many pole self energy
  !  NRPts      - Number of points to use inside atom.
  !               Other points are linearly interpolated.
  !  WpCorr     - Array of frequencies for many pole self energy.
  !  rsTmp      - Temp var for rs
  !  WpTmp      - Temp var for Wp
  !  AmpFac     - g_i (pole strengths)
  !  Gamma      - pole broadening
  !  RsInt      - Rs in the intersitial
  !  DRs        - RsCore - RsInt
  !  delrHL     - Re[delta sigma] for many pole self energy
  !  deliHL     - Im[delta sigma] for many pole self energy
  !  Rs1(NRPts) - Array of Rs points for interpolation

  CHARACTER*512 slog      
  LOGICAL csig
  INTEGER NRPts, NPoles
  PARAMETER (tol=0.0004)
  PARAMETER (NRPts=10)
  COMPLEX*16 :: delta, deltav, ZRnrm, cmu, SigF(NRPts), deltaHL(NRPts)
  DOUBLE PRECISION WpCorr(MxPole), rsTmp, WpTmp, AmpFac(MxPole), Gamma(MxPole)
  DOUBLE PRECISION RsInt, DRs, delrHL(NRPts), deliHL(NRPts), Rs1(NRPts), EGap, RsMin, RsMax
  CHARACTER(LEN=10) ColumnLabels(20)
  SAVE SigF
  ! Josh END

  !  First calculate vxc to correct the local momentum dispersion
  !  relation, delta = vxc(e,k) - vxc(mu,k), and
  !            p^2 = k^2 -mu + kf^2 - delta.
  !  In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
  !                       = vcoul(r) + vxcgs(r) + delta(e,r).


  ! at jri potential is smooth continuation of potential to r(jri)
  ! at this point potential jumps to interstitial value at jri+1
  csig = .FALSE.
  jri1 = jri + 1
  nmax = 1
  nul = 0
  ibp = index / 10
  ixc = MOD(index, 10)
  ixcTmp = ixc
  ZRNrm = 1.d0
  IF ((iPl.GT.0).AND.(ixc.EQ.0)) THEN
    csig = .TRUE.
  ENDIF

  ! IF (ixc .EQ. 2 .OR. DBLE(em).LE.xmu) THEN
  !    DO i = 1, jri1
  !       v(i) = vtot(i)
  !       vval(i) = vvalgs(i)
  !    ENDDO
  !    !  Ground state exchange, no self energy calculation
  !    GOTO 888
  ! ENDIF

  ! TTS - Apply self-energy calculation below Fermi level at finite temperature
  em = em_in
  IF (index.EQ.2) THEN
    DO i = 1, jri1
      v(i) = vtot(i)
      vval(i) = vvalgs(i)
    ENDDO
    !  Ground state exchange, no self energy calculation
    GOTO 888
  ELSE!IF ((index.EQ.6).OR.(index.EQ.7)) THEN
    ! For finite temperature self-energies
    IF (DBLE(em_in).LE.xmu) THEN
      ! Use self energy at xmu for E < xmu
      IF (xmu.LT.0.0d0) THEN
        em = DCMPLX(xmu*0.9999d0, DIMAG(em_in))
      ELSE
        em = DCMPLX(xmu*1.0001d0, DIMAG(em_in))
      ENDIF
    ENDIF
  ENDIF
  IF (index.EQ.-2) index = 6 ! JJK - calculate self-energy if temperature dependent ground state

  !  Josh - Added CSigma to calculate HL Sigma with broadening and such.
  !  Calculate Rs at the core and interstitial densities.
  IF (densty(jri1).LE.0) THEN
    RsInt =10
  ELSE
    RsInt = (3 / (4*pi*densty(jri1))) ** third
  ENDIF
  IF (densty(1).LE.0) THEN
    rscore =101.d0
  ELSE
    rscore = (3 / (4*pi*densty(1))) ** third
  ENDIF
  IF (MAXVAL(densty(1:jri1)).LE.0) THEN
    RsMin =RsInt*1.d-3
  ELSE
    RsMin = (3 / (4*pi*MAXVAL(densty(1:jri1)))) ** third
  ENDIF
  IF (MINVAL(densty(1:jri1)).LE.0) THEN
    RsMax=RsInt*2.d0
  ELSE
    RsMax = (3 / (4*pi*MINVAL(densty(1:jri1)))) ** third
  ENDIF

  !RsMax = RsInt

  DRs = (RsMax-RsMin)/(NRPts-2)
  IF (iPl.EQ.2) THEN
    Rs1(1) = RsMin
    DO i= 2, NRPts
        IF (RsMax.GT.RsInt) THEN
          IF (Rs1(i-1).LT.RsInt) THEN
              Rs1(i) = RsMin + DBLE(i-1)*DRs
              IF (Rs1(i).GT.RsInt) Rs1(i) = RsInt
          ELSE  
              IF (i.GT.2) THEN
                Rs1(i) = RsMin + DBLE(i-2)*DRs
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
  !  Now calculate delta sigma as a function of Rs and energy
  IF (csig) THEN
    ! Count the number of poles
    DO i = 1, MxPole
      IF (WpCorr(i).LT.-1.d0) THEN
        NPoles = i - 1
        EXIT
      ENDIF
    ENDDO
    ! Self energy at the fermi level is calculated once only since it is independent of
    ! energy.
    IF (ifirst.EQ.0) THEN            
      cmu = xmu*1.00001d0
      DO i = 1, NRPts
        SigF(i) = 0.d0
        !Rs1(i)=RsMin+DBLE(i-1)*DRs
        
        ! If iPl = 1, use r independent sigma with parameters calculated
        ! from the interstitial density, i.e. bulk self-energy
        ! If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RInt)]
        IF ((iPl.EQ.2).OR.(i.EQ.NRPts)) THEN
          ! The next line will not work with ipl = 2. Get rid of ipl = 2 in future.
          WpCorr(:) = WpCorr(:)*SQRT(3.d0/Rs1(i)**3)
          CALL CSigZ(cmu,xmu,Rs1(i),SigF(i),ZRnrm,WpCorr,Gamma,AmpFac,EGap,NPoles,.TRUE.,.FALSE.)
          WpCorr(:) = WpCorr(:)/SQRT(3.d0/Rs1(i)**3)
        ENDIF
      ENDDO
      IF (iPl.EQ.1) THEN
        SigF(1:NRPts-1) = SigF(NRPts)
      ENDIF          
    ENDIF

    DO i = 1, NRPts
      deltaHL(i) = 0.d0
      !Rs1(i)=RsMin+DBLE(i-1)*DRs
      ! If iPl = 1, use r independent sigma with parameters calculated
      ! from the interstitial density, i.e. bulk self-energy
      ! If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RInt)]
      IF ((iPl.EQ.2).OR.(i.EQ.NRPts)) THEN
        WpCorr(:) = WpCorr(:)*SQRT(3.d0/Rs1(i)**3)
        CALL CSigZ(em,xmu,Rs1(i),deltaHL(i),ZRnrm,WpCorr,Gamma,AmpFac,EGap,NPoles,.TRUE.,.FALSE.)
        WpCorr(:) = WpCorr(:)/SQRT(3.d0/Rs1(i)**3)
        deltaHL(i) = ZRnrm*(deltaHL(i) - SigF(i))
      ENDIF
    ENDDO
    IF (iPl.EQ.1) THEN
      deltaHL(1:NRPts-1) = deltaHL(NRPts)
    ENDIF
    delrHL(:) = DBLE(deltaHL(:))
    deliHL(:) = DIMAG(deltaHL(:))
  ENDIF

  !  END Josh
  !  Add the self energy correction
  DO i = jri1, 1, -1
    niter = 0
    IF (densty(i).LE.0) THEN
      rs = 10
    ELSE
      rs = (3 / (4*pi*densty(i))) ** third
    ENDIF
    IF(ixc.EQ.-2) rs = RsInt ! DEBUG - JK use only constant imaginary part at the interstitial for 
                             !         temperature dependent ground state calculations.
    !   write(22,*) 1.d0*exp(dble(i)*0.01), densty(i)         
    !  Josh - If csigma is turned on, interpolate onto rs.
    !  Then skip to 15 (skip other calculations and self
    !  consistency)
    ! Test with constant SE.
    IF (.FALSE.) THEN
      delrHL(:) = 0.d0/hart
      deliHL(:) = -5.d0/hart
    ENDIF
    
    IF (csig) THEN
      IF (iPl.EQ.2) THEN
        IF ((rs.LT.RsMin).OR.(rs.GT.RsMax)) THEN
          delrHL = 0.d0
        ELSE
          CALL terp(Rs1, delrHL, NRPts, 1, rs, delr)
          CALL terp(Rs1, deliHL, NRPts, 1, rs, deli)
        END IF
      ELSE
        delr = delrHL(1)
        deli = deliHL(1)
      END IF
      GOTO 15
    ENDIF
    !  END Josh
    
    ! xf = 1.9191.../rs
    xf = fa / rs
    rsm = rs / (1+dmag(i))**third
    xfm = fa / rsm

    ! IF (ixc.EQ.5) THEN
    !    IF ( denval(i).GT.0.00001) THEN
    !       rsval = (3 / (4*pi*denval(i))) ** third
    !       IF (rsval.GT.10.0) rsval=10.0
    !    ELSE
    !       rsval = 10.0
    !    ENDIF
    !    xfval = fa / rsval
    ! ELSEIF (ixc.GE.6) THEN
    !    IF (densty(i).LE.denval(i)) THEN
    !       rscore = 101.0
    !    ELSE
    !       rscore = (3 / (4*pi*(densty(i)-denval(i)))) ** third
    !    ENDIF
    ! ENDIF

    IF ((index.EQ.5).OR.(index.EQ.15)) THEN
      IF ( denval(i).GT.0.00001) THEN
        rsval = (3 / (4*pi*denval(i))) ** third
        IF (rsval.GT.10.0) rsval=10.0
      ELSE
        rsval = 10.0
      ENDIF
      xfval = fa / rsval
    ELSEIF (index.EQ.9) THEN
      IF (densty(i).LE.denval(i)) THEN
        rscore = 101.0
      ELSE
        rscore = (3 / (4*pi*(densty(i)-denval(i)))) ** third
      ENDIF
    ENDIF

    ! IF (ifirst.EQ.0) THEN
    !    !  vxc_mu indep of energy, calc only once
    !    !  Calculate vxc at fermi level e = mu, j.m. 1/12/89
    !    xk = xf * 1.00001
    !    gsrel(i) = 1.0d0
    !    IF (ixc.LT.5) THEN
    !       CALL sigma(ixc, ibp,rs,rscore,xk,vxcrmu(i),vxcimu(i),.TRUE.)
    !       IF (index.EQ.0) THEN
    !       !  do not need 4 following lines for gs difference in potential
    !          ! xmag = 1.0d0+ dmag(i)
    !          ! call vbh(rs,xmag,v1)
    !          ! call vbh(rs, 1.0d0,v0)
    !          ! if (v0 .ne. 0) gsrel(i) = v1/v0
    !       ENDIF
    !    ELSE
    !       CALL sigma(nul,ibp, rs, rscore,xk,vxcrmu(i),vxcimu(i),.TRUE.)
    !    ENDIF
    !    IF (ixc.EQ.5) THEN
    !       xkpp = xfval * 1.00001
    !       CALL sigma(ixc, ibp, rsval, rscore, xkpp, vvxcrm(i),vvxcim(i),.FALSE.)
    !       IF (ixc.EQ.5 .AND. i.EQ.jri1) THEN
    !          vvxcrm(jri1) =  vxcrmu(jri1)
    !          vvxcim(jri1) =  vxcimu(jri1)
    !       ENDIF
    !    ELSEIF (ixc.GE.6) THEN
    !       CALL sigma(ixc, ibp, rs, rscore, xk, vvxcrm(i), vvxcim(i),.FALSE.)
    !       IF (ixc.EQ.6.AND.i.EQ.jri1) THEN
    !          vvxcrm(jri1) =  vxcrmu(jri1)
    !          vvxcim(jri1) =  vxcimu(jri1)
    !       ENDIF
    !    ELSE
    !       vvxcrm(i) = 0.0d0
    !       vvxcim(i) = 0.0d0
    !    ENDIF
    ! ENDIF

    IF (ifirst.EQ.0) THEN
      !  vxc_mu indep of energy, calc only once
      !  Calculate vxc at fermi level e = mu, j.m. 1/12/89
      xk = xf * 1.00001
      gsrel(i) = 1.0d0 
      IF ((index.EQ.0).OR.(index.EQ.1).OR.(index.EQ.2).OR.(index.EQ.3).OR.(index.EQ.6).OR.(index.EQ.7) &
          .OR.(index.EQ.11).OR.(index.EQ.13)) THEN
        CALL sigma(index, ibp,rs, rscore, xk, vxcrmu(i), vxcimu(i), xmu, .TRUE.)
        ! Only get the real part
      ELSE
        ! Use index = 0 -> HL
        ! index = 5, 15, 9
        CALL sigma(nul, ibp, rs, rscore, xk, vxcrmu(i), vxcimu(i), xmu, .FALSE.)
      ENDIF
  
      IF ((index.EQ.5).OR.(index.EQ.15)) THEN
        xkpp = xfval * 1.00001
        CALL sigma(index, ibp, rsval, rscore, xkpp, vvxcrm(i),vvxcim(i), xmu, .FALSE.)
        IF (((index.EQ.5).OR.(index.EQ.15)).AND. i.EQ.jri1) THEN
          vvxcrm(jri1) =  vxcrmu(jri1)
          vvxcim(jri1) =  vxcimu(jri1)
        ENDIF
      ELSEIF (index.EQ.9) THEN
        ! index = 9
        CALL sigma(index, ibp, rs, rscore, xk, vvxcrm(i), vvxcim(i), xmu, .FALSE.)
        IF (index.EQ.9.AND.i.EQ.jri1) THEN
          vvxcrm(jri1) =  vxcrmu(jri1)
          vvxcim(jri1) =  vxcimu(jri1)
        ENDIF
      ELSE
        ! index = 0, 1, 11, 2, 3, 13, 6, 7
        vvxcrm(i) = 0.0d0
        vvxcim(i) = 0.0d0
      ENDIF
    ENDIF

    !  xk2 is the local momentum squared, p^2 = k^2 - 2*mu + kf^2,
    !  k^2 represents energy measured from vacuum.
    !  See formula 2.15 in Lee and Beni's paper with the last 2
    !  terms neglected.  Phys Rev B 15 (6) 2862 (1977)
    xk2 = 2 * (DBLE(em) - xmu) + xf**2
    xk = SQRT(xk2)
    xkm2 = 2 * (DBLE(em) - xmu) + xfm**2
    !  quick fix
    IF (xkm2.LT.0) xkm2 = xk2
    xkm = SQRT(xkm2)

    !  find \delta_1
    ! IF (ixc .lt. 5) THEN
    !    CALL sigma(ixc, ibp, rs, rscore, xk, vxcr, vxci,.FALSE.)
    ! ELSE
    !    CALL sigma(nul, ibp, rs, rscore, xk, vxcr, vxci,.FALSE.)
    ! ENDIF
    IF ((index.EQ.0).OR.(index.EQ.1).OR.(index.EQ.2).OR.(index.EQ.3).OR.(index.EQ.6).OR.(index.EQ.7) &
          .OR.(index.EQ.11).OR.(index.EQ.13)) THEN
      CALL sigma(index, ibp, rs, rscore, xk, vxcr, vxci, xmu, .FALSE.)
    ELSE
      ! index = 5, 15, 9, (4 ?)
      CALL sigma(nul, ibp, rs, rscore, xk, vxcr, vxci, xmu, .FALSE.)
    ENDIF
    del1r = gsrel(i) * (vxcr - vxcrmu(i))

    !  Correct local momentum according to the formula
    !  p^2 = k^2 - 2*mu + kf^2 - 2*delta.  Note that imag part
    !  of delta is ignored, since xk2 is a real quantity.

    !  find xk(em) by iterative solution of dyson equation
    50 CONTINUE
    
    xk2 = 2*(DBLE(em) - xmu - del1r) + xf**2
    IF (xk2.LT.0) THEN
      WRITE(slog,'(1pe13.5, 3i8, a)') xk2, i, ie, iph, ' xk2, i, ie, iph'
      CALL wlog(slog)
      CALL wlog(' em, xf**2, xmu, delta')
      WRITE(slog,'(1p, 5e13.5)') DBLE(em), xf**2, xmu, del1r
      CALL wlog(slog)
      CALL par_stop('XCPOT-2')
    ENDIF
    xk = SQRT(xk2)

    !  calculate \delta_2 and \delta_v,2 with the corrected local momentum
    CALL sigma(index, ibp, rs, rscore, xk, vxcr, vxci, xmu, .FALSE.)

    !  delta corrected calculated with new local momentum
    delr = gsrel(i) * (vxcr - vxcrmu(i))
    deli = vxci - vxcimu(i)

    IF (((index.EQ.5).OR.(index.EQ.15).OR.(index.EQ.9)).AND.(i.EQ.jri1).AND.(xk.GT.xf)) THEN
      delvr = delr
      delvi = deli
    ENDIF

    IF (niter.LT.nmax) THEN
      del1r=delr
      niter=niter+1
      GOTO 50
    ENDIF

    ! IF (ixc .GE. 5 .AND. i.LT.jri1 .AND. xk.GT.xf) THEN
    !    IF (ixc.EQ.5) THEN
    !       xkpp = SQRT(xk**2-xf**2+xfval**2)
    !       CALL sigma(ixc, ibp, rsval,rscore,xkpp,vxcvr,vxcvi,.FALSE.)
    !    ELSE
    !       CALL sigma(ixc, ibp, rs, rscore, xk, vxcvr, vxcvi,.FALSE.)
    !    ENDIF
    !    delvr = vxcvr-vvxcrm(i)
    !    delvi = vxcvi-vvxcim(i)
    ! ENDIF
    IF ( ((index.EQ.5).OR.(index.EQ.15).OR.(index.EQ.9)).AND.(i.LT.jri1).AND.(xk.GT.xf) ) THEN
      IF ((index.EQ.5).OR.(index.EQ.15)) THEN
        xkpp = SQRT(xk**2-xf**2+xfval**2)
        CALL sigma(index, ibp, rsval, rscore, xkpp, vxcvr, vxcvi, xmu, .FALSE.)
      ELSE
        ! index = 9
        CALL sigma(index, ibp, rs, rscore, xk, vxcvr, vxcvi, xmu, .FALSE.)
      ENDIF
      delvr = vxcvr-vvxcrm(i)
      delvi = vxcvi-vvxcim(i)
    ENDIF

    !  Josh - Skip SC loop if CSigma is called. CSigma calculates self consistently.
    15 CONTINUE
    
    delta = DCMPLX(delr, deli)

    !  Josh - write out delta sigma at interstitial level to sigma.dat.
    !  TTS  - Replaced em with em_in
    IF (i.EQ.jri1) THEN
        !         write(45,'(X,20e14.6)') (DBLE(em) - xmu)*hart, delr*hart,   &
        !  &                        deli*hart, DBLE(ZRnrm), DIMAG(ZRnrm),     &
        !  &                        SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2),     &
        !  &                        ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),          &
        !  &                        SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr
        ColumnLabels(:) = ' '
        ColumnLabels(1) = 'E-Mu'
        ColumnLabels(2) = 'Re[Sigma]'
        ColumnLabels(3) = 'Im[Sigma]'
        ColumnLabels(4) = 'Re[Z]'
        ColumnLabels(5) = 'Im[Z]'
        ColumnLabels(6) = '|Z|'
        ColumnLabels(7) = 'Phase[Z]'
        ColumnLabels(8) = 'IMFP'
        CALL WriteData(TRIM(fname)//'.dat',                          & ! Specify file name.
              & Double1 = (DBLE(em_in) - xmu)*hart,                  & ! Specify 1st col.
              & DComplex2 = delta*hart,                              & ! Specify 2nd col.
              & DComplex3 = ZRnrm,                                   & ! Specify 3nd col.
              & Double4 = SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2),      & ! Specify 4rd col.
              & Double5 = ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),           & ! 5th col.
              & Double6 = SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr,    & ! 6th col
              & ColumnLabels = ColumnLabels)
    ENDIF
    !  Josh END

    ! IF (ixc.EQ.5) delta = dcmplx(delr,delvi)
    IF ((index.EQ.5).OR.(index.EQ.15)) delta = DCMPLX(delr, delvi)

    v(i) = vtot(i) + delta

    IF ((index.EQ.5).OR.(index.EQ.15).OR.(index.EQ.9)) THEN
      deltav = DCMPLX(delvr, delvi)
      vval(i) = vvalgs(i) + deltav
    ENDIF
  ENDDO

  ifirst = 1

  !  Reference the potential with respect to mt potential, ie,
  !  first interstitial point.  v(jri1) = 0

  !  Note that the reference does not contain the core hole lifetime
  !  since the total atomic potential should have it. However in the
  !  perturbation  deltav = v - vmt it cancels out.
  !  ( deltav = vat - igamma - (vatmt-igamma) ).

888  IF (ifirst.eq.0) THEN

ENDIF
  
eref = v(jri1)

DO i = 1, jri1
  v(i) = v(i) - eref
ENDDO

IF ((index.EQ.5).OR.(index.EQ.15).OR.(index.EQ.9)) THEN
  DO i = 1, jri1
    vval(i) = vval(i) - eref
  ENDDO
ELSE
  ! index = 0, 1, 11, 2, 3, 13, 6, 7
  DO i = 1, jri1
    vval(i) = v(i)
  ENDDO
ENDIF

!  Real self energy, zero imag part
IF (lreal.GT.0) THEN
  DO i = 1, jri1
      v(i) = DBLE(v(i))
      ! IF (ixc.GT.4)  vval(i) = DBLE(vval(i))
      IF ((index.EQ.5).OR.(index.EQ.9).OR.(index.EQ.15)) vval(i) = DBLE(vval(i))
  ENDDO
  eref = DBLE(eref)
ENDIF
! We need to restore index to original value
IF (ixc.EQ.-2) index = -2
RETURN
END SUBROUTINE xcpot
  
SUBROUTINE sigma(ixc, ibp, rs, rscore, xk, vr, vi, xmu, realonly)
  ! USE m_COHSEX, ONLY: gw_exact
  USE potential_inp, ONLY: scf_temperature
  USE constants, ONLY: hart
  IMPLICIT NONE
  ! Computes the self-energy given radius rs and momentum xk
  ! ixc  - index
  ! ibp  - Deprecated
  ! rs   - Radius
  ! rscore -  Radius at core
  ! vr   - Real part in a.u.
  ! vi   - Imaginary part in a.u.
  ! realonly - Return only the real part of self-energy 
  INTEGER, INTENT(IN) :: ixc, ibp
  REAL(8), INTENT(IN) :: rs, rscore, xk, xmu
  REAL(8), INTENT(OUT) :: vr, vi
  LOGICAL, INTENT(IN) :: realonly

  INTEGER :: icusp
  REAL(8) :: temperature, vrp

  SELECT CASE(ixc)
    CASE(0)
      ! Hedin-Lundqvist + constant imaginary
      CALL rhl(rs, xk, vr, vi)
    CASE(1)
      ! Dirac-Hara + constant imaginary
      vi = 0.d0
      CALL edp(rs, xk, vr)
    CASE(2)
      ! Ground state + const real & imag part
      ! This seems to be redundant since we never call sigma
      vi = 0.d0
      vr = 0.d0
      ! temperature = scf_temperature / hart
      ! IF (realonly) THEN
      !   vi = 0.d0
      ! ELSE
      !   CALL imgw_v2(rs, xk, temperature, vi)
      ! ENDIF
    CASE(3)
      ! Dirac-Hara + HL imag part + constant imaginary
      CALL edp(rs, xk, vr)
      CALL imhl(rs, xk, vi, icusp)
    CASE(5)
      ! Partially nonlocal: Dirac-Fock for core + HL for valence electrons + constant imaginary
      ! TODO: Why is this the same as HL (case 0) ? Update documentation ?
      CALL rhl(rs, xk, vr, vi)
    CASE(6)
      ! Finite temperature RPA GW
      temperature = scf_temperature / hart
      CALL rgw(rs, xk, temperature, vr)
      IF (realonly) THEN
        vi = 0.d0
      ELSE
        CALL imgw_v2(rs, xk, temperature, vi)
      ENDIF
    CASE(7)
      ! Temperature independent RPA GW (For testing purposes)
      temperature = 0.d0
      CALL rgw(rs, xk, temperature, vr)
      IF (realonly) THEN
        vi = 0.d0
      ELSE
        CALL imgw_v2(rs, xk, temperature, vi)
      ENDIF
    CASE(9)
      ! TODO: What is this ?
      ! This is the original ixc > 6 
      ! IF (ixc.GE.6) THEN
      CALL edp(rscore, xk, vrp)
      vr = vr - vrp
      ! ENDIF
    CASE(10)
      ! Same as ixc=0 with broadened plasmon HL selfenergy
      ! using interpolation for the real and imaginary part.
      CALL rhlbp(rs, xk, vr, vi)
    CASE(13)
      ! Same as ixc=3 with broadened plasmon HL selfenergy
      ! using interpolation for the real and imaginary part.
      CALL edp(rs, xk, vr)
      CALL imhl(rs, xk, vi, icusp)
    CASE(15)
      ! Same as ixc=5 with broadened plasmon HL selfenergy
      ! using interpolation for the real and imaginary part.
      CALL rhlbp (rs, xk, vr, vi)
    CASE DEFAULT
      vr = 0.d0
      vi = 0.d0
  END SELECT
  ! original script
  ! if ((ixc.eq.0 .or. ixc.ge.5) .and. ibp .eq. 0) then
  !    call rhl (rs, xk, vr, vi)
  ! elseif ((ixc.eq.0.or. ixc.ge.5) .and. ibp .eq. 1) then
  !    call rhlbp (rs, xk, vr, vi)
  ! elseif (ixc .eq. 1) then
  !    vi = 0
  !    call edp(rs,xk,vr)
  ! elseif (ixc .eq. 3) then
  !    call edp(rs,xk,vr)
  !    call imhl (rs,xk,vi,icusp)
  ! endif

  ! if (ixc .ge. 6) then
  !    call edp(rscore,xk,vrp)
  !    vr = vr - vrp
  ! endif

  ! call gw_exact(rs, xk, vr, vi)
  RETURN
END SUBROUTINE sigma
  
  
