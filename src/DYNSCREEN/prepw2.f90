!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: prep.f90,v $:
! $Revision: 1.14 $
! $Author: hebhop $
! $Date: 2012/01/07 00:55:48 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     PREP
!=======================================================================
subroutine prepw2 (vr0, nrx, ri, dx, x0, ilast)
!KJ removed unused ibounc
  use DimsMod, only: natx, nrptx, nphx=>nphu, nex, lx, ltot
  use constants
  use atoms_inp
  use ldos_inp
  use screenw_inp
  use crpa_inp
  use par
  implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  integer, intent(in) :: nrx
  double precision, intent(in) :: vr0, dx, x0, ri(nrptx) 
  integer, intent(out) :: ilast
  logical :: CRPA
  !     W, the result of calculation
  double precision :: vch(nrx) 
  double precision fxc(nrx)
  integer isc
  integer, parameter :: mxsc = 3 
  complex*16 wscrn(nrx), wscrn_SC(mxsc,nrx), wscrn0(nrx)

  !     =================================================================
  !     function getiat
  integer getiat
  !     -----------------------------------------------------------------
  real srat(3,natx)
  !     pot.bin --------------------------------------------------------
  integer ntitle
  character*80 title(nheadx), flname
  integer nohole, ihole, inters, iafolp, jumprm, iunf
  double precision rnrmav, xmu, vint, rhoint
  double precision emu, s02, erelax, wp, ecv, rs, xf, qtotel, totvol
  integer imt(0:nphx), inrm(0:nphx)
  double precision rmt(0:nphx), rnrm(0:nphx)
  double precision folp(0:nphx),folpx(0:nphx)
  double precision dgc0(251), dpc0(251)
  double precision dgc(251, 41, 0:nphx+1), dpc(251, 41, 0:nphx+1)
  double precision adgc(10, 41, 0:nphx+1), adpc(10, 41, 0:nphx+1)
  double precision edens(251, 0:nphx), vclap(251, 0:nphx)
  double precision vtot(251, 0:nphx), edenvl(251, 0:nphx)
  double precision vvalgs(251, 0:nphx), dmag(251, 0:nphx)
  double precision xnval(41,0:nphx)
  double precision eorb(41), dx0
  complex*16 chi0_1(nrx,nrx), chi0_2(nrx,nrx), chi0(nrx,nrx)
  integer kappa(41)
  double precision qnrm(0:nphx)

  real*8, allocatable :: xnmues(:,:)

  integer iorb(-5:4,0:nphx)
  integer iz(0:nphx)
  double precision xion(0:nphx)
  double precision xnatph(0:nphx)
  !     -----------------------------------------------------------------
  double precision vtotph(nrptx,0:nphx), vvalph(nrptx,0:nphx)
  double precision rhoph(nrptx), rhphvl(nrptx), dmagx(nrptx)
  double precision dgcx(nrptx), dpcx(nrptx), rhoc(nrptx)
  double precision dgcn(nrptx,41,0:nphx), dpcn(nrptx,41,0:nphx)
  integer iph, jnew
  double precision edge, vjump

  integer jri, jri1, ne, ne1, ne3, i, j, k, iholetmp
  integer inclus, lmax(0:nphx)
  complex*16 ph1(nex, ltot+1, 0:nphx)
  complex*16 ph2(nex, ltot+1, 0:nphx)
  complex*16 ph3(nex, ltot+1, 0:nphx)
  complex*16 ph4(nex, ltot+1, 0:nphx)
  double precision de, enext
  complex*16 eref(nex), em(nex), em2(nex), em3(nex), em4(nex)
  integer ne2
  double precision domega, omega, omega_max, deltaE, emin2
  integer iomega, iomega0, nomega, ir, ie, ixc0
  COMPLEX*16 gtrl1(0:lx, nex), gtrl2(0:lx, nex), gtrl3(0:lx, nex), gtrl4(0:lx, nex)
  COMPLEX*16 ggmmp1(0:lx,nex,2*lx+1,2*lx+1), ggmmp2(0:lx,nex,2*lx+1,2*lx+1), ggmmp3(0:lx,nex,2*lx+1,2*lx+1), ggmmp4(0:lx,nex,2*lx+1,2*lx+1)
  COMPLEX*16 Kmat(nrx,nrx), Amat(nrx,nrx), Amat2(nrx,nrx), Qmatm1(nrx)
  COMPLEX*16 beta0(2), beta, beta2, beta3, beta4, beta_SC(mxsc), beta_tilde_SC(mxsc), alpha, alpha0
  COMPLEX*16 beta_decomp(nrx), betaL, betaR, beta_sum, beta_Fluct
  COMPLEX*16 norm, norm_max
  COMPLEX*16 EigVals(nrx), EigVecsL(nrx,nrx), EigVecsR(nrx,nrx), Work(2*nrx)
  DOUBLE PRECISION RWork(2*nrx)
  integer info, lda, ldb, ipiv(nrx), nrhs, ilast0, iproc, iprint, ifirst
  character, parameter :: trans = 'N'

  LOGICAL FreeElectronTest, DiagonalizeEps
  DOUBLE PRECISION kf

  DiagonalizeEps = .TRUE.
  ! Allocate variables
  allocate(xnmues(0:lx,0:nphx))

  FreeElectronTest = .FALSE.
  ifirst = 1
  iprint = 0
  ! Open files
  if(master) then
     OPEN(17,FILE='beta0.dat',STATUS='REPLACE')
     OPEN(18,FILE='beta.dat',STATUS='REPLACE')
     OPEN(20,FILE='betaSC.dat',STATUS='REPLACE')
     OPEN(21,FILE='eigvals.dat',STATUS='REPLACE')
     OPEN(23,FILE='beta_sum.dat',STATUS='REPLACE')
     OPEN(19,FILE='wscrnw.dat',STATUS='REPLACE')

     iprint = 1
  end if
  OPEN(22,FILE='eigvecs.dat',STATUS='REPLACE')
  !     read pot.bin
  call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,            &
       &                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
       &                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
       &                  dgc0, dpc0, dgc, dpc, adgc, adpc,               &
       &                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
       &                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
       &                  inters, totvol, iafolp, xion, iunf, iz, jumprm)

  dx0 = 0.05
!KJ here, it is necessary to be careful ...  If pot ran in k-space, rfms 
  rfms2  = real(ScreenI%rfms)
  rdirec = real(ScreenI%rfms*2)

  edge = xmu - vr0
  emu  = emu - vr0
  ! xmu for rutile
  !xmu = -15.d0/hart
  IF(ScreenI%xmu.LT.0.d0) xmu = ScreenI%xmu
  IF(ScreenI%xmu.GT.1.d-1) FreeElectronTest = .TRUE.

  ! Read these from input in future.
  !domega = 0.33333d0/hart
  nomega = ScreenI%nomega
  omega_max = ScreenI%omega_max
  domega = omega_max/(nomega-1)
  ne2    = ScreenI%ne2
  ilast0 = ScreenI%nrptx0
  IF(ilast0.LT.0) ilast0=MIN(INT((LOG(-ilast0*rnrm(0))+8.8d0)/dx)+1,nrx)
  IF(ilast0.GT.nrx) THEN
     IF(master) PRINT*, 'nrptx0 > nrx, resetting'
     ilast0 = nrx
  END IF
  IF(ScreenI%maxl.LT.lmaxph(0)+1) THEN
     IF(master) PRINT*, 'maxl < lmaxph+1, resetting to lmaxph+1' 
     ScreenI%maxl = lmaxph(0) + 1
  END IF
  IF(master) THEN
     PRINT*, 'ilast=', ilast0
     PRINT*, 'ne2=', ne2
     PRINT*, 'nomega=', nomega
     PRINT*, 'omega_max=', omega_max*hart
     PRINT*, 'gam1=', ScreenI%gam1*hart
     PRINT*, 'maxl=', ScreenI%maxl
     PRINT*, 'lmaxph(0)', lmaxph(0)
     PRINT*, 'xmu=', xmu*hart
     PRINT*, 'rmax = ', ri(ilast0)
     PRINT*, 'Number of processors:', numprocs
     PRINT*, 'Free Electron Test?', FreeElectronTest
  END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate omega independent quantities.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Free electron gas - rs = 4.0
  IF(FreeElectronTest) THEN
     rs = ScreenI%xmu*hart
     IF(master) PRINT*, 'rs = ', rs
     kf = fa/rs
     xmu = kf**2/2.d0
     IF(master) PRINT*, 'xmu =', xmu
     ScreenI%emin = -2*xmu 
  END IF

  !PRINT*, 'DEBUG: 1' 
  ! Set omega independent energy grids
  call setegi(xmu+ScreenI%emin, xmu, ScreenI%eimax, ScreenI%gam1, ScreenI%ner, ScreenI%nei, em, ne)
  !PRINT*, em(1)*hart, em(20)*hart, em(ne-20)*hart, em(ne)*hart

  !     transform to single precision
  srat=real(rat)

  
  iph = 0 
  !PRINT*, 'DEBUG: 2' 
  call yprep(iph, nat, inclus, nph, iphat, rfms2, srat)

  ! Calculate phase shifts 
  do iph = 0, nph
     IF(iph.EQ.0) THEN ! JK Fixed ihole bug 10/2010
        iholetmp = ihole
     ELSE
        iholetmp = 0
     END IF
     !print *, '     get phase shift for potential type ', iph
     lmax(iph) = lx

     !PRINT*, 'DEBUG: 3' 
     call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),vint,&
          &           rhoint,dx0,dx,jumprm,vjump,ri,vtotph(1,iph),rhoph,dmagx)
     call fixdsx (iph, dx0, dx, dgc, dpc, dgcn(1,1,iph),dpcn(1,1,iph))

     jri     = getiat(x0, dx, rmt(iph)) + 1
     jri1    = jri + 1
     eref(1) = vtotph(jri1,iph)
     do i = 1, jri1
        vtotph(i,iph) = vtotph(i,iph) - eref(1)
     end do

     if (ixc .ge. 5) then
        do i = 1, jri1
           vvalph(i,iph) = vvalph(i,iph) - eref(1)
        end do
     else
        do i = 1, jri1
           vvalph(i,iph) = vtotph(i,iph)
        end do
     endif

     !PRINT*, 'DEBUG: 4' 
     call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne, em, ixc,       &
          &             vtotph(1,iph), vvalph(1,iph),                        &
          &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
          &             adgc(1,1,iph), adpc(1,1,iph),                        &
          &             iz(iph), xion(iph), iunf,                            &
          &             iholetmp, MAXVAL(lmaxph), xnval(1,iph), ph1(1,1,iph), iph) !KJ iph
     !PRINT*, 'DEBUG: 5' 
  end do
  !PRINT*, 'DEBUG: 6' 
  dgcx(:) = dgcn(:,ihole,0)
  dpcx(:) = dpcn(:,ihole,0)
  rhoc = 0.d0
  DO i=1, nrptx
     rhoc(i) = (dgcx(i)**2 + dpcx(i)**2)*dx*ri(i)
  END DO
  IF(ilast0.GT.nrptx) PRINT*, 'WARNING: ilast > nrptx in prepw2'
  ! Calculate vch
  ! definition of vch and K can be moved outside of frequency loop.
  vch(:) = 0.d0
  wscrn(:) = 0.d0

  do i = 1, ilast0
     vch(i)  = rhoc(i)
  end do
  do i = 1, ilast0
     wscrn(i) = 0.0d0
     do j = 1, i
        wscrn(i) = wscrn(i) + vch(j)
     end do
     wscrn(i) = wscrn(i)/ri(i)
     do j = i + 1, ilast0
        wscrn(i) = wscrn(i) + vch(j)/ri(j)
     end do
  end do
  do i = 1, ilast0
     vch(i) = wscrn(i)
  end do
  ! Calculate beta_SC(1)
  beta_SC(1) = 0.d0
  do i = 1, ilast0
     beta_SC(1) = beta_SC(1) + rhoc(i)*vch(i)
  end do

  ! Calculate K_00 (for now K = 4pi/r>)
  Kmat(:,:) = 0.d0
  do i = 1, ilast0
     do j = i, ilast0
        Kmat(i,j) =  4.0d0*pi * 1/ri(j)
        Kmat(j,i) = Kmat(i,j)
     end do
  end do
  !KMat = 0.d0! Test

  !=================================================================
  !     LDA Exchange Correlation (actually this probably doesn't work)
  !=================================================================  
  if (ScreenI%lfxc .gt. 0) then
     call ldafxc(ilast, ri, rhoph, ScreenI%lfxc, fxc)
     call wlog('Using TDLDA kernel.')
     do i = 1, ilast
        Kmat(i,i) = Kmat(i,i) + 4.0d0*pi * fxc(i)
     end do
  end if

  ! Perform FMS for these first 2 grids and get gtrl for each
  call getgtrl(em, ne, eref(1), ph1, lmaxph, 0, nph, lfms2, inclus, rdirec, toler1, toler2, iprint, gtrl1) 
  !PRINT*, 'gtrl', gtrl1(0,1)
  !STOP
  
  ! Loop over omega
  iomega0 = 0
  DO iomega = 0, nomega
  IF(MOD(iomega,numprocs) .NE. this_process) CYCLE
  omega = domega*iomega
  if(master) PRINT*, 'Calculating beta(omega) for omega = ', omega, 'iproc = ', this_process
  if(master) PRINT*, 'iomega, nomega = ', iomega, nomega
  !ilast0 = ScreenI%nrptx0

  IF(FreeElectronTest) THEN
     deltaE = omega/(ne2-1)
     emin2 = xmu - omega
  ELSE
     deltaE = omega/(ne2-1)
     emin2 = xmu - omega
  END IF
  DO ie = 1, ne2
     em2(ie) = emin2 + (ie-1)*deltaE + ScreenI%gam1*coni
  END DO
  em3(1:ne) = em(1:ne) - omega
  em4(1:ne2) = em2(1:ne2) + omega
  ! Calculate phase shifts for frequency dependent grid.
  do iph = 0, nph
     IF(iph.EQ.0) THEN ! JK Fixed ihole bug 10/2010
        iholetmp = ihole
     ELSE
        iholetmp = 0
     END IF
     !print *, '     get phase shift for potential type ', iph
     lmax(iph) = lx

     call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),vint,&
          &           rhoint,dx0,dx,jumprm,vjump,ri,vtotph(1,iph),rhoph,dmagx)
     call fixdsx (iph, dx0, dx, dgc, dpc, dgcn(1,1,iph),dpcn(1,1,iph))

     jri     = getiat(x0, dx, rmt(iph)) + 1
     jri1    = jri + 1
     eref(1) = vtotph(jri1,iph)
     do i = 1, jri1
        vtotph(i,iph) = vtotph(i,iph) - eref(1)
     end do

     if (ixc .ge. 5) then
        do i = 1, jri1
           vvalph(i,iph) = vvalph(i,iph) - eref(1)
        end do
     else
        do i = 1, jri1
           vvalph(i,iph) = vtotph(i,iph)
        end do
     endif

     call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne2, em2, ixc,       &
          &             vtotph(1,iph), vvalph(1,iph),                        &
          &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
          &             adgc(1,1,iph), adpc(1,1,iph),                        &
          &             iz(iph), xion(iph), iunf,                            &
          &             iholetmp, MAXVAL(lmaxph), xnval(1,iph), ph2(1,1,iph), iph) !KJ iph
     call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne, em3, ixc,       &
          &             vtotph(1,iph), vvalph(1,iph),                        &
          &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
          &             adgc(1,1,iph), adpc(1,1,iph),                        &
          &             iz(iph), xion(iph), iunf,                            &
          &             iholetmp, MAXVAL(lmaxph), xnval(1,iph), ph3(1,1,iph), iph) !KJ iph
     call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne2, em4, ixc,       &
          &             vtotph(1,iph), vvalph(1,iph),                        &
          &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
          &             adgc(1,1,iph), adpc(1,1,iph),                        &
          &             iz(iph), xion(iph), iunf,                            &
          &             iholetmp, MAXVAL(lmaxph), xnval(1,iph), ph4(1,1,iph), iph) !KJ iph
  end do
  call getgtrl(em2, ne2, eref(1), ph2, lmaxph, 0, nph, lfms2, inclus, rdirec, toler1, toler2, iprint, gtrl2) 
  call getgtrl(em3, ne, eref(1), ph3, lmaxph, 0, nph, lfms2, inclus, rdirec, toler1, toler2, iprint, gtrl3) 
  call getgtrl(em4, ne2, eref(1), ph4, lmaxph, 0, nph, lfms2, inclus, rdirec, toler1, toler2, iprint, gtrl4) 

  ! Calculate chi01 and chi02 separately.
  ixc0 = 0
  IF(master) PRINT*, 'Calculating Chi0_1.'
  call calculate_chi01( inclus, master, FreeElectronTest, nrx,                              &
     &                   ScreenI%maxl, lmaxph(0), ScreenI%irrh, ScreenI%iend, ScreenI%nrptx0,  &
     &                   ihole, xmu, adgc, adpc, iz(0),             &
     &                   xion(0), iunf, xnval,                      &
     &                   ri, dx, x0, rmt(0), rnrm(0), em, ne, omega, ixc0, &
     &                   iph, vtotph, vvalph, eref(1), dgcn, dpcn,    &
     &                   ilast0, gtrl1, gtrl3, chi0_1)

  !PRINT*, 'Before calc2'
  IF(master) PRINT*, 'Calculating Chi0_2.'
  call calculate_chi02( inclus, master, FreeElectronTest, nrx,                                &
     &                   ScreenI%maxl, lmaxph(0), ScreenI%irrh, ScreenI%iend, ScreenI%nrptx0,  &
     &                   ihole, xmu, adgc, adpc, iz(0),             &
     &                   xion(0), iunf, xnval,                      &
     &                   ri, dx, x0, rmt(0), rnrm(0), em2, ne2, omega, ixc0, &
     &                   iph, vtotph, vvalph, eref(1), dgcn, dpcn,    &
     &                   ilast0, gtrl2, gtrl4, chi0_2)
   
  DO i = 1, ilast0
     DO j = 1, ilast0
        chi0(i,j) = (DIMAG(chi0_1(i,j)) + chi0_2(i,j))*EXP(-0.d0*abs(ri(i)-ri(j))/40.d0)
        chi0(j,i) = chi0(i,j)
        chi0_1(j,i) = chi0_1(i,j)
        chi0_2(j,i) = chi0_2(i,j)
     END DO
  END DO
  
  ! Calculate beta0 
  beta0(:) = 0.d0
  do i = 1, ilast0
     do j = 1, ilast0
        beta0(1) = beta0(1) + vch(i)*chi0_1(i,j)*vch(j)
        beta0(2) = beta0(2) + vch(i)*chi0(i,j)*vch(j)
     end do
  end do
  beta0 = beta0*4
 

  ! Calculate 1 - K*chi0
  do i = 1, ilast0
     do j = 1, ilast0
        if (i .eq. j) then
           Amat(i,j) = 1.0d0
        else
           Amat(i,j) = 0.0d0
        end if
        do k = 1, ilast0
           Amat(i,j) = Amat(i,j)-Kmat(i,k)*chi0(k,j)
        end do
     end do
     Qmatm1(i) = 1.d0/Amat(i,i)
  end do
  Amat2 = 0.d0
  do i = 1, ilast0
     do j = 1, ilast0
        do k = 1, ilast0
           Amat2(i,j) = Amat2(i,j)+Kmat(i,k)*chi0(k,j)
        end do
     end do
  end do

  ! Matrix inverse to get w_ch
  IF(master) call wlog('Computing (1 - K Chi0)^-1 v_ch')

  IF(.FALSE.) GOTO 105 ! Skip
  IF(DiagonalizeEps) THEN
     Amat2 = Amat 
     DO i = 1, ilast0
        Amat2(i,i) = Amat(i,i) - 1.d0
     END DO
     !IF(.FALSE.) THEN
        !Amat2(1:4,1:4) = TRANSPOSE(reshape((/ ( 5.0, 9.0), ( 5.0,  5.0), (-6.0, -6.0), (-7.0, -7.0), &
        !        &                    ( 3.0, 3.0), ( 6.0, 10.0), (-5.0, -5.0), (-6.0, -6.0), &
        !        &                    ( 2.0, 2.0), ( 3.0,  3.0), (-1.0,  3.0), (-5.0, -5.0), &
        !        &                    ( 1.0, 1.0), ( 2.0,  2.0), (-3.0, -3.0), ( 0.0,  4.0) /), shape(Amat2(1:4,1:4))))
        !Amat = Amat2
     !END IF

     ! Diagonalize 1-Kchi0 to get fluctuation potentials. 
     CALL ZGEEV( 'V', 'V', ilast0, AMat2, nrx, EigVals, EigVecsL, nrx, EigVecsR, nrx, &
     &                  WORK, 2*nrx, RWORK, info )
     IF(info.GT.0) STOP 'QR failed to converge in ZGEEV.'
     ! Normalize eigenvectors CONJG(L_i)*R_j = \delta_ij
     norm_max = 1.d10
     DO i = 1, ilast0
        DO ir = 1, ilast0
           norm = norm + CONJG(EigVecsL(ir,i))*EigVecsR(ir,i)
        END DO      
        norm_max = MIN(ABS(norm_max),ABS(norm))
     END DO
     PRINT*, 'norm_max=', norm_max
     DO i = 1, ilast0
        DO j = 1, ilast0
           norm = 0.d0
           DO ir = 1, ilast0
              norm = norm + CONJG(EigVecsL(ir,i))*EigVecsR(ir,j)
           END DO
           !WRITE(222,*) omega, EigVecsR(i,j), EigVecsL(i,j)
           IF(i.NE.j) THEN
              IF(ABS(norm)/ABS(norm_max).GT.1.d-3) THEN
                 PRINT*, 'norm', i, j, ABS(norm), ABS(norm/norm_max)
                 STOP 'Eigenvectors not orthogonal from zgeev.'
              END IF
           ELSE
              EigVecsR(1:ilast0,i) = EigVecsR(1:ilast0,i)/norm
           END IF
        END DO
     END DO
     EigVals = EigVals + 1.d0
     ! Check if we can reproduce original matrix.
     !DO i = 1, ilast0
     !   DO j = 1, ilast0
     !      Amat2(i,j) = 0.d0
     !      DO ir = 1, ilast0
     !         Amat2(i,j) = Amat2(i,j) + CONJG(EigVecsL(j,ir))*(EigVals(ir))*EigVecsR(i,ir)
     !      END DO
           !WRITE(221,*) i,j,omega, Amat2(i,j), Amat(i,j)
     !   END DO
        !WRITE(223,*) omega, EigVals(i)
     !END DO
  END IF
  !IF(omega.GT.0.d0) STOP
  nrhs = 1
  !     size(Amat,1)
  lda  = nrx
  !     size(wscrn)
  ldb  = nrx
  !     actual size of nxn matrix is ilast
  wscrn = vch
  call zgetrf(ilast0, ilast0, Amat, lda, ipiv, info)

  if (info .eq. 0) then
     call zgetrs(trans, ilast0, nrhs, Amat, lda, ipiv, wscrn, ldb, info)

     !if (info .ne. 0) print *, info1
  else
     call wlog('The factor U is singular')
     call par_stop('Stopping now.')
  end if

  ! Iterate to self-consistency (integral equation for w).
  wscrn_SC(1,:) = vch(:)
  ! Calculate beta, beta_tilde for first iteration.
  beta_tilde_SC(2:mxsc) = 0.d0
  ifirst = 0
  DO isc = 2, mxsc-1
     wscrn_SC(mxsc,:) = 0.d0
     ! Calculate beta_tilde(
     DO i = 1, ilast0
        wscrn_SC(mxsc,i) = vch(i)
        DO j = 1, ilast0
           wscrn_SC(mxsc,i) = wscrn_SC(mxsc,i) + Amat2(i,j)*wscrn_SC(isc-1,j) 
        END DO
     END DO
     !DO i = 1, ilast0
     !   wscrn_SC(isc,i) = vch(i)
     !   DO j = 1, ilast0
     !      wscrn_SC(isc,i) = wscrn_SC(isc,i) + Amat2(i,j)*wscrn_SC(mxsc,j) 
     !   END DO
     !END DO
           
     ! wscrn_SC(:,mxsc) now holds wscrn_{i-1}^{out}, while wscrn_SC(:,isc-1) holds wscrn_{i-1}^{in}, and wscrn_SC(:,isc) 
     ! Calculate beta_in
     !beta_SC(isc-1:isc+1) = 0.d0
     !DO i = 1, ilast0
     !   beta_SC(isc-1) = beta_SC(isc-1) + rhoc(i)*wscrn_SC(i,isc-1)
     !   beta_SC(isc) = beta_SC(isc) + rhoc(i)*wscrn_SC(i,mxsc)
     !   beta_SC(isc+1) = beta_SC(isc+1) + rhoc(i)*wscrn_SC(i,isc)
     !END DO
     ! wscrn_SC(:,mxsc) now holds wscrn_{i-1}^{out}, while wscrn_SC(:,isc-1) holds wscrn_{i-1}^{in}. Calculate optimal alpha.
     !alpha = 0.0d0
     !DO i = 1, ilast0
     !   alpha0 = (wscrn_SC(isc,i) - 2.d0*wscrn_SC(mxsc,i) - wscrn_SC(isc-1,i))
     !   IF(ABS(wscrn_SC(mxsc,i)).GT.1.d-4) alpha0 = alpha0/wscrn_SC(mxsc,i)
     !   IF(ABS(alpha0).GT.1.d-2.AND.ifirst.EQ.0) THEN
     !      alpha0 = (wscrn_SC(isc-1,i) - wscrn_SC(mxsc,i))/(wscrn_SC(isc,i) - 2.d0*wscrn_SC(mxsc,i) - wscrn_SC(isc-1,i))
     !   ELSE
     !      alpha0 = 0.01d0
     !      !ifirst = 1
     !   END IF 
     !   alpha = alpha + alpha0
     !END DO
     !alpha = alpha/ilast0
     !IF(ABS(alpha).GT.1.d0) alpha = alpha/ABS(alpha)
     !alpha = 0.01d0
     !IF(ABS(alpha).LT.0.001d0) alpha = 0.001d0
     !PRINT*, omega, isc, DBLE(alpha), DIMAG(alpha)
     !IF(ABS(beta_SC(isc+1)-beta_SC(isc)).GT.1.d-8.AND.ABS(beta_SC(isc)-beta_SC(isc-1)).GT.1.d-8) THEN 
     !   alpha = -(beta_SC(isc) - beta_SC(isc-1))/(beta_SC(isc+1) - beta_SC(isc))
     !ELSE
     !   alpha = 0.2d0
     !END IF
     alpha = 0.001d0
     !IF(ABS(wscrn_SC(mxsc,1) - wscrn_SC(isc-1,1)).GT.0.1d0) THEN
     !       alpha = 0.1d0/(wscrn_SC(mxsc,1) - wscrn_SC(isc-1,1))
     !ELSE
     !       alpha = 0.1d0
     !END IF
     ! Use Aitken's method to update
     wscrn_SC(isc,:) = alpha*wscrn_SC(mxsc,:) + (1.d0-alpha)*wscrn_SC(isc-1,:)
     IF(isc.GT.2) THEN
        DO i = 1, ilast0
           IF(ABS(wscrn_SC(isc,i) - 2.d0*wscrn_SC(isc-1,i) + wscrn_SC(isc-2,i)).GT.1.d-3.AND.ABS(wscrn_SC(isc,i)-wscrn_SC(isc-1,i)).LT.0.1d0) THEN
              wscrn_SC(isc,i) = wscrn_SC(isc,i) - (wscrn_SC(isc,i) - wscrn_SC(isc-1,i))**2/(wscrn_SC(isc,i) - 2.d0*wscrn_SC(isc-1,i) + wscrn_SC(isc-2,i))
           END IF
        END DO
     END IF
     beta_SC(1) = 0.d0
     do i = 1, ilast0
        beta_SC(1) = beta_SC(1) + rhoc(i)*wscrn_SC(isc,i)
     end do
     !PRINT*, this_process, isc
     !IF(master) WRITE(135,'(1e20.10, I4, 20e20.10)') omega, isc, DBLE(alpha), DIMAG(alpha), DBLE(beta_SC(isc-1)), DIMAG(beta_SC(isc-1)), DBLE(beta_SC(isc)), DIMAG(beta_SC(isc)), DBLE(beta_SC(isc+1)), DIMAG(beta_SC(isc+1))
  END DO
  !IF(master) WRITE(135,*)
105 CONTINUE
  IF(iomega0.EQ.0) THEN
     !IF(master) THEN
        wscrn0 = wscrn
        ! Send wscrn0 to other processes
        CALL par_bcast_dc(wscrn0, nrx, 0) 
     !ELSE
     !   CALL par_recv_dc(wscrn0, nrx, 0, 0)
     !END IF
  END IF
  ifirst = 0
  beta = 0.d0
  beta2 = 0.d0
  beta3 = 0.d0
  beta4 = 0.d0
  !beta0 = 0.d0
  do i = 1, ilast0
     do j = 1, ilast0
        beta2 = beta2 + wscrn(i)*chi0(i,j)*CONJG(wscrn(j))  !*dx**2*ri(i)*ri(j)
        beta4 = beta4 + wscrn0(i)*chi0(i,j)*CONJG(wscrn0(j))
        beta = beta + vch(i)*chi0(i,j)*wscrn(j)  !*dx**2*ri(i)*ri(j)
        !beta0 = beta0 + vch(i)*chi0(i,j)*vch(j)  !*dx**2*ri(i)*ri(j)
        !beta_SC(isc) = beta_SC(isc) + vch(i)*chi0(i,j)*vch(j)
     end do
     beta3 = beta3 + (dpcx(i)**2 + dgcx(i)**2)*wscrn(i)*ri(i)*dx
  end do
  beta = beta*4
  beta2 = beta2*4
  beta3 = beta3/pi
  beta4 = beta4*4
  beta_SC(1) = beta_SC(1)*4
  IF(DiagonalizeEps) THEN
     beta_sum = 0.d0
     DO i = 1, ilast0
        beta_decomp(i) = 0.d0
        betaL = 0.d0
        betaR = 0.d0
        DO ir = 1, ilast0
           betaL = betaL + (dpcx(ir)**2 + dgcx(ir)**2)*EigVecsR(ir,i)*ri(ir)*dx
           betaR = betaR + vch(ir)*CONJG(EigVecsL(ir,i))
        END DO

        beta_decomp(i) = betaL*betaR/EigVals(i)/pi
        beta_sum = beta_sum + beta_decomp(i)
     END DO
     PRINT*, 'iomega0 = ', iomega0
     !IF(iomega.EQ.38) THEN
        ! Integrate rho_c*  fluctuation potential V_s
        betaL = 0.d0
        betaR = 0.d0
        DO ir = 1, ilast0
           betaL = betaL + (dpcx(ir)**2 + dgcx(ir)**2)*EigVecsR(ir,1)*ri(ir)*dx
           betaR = betaR + vch(ir)*CONJG(EigVecsL(ir,1))
           WRITE(22,'(20e20.10)') ri(i), DBLE(EigVecsR(i,1)), DIMAG(EigVecsR(i,1)), DBLE(EigVecsL(i,1)), DIMAG(EigVecsL(i,1))
        END DO
        beta_Fluct = betaR*betaL
        PRINT*, 'omega = ', omega
        PRINT*, 'Fluctuation amplitude of beta for rutile: ', beta_Fluct
        PRINT*, 'Plasmon frequency of beta for rutile: 0.535929 - 0.031 I'
     !END IF

  END IF
  iomega0 = iomega0 + 1

  IF(master) THEN
     do i = 1, ilast0
        WRITE(19,'(5e20.10)') ri(i), omega, DBLE(wscrn(i)), DIMAG(wscrn(i)), vch(i)
     end do
     WRITE(19,*)
     FLUSH(19)

     WRITE(18,'(6e20.10)') omega, DBLE(beta), DIMAG(beta), DIMAG(beta2), DIMAG(beta3), DIMAG(beta4)
     FLUSH(18)
     WRITE(20,'(40e20.10)') omega, DBLE(beta_SC(1)), DIMAG(beta_SC(1))
     FLUSH(20)
     WRITE(17,'(5e20.10)') omega, DBLE(beta0(1)), DIMAG(beta0(1)), DBLE(beta0(2)), DIMAG(beta0(2))
     FLUSH(17)
     IF(DiagonalizeEps) THEN
        !WRITE(22,'(20e20.10)') ri(i), omega, DBLE(EigVecsR(i,1)), DIMAG(EigVecsR(i,1)), DBLE(EigVecsL(i,1)), DIMAG(EigVecsL(i,1))

        DO i = 1, ilast0
           IF(omega.EQ.0.d0) THEN
              WRITE(flname, '(A5,I4.4)') 'beta_', i
              OPEN(i+23,FILE=flname,STATUS='REPLACE')
           END IF
           WRITE(21,'(1e20.10, I5, 2e20.10)') omega, i, DBLE(EigVals(i)), DIMAG(EigVals(i))
           DO j=1, ilast0
           END DO
           WRITE(i+23,'(5e20.10)') omega, DBLE(beta_decomp(i)), DIMAG(beta_decomp(i)), DBLE(beta_decomp(1)*EigVals(1)/EigVals(i)), DIMAG(beta_decomp(1)*EigVals(1)/EigVals(i))
           FLUSH(i+23)
        END DO
        WRITE(23,*) omega, DBLE(beta_sum), DIMAG(beta_sum)
        WRITE(21,*)
        FLUSH(21)
     END IF
     DO iproc = 1, numprocs-1
        IF(iomega0.GT.nomega) EXIT
        CALL par_recv_dc(beta, 1, iproc, iproc)
        CALL par_recv_dc(beta2, 1, iproc, iproc)
        CALL par_recv_dc(beta3, 1, iproc, iproc)
        CALL par_recv_dc(beta4, 1, iproc, iproc)
        CALL par_recv_dc(beta0, 2, iproc, iproc)
        CALL par_recv_dc(beta_SC(1), 1, iproc, iproc)
        CALL par_recv_dc(omega, 1, iproc, iproc)
        CALL par_recv_dc(wscrn, nrx, iproc, iproc)
        IF(DiagonalizeEps) THEN
           CALL par_recv_dc(EigVals, nrx, iproc, iproc)
           CALL par_recv_dc(beta_decomp, nrx, iproc, iproc)
           CALL par_recv_dc(beta_sum, 1, iproc, iproc)
           !CALL par_recv_dc(EigVecsR(1,1), nrx**2, iproc, iproc)
           !CALL par_recv_dc(EigVecsL(1,1), nrx**2, iproc, iproc)
           DO i = 1, ilast0
              WRITE(i+23,'(5e20.10)') omega, DBLE(beta_decomp(i)), DIMAG(beta_decomp(i)), DBLE(beta_decomp(1)*EigVals(1)/EigVals(i)), DIMAG(beta_decomp(1)*EigVals(1)/EigVals(i))
              FLUSH(i+23)
           END DO
           WRITE(23,*) omega, DBLE(beta_sum), DIMAG(beta_sum)
        END IF
           
        do i = 1, ilast0
           WRITE(19,'(5e20.10)') ri(i), omega, DBLE(wscrn(i)), DIMAG(wscrn(i)), vch(i)
        end do
        WRITE(19,*)
        FLUSH(19)

        WRITE(18,'(6e20.10)') omega, DBLE(beta), DIMAG(beta), DIMAG(beta2), DIMAG(beta3), DIMAG(beta4)
        FLUSH(18)
        WRITE(20,'(40e20.10)') omega, DBLE(beta_SC(1)), DIMAG(beta_SC(1))
        FLUSH(20)
        WRITE(17,'(5e20.10)') omega, DBLE(beta0(1)), DIMAG(beta0(1)), DBLE(beta0(2)), DIMAG(beta0(2))
        FLUSH(17)
        IF(DiagonalizeEps) THEN
           DO i = 1, ilast0
              WRITE(21,'(1e20.10, I5, 2e20.10)') omega, i, DBLE(EigVals(i)), DIMAG(EigVals(i))
           END DO
           WRITE(21,*)
           FLUSH(21)
        END IF
        iomega0 = iomega0 + 1
     END DO
  ELSE
     call par_send_dc(beta, 1, 0, this_process)
     CALL par_send_dc(beta2, 1, 0, this_process)
     CALL par_send_dc(beta3, 1, 0, this_process)
     CALL par_send_dc(beta4, 1, 0, this_process)
     CALL par_send_dc(beta0, 2, 0, this_process)
     CALL par_send_dc(beta_SC(1), 1, 0, this_process)
     CALL par_send_dc(omega, 1, 0, this_process)
     CALL par_send_dc(wscrn, nrx, 0, this_process)
     IF(DiagonalizeEps) THEN
        CALL par_send_dc(EigVals, nrx, 0, this_process)
        CALL par_send_dc(beta_decomp, nrx, 0, this_process)
        CALL par_send_dc(beta_sum, 1, 0, this_process)
        !CALL par_send_dc(EigVecsR(1,1), nrx**2, 0, this_process)
        !CALL par_send_dc(EigVecsL(1,1), nrx**2, 0, this_process)
     END IF
  END IF
100 CONTINUE
  END DO ! end of loop over omega
  call par_barrier
  deallocate(xnmues)
  ilast = ilast0

  return
end subroutine prepw2
!=======================================================================
!     END PREP
!=======================================================================

subroutine getgtrl(em, ne, eref, ph, lmaxph, iph, nph, lfms2, inclus, rdirec, toler1, toler2, ipr, gtrl) 
  use DimsMod, only: lx, nrptx, natx, nphx=>nphu, nspx=>nspu, ltot, nex
  use constants
  IMPLICIT NONE
  ! Input
  COMPLEX*16 em(nex), eref ! energy grid, energy reference
  complex*16 ph(nex, ltot+1, 0:nphx) ! phase shifts
  integer lmaxph(0:nphx) ! max l for each unique potential
  INTEGER nph, lfms2, inclus, ne, iph, ipr ! ipr prints message to screen
  REAL rdirec,toler1,toler2

  ! Output
  COMPLEX*16 gtrl(0:lx, nex)
  COMPLEX*16 ggmmp(0:lx,nex,2*lx+1,2*lx+1)

  ! Local
  COMPLEX*16 p2, ck
  INTEGER iverb, minv, ispin, nsp 
  complex, allocatable :: gg(:,:,:), gtr(:,:,:), xphase(:,:,:)
  logical, allocatable :: lcalc(:)

  !Iterators
  INTEGER i, j, ill, ipp, il, im, imp, ix, ie
  REAL rpart, aipart
  COMPLEX cks(nspx)

  ! Allocate variables
  allocate(gg(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx))
  allocate(gtr(0:lx, 0:nphx, nex))
  allocate(xphase(nspx, -lx:lx, 0:nphx))
  allocate(lcalc(0:lx))
  if (inclus .le. 1) return
  ! Loop over energy grid.
  DO ie = 1, ne
     
     p2 = em(ie) - eref
     ck   = sqrt(2*p2 + (p2*alphfs)**2)
     rpart  = real( dble(ck))
     aipart = real(dimag(ck))
     cks(1) = cmplx(rpart, aipart)
        
     !=============================================================
     !         FMS
     !=============================================================
     do ipp = 0, nph
        do ill = -lmaxph(ipp), lmaxph(ipp)
           rpart  = dble( ph(ie, abs(ill)+1, ipp))
           aipart = dimag(ph(ie, abs(ill)+1, ipp))
           xphase(1, ill, ipp) = cmplx(rpart, aipart)
        end do
     end do
     
     iverb = 0
     nsp = 1
     ispin = 0
     do i = 0, lx
        lcalc(i) = .true.
     end do
     minv = 0
     
     do i = 0, lx
        do j = 0, nphx
           gtr(i,j,ie) = 0.0d0
        end do
     end do

     !WRITE(115,*) 'lmaxph=',lmaxph
     !WRITE(115,*) 'xphase=',xphase
     !WRITE(115,*) 'ie=',ie
     !WRITE(115,*) iverb,minv,rdirec,toler1,toler2
     !WRITE(115,*) lcalc
     !STOP
     IF(ipr.EQ.1.and.MOD(ie,10).EQ.0) PRINT*, 'FMS for energy point', ie, 'of', ne
     call fms(lfms2, nsp, ispin, inclus, nph, cks,lmaxph, &
          & xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc, gg)

     do il = 0, lmaxph(iph)
        ix = il**2
        do im = 1, 2*il+1
           gtr(il,iph,ie) = gtr(il,iph,ie) + gg(ix+im,ix+im,iph)
           do imp = 1, 2*il+1
              ggmmp(il, ie, im, imp) = gg(ix+im, ix+imp, iph) 
           end do
        end do
        
        gtrl(il,ie) = ( dble(real( gtr(il,iph,ie))) &
             & + coni*dble(aimag(gtr(il,iph,ie))) )&
             & * exp(2*coni*ph(ie,il+1,iph))/(2*il+1.0d0)
     end do
     !PRINT*, gtrl(0,1)
     !STOP
  end do 
end subroutine getgtrl

!subroutine get_fluct_pot(nrx, VecL,VecR,lambda,omega,V_s,omega_s)
!  INTEGER nrx
!  COMPLEX*16 VecL(nrx,nrx,2), VecR(nrx,nrx,2) ! EigenVectors of epsilon
!  COMPLEX*16 lambda(nrx,2) ! Eigenvalues corresponding to above eigenvectors
!  DOUBLE PRECISION omega ! energy at which epsilon was diagonalized 
  
