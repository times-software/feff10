! This program will calculate the RIXS spectrum.
! Need the following ingredients, which for now will be read from file.
! 1 and 2 denote different core-hole potentials.
! rho1, rhoc1, rho2, rhoc2
! rl1, rl2
! vch1, vch2
! rkk1
! xsnorm1
! rkk and xsnorm will be interpolated onto the energy grid defined by rho.
PROGRAM RIXS
  USE dimsmod, only:init_dimensions,nex, nspx=>nspu, nphx=>nphu, lx, ltot, nrptx
  use errorfile
  USE IOMOD
  USE par
  use constants
  use fms_inp,only: minv,rclust=>rfms2,lmaxph,ipr3 !KJ this and following from fmstot
  use global_inp,only: ipol,ispin,le2,angks,ptz,do_nrixs
  use atoms_inp,only: nat,iphat,ratdbl=>rat
  use nrixs_inp,only: jmax,kfinmax,jinit
  use eels_inp,only: ipmin,ipmax,ipstep,eels
  use rixs_inp
  use AtomicPotIO
  use xmuio 
  IMPLICIT NONE

  INTEGER ne, ne1, ne3, nph, ihole, ik0, lmaxp1
  REAL(8) rnrmav, xmu, edge, emu, emu2
  complex*16, allocatable  :: em(:), eref(:,:)
  integer NTotal, iTotal
  integer, allocatable :: iz(:)
  integer, allocatable :: lmax(:,:)
  character*6, allocatable  :: potlbl(:)
  complex*16, allocatable :: ph(:,:,:,:)
  complex*16, allocatable :: rkk(:,:,:)
  COMPLEX*16, allocatable :: gg(:,:,:,:), gg2(:,:,:,:)
  !REAL(8), allocatable :: rho(:,:), rhoc(:,:)

  INTEGER ne_2, ne1_2, ne3_2, nph_2, ihole_2, ik0_2, lmaxp1_2, itmp(10), ihl(10)
  REAL(8) rnrmav_2, xmu_2, edge_2
  complex*16, allocatable  :: em_2(:), eref_2(:,:)
  integer, allocatable :: iz_2(:)
  integer, allocatable :: lmax_2(:,:)
  character*6, allocatable  :: potlbl_2(:)
  complex*16, allocatable :: ph_2(:,:,:,:)
  complex*16, allocatable :: rkk_2(:,:,:)
  !REAL(8), allocatable :: rho_2(:,:), rhoc_2(:,:)

  complex*16, allocatable :: bmat(:,:,:,:,:,:,:) !KJ added last index
  complex*16, allocatable :: bmat0(:,:,:,:,:,:)  !KJ new variable 1-06      
  complex*16, allocatable :: TLb(:,:,:,:), totTLb(:)
  complex*16, allocatable :: xsect_rxs(:,:,:,:), xsect_tmp(:,:,:)
  complex*16 xInt(2)
  complex*16, allocatable :: Rl(:,:,:), Rl_2(:,:,:), Sigma(:)
  complex*16, allocatable :: ctmp(:), k2(:), k1(:)
  complex*16 phx, rInt(nrptx)

  real*8 vch1(nrptx), vch2(nrptx), DeltaV(nrptx), ri(nrptx), rmt, dtmp
  real*8 wall_start, wall_end, gam_ch, w, w0, wmin,wmax, dw, w1, &
       & w2, width, rem(nex), remtmp(nex), xsnorm(nex), &
       & dx1, dx2, x0_1, x0_2, tlRe(nex), tlIm(nex), xtmp(nex), xasEI(nex,2), xasEF(nex,2), &
       & xes(nex,2), wXES(10000), deltaMin, delta
  integer :: ios, iE1, iE2, ikap, kap1, kinit, linit, iw, ip, ind, is1, imu, jri1, jri2
  integer :: knd(8), lnd(8), lll, ir, lllmax, llltmp, imt, jmt, iPercent, m1, m2, iEdge, nEdge, i1, i2, L1, L2, iph, ix, nSig, nw
  character(300) fname
  logical ltrace, ReadSigma, SkipCalc, MBConv, ReadPoles
  COMPLEX*16 F0, F1, bmatavg
  real(8) gam_exp(2), Edge1(2), EdgeSplit(10), EdgeAmp(10), gam_tot, gam_edge(10), EMin(2), EMax(2), EEMin, EEMax, ET, &
       & EE, EIMin, EIMax, EI, muXES(10000)
  character c
  CHARACTER(512) message
  CHARACTER(LEN=30) rlfiles(2), ggfiles(2), phasefiles(2), wscrnfiles(2), xsectfile
  COMPLEX*16 KKInt
  REAL(8), PARAMETER :: xInfinity = 1.d30
  REAL(8) IntDoubleLorentz
  EXTERNAL KKInt, IntDoubleLorentz

  call par_begin
  if (worker) go to 400
  if (master) call OpenErrorfileAtLaunch('rixs')
  MBConv = .TRUE.
  ! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0

  ReadSigma = .TRUE.
!  SkipCalc = .FALSE.
  CALL init_dimensions

  allocate(iz(0:nphx))

  CALL reafms

  CALL rixs_read
  CALL rixs_set(gam_ch, gam_exp, EMin, EMax, xmu, ReadPoles, SkipCalc, MBConv, ReadSigma)
  !PRINT*, 'lx, ipmin, ipmax'
  !PRINT*, lx, ipmin, ipmax
  rlfiles(1) = 'rl_1.dat' !TRIM(ADJUSTL(RixsI%Edges(1))) // '/rl.dat'
  rlfiles(2) = 'rl_2.dat' !TRIM(ADJUSTL(RixsI%Edges(2))) // '/rl.dat'
  ggfiles(1) = 'gg_1.bin' !TRIM(ADJUSTL(RixsI%Edges(1))) // '/gg.bin'
  ggfiles(2) = 'gg_2.bin' !TRIM(ADJUSTL(RixsI%Edges(2))) // '/gg.bin'
  phasefiles(1) = 'phase_1.bin' !TRIM(ADJUSTL(RixsI%Edges(1))) // '/phase.bin'
  phasefiles(2) = 'phase_2.bin' !TRIM(ADJUSTL(RixsI%Edges(2))) // '/phase.bin'
  !PRINT*, phasefiles(1)
  wscrnfiles(1) = 'wscrn_1.dat' !TRIM(ADJUSTL(RixsI%Edges(1))) // '/wscrn.dat'
  wscrnfiles(2) = 'wscrn_2.dat' !TRIM(ADJUSTL(RixsI%Edges(2))) // '/wscrn.dat'

  CALL ReadData(rlfiles(1), Double1 = rmt, Int2 = lllmax, Int3 = jri1)
  CALL ReadData(rlfiles(1), Double1 = dx1, Double2 = x0_1)

  CALL ReadData(rlfiles(2), Double1 = rmt, Int2 = lllmax, Int3 = jri2)
  CALL ReadData(rlfiles(2), Double1 = dx2, Double2 = x0_2)

  ! Initialize clock
  call seconds(wall_start)
  wall_comm = 0.0

  ! Open the log file, unit 11.  See subroutine wlog.
  if (master) then
     open (unit=11, file='logrixs.dat', status='unknown', iostat=ios)
     call chopen (ios, 'logrixs.dat', 'rixs')
  else
     par_type = 2
  endif


  ! Allocate local variables
  allocate(em(nex), eref(nex, nspx))
  allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(rkk(nex,8,nspx))
  allocate(lmax(nex, 0:nphx))
  allocate(potlbl(0:nphx))
  !allocate(rho(nex,4),rhoc(nex,4))

  allocate(ph_2(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(lmax_2(nex, 0:nphx))
  allocate(em_2(nex), eref_2(nex, nspx))
  allocate(potlbl_2(0:nphx))
  allocate(rkk_2(nex,8,nspx))
  !allocate(rho_2(nex,lx+1),rhoc_2(nex,lx+1))
  allocate(iz_2(0:nphx))

  ! Read phase.bin first to get ne to minimize memory usage
  ! phase.bin for both core holes (get rkk).
  call rdxsphrxs(ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,         &
            ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1, phasefiles(1))
  call rdxsphrxs(ne_2, ne1_2, ne3_2, nph_2, ihole_2, rnrmav_2, xmu_2, edge_2,         &
            ik0_2, em_2, eref_2, iz_2, potlbl_2, ph_2, rkk_2, lmax_2, lmaxp1_2, phasefiles(2))
  allocate(Rl(jri1,-lx:lx,ne), Rl_2(jri2,-lx:lx,ne))
  allocate(bmat(-lx:lx,0:1,8, -lx:lx,0:1,8,ipmin:ipmax))
  allocate(bmat0(-lx:lx,0:1,8,-lx:lx,0:1,8))
  allocate(gg(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx,ne),gg2(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx,ne))
  allocate(TLb(ne,ne,nspx*(lx+1)**2,8),totTLb(ne))
  allocate(xsect_rxs(ne,ne,2,10), xsect_tmp(ne,ne,2))
  allocate(ctmp(ne),k2(ne),k1(ne))
  allocate(Sigma(ne))
  Sigma(:) = 0.d0
  
  ! Read data
  !PRINT*, ''
  CALL wlog(' ')
  CALL wlog('Reading data.')

  ! ! phase.bin for both core holes (get rkk).
  ! call rdxsphrxs(ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,         &
  !           ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1, phasefiles(1))
  ! call rdxsphrxs(ne_2, ne1_2, ne3_2, nph_2, ihole_2, rnrmav_2, xmu_2, edge_2,         &
  !           ik0_2, em_2, eref_2, iz_2, potlbl_2, ph_2, rkk_2, lmax_2, lmaxp1_2, phasefiles(2))
  rem(1:ne1) = DBLE(em(1:ne1))
!  PRINT*, 'Enter mu:'
!  READ*, xmu
!  xmu = xmu/hart
  !PRINT*, 'Read edge energies from poles.dat? 0) no. 1) yes.'
  EdgeSplit(:) = 0.d0
  !CALL ReadArrayData('edge.dat', Int1 = itmp(1:10), Int2 = ihl(1:10), Double3 = EdgeSplit(1:10), NumElements = nEdge)
!  gam_ch = gam_Edge(1)
!  gam_Edge(1:9) = gam_edge(2:10)
!  READ*, EdgeSplit(1)
  IF(ReadPoles) THEN
     CALL wlog('Reading data from edges.dat')
     CALL ReadArrayData('edges.dat', Double1 = EdgeSplit, Double2 = EdgeAmp, Double3 = gam_Edge, NumElements = nEdge)
     CALL CloseFl('edges.dat')
     IF(EdgeAmp(nEdge).LT.0.d0) EdgeAmp(nEdge) = 1
     Edge1(1) = EdgeSplit(1)
     EdgeSplit(1:9) = EdgeSplit(2:10)
     gam_ch = gam_Edge(1) 
     gam_Edge(1:9) = gam_edge(2:10)
     EdgeAmp(1:9) = EdgeAmp(2:10)
     nEdge = nEdge - 1
     DO iEdge = 1, INT(nEdge/2)
        xtmp(1) = EdgeSplit(iEdge)
        EdgeSplit(iEdge) = EdgeSplit(nEdge - iEdge + 1)
        EdgeSplit(nEdge-iEdge + 1) = xtmp(1) 
        xtmp(1) = gam_Edge(iEdge) 
        gam_Edge(iEdge) = gam_Edge(nEdge - iEdge + 1)
        gam_Edge(nEdge - iEdge + 1) = xtmp(1)
        xtmp(1) = EdgeAmp(iEdge) 
        EdgeAmp(iEdge) = EdgeAmp(nEdge - iEdge + 1)
        EdgeAmp(nEdge - iEdge + 1) = xtmp(1)
     END DO
     DO iEdge = 1, nEdge
        IF(EdgeSplit(iEdge).LE.0.d0) EdgeSplit(iEdge) = -xmu
     END DO
     Edge1(2) = EdgeSplit(1)
     Edge1(1:nEdge) = Edge1(1:nEdge) - xmu
     EdgeSplit(:) = EdgeSplit(:) - Edge1(2)
     WRITE(message,'(a,I3)') 'Number of poles: ', nEdge
     DO iEdge = 1, nEdge
       CALL wlog('Pole data: iEdge, Edge energy, Edge split, Amp, Gamma')
       WRITE(message,'(I3,5F10.5)') iEdge, Edge1(iEdge), EdgeSplit(iEdge), EdgeAmp(iEdge), gam_edge(iEdge)
       CALL wlog(message)
     END DO
  ELSE
     nEdge = 1
     EdgeSplit(1) = 0.d0
     EdgeAmp(1) = 1.d0
  END IF
  !STOP 
  !READ*, EdgeSplit(2), EdgeAmp(1), gam_Edge(1), EdgeAmp(2), gam_Edge(2)
  !EdgeSplit(1) = 0.d0
  !EdgeSplit(:) = EdgeSplit(:)/hart
  !gam_Edge(:) = gam_Edge(:)/hart
  !Edge1 = Edge1/hart
  !EdgeAmp(1) = 4.d0
  !EdgeAmp(2) = 2.d0
  !IF(EdgeSplit(2).GT.0.d0) nEdge = 2
  IF(RixsI%xmu.GT.-1.d2) xmu = RixsI%xmu
  WRITE(message,'(a,f10.5)') 'xmu = ', xmu
  CALL wlog(message)
  IF(SkipCalc) GOTO 80
  ! Read gg
  iph = 0
  CALL wlog('Reading gg_1.bin')
  DO iE1 = 1, ne1
     CALL Read2D(ggfiles(1), gg(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,iph,iE1), L1, L2)
     ! Check that bounds are correct.
     IF(L1.ne.nspx*(lmaxph(iph)+1)**2.or.L2.ne.nspx*(lmaxph(iph)+1)**2) CALL Error('Error when reading gg.bin')
  END DO
  CALL CloseFl(ggfiles(1))
  CALL wlog('Reading gg_2.bin')
  DO iE1 = 1, ne1
     CALL Read2D(ggfiles(2), gg2(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,iph,iE1), L1, L2)
     ! Check that bounds are correct.
     IF(L1.ne.nspx*(lmaxph(iph)+1)**2.or.L2.ne.nspx*(lmaxph(iph)+1)**2) CALL Error('Error when reading gg.bin')
  END DO
  CALL CloseFl(ggfiles(2))
  CALL wlog('Finished reading gg files.')

  ! Transform to real spherical harmonics and take imaginary parts
  !CALL ToRealYlm(gg,lmaxph(iph),nspx,iph,ne1)
  !CALL ToRealYlm(gg2,lmaxph(iph),nspx,iph,ne1)
  IF(.FALSE.) THEN ! Test gg tranformation
     OPEN(UNIT = 14,FILE = 'rhoLL.dat', STATUS = 'REPLACE')
     DO iE1 = 1, ne1
        WRITE(14, '(100E20.10)') DBLE(em(iE1)), &
         ((DIMAG(gg(L1**2+L1+m1+1,L1**2+L1+m1+1,iph,iE1)* &
          EXP(2.d0*coni*ph(iE1, L1, 1, 0))), m1 = -L1, L1), L1 = 0, lx)
     END DO
     CLOSE(14)
  END IF
              
  ! xsect.dat
  xsectfile = 'xsect_2.dat' !TRIM(ADJUSTL(RixsI%Edges(2))) // '/xsect.dat'
  CALL ReadArrayData(xsectfile, Double1 = remtmp, Double2 = remtmp, Double3 = xsnorm)
  remtmp(:) = DBLE(em(:))

  ! core hole potentials
  CALL ReadArrayData(wscrnfiles(1), Double1 = DeltaV, Double2 = vch1)
  IF(TRIM(ADJUSTL(RixsI%Edges(2))).NE.'VAL') THEN
     CALL ReadArrayData(wscrnfiles(2), Double1 = DeltaV, Double2 = vch2)
  ELSE
     vch2 = 0.d0
  END IF
  DeltaV = -vch1 + vch2
  IF(.FALSE.) DeltaV = 1.d0
   
  ! Make radial grid.
  !dx = 0.05d0
  DO ir = 1, nrptx
     ri(ir) = exp(-x0_1+dx1*(ir-1))
  END DO
  imt = (log(rmt) + x0_1) / dx1  +  1
  jmt = imt+1
  !vch1(:) = -10.d0
  !vch2(:) = -5.d0
  
  ! Enter energy loop. E2 is outer energy related to energy loss. E1 is inner energy related to incoming photon.
  k2(1:ne) = 2.d0*(rem(1:ne) - DBLE(eref_2(1,nspx)))
  DO iE2 = 1, ne1
     k2(iE2) = SQRT(k2(iE2))
  END DO
  DO iE1 = 1, ne1
     k1(iE1) = SQRT(2.d0*(rem(iE1)-eref(1,nspx)))
  END DO
  
  call setkap (ihole, kinit, linit)

  !KJ 1-06  I added the do-loop around the call to bcoef for ELNES calculation.
  ltrace = .FALSE.
  do ip=ipmin,ipmax,ipstep
     if (eels.eq.1) call iniptz(ptz,ip,2)  !KJ Only change ptz for ELNES
     call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, knd, lnd, bmat0)
     bmat(:,:,:,:,:,:,ip)= bmat0(:,:,:,:,:,:)
  enddo
  WRITE(66,'(e20.10)') bmat
  ! Read in Rl and Rl2
  CALL wlog(' ')
  CALL wlog('Reading Rl_1.')
  DO iE1 = 1, ne1
     DO lll = 0, lllmax
        CALL ReadData(rlfiles(1), DComplex1 = em(iE1), Int2 = llltmp, DComplex3 = phx)
        CALL ReadArrayData(rlfiles(1), DComplex1 = Rl(:,llltmp,iE1))
        !ph(iE1,lll,1,0) = phx
     END DO
     !IF(MOD(iE1,30).EQ.0) PRINT '(I3,A4,I3)', iE1, ' of ', ne1
  END DO
  CALL CloseFl(rlfiles(1))
  CALL wlog(' ')
  CALL wlog('Reading Rl_2.')
  DO iE2 = 1, ne1
     DO lll = 0, lllmax
        CALL ReadData(rlfiles(2), DComplex1 = em_2(iE2), Int2 = llltmp, DComplex3 = phx)
        CALL ReadArrayData(rlfiles(2), DComplex1 = Rl_2(:,llltmp,iE2))
        !ph_2(iE2,lll,1,0) = phx
     END DO
     !IF(MOD(iE2,30).EQ.0) PRINT '(I3,A4,I3)', iE2, ' of ', ne1
  END DO  
  CALL CloseFl(rlfiles(2))
  ! Read self-energy from mpse.dat
  remtmp(:) = -1.d10
  IF(ReadSigma) THEN
     CALL wlog('Reading mpse.dat')
     CALL ReadArrayData('mpse.dat', Double1 = remtmp, DComplex2 = Sigma)
     CALL CloseFl('mpse.dat')
     remtmp(:) = remtmp(:)/hart
     Sigma(:) = Sigma(:)/hart
     DO iE1 = 1, ne
        IF(remtmp(iE1).LT.-1.d5) EXIT
     END DO
     nSig = iE1 - 1
     ctmp(:) = -1.d20
     DO iE1 = 1, ne1
        IF((rem(iE1).GE.xmu).AND.(rem(iE1).LE.remtmp(nSig))) THEN
           CALL terpc(remtmp,Sigma,nSig,1, rem(iE1)-xmu, ctmp(iE1))
        ELSEIF((rem(iE1).LT.xmu)) THEN
           ctmp(iE1) = Sigma(1)
        ELSE
           ctmp(iE1) = Sigma(nSig)
        END IF
     END DO
     Sigma(:) = ctmp(:)
     CALL wlog('Done reading mpse.dat')
  END IF
  !xmu = -14.d0/hart
  ! Outer loop over L (kappa)
  CALL wlog('Forming T.')
  Tlb(:,:,:,:) = 0.d0
  DO ind = 1, 8
     ix = lnd(ind)**2 + lnd(ind) + 1
     IF(lnd(ind).GE.0) THEN
        DO is1 = 1, nspx
           DO iE2 = 1, ne1
              ! For each E2, calculate TLb
              DO iE1 = 1, ne1
                 IF((rem(iE1).GE.xmu).AND.rem(iE2).GE.xmu) THEN
                    totTlb(1:) = 0.d0
                    rInt(1:jmt) = Rl(1:jmt,lnd(ind),iE1)*DeltaV(1:jmt)*Rl_2(1:jmt,lnd(ind),iE2)
                    CALL csomm2(ri,rInt,dx1, totTlb(1), rmt, jmt)
                    !rInt(1:251) = Rl(1:251,lnd(ind),iE1)*Rl_2(1:251,lnd(ind),iE2)                  
                    !CALL csomm2(ri,rInt,dx, totTlb(2), rmt, jmt)
                    IF(ind.gt.3) THEN
                       WRITE(17,'(6f20.10)') rem(iE1), rem(iE2), &
                         DBLE(totTlb(1)), DIMAG(totTlb(1)), DBLE(totTlb(2)), &
                         DIMAG(totTlb(2))
                    END IF
                    ! Add matrix elements.
                    DO m1 = -lnd(ind), lnd(ind)
                       ! Average bmat
                       !
                       Tlb(iE1,iE2,ix + m1,ind)  = DBLE(totTLb(1))*ABS(rkk(iE1,ind,is1) *EXP(-coni*ph(iE1, -lnd(ind), 1, 0)))* & 
                            & (1.d0 + DIMAG(gg(ix + m1,ix + m1,iph,iE1)*EXP(2.d0*coni*ph(iE1, -lnd(ind), 1, 0))))* &
                            & SQRT(xsnorm(iE1))
                       !Tlb(iE1,iE2,ix + m1,ind) = (Tlb(iE1,iE2,ix + m1,ind) - DBLE(rkk(iE2,ind,is1)*EXP(-coni*ph_2(iE2, lnd(ind), 1, 0)))*SQRT(xsnorm(iE2)))*sqrt(ABS(bmat(m1,is1-1,ind,m1,is1-1,ind,ipmin)))
                    END DO
                    !Tlb(iE1,iE2,ind) = 1.d0
                 END IF
              END DO
              IF(ind.GT.3) WRITE(17,*)
           END DO
        END DO
     END IF
  END DO
  !IF(.TRUE.) ne1 = 50
    !Test
  IF(.FALSE.) THEN
     ne1 = 200
     xmu = 0.d0
     DO iE1 = 1, ne1
        rem(iE1) = 0.1d0*DBLE(iE1) 
        Tlb(iE1,:,:,:) = 1.5d0/((rem(iE1) - 10.d0)**2 + (1.5d0)**2)
     END DO
     !Tlb(:24,:,:,:) = 0.d0
     gg2 = 0.d0
  END IF
  open(13, FILE = 'tlb0.dat', STATUS = 'REPLACE')
  DO iE2 = 1, ne1
     DO iE1 = 1, ne1
        WRITE(13,'(20f20.10)') DBLE(rem(iE1)), DBLE(rem(iE2)), &
         (DBLE(Tlb(iE1,iE2, lnd(ind)**2 + lnd(ind) + 1, ind)), &
         DIMAG(Tlb(iE1,iE2,lnd(ind)**2 + lnd(ind) + 1, ind)), ind = 1, 8)
     END DO
     WRITE(13,*)
  END DO
  CLOSE(13)
  IF(.FALSE.) STOP ! debugging
  ! Now do convolution to get edge right.
  CALL wlog('Performing convolution of T.')
  NTotal = 0
  DO ind = 1, 8
        IF(lnd(ind).GE.0) NTotal = NTotal + nspx*ne1**2*(ne1-1)*(2*lnd(ind)+1)
  END DO
  iTotal = 1
  iPercent = 0

  DO ind = 1, 8
     ix = lnd(ind)**2 + lnd(ind) + 1
     IF(lnd(ind).GE.0) THEN
        DO is1 = 1, nspx
           DO m1 = -lnd(ind), lnd(ind)              
              DO iE2 = 1, ne1         
              ! DO iE2 = 1, 1
                 ctmp(1) = 0.d0
                 IF(rem(iE2).GE.xmu) THEN
                    ctmp(1) =  ABS(rkk_2(iE2,ind,is1)*EXP(-coni*ph_2(iE2, -lnd(ind), 1, 0)))*SQRT(xsnorm(iE2))
                 END IF
                 DO iE1 = 1, ne1
                    totTLb(iE1) = 0.d0
                    DO iw = 1, ne1 - 1
                       gam_tot = gam_ch !+ ABS(DIMAG(Sigma(iE1)))
!KJ gives KKInt errors when xmu=rem(iw+1)                       IF(rem(iw+1).GE.xmu) THEN 
                       IF(rem(iw+1).GT.xmu) THEN 
!                           print*,'iE2,iE1,iw,rem(iw+1),rem(iw),xmu',iE2,iE1,iw,rem(iw+1),rem(iw),xmu
                          IF(rem(iw).LT.xmu) THEN
                             F1 = Tlb(iw + 2, iE2, ix + m1, ind) - Tlb(iw + 1, iE2, ix + m1, ind)
                             F1 = F1/(rem(iw+2) - rem(iw+1))
                             F0 = Tlb(iw + 2, iE2, ix + m1, ind) - F1*rem(iw+2)
                             totTlb(iE1) = totTlb(iE1) + KKInt(F1,F0,xmu,rem(iw+1),gam_tot,rem(iE1)*1.00001d0)
                          ELSE
                             F1 = Tlb(iw + 1, iE2, ix + m1, ind) - Tlb(iw, iE2, ix + m1, ind)
                             F1 = F1/(rem(iw+1) - rem(iw))
                             F0 = Tlb(iw + 1, iE2, ix + m1, ind) - F1*rem(iw+1)
                             totTlb(iE1) = totTlb(iE1) + KKInt(F1,F0,rem(iw),rem(iw+1),gam_tot,rem(iE1)*1.00001d0)
                          END IF
                       END IF
                       iTotal = iTotal + 1
                       IF((DBLE(iTotal)/DBLE(NTotal)).GT.0.1d0) THEN
                          iTotal = 1
                          iPercent = iPercent+10
                          WRITE(message,'(I4,a)') iPercent, '%'
                          CALL wlog(message)
                       END IF
                    END DO
                    ! Add last point. Model Tlb(E) = a/E**2 for point beyond last.                                                                                            
                    F0 = Tlb(ne1,iE2,ix+m1,ind)*rem(ne1)**2
                    totTlb(iE1) = totTlb(iE1) + (F0*((0,-1)*gam_tot + rem(iE1) +   &
                         &   rem(ne1)*Log(((0,1)*gam_tot - rem(iE1) + rem(ne1))/rem(ne1))))/  &
                         &   ((gam_tot + (0,1)*rem(iE1))**2*rem(ne1))
                    
                    totTlb(iE1) = (-totTlb(iE1)/pi*sqrt(k2(iE2)/pi) + ctmp(1))
                    ! Test with larger tlb
                    !totTlb(iE1) = (-1000*totTlb(iE1)/pi*sqrt(k2(iE2)/pi) + ctmp(1))
                   
                    ! Test with no deltav, i.e., T = d.
                    IF(.FALSE.) totTlb(iE1) = ctmp(1)
                 END DO 
                 !TLb(:,iE2,ix + m1, ind) = -coni*DBLE(totTLb(:)) + DIMAG(totTLb)/pi
                 TLb(:,iE2,ix + m1, ind) = - totTlb(:)
              END DO
              Tlb(:,:,ix + m1, ind) = Tlb(:,:,ix + m1, ind)*sqrt(ABS(bmat(m1,is1-1,ind,m1,is1-1,ind,ipmin)))
           END DO
        END DO
     END IF
  END DO
10 CONTINUE
  open(13, FILE = 'tlb1.dat', STATUS = 'REPLACE')
  DO iE2 = 1, ne1
  !DO iE2 = 1, 1
     DO iE1 = 1, ne1
        WRITE(13,'(20f20.10)') DBLE(rem(iE1)), DBLE(rem(iE2)), &
          (DBLE(Tlb(iE1,iE2, lnd(ind)**2 + lnd(ind) + 1, ind)), &
          DIMAG(Tlb(iE1,iE2,lnd(ind)**2 + lnd(ind) + 1, ind)), ind = 1, 8)
     END DO
     WRITE(13,*)
  END DO
  CLOSE(13)
  xsect_rxs(:,:,:,1) = 0.d0
  DO iE1 = 1, ne1
     DO iE2 = 1, ne1
        DO ind = 1, 8
           ix = lnd(ind)**2 + lnd(ind) + 1
           DO is1 = 1, nspx
              IF(lnd(ind).GE.0) THEN
                 IF(ind.LE.3) THEN
                    DO m1 = -lnd(ind), lnd(ind)
                       DO m2 = m1, m1 !-lnd(ind), lnd(ind)
                          xsect_rxs(iE2,iE1,1,1) = xsect_rxs(iE2,iE1,1,1) +  &
                            & ABS(TLb(iE1,iE2,ix + m2, ind))**2 *&
                            & (1.d0 + DIMAG(gg2(ix + m1,ix + m2,iph,iE2)*EXP(2.d0*coni*ph_2(iE2, -lnd(ind), 1, 0))))
                       END DO
                    END DO
                 ELSE
                    DO m1 = -lnd(ind), lnd(ind)
                       xsect_rxs(iE2,iE1,2,1) = xsect_rxs(iE2,iE1,2,1) + ABS(TLb(iE1,iE2,ix + m1, ind))**2 *&
                            & (1.d0 + DIMAG(gg2(ix + m1,ix + m1,iph,iE2)*EXP(2.d0*coni*ph_2(iE2, -lnd(ind), 1, 0))))
                    END DO
                 END IF
              END IF
           END DO
        END DO
     END DO
  END DO
  open(14, FILE = 'rixs0.dat', STATUS = 'REPLACE')
  open(15, FILE = 'xas0.dat', STATUS = 'REPLACE')
  !GOTO 35
  DO iE1 = 1, ne1
     DO iE2 = 1, ne1
        WRITE(14,'(20E30.15)') rem(iE2), rem(iE1), xsect_rxs(iE1,iE2,1,1), xsect_rxs(iE1,iE2,2,1)
     END DO
     WRITE(15,'(20E30.15)') rem(iE1), xsect_rxs(iE1,iE1,1,1), xsect_rxs(iE1,iE1,2,1)
     WRITE(14,*)
  END DO
  CLOSE(14)
  CLOSE(15)

  IF(.FALSE.) STOP
  ! Debug
  IF(.FALSE.) THEN
     DO iE1 = 1, 100
        DO iE2 = 1, 100
           xsect_rxs(iE2,iE1,:,:) = 1.d0 + 0.01d0/((rem(iE1)**2 + 0.01**2))
        END DO
     END DO
     xsect_rxs(:50,:,:,:) = 0.d0
     xsect_rxs(:,:50,:,:) = 0.d0
     ne1 = 100
  END IF
  CALL wlog('Performing convolution of xsect_rxs.')
  WRITE(message,'(I3,a)') nEdge, ' poles.'
  CALL wlog(message)
  !rem(:) = DBLE(em(:)-eref(:,nspx))
  ! Now perform the integral for broadening.
  NTotal = ne1*ne1*2
  iTotal = 1
  DO iEdge = 2, nEdge
     xsect_rxs(:,:,:,iEdge) = xsect_rxs(:,:,:,1)
  END DO
  DO iEdge = 1, nEdge
  WRITE(message,'(a,I3,3f10.5)') 'Edge: ', iEdge, EdgeSplit(iEdge), EdgeAmp(iEdge), gam_Edge(iEdge)
  CALL wlog(message)

  ! Now perform convolution over incident energy
  xsect_tmp(:,:,:) = 0.d0
  NTotal = ne1*ne1*2
  iTotal = 0
  iPercent = 0
  WRITE(message,'(a,2f10.5)') 'gam_exp = ', gam_exp(1), gam_exp(2)
  CALL wlog('Convolving over incident energy.')
  deltaMin = 100.d0
  DO iE1 = 2, ne1
     IF((rem(iE1) - rem(iE1-1)).LT.deltaMin) deltaMin = rem(iE1) - rem(iE1-1)
  END DO
  ctmp(1) = 0.d0
  DO ind = 1, 2
     DO iE1 = 1, ne1
        IF(DBLE(iTotal)/DBLE(NTotal)*100.d0.GT.10.d0) THEN
           iPercent = iPercent + 10
           iTotal = 0
           WRITE(message,'(I3,a1)') iPercent, '%'
           CALL wlog(message)
        END IF        
        !WRITE(33,*) rem(iE1), ABS(DIMAG(Sigma(iE1)))
        DO iE2 = 1, ne1
           iTotal = iTotal + 1
           totTLb(iE2) = 0.d0
           xInt(:) = 0.d0
           ctmp(:) = 0.d0
           !GOTO 15
           IF(.TRUE.) THEN ! Use interpolation
              !wmin = rem(1) - MIN(rem(iE1),rem(iE2))
              !wmax = rem(ne1) - MAX(rem(iE1),rem(iE2))
              ! wmin = xmu 
!               wmin = MAX(rem(1)+rem(iE2)-rem(iE1),wmin)
!               !wmin = rem(1) - rem(iE1)
!               !wmin = MAX(wmin,xmu-rem(iE2))
!               wmin = MAX(wmin,rem(iE2)-10.d0*(gam_tot+ABS(DIMAG(Sigma(iE1)))))
!               wmax = rem(ne1) !+ 10.d0*(gam_tot+ABS(DIMAG(Sigma(iE1))))
!               wmax = MIN(rem(ne1)-rem(iE1)+rem(iE2),wmax)
!               wmax = MIN(wmax,rem(iE2)+10.d0*(gam_tot+ABS(DIMAG(Sigma(iE1)))))
!               dw = (wmax - wmin)/400.d0
!               wmax = wmax - dw
              !IF(wmax.LT.wmin) CYCLE
              ! if delE > 0, greatest point given by w = Emax - delE, 
              ! least point given by w = Emin
              ! if delE < 0, greatest point given by w = Emax,
              ! least point given by w = Emin - delE
              IF((rem(iE1) - rem(iE2)).GT.0) THEN
                 xtmp(1) = rem(ne1) - rem(iE1) + rem(iE2)
                 xtmp(2) = rem(1)
              ELSE
                 xtmp(1) = rem(ne1)
                 xtmp(2) = rem(1) - rem(iE1) + rem(iE2)
              END IF
              w1 = MIN(xtmp(2),xtmp(1))                 
              wmin = w1
              wmax = MAX(xtmp(2),xtmp(1))

              IF(iE1.GT.1) THEN
                 xInt(1) = rem(iE1) - rem(iE1-1)
              ELSE
                 xInt(1) = rem(iE1+1) - rem(iE1)
              END IF
              IF(iE2.GT.1) THEN
                 xInt(2) = rem(iE2) - rem(iE2-1)
              ELSE
                 xInt(2) = rem(iE2+1) - rem(iE2)
              END IF
              delta = MIN(DBLE(xInt(1)),DBLE(xInt(2)))
              !dw = MIN(wmax - wmin)/nw,deltaMin)
              ! Try dw = deltaMin
              IF(.TRUE.) THEN
                 dw = MAX((delta + deltaMin)/4.d0,deltaMin)
                 nw = INT((wmax - wmin)/dw)
                 dw = (wmax - wmin)/DBLE(nw)
              ELSE
                 nw = INT(5.d0*(wmax - wmin)/gam_tot)
                 dw = (wmax - wmin)/DBLE(nw)
              END IF
              
              IF(wmin.GT.xmu) THEN
                 CALL terpc(rem,Sigma,ne1,1, wmin, ctmp(1))
              ELSE
                 ctmp(1) = Sigma(1)
              END IF        
              ctmp(1) = 0.d0
              CALL BLInterp2D(rem,rem,xsect_rxs(:,:,ind,iEdge),ne1,ne1,ne,w1, rem(iE1) + w1 - rem(iE2), xInt(1))
              IF(.FALSE.) THEN
                 xInt(1) = xInt(1)*(gam_tot+ABS(DIMAG(ctmp(1))))/ &
                ((wmin-rem(iE2)+DBLE(ctmp(1)))**2 + (gam_tot+ABS(DIMAG(ctmp(1))))**2)/pi
              END IF
              w0 = rem(iE2) - DBLE(ctmp(1))
              gam_tot = gam_exp(1) + ABS(DIMAG(ctmp(1)))
              DO iw = 1,  nw
                 w = wmin + dw*DBLE(iw - 1)
                 w1 = MIN(MAX(w+dw,xtmp(2)),xtmp(1)) 
                 !IF((rem(iE2)+w+dw).LT.xmu) CYCLE
                 !PRINT*, iE1, iE2 iw
                 IF(w.GE.xmu) THEN
                    CALL terpc(rem,Sigma,ne1,1,rem(iE1), ctmp(1))                   
                 ELSE
                    ctmp(1) = 0.d0
                 END IF
                 CALL terpc(rem,Sigma,ne1,1,rem(iE1) + w1 - rem(iE2), ctmp(2))
                 CALL terpc(rem,Sigma,ne1,1,rem(iE2), ctmp(2))
                 gam_tot = gam_exp(1) + ABS(DIMAG(ctmp(1)+ctmp(2)))/2.d0
                 ctmp(1) = (ctmp(1)+ctmp(2))/2.d0
                 CALL BLInterp2D(rem,rem,xsect_rxs(:,:,ind,iEdge),ne1,ne1,ne,w1, rem(iE1) + w1 - rem(iE2), xInt(2))
                 IF(.TRUE.) THEN
                    ! Use linearized xsect.
                    F1 = (xInt(2) - xInt(1))/dw
                    F0 = xInt(2) - F1*w1
                    xInt(1) = xInt(2)
                    xInt(2) = -2.d0*(F0 + F1*w0)*( ATAN((w-w0)/gam_tot) - ATAN((w+dw-w0)/gam_tot) )
                    xInt(2) = xInt(2) + F1*gam_tot*LOG( ((w+dw-w0)**2 + gam_tot**2)/((w-w0)**2 + gam_tot**2) )
                    xInt(2) = xInt(2)*0.5d0/pi
                    xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + xInt(2)
                 ELSE
                    xInt(2) = xInt(2)*gam_tot/((w+dw-rem(iE2)+DBLE(ctmp(1)))**2 + gam_tot**2)/pi
                    xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + 0.5d0*(xInt(2) + xInt(1))*dw
                    xInt(1) = xInt(2)
                 END IF
              END DO
              ! Add endpoint
              CALL BLInterp2D(rem,rem,xsect_rxs(:,:,ind,iEdge),ne1,ne1,ne,w1, rem(iE1) + wmax - rem(iE2), xInt(2))
              xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + xInt(2)*(0.5d0 - ATAN((wmax-rem(iE2)+DBLE(ctmp(1)))/gam_tot)/pi)
           ELSE
              DO iw = 1, ne1 - 1
                 xInt(1) = xInt(2)
                 w = rem(iw)
                 dw = rem(iw+1) - w
                 ! w2 = E2-E2+w
                 w2 = rem(iE2)-rem(iE1)+w
                 IF(w2.LT.rem(1)) THEN
                    xInt(2) = xsect_rxs(1,iw,ind,iEdge)
                 ELSE IF(w2.GT.rem(ne1)) THEN
                    xInt(2) = xsect_rxs(ne1,iw,ind,iEdge)
                 ELSE
                    CALL terpc(rem,xsect_rxs(:,iw,ind,iEdge),ne1,1,w2,xInt(2))
                 END IF
                 !xInt(2) = xInt(2)*gam_tot/pi*1.d0/((w2+DBLE(Sigma(iE1)))**2 + gam_tot**2)
                 !xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) +
                 ! 0.5d0*(xInt(1) + xInt(2))*dw
                 xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) - &
                      &  xInt(2)/pi*(ATAN2(rem(iE1) - w - dw            &
                      & -DBLE(Sigma(iE1)),gam_tot) - ATAN2(rem(iE1)-w       &
                      & -DBLE(Sigma(iE1)),gam_tot))
              END DO
15            CONTINUE
              ! Add end points. 
              ! Bottom of grid.
              !           w = rem(1)
              !           w2 = rem(iE2)-rem(iE1) + w
              !           IF(w2.LT.rem(1)) THEN
              !              xInt(2) = xsect_rxs(1,1,ind,iEdge)
              !           ELSE IF(w2.GT.rem(ne1)) THEN
              !              xInt(2) = xsect_rxs(1,ne1,ind,iEdge)
              !           ELSE
              !              CALL terpc(rem,xsect_rxs(:,1,ind,iEdge),ne1,1,w2,xInt(2))
              !           END IF           
              !           xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + &
              !                &  xInt(2)*(ATAN2(rem(iE1)-w            &
              !                & -DBLE(Sigma(iE1)),gam_tot)/pi + 1.d0/2.d0) 
              !
              !           ! Top of grid
              !           w = rem(ne1)
              !           w2 = rem(iE2)-rem(iE1) + w
              !           IF(w2.LT.rem(1)) THEN
              !              xInt(2) = xsect_rxs(1,ne1,ind,iEdge)
              !           ELSE IF(w2.GT.rem(ne1)) THEN
              !              xInt(2) = xsect_rxs(ne1,ne1,ind,iEdge)
              !           ELSE
              !              CALL terpc(rem,xsect_rxs(:,ne1,ind,iEdge),ne1,1,w2,xInt(2))
              !           END IF      
              !           xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + &
              !                &  xInt(2)*(1.d0/2.d0-ATAN2(rem(iE1)-w            &
              !                & -DBLE(Sigma(iE1)),gam_tot)/pi) 
           END IF
        END DO
     END DO
  END DO
  xsect_rxs(:,:,:,iEdge) = xsect_tmp(:,:,:)
  
  open(14, FILE = 'rixs1.dat', STATUS = 'REPLACE')
  open(15, FILE = 'xas1.dat', STATUS = 'REPLACE')
  !GOTO 35
  DO iE1 = 1, ne1
     DO iE2 = 1, ne1
        WRITE(14,'(20E30.15)') rem(iE2), rem(iE1), xsect_rxs(iE1,iE2,1,1), xsect_rxs(iE1,iE2,2,1)
     END DO
     WRITE(15,'(20E30.15)') rem(iE1), xsect_rxs(iE1,iE1,1,1), xsect_rxs(iE1,iE1,2,1)
     WRITE(14,*)
  END DO
  CLOSE(14)
  CLOSE(15)
  xtmp(1) = gam_ch !+ ABS(DIMAG(Sigma(iE1)))
  iPercent = 0
  DO ind = 1, 2
     DO iE1 = 1, ne1
        DO iE2 = 1, ne1
           gam_tot = gam_exp(2) + gam_Edge(iEdge) !+ ABS(DIMAG(Sigma(iE2)))
           IF(gam_tot.LT.0.01d0/hart) THEN
           !IF(.TRUE.) THEN
              ! Skip convolution in E2 direction. 
              totTlb(iE2) = xsect_rxs(iE2,iE1,ind,iEdge)/((rem(iE1)-rem(iE2))**2 + gam_ch**2)
              GOTO 55
           END IF

           totTLb(iE2) = 0.d0
           DO iw = 1, ne1-1
              ! Find linear interpolation
              xtmp(1) = (xsect_rxs(iw+1,iE1,ind,iEdge) - xsect_rxs(iw,iE1,ind,iEdge))/(rem(iw+1)-rem(iw))
              xtmp(2) = xsect_rxs(iw,iE1,ind,iEdge) - xtmp(1)*rem(iw)
              totTLb(iE2) = totTLb(iE2) + (IntDoubleLorentz(rem(iE1),rem(iE2),gam_ch,gam_tot,xtmp(2),xtmp(1),rem(iw+1)*1.001d0,1) &
                 & - IntDoubleLorentz(rem(iE1),rem(iE2),gam_ch,gam_tot,xtmp(2),xtmp(1),rem(iw)*1.001d0,1)) 
           END DO
           ! end point correction
           xtmp(1) = xsect_rxs(ne1,iE1,ind,iEdge)
           xtmp(2) = 0.d0
           totTLb(iE2) =  totTLb(iE2) + (IntDoubleLorentz(rem(iE1),rem(iE2),gam_ch,gam_tot,xtmp(1),xtmp(2),xInfinity,-1) - &
              &  IntDoubleLorentz(rem(iE1),rem(iE2),gam_ch,gam_tot,xtmp(1),xtmp(2),rem(ne1)*1.001d0,1)) 
           GOTO 55

! Note by FDV
! According to Solaris Studio the statements below can not be reached.

           ! width set integral based on width*gam
           width = 10.d0
           ! nw set number of points per gam_min
           nw = 10

           w0 = MIN(rem(iE1)-width*xtmp(1),rem(iE2)-width*gam_tot)
           wmax = MAX(rem(iE1)+width*xtmp(1),rem(iE2)+width*gam_tot)
           dw = MIN(gam_tot,xtmp(1))/DBLE(nw)
           nw = FLOOR((wmax-w0)/dw)
           xInt(:) = 0.d0
           DO iw = 1, nw + 1
              xInt(1) = xInt(2)
              w = w0 + DBLE(iw-1)*dw
              !w = rem(iE2) - rem(ne1) - 10.0*gam_exp(2) + (iw-1)*dw
              IF(w.GT.rem(ne1)) THEN
                 xInt(2) = xsect_rxs(ne1,iE1,ind,iEdge)
                 !ctmp(1) = Sigma(ne1)
              ELSE IF(w.LT.rem(1)) THEN
                 xInt(2) = xsect_rxs(1,iE1,ind,iEdge)
                 !ctmp(1) = Sigma(1)
              ELSE
                 CALL terpc(rem,xsect_rxs(:,iE1,ind,iEdge),ne1,1,w,xInt(2)) !- DBLE(Sigma(iE2)),xInt(2))
              END IF
              !xInt(2) = (gam_exp(2)+ABS(DIMAG(Sigma(iE2))))/pi*xInt(2)*1.d0/(w**2 + (gam_exp(2)+ABS(DIMAG(Sigma(iE2))))**2)/((rem(iE1)-DBLE(Sigma(iE1))-rem(iE2)+DBLE(Sigma(iE2))+w)**2 + (gam_ch+ABS(DIMAG(Sigma(iE1))))**2)
              !xInt(2) = (gam_tot + ABS(DIMAG(Sigma(iE1)-ctmp(1))))/pi*xInt(2)*1.d0/(w**2 + gam_tot**2)/((rem(iE1)-DBLE(Sigma(iE1))-rem(iE2)+DBLE(ctmp(1))+w)**2 + (gam_ch+ABS(DIMAG(Sigma(iE1)-ctmp(1))))**2)
              !xInt(2) = gam_tot/pi*xInt(2)*1.d0/(w**2 + gam_tot**2)/((rem(iE1)+ DBLE(Sigma(iE2) - DBLE(Sigma(iE1) - rem(iE2)+w)**2 + (gam_ch+DBLE()**2) !+ABS(DIMAG(Sigma(iE1))))**2)
              xInt(2) = xInt(2) * 1.d0/((rem(iE1) - w)**2 + xtmp(1)**2) * gam_tot/pi/((rem(iE2) - w)**2 + gam_tot**2)
              totTLb(iE2) = totTLb(iE2) + 0.5d0*(xInt(2) + xInt(1))*dw
           END DO
55         CONTINUE
           iTotal = iTotal + 1
           IF(DBLE(iTotal)/DBLE(NTotal).GT.0.1d0) THEN
              iTotal = 1
              iPercent = iPercent+10
              WRITE(message,'(I3,a)') iPercent, '%'
              CALL wlog(message)
           END IF
        END DO
        xsect_rxs(:,iE1,ind,iEdge) = totTLb(:)*EdgeAmp(iEdge)
     END DO
  END DO
  xsect_tmp(:,:,:) = xsect_rxs(:,:,:,iEdge)
  !GOTO 25
                 
              
              
  
35 CONTINUE
  END DO

  ! Now add all edges if more than one.
  xsect_tmp(:,:,:) = 0.d0
  IF(nEdge.GT.0) THEN
     DO ind = 1, 2
        DO iE1 = 1, ne1
           DO iEdge = 1, nEdge   
           DO iE2 = 1, ne1
              IF((rem(iE2)-EdgeSplit(iEdge)).GT.rem(1)) THEN
                 CALL terpc(rem,xsect_rxs(:,iE1,ind,iEdge),ne1,1,rem(iE2) - EdgeSplit(iEdge),xInt(1))
              ELSE
                 xInt(1) = xsect_rxs(1,iE1,ind,iEdge)
              END IF
              xsect_tmp(iE2,iE1,ind) = xsect_tmp(iE2,iE1,ind) + xInt(1)
           END DO
           END DO
        END DO
     END DO
  END IF
80 CONTINUE
  IF(SkipCalc) THEN
     ! Read data from rixsET.dat.
     open(14, FILE = 'rixsET.dat', STATUS = 'OLD')
     DO iE1 = 1, ne1
        DO iE2 = 1, ne1
           READ(14,*,END = 85) rem(iE2), remtmp(iE1), xtmp(1), xtmp(2), xtmp(3), xtmp(4)
           xsect_tmp(iE1,iE2,1) = xtmp(1)
           xsect_tmp(iE1,iE2,2) = xtmp(3)
        END DO
        READ(14,*,END = 85)
     END DO
85   CONTINUE
     ne1 = iE1 - 1
     !Edge1 = Edge1 - rem(ne1)
     !Edge1 = Edge1/hart
     rem = rem/hart - Edge1(1)
     CLOSE(14)
  END IF
  CALL wlog('Writing results.')
  ! Integrate over E1 and E2 to get line spectra.
  xasEI = 0.d0
  xasEF = 0.d0
  xes = 0.d0
  DO iE1 = 1, ne1
     DO iE2 = 2, ne1
        IF(rem(iE2).GT.EMax(1)) EXIT
        IF(rem(iE2).GT.EMin(1)) THEN
           xasEI(iE1,1) = xasEI(iE1,1) + 0.5d0*(xsect_tmp(iE2, iE1,1) + xsect_tmp(iE2-1,iE1,1))*(rem(iE2)-rem(iE2-1))
           xasEI(iE1,2) = xasEI(iE1,2) + 0.5d0*(xsect_tmp(iE2, iE1,2) + xsect_tmp(iE2-1,iE1,2))*(rem(iE2)-rem(iE2-1))
        END IF
     END DO
  END DO
  DO iE1 = 1, ne1
     DO iE2 = 2, ne1
        IF(rem(iE2).GT.EMax(2)) EXIT
        IF(rem(iE2).GT.EMin(2)) THEN
           xasEF(iE1,1) = xasEF(iE1,1) + 0.5d0*(xsect_tmp(iE1, iE2,1) + xsect_tmp(iE1,iE2-1,1))*(rem(iE2)-rem(iE2-1))
           xasEF(iE1,2) = xasEF(iE1,2) + 0.5d0*(xsect_tmp(iE1, iE2,2) + xsect_tmp(iE1,iE2-1,2))*(rem(iE2)-rem(iE2-1))
        END IF
     END DO
  END DO
  open(14, FILE = 'rixsET.dat', STATUS = 'REPLACE')
  open(15, FILE = 'herfd.dat', STATUS = 'REPLACE')
  open(16, FILE = 'rixsEE.dat', STATUS = 'REPLACE')
  open(17, FILE = 'xasEI.dat', STATUS = 'REPLACE')
  open(18, FILE = 'xasEF.dat', STATUS = 'REPLACE')
  !open(19, FILE = 'rixsEI.dat', STATUS = 'REPLACE')
  DO iE1 = 1, ne1
     WRITE(17,'(20E30.15)') (rem(ie1) + Edge1(1))*hart, xasEI(iE1,1), xasEI(iE1,2)
  END DO
  DO iE1 = 1, ne1
     WRITE(18,'(20E30.15)') (rem(ie1) + Edge1(1) - Edge1(2))*hart, xasEF(iE1,1), xasEF(iE1,2)
  END DO
  DO iE1 = 1, ne1
     DO iE2 = 1, ne1
        WRITE(14,'(20E30.15)') (rem(iE2)+Edge1(1))*hart, &
          (rem(iE1)+Edge1(2))*hart, xsect_tmp(iE1,iE2,1), xsect_tmp(iE1,iE2,2)
     END DO
     WRITE(15,'(20E30.15)') (rem(iE1)+Edge1(1))*hart, xsect_tmp(iE1,iE1,1), xsect_tmp(iE1,iE1,2)
     WRITE(14,*)
  END DO
  !DO iE1 = 1, ne1
  !   DO iE2 = 1, ne1
  !      WRITE(19,'(20E30.15)') (rem(iE2)+Edge1(1)-Edge1(2))*hart, (rem(iE1)+Edge1(1))*hart, xsect_tmp(iE2,iE1,1), xsect_tmp(iE2,iE1,2)
  !   END DO
  !   WRITE(19,*)
  !END DO
  ! Interpolate onto a new grid, constant in incident and emission energies, so
  ! that one can plot HERFD, or emission at different energies.
  EEMin = rem(1) - rem(ne1)
  EEMax = rem(ne1) - rem(1)
  EIMin = rem(1)
  EIMax = rem(ne1)
  ! Loop over incident energy.
  DO iE2 = 1, ne1
     ! Loop over emission energy.
     EE = (EEMax - EEMin)/DBLE(ne1-1)*DBLE(iE2-1) + EEMin
     DO iE1 = 1, ne1
        EI = (EIMax - EIMin)/DBLE(ne1-1)*DBLE(iE1-1) + EIMin
        ! Get energy transfer.
        ET = EI - EE
        ! If out of range, don't extrapolate. Just output zero.
        IF((ET.LT.rem(1)).OR.(ET.GT.rem(ne1))) THEN 
           WRITE(16,'(20E30.15)') (EI + Edge1(1))*hart, (EE + Edge1(1) - Edge1(2))*hart, 0.d0, 0.d0
        ELSE
           CALL BLInterp2D(rem,rem,xsect_tmp(:,:,1),ne1,ne1,ne,EI,ET,xInt(1))
           CALL BLInterp2D(rem,rem,xsect_tmp(:,:,2),ne1,ne1,ne,EI,ET,xInt(2))
           WRITE(16,'(20E30.15)') (EI + Edge1(1))*hart, (EE + Edge1(1) - Edge1(2))*hart, DBLE(xInt(1)), DBLE(xInt(2))
        END IF
     END DO
     WRITE(16,*)
  END DO

  !DO iE1 = 1, ne1
  !   DO iE2 = 1, ne1 - iE1 + 1
  !      WRITE(16,'(20E30.15)') rem(iE2)*hart, (rem(iE2) - rem(iE2+iE1-1))*hart, xsect_tmp(iE2+iE1-1,iE2,1), xsect_tmp(iE2+iE1-1,iE2,2)
  !   END DO
  !   DO iE2 = ne1 - iE1 + 2, ne1
  !      WRITE(16,'(20E30.15)') rem(iE2)*hart, (rem(iE2) - rem(iE2+iE1-1))*hart, 0.d0, 0.d0
  !   END DO
  !   WRITE(16,*)
  !END DO
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
  CLOSE(19)
  IF(MBConv) THEN ! Read spectral function and convolve.
     xsect_rxs(:,:,:,:) = 0.d0 
     CALL ReadXMU('XES/xmu.dat', nw, xe = wXES, mu = muXES)
     !CALL ReadArrayData('XES/xmu.dat',Double1 = wXES, Double2 = muXES, NumElements = nw)
     wXES(:) = -(wXES(:)/hart - xmu)
     DO iw = 1, nw
        IF(iw.GE.(nw+1-iw)) EXIT
        xtmp(1) = wXES(iw)
        wXES(iw) = wXES(nw+1-iw)
        wXES(nw+1-iw) = xtmp(1)
        xtmp(1) = muXES(iw)
        muXES(iw) = muXES(nw+1-iw)
        muXES(nw+1-iw) = xtmp(1)
     END DO
     !muXES(:) = muXES(:)
     DO ind = 1, 2
        DO iE1 = 1, ne1
           !PRINT*, iE1, 'of', 2*ne1
           DO iE2 = 1, ne1
              ctmp(1) = 0.d0
              !IF(rem(iE2).LT.xmu) THEN
              !   xsect_rxs(iE2,iE1,ind,1) = xsect_tmp(iE2,iE1,ind)
              !   CYCLE
              !END IF
              DO iw = 2, nw
                 IF((rem(iE2) - wXES(iw)).LT.rem(1)) THEN
                    xInt(1) = 0.d0
                 ELSEIF((rem(iE2)-wXES(iw)).GT.rem(ne1)) THEN
                    xInt(1) = xsect_tmp(ne1,iE1,ind)
                 ELSE
                    CALL terpc(rem,xsect_tmp(:,iE1,ind),ne1,1,rem(iE2) - wXES(iw),xInt(1))
                 END IF
                 IF((rem(iE2) - wXES(iw-1)).LT.xmu) THEN
                    xInt(2) = 0.d0
                 ELSEIF((rem(iE2)-wXES(iw-1)).GT.rem(ne1)) THEN
                    xInt(2) = xsect_tmp(ne1,iE1,ind)
                 ELSE
                    CALL terpc(rem,xsect_tmp(:,iE1,ind),ne1,1,rem(iE2) - wXES(iw-1),xInt(2))
                 END IF
                 xsect_rxs(iE2,iE1,ind,1) = xsect_rxs(iE2,iE1,ind,1) + &
                   0.5d0*(xInt(1)*muXES(iw)+xInt(2)*muXES(iw-1))*(wXES(iw)-wXES(iw-1))
              END DO
           END DO
        END DO
     END DO
     xsect_tmp(:,:,:) = xsect_rxs(:,:,:,1)
  END IF
  IF(.FALSE.) THEN ! Read spectral function and convolve.
     CALL ReadArrayData('muXES.dat',Double1 = wXES, Double2 = muXES, NumElements = nw)
     wXES(:) = wXES(:)/hart
     muXES(:) = muXES(:)*hart
     DO ind = 1, 2
        DO iE1 = 1, ne1
           !PRINT*, iE1, 'of', 2*ne1           
           DO iE2 = 1, ne1
              ctmp(1) = 0.d0
              IF(rem(iE2).LT.xmu) THEN
                 xsect_rxs(iE2,iE1,ind,1) = xsect_tmp(iE2,iE1,ind)
                 CYCLE
              END IF
              DO iw = 2, nw
                 IF(wXES(iw).LT.(xmu-rem(iE2))) CYCLE
                 ctmp(1) = ctmp(1) + 0.5d0*(muXES(iw) + muXES(iw-1))*(wXES(iw) - wXES(iw-1))
                 IF((rem(iE2) + wXES(iw)).LT.xmu) THEN
                    xInt(1) = 0.d0
                 ELSEIF((rem(iE2)+wXES(iw)).GT.rem(ne1)) THEN
                    xInt(1) = xsect_tmp(ne1,iE1,ind)
                 ELSE
                    CALL terpc(rem,xsect_tmp(:,iE1,ind),ne1,1,rem(iE2) + wXES(iw),xInt(1))
                 END IF
                 IF((rem(iE2) + wXES(iw-1)).LT.xmu) THEN
                    xInt(2) = 0.d0
                 ELSEIF((rem(iE2)+wXES(iw-1)).GT.rem(ne1)) THEN
                    xInt(1) = xsect_tmp(ne1,iE1,ind)
                 ELSE
                    CALL terpc(rem,xsect_tmp(:,iE1,ind),ne1,1,rem(iE2) + wXES(iw-1),xInt(2))
                 END IF
                 xsect_rxs(iE2,iE1,ind,1) = xsect_rxs(iE2,iE1,ind,1) + &
                   0.5d0*(xInt(1)*muXES(iw)+xInt(2)*muXES(iw-1))*(wXES(iw)-wXES(iw-1))
              END DO
              xsect_rxs(iE2,iE1,ind,1) = xsect_rxs(iE2,iE1,ind,1)/ctmp(1)
           END DO
        END DO
     END DO
     xsect_tmp(:,:,:) = xsect_rxs(:,:,:,1)
  END IF

              
  ! Integrate over E1 and E2 to get line spectra.
  xasEI = 0.d0
  xasEF = 0.d0
  xes = 0.d0
  DO iE1 = 1, ne1
     DO iE2 = 2, ne1
        IF(rem(iE2).GT.EMax(1)) EXIT
        IF(rem(iE2).GT.EMin(1)) THEN
           xasEI(iE1,1) = xasEI(iE1,1) + 0.5d0*(xsect_tmp(iE2, iE1,1) + xsect_tmp(iE2-1,iE1,1))*(rem(iE2)-rem(iE2-1))
           xasEI(iE1,2) = xasEI(iE1,2) + 0.5d0*(xsect_tmp(iE2, iE1,2) + xsect_tmp(iE2-1,iE1,2))*(rem(iE2)-rem(iE2-1))
        END IF
     END DO
  END DO
  DO iE1 = 1, ne1
     DO iE2 = 2, ne1
        IF(rem(iE2).GT.EMax(2)) EXIT
        IF(rem(iE2).GT.EMin(2)) THEN
           xasEF(iE1,1) = xasEF(iE1,1) + 0.5d0*(xsect_tmp(iE1, iE2,1) + xsect_tmp(iE1,iE2-1,1))*(rem(iE2)-rem(iE2-1))
           xasEF(iE1,2) = xasEF(iE1,2) + 0.5d0*(xsect_tmp(iE1, iE2,2) + xsect_tmp(iE1,iE2-1,2))*(rem(iE2)-rem(iE2-1))
        END IF
     END DO
  END DO
  !DO iE1 = 1, ne1
  !   DO iE2 = 2, ne1
  !      xes(iE1) = xasEI(iE1) + 0.5d0*(xsect_tmp(iE1, iE2) + xsect_tmp(iE1,iE2-1))*(rem(iE2)-rem(iE2-1))
  !   END DO
  !END DO
  CALL wlog('Writing results.')
  open(14, FILE = 'rixsET-sat.dat', STATUS = 'REPLACE')
  open(15, FILE = 'herfd-sat.dat', STATUS = 'REPLACE')
  open(16, FILE = 'rixsEE-sat.dat', STATUS = 'REPLACE')
  open(17, FILE = 'xasEI-sat.dat', STATUS = 'REPLACE')
  open(18, FILE = 'xasEF-sat.dat', STATUS = 'REPLACE')
  !open(19, FILE = 'rixsEI-sat.dat', STATUS = 'REPLACE')
  DO iE1 = 1, ne1
     WRITE(17,'(20E30.15)') rem(ie1)*hart, xasEI(iE1,1), xasEI(iE1,2)
  END DO
  DO iE1 = 1, ne1
     WRITE(18,'(20E30.15)') rem(ie1)*hart, xasEF(iE1,1), xasEF(iE1,2)
  END DO
  DO iE1 = 1, ne1
     DO iE2 = 1, ne1
        WRITE(14,'(20E30.15)') (rem(iE2)+Edge1(1))*hart, (rem(iE1)+Edge1(2))*hart, xsect_tmp(iE1,iE2,1), xsect_tmp(iE1,iE2,2)        
     END DO
     WRITE(15,'(20E30.15)') (rem(iE1)+Edge1(1))*hart, xsect_tmp(iE1,iE1,1), xsect_tmp(iE1,iE1,2)        
     WRITE(14,*)   
  END DO
  !DO iE1 = 1, ne1
  !   DO iE2 = 1, ne1
  !      WRITE(19,'(20E30.15)') (rem(iE2)+Edge1(1))*hart, (rem(iE1)+Edge1(2))*hart, xsect_tmp(iE2,iE1,1), xsect_tmp(iE2,iE1,2)        
  !   END DO
  !   WRITE(19,*)   
  !END DO  
  ! Interpolate onto a new grid, constant in incident and emission energies, so
  ! that one can plot HERFD, or emission at different energies.
  EEMin = rem(1) - rem(ne1)
  EEMax = rem(ne1) - rem(1)
  EIMin = rem(1)
  EIMax = rem(ne1)
  ! Loop over incident energy.
  DO iE2 = 1, ne1
     ! Loop over emission energy.
     EE = (EEMax - EEMin)/DBLE(ne1-1)*DBLE(iE2-1) + EEMin
     DO iE1 = 1, ne1
        EI = (EIMax - EIMin)/DBLE(ne1-1)*DBLE(iE1-1) + EIMin
        ! Get energy transfer.
        ET = EI - EE
        ! If out of range, don't extrapolate. Just output zero.
        IF((ET.LT.rem(1)).OR.(ET.GT.rem(ne1))) THEN 
           WRITE(16,'(20E30.15)') (EI + Edge1(1))*hart, (EE + Edge1(1) - Edge1(2))*hart, 0.d0, 0.d0
        ELSE
           CALL BLInterp2D(rem,rem,xsect_tmp(:,:,1),ne1,ne1,ne,EI,ET,xInt(1))
           CALL BLInterp2D(rem,rem,xsect_tmp(:,:,2),ne1,ne1,ne,EI,ET,xInt(2))
           WRITE(16,'(20E30.15)') (EI + Edge1(1))*hart, (EE + Edge1(1) - Edge1(2))*hart, DBLE(xInt(1)), DBLE(xInt(2))
        END IF
     END DO
     WRITE(16,*)
  END DO
 
  !DO iE1 = 1, ne1
  !   DO iE2 = 1, ne1 - iE1 + 1
  !      WRITE(16,'(20E30.15)') rem(iE2)*hart, (rem(iE2) - rem(iE2+iE1-1))*hart, xsect_tmp(iE2+iE1-1,iE2,1), xsect_tmp(iE2+iE1-1,iE2,2)
  !   END DO
  !   DO iE2 = ne1 - iE1 + 2, ne1
  !      WRITE(16,'(20E30.15)') rem(iE2)*hart, (rem(iE2) - rem(iE2+iE1-1))*hart, 0.d0, 0.d0
  !   END DO
  !   WRITE(16,*)
  !END DO



  400 call par_barrier
      call par_end
  if(master)call WipeErrorfileAtFinish
  stop

END PROGRAM RIXS
