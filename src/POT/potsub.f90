!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: potsub.f90,v $:
! $Revision: 1.19 $
! $Author: jorissen $
! $Date: 2012/10/23 20:08:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pot !KJ put everything in modules 7-09

  USE AtomicPotIO
  use DimsMod, only: nphx=>nphu, lx, nrptx, nex
  use par
  use constants
  use potential_inp, ecv0=>ecv, folp0=>folp !KJ 7-09 this "renaming" used to happen through calling arguments having different local name
  use atoms_inp !KJ
  use workstrfacs2,only: eta,eta0 !KJ
  use controls,only: ispace !KJ
  use m_thermal_scf, only: thscf_main, thscf_init, thscf_deinit
  use m_sommerfeld_scf, only: sommerfeld_scf_main, sommerfeld_scf_init, sommerfeld_scf_deinit
  !     Cluster code -- multiple shell single scattering version of FEFF
  !     This program (or subroutine) calculates potentials
  !     for unique potentials specifed by atoms and overlap cards.

  implicit none !double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  !     Notes:
  !        nat    number of atoms in problem
  !        nph    number of unique potentials
  !        ihole  hole code of absorbing atom
  !        iph=0 for central atom

  !     ATOM output
  !     Note that ATOM output is dimensioned 251, all other r grid data is set to nrptx, currently 250
  !     rho(251,0:nphx+1)     -   density*4*pi
  real*8 rho(251,0:nphx+1)
  !     vcoul(251,0:nphx+1)   -   coulomb potential
  real*8 vcoul(251,0:nphx+1)
  real*8 dr(251), drho(251), dvcoul(251)

  !     Overlap calculation results
  !     overlapped density*4*pi
  real*8 edens(251,0:nphx), edensTmp(251,0:nphx), edensTrack(251,0:100)
  !     overlapped coul pot
  real*8 vclap(251,0:nphx), vclapp(251,0:nphx)
  !     overlapped total potential
  real*8 vtot (251,0:nphx), vtotTmp(251,0:nphx)

  real*8 ratTmp(3,natx)
  !     folp(0:nphx)  - overlap factor for rmt calculation
  real*8 folp(0:nphx), folpx(0:nphx)

  !     Muffin tin calculation results
  !     r mesh index just inside rmt
  integer imt(0:nphx)
  !     r mesh index just inside rnorman
  integer inrm(0:nphx)
  !     muffin tin radius
  real*8 rmt(0:nphx)
  !     norman radius
  real*8 rnrm(0:nphx), qnrm(0:nphx), qold(0:nphx)

  logical ok
  real*8, allocatable, dimension(:,:) :: xnmues, xnvmu, xnmues_old
  complex,ALLOCATABLE :: gtr(:)
  complex*16,ALLOCATABLE ::  xrhoce(:)
  complex*16,ALLOCATABLE ::  xrhole(:)
  complex*16,ALLOCATABLE ::  yrhoce(:)
  complex*16,ALLOCATABLE ::  yrhole(:)
  !     need irregular solution for complex potential. fix later
  real*8 dgc0(251), dpc0(251)

  !     additional data needed for relativistic version
  real*8 dgc(251,41,0:nphx+1), dpc(251,41,0:nphx+1)
  real*8 adgc(10,41,0:nphx+1), adpc(10,41,0:nphx+1)
  real*8 rhoval(251,0:nphx+1), edenvl(251,0:nphx)
  real*8 vvalgs (251,0:nphx)

  real*8 ri(nrptx)
  real*8 dmag(251,0:nphx+1)
  real*8 xnval(41,0:nphx+1), eorb(41,0:nphx+1)
  integer kappa(41,0:nphx+1), iorb(-5:4,0:nphx+1), norb(0:nphx+1)
  logical lpass, EmptyCell(nph)
  character*512 slog
  real*8 qtotel,xmu,wp,dx,x0,xntot,xnvmup,sum,e_chsh,ecv,s02,rhoint,vint,rs,xf,xmunew,rnrmav,ecv_save
  integer i,iph,il,ir,ip,iscmt,npr,ios,j,iCenter,nfree,idmag,ilpass, iii,lll
  !     Josh use nhtmp to save nohole value
  integer nhtmp, nstarts !nstarts : how often have we attempted scf
  integer, parameter :: nstartsmax = 3
  integer nscmt_min !KJ minimum number of scf iterations (3 or 2 depending on circumstances) 
  ! Dimension emu and erelax to hold all possible edge values.
  real*8 emu, erelax, vTmp
  real*8 chargedistance, partialchargedistance !KJ test of SCF convergence
  integer,parameter :: negx=80
  real*8 scfdos(negx,0:lx,0:nphx),edos(negx)
  real*8 dosTrack(negx,0:100)
  logical, parameter :: track_dos_convergence = .false.
  ! Debug: Fer
  ! correorb input
  integer          sh_iz, sh_ihole, sh_jri, sh_kappa(41)
  real*8    	   sh_rmt, sh_dx, sh_p2f, sh_edge, sh_ri(nrptx), &
  sh_dgcn(nrptx,41), sh_dpcn(nrptx,41), &
  sh_adgc(10,41), sh_adpc(10,41), &
  sh_eorb(41), rfms1_end
  complex*16       sh_vxc(nrptx)
  ! correorb output
  integer          sh_neg(41), sh_norbp, iramp
  real*8    	   sh_eng(nex,41), sh_rhoj(nex,41)

  ! Added by FDV
  real*8 Q_Tot

  10 format (4x, a, i5)

  ! Allocate arrays.
  allocate(xnmues(0:lx,0:nphx), xnvmu(0:lx,0:nphx+1), xnmues_old(0:lx,0:nphx))
  ALLOCATE(gtr((lx + 1) * (nphx + 1) * numprocs))
  ALLOCATE(xrhoce((lx + 1) * (nphx + 1) * numprocs))
  ALLOCATE(xrhole((lx + 1) * (nphx + 1) * numprocs))
  ALLOCATE(yrhoce(251 * (nphx + 1) * numprocs ))
  ALLOCATE(yrhole(251 * (lx + 1) * (nphx + 1) * numprocs ))

  ! Allocate some more arrays:

  ! Set final scf radius if ramping requested.
  IF(ramp_scf) THEN
     rfms1_end = rfms1
     rfms1 = rfms1_start
     IF(rfms1.GE.rfms1_end) THEN
        rfms1 = rfms1_end
        nramp = 0
     END IF
     iramp = 1
  END IF
     
  !     Josh - for now if nohole=2 reset to 0 so that regular nohole potential is used
  nhtmp = nohole
  if (nohole.eq.2) then
     ! nohole = 2 means rpa core-hole. For XAS, add screened core-hole to nohole calculation
     ! while for XES, subtract screened core-hole to core-hole calculation.
     if(ispec.EQ.2) then
        nohole = -1
     else
        nohole = 0
     end if
  end if
  !     Josh

  if (StartFromFile) then
    nscmt_min=2
   else
     nscmt_min=3
   endif
 
 
   !     variables ecv0 and folp0 serve as input only; do not change them
   !     since it will change file feff.ior content
   !     ecv and folp are passed through pot.bin to next modules.
   ecv = ecv0
   folp(0:nph) = folp0(0:nph)
   ok=.true. !KJ
   if (master) open(28,file='convergence.scf',access='append',status='unknown',form='formatted') !KJ
   if (master) write(28,*) '# it. E_fermi(eV)  Charge Distance  Partial Chg. D.  Convergence'
   if (master) open(29,file='convergence.scf.fine',access='append',status='unknown',form='formatted') !KJ
 
   kappa(:,:) = 0
 
   call inipot (dgc, dpc, edenvl, vvalgs, xnmues) !simply initializes these arrays to 0 LOL
 
   !     increase the length of hydrogen bonds for potential only
   call moveh (nat, iphat, iz, rat)


  nfree = 1
  do i=0,nph
    if (abs(xion(i)) .gt. 1.d-3) nfree = 2
  enddo

  CALL ReadAtomicPots(nph, iz(0:nph), ihole, rho, dmag(:,0:nph+1), rhoval, vcoul, dgc0,  &
  & dpc0, dgc(:,:,0:nph+1), dpc(:,:,0:nph+1), adgc(:,:,0:nph+1), adpc(:,:,0:nph+1), &
  & erelax, emu, xnvmu, xnval(:,0:nph+1), norb, eorb, drho, dvcoul, iphat,    &
  & ratTmp, iatph(0:nph), novr(0:nph), iphovr, nnovr, rovr, nat, edens, &
  & edenvl, vclap,  rnrm(0:nph), kappa(:,0:nph+1), iorb(:,0:nph+1), s02)

  edensTrack(:,0)=edens(:,0) !atomic overlap for ipot=0 goes into iteration 0 of density tracker

  DO iph = 1, nph
    IF(iz(iph).eq.0) THEN
      EmptyCell(iph) = .TRUE.
      !            folp(iph) = 1.d0
    ELSE
      EmptyCell(iph) = .FALSE.
    END IF
  END DO

  !  end of free atom calculations (might be done twice if ION used)
  !       Debugging - read values from old version that matches feff84.
  !       OPEN(unit=39,file='fort.39',status='old')
  !       READ(39,*) nph, iz, rho, dmag, rhoval, vcoul, dgc0, dpc0, &
  !            &     dgc, dpc, adgc, adpc, erelax(ihole), emu(ihole), xnvmu, xnval, norb,     &
  !           & eorb, drho, dvcoul, iphat, iatph, novr, iphovr, nnovr, &
  !           &     rovr, nat, edens, edenvl, vclap, rnrm, kappa, iorb

  !c new patch
  !     itest = 1
  !     if (itest.eq.1) then
  !c      use orbitals with core-hole for initial orbitals
  !c      orthogonaliztion problem for NRIXS calculations
  !       do i = 1, 251
  !         dgc0(i) = dgc(i,0)
  !         dpc0(i) = dpc(i,0)
  !       enddo
  !     endif
  !c end new patch

  !     Find total charges for istprm
  !     qtotel - total number of e in a cluster
  qtotel = 0
  do iph = 0,nph
    qtotel = qtotel + (iz(iph)-xion(iph)) * xnatph(iph)
  enddo
  !     photoelectron moves out of the system
  !     do not remove now since we are putting screening electron back


  !     Find muffin tin radii, add gsxc to potentials, and find interstitial parameters
  if(master)call wlog('Muffin tin radii and interstitial parameters [bohr]:')

  rmt(0) = -1
  xmu = 100.d0
  if (iafolp.ge.0) then
    do iph=0,nph
      folpx(iph) = folp(iph)
      folp(iph) = 1
    enddo
  endif

  ! Debug: FDV
  !       write(6,*) ' Entering afolp'
  !       do iph=0,nph
  !         write(6,fmt='(a,i4,3f16.10)'), &
  !         'rmt, folpx, folp: ', iph, rmt(iph)*bohr, folpx(iph), folp(iph)
  !       end do

  idmag = 0

  call istprm (nph, nat, iphat, rat, iatph, xnatph,                 &
  &            novr, iphovr, nnovr, rovr, folp, folpx, iafolp,       &
  &            edens, edenvl, idmag,                                 &
  &            dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,      &
  &            ixc, rhoint,vint, rs, xf, xmu, xmunew,                &
  &            rnrmav, qtotel, inters, totvol) !KJ ,EmptyCell)
  xmu = xmunew

  ! Debug: FDV
  !       write(6,*) ' Entering afolp'
  !       do iph=0,nph
  !         write(6,fmt='(a,i4,3f16.10)'), &
  !         'rmt, folpx, folp: ', iph, rmt(iph)*bohr, folpx(iph), folp(iph)
  !       end do

  !     Automatic max reasonable overlap
  if (iafolp .ge. 0)  then
    call afolp (nph, nat, iphat, rat, iatph, xnatph,               &
    &               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
    &               edens, edenvl,                                     &
    &               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
    &               ixc, rhoint,vint, rs, xf, xmu, xmunew,             &
    &               rnrmav, qtotel, inters, totvol)
    xmu =xmunew
    ! Debug: FDV
    !     do iph=0,nph
    !       write(6,fmt='(a,i4,f16.10)'), 'rmt: ', iph, rmt(iph)*bohr
    !     end do
    !     stop

  endif

  !     wp is plasmon frequency in hart
  wp = sqrt(12.*rs/fa**4) * xf**2 / 2.d0

  !     Phase shift calculation
  !     Atom r grid
  dx = 0.05d0
  x0 = 8.8d0


  !     Find self-consistent muffin-tin potential.
  qnrm(:)=0
  qold(:)=0
  xnmues_old(:,:)=0.d0
  nstarts=1

  100 continue

  !KJ 12-2010 modified following loop to prevent iteration being repeated ad infinitum.
  if (nscmt.gt.0 .or. (ispec.ne.0 .and. ispec.lt.4)) then
    ecv_save=ecv
    call corval(ecv, xnvmu, eorb, norb, xnval, kappa, rgrd, nohole, &
    nph, edens, edenvl, vtot, vvalgs, rmt, rnrm, ixc, rhoint, vint, jumprm,  &
    x0, ri, dx, xion, iunf, iz, adgc, adpc, dgc, dpc, ihole, lmaxsc)
    if(.not.ok) then
      nstarts=nstarts+1
      if((abs(ecv-ecv_save).lt.0.05) .and. (nstarts.eq.2)) nstarts=3
      if((abs(ecv-ecv_save).lt.0.05) .or. (nstarts.eq.3)) ca1=max(ca1/5.d0,0.01)   ! 0.05Ha = 1.3eV
      !KJ         2 alternative attempts are allowed: one for a change in e_cv, and one for reduced mixing factor.
      if(nstarts.gt.nstartsmax) then
        call wlog('Bad counts found repeatedly - pot will exit now.')
        goto 400
      endif
    endif
  endif

  !      CALL WriteExternalPot(vtot, vint, edens, rhoint, rat(:,1:nat), xmu, imt, rmt)
  IF(ExternalPot) THEN
    CALL ReadExternalPot(vtot(:,0:nat-1), vint, edens(:,0:nat-1), rhoint, rat(1:3,1:nat), xmu, imt(1:nat), rmt(1:nat))
    !         IF(vint.gt.-0.05d0) vint = -0.05d0
  END IF
  !      vtot(:,5:8) = 0.d0
  !      edens(:,5:8) = 0.d0

  if (StartFromFile .and. nstarts.eq.1) then
    !KJ forget about the atomic overlap.  Open an existing pot.bin instead and use the potentials found there
    !KJ as starting point for this scf cycle.
    if(master)call wlog('Taking initial potentials from pre-calculated file.')
    if(master)call wlog(':WARNING If the old and new calculation were not precisely matched, all bets are off.')
    !KJ It may be desirable to overwrite only select arrays instead of using the entire file pot.bin .
    !KJ Need to discuss this.
    call importpot(vtot,vint,edens,rhoint,xmu)
    !         call importpot(rnrmav, xmu, vint, rhoint, emu, erelax, ecv,rs,xf,        &
    !                       edens, vclap, vtot, edenvl, vvalgs, qnrm, xnmues, inters, totvol)
  endif
  !           CALL ReadAtomicPots(nph, iz(0:nph), ihole, rho, dmag(:,0:nph+1), rhoval, vcoul, dgc0,  &
  !      & dpc0, dgc(:,:,0:nph+1), dpc(:,:,0:nph+1), adgc(:,:,0:nph+1), adpc(:,:,0:nph+1), &
  !      & erelax, emu, xnvmu, xnval(:,0:nph+1), norb, eorb, drho, dvcoul, iphat,    &
  !      & ratTmp, iatph(0:nph), novr(0:nph), iphovr, nnovr, rovr, nat, edens, &
  !      & edenvl, vclap,  rnrm(0:nph), kappa(:,0:nph+1), iorb(:,0:nph+1), s02)

  !     find a total number of valence electrons
  !     xntot - required number of valence electrons below fermi level
  !     xnvmu(iph) = xnvmu(iph)-xion(iph)
  !     xnvmu - number of valence electron within norman sphere
  xntot=0.0d0
  do iph=0,nph
    xnvmup = 0
    do i = 0,lx
      xnvmup = xnvmup + xnvmu(i,iph)
    enddo
    !    	 x35 and earlier   xntot = xntot + xnatph(iph)*(xnvmup+xion(iph))
    xntot = xntot + xnatph(iph) * xnvmup
    ! Debug: FDV
    !        write(6,fmt='(i3,f12.3)') iph, xnatph(iph) * xnvmup
  enddo

  ! Debug: FDV
  ! Testing a total charge change in the SCF
  !     read(956,*) Q_Tot
  !     xntot = xntot - Q_Tot

  !write(*,*) 'xntot is: ',xntot
  !write(*,*) 'xnvmu(:,0)',xnvmu(:,0)
  !write(*,*) 'xnvmu(:,1)',xnvmu(:,1)

  !     need to update vxcval in case the core-valence separation was
  !     made in subroutine corval. Need vxcval only for nonlocal exchange.
  if (mod(ixc,10).ge.5) then
    call  istprm (nph, nat, iphat, rat, iatph, xnatph,             &
    &               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
    &               edens, edenvl, idmag,                              &
    &               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
    &               ixc, rhoint,vint, rs, xf, xmu, xmunew,             &
    &               rnrmav, qtotel, inters, totvol)
    xmunew = xmu
  endif
195 CONTINUE
  if(master) then
    write(slog,131) ecv*hart
    call wlog(slog)
    write(slog,130) xmu*hart
    call wlog(slog)
    !         write(slog,132) ca1
    !         call wlog(slog)
    130    format('Initial Fermi level:              mu= ',f9.3,' eV')
    131    format('Core-valence separation energy:  ecv= ',f9.3,' eV')
    132    format(' ca1= ',f9.3)
  endif
  dosTrack=0.d0

  !     do first nmix iterations with mixing scheme. Need for f-elements.
  140 nmix=nmix-1

  !     number of processors for parallel execution
  npr = numprocs
  if (nscmt.le.0) goto 210
  !**************** BEGIN OF SCF CYCLE *******************
  do 200 iscmt =1,nscmt
    !        need to store coulomb potential
    vclapp(:,:)=vclap(:,:)

    if( (ispace.eq.0) .and. (iscmt.gt.1) .and. (abs(eta-eta0).gt.0.1) ) eta=eta0
    xnmues_old=xnmues

    if (scf_temperature.gt.0) then
      ! call thermal SCF routine. the init and de-init will eventually be moved
      ! outside of this loop (need to be careful with gotos to avoid leaks...)
      if (iscfth.EQ.1) then
          ! Sommerfeld expansion
          call sommerfeld_scf_init(ecv, xmunew, iscmt)
          call sommerfeld_scf_main(iscmt, ecv, vclap, edens, edenvl, vtot, vvalgs, &
          rmt, rnrm, qnrm, rhoint, vint, xmunew, xntot, xnvmu, xnval, x0, &
          ri, dx, adgc, adpc, dgc, dpc, rhoval, xnmues, ok, rgrd)
          call sommerfeld_scf_deinit
      elseif (iscfth.EQ.2) then
          ! m_thermal_scf
          call thscf_init(ecv, xmunew, iscmt)
          call thscf_main(iscmt, ecv, vclap, edens, edenvl, vtot, vvalgs, &
          rmt, rnrm, qnrm, rhoint, vint, xmunew, xntot, xnvmu, xnval, x0, &
          ri, dx, adgc, adpc, dgc, dpc, rhoval, xnmues, ok, rgrd)
          call thscf_deinit
      endif
    else
      IF(master) PRINT*, 'Running SCF with rscf = ', rfms1*bohr, iramp, nramp+1
      if (npr.le.1) then
        PRINT*, "Zero temperature single thread"
        call scmt(iscmt, ecv, nph, nat, vclap, edens, edenvl, vtot, vvalgs, rmt, rnrm, qnrm,   &
        &                  ixc, rhoint, vint, xmunew, jumprm, xntot, xnvmu, xnval,        &
        &                  x0, ri, dx, xnatph, xion, iunf, iz, adgc, adpc, dgc,dpc, ihole,        &
        &                  rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok, rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1 & !)
        ,edos,scfdos)
      else
        call scmtmp (npr,  iscmt, ecv, nph, nat, vclap, edens, edenvl, vtot, vvalgs, rmt, rnrm, qnrm,            &
        &                  ixc, rhoint, vint, xmunew, jumprm, xntot, xnvmu, xnval,     &
        &                  x0, ri, dx, xnatph, xion, iunf, iz, adgc, adpc, dgc,dpc, ihole,   &
        &                  rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok, rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1,    &
        &                  gtr, xrhole, xrhoce, yrhole, yrhoce )
      endif
    endif

    if (track_dos_convergence) then
      edensTrack(:,iscmt)=edens(:,0)
      do lll=0,lx
        do iii=1,negx
          dosTrack(iii,iscmt)=dosTrack(iii,iscmt)+scfdos(iii,lll,0)
        enddo
      enddo
      open(49,file="ldos0_history.dat")
      do i=1,negx !251
        write(49,'(101f20.5)') edos(i),dosTrack(i,0:iscmt)
      enddo
      close(49)
    endif


    if (.not. ok) goto 100
    !        if need to change core-valence separation then start scmt loop all over again

    !        write out Fermi level and charge transfers  and do tests of self-consistency
    lpass = .true.
    chargedistance=0.d0
    partialchargedistance=0.d0
    if (iscmt.lt.nscmt .and. iscmt.le.nscmt_min) lpass =.false. !have to do at least nscmt_min iterations
    write (slog,150)   xmunew*hart
    150    format ('New Fermi level:    mu= ',f9.3,' eV')
    !if(master)call wlog(slog)
    ! CONVERGENCE TEST 1 : the FERMI LEVEL must be converged to tolmu
    if (abs (xmunew - xmu) .gt. tolmu) lpass = .false.
    !if(master)call wlog(' Charge transfer:  type  charge ')
    if(master) write(29,*) ' Charge transfer:  type  charge '
    do iph=0,nph
      write (slog,180) iph, -qnrm(iph) + xion(iph)
      !call wlog(slog)
      write(29,*) trim(slog)
      ! CONVERGENCE TEST 2 : the TOTAL CHARGE ON ATOM IPH must be converged to tolq
      if (abs(qnrm(iph)-qold(iph)).gt.tolq) lpass = .false.
      chargedistance=max(chargedistance,abs(qnrm(iph)-qold(iph)))
      qold(iph) = qnrm(iph)
      !           check self-consistency of charges
      ! CONVERGENCE TEST 3 : the TOTAL VALENCE ELECTRON CHARGE ON ATOM IPH must equal the formal number of electrons to within tolsum
      sum = -qnrm(iph)
      do il=0,lx
        sum = sum + xnmues(il,iph) - xnvmu(il,iph)
        ! CONVERGENCE TEST 4 : the LOCAL PARTIAL VALENCE CHARGE must be converged to tolqp
        if (abs(xnmues(il,iph)-xnmues_old(il,iph)) .gt. tolqp) lpass = .false.
        partialchargedistance=max(partialchargedistance,abs(xnmues(il,iph)-xnmues_old(il,iph)))
      enddo
      !if(master) write(*,'(i3,x,20f12.5)') iph,qnrm(iph),qold(iph),abs(qnrm(iph)-qold(iph)),tolq
      !do il=0,lx
      !   if(master) write(*,'(a3,i3,x,20f12.5)') '   ',il,xnmues(il,iph),xnvmu(il,iph),xnmues(il,iph)-xnvmu(il,iph),sum
      !enddo
      if (abs(sum).gt.tolsum) lpass = .false.
    enddo
    180    format('     ',i3, 2f9.3)
    !KJ:     provide easily trackable information on SCF convergence :
    if (iscmt.eq.1 .and. master) then
      write(slog,'(i3,2x,f10.3,5x,f10.3,5x,f10.4,5x,i1)') 0,     xmu*hart,    0.d0,  0.d0,         0
      write(28,*) trim(slog)
      write(29,*) trim(slog)
    endif
    if (lpass) ilpass=1
    if (.not.lpass) ilpass=0
    write(slog,'(i3,2x,f10.3,5x,f10.4,5x,f10.4,5x,i1)') iscmt, xmunew*hart, chargedistance, partialchargedistance,ilpass
    if(master) write(28,*) trim(slog)
    if(master) write(29,*) trim(slog)

    write(slog,'("New Fermi level:    mu=",f8.3," eV", &
    "  Charge distance=",f8.4," (partial c.d.=",f8.4,")")') &
    xmunew*hart, chargedistance, partialchargedistance
    if(master) call wlog(slog)

    xmu = xmunew
    !:KJ
    !        recalculate core density (edens) here. fix later. ala
    !        call scfdat
    !        for now use the old core density
    if (iscmt.eq.nscmt .or. lpass) then
      !           restore  total density from previous iteration
      do ip=0,nph
        do ir=1,251
          !                need total density for istprm
          edens(ir,ip) = edens(ir,ip)-rhoval(ir,ip)+edenvl(ir,ip)
          vclap(ir,ip) = vclapp(ir,ip)
        enddo
        !             remember the reported charge transfer
        qnrm(ip) = -qnrm(ip) + xion(ip)
      enddo
      !!           exit self-consistency loop
      !            goto 210
    else
      !           update valence density
      edenvl(:,0:nph)=rhoval(:,0:nph)
      !           need total density for istprm
    endif
    if(lpass .or. (nscmt.eq.0)) then
      goto 210
    elseif(iscmt.eq.nscmt) then
      goto 200
    endif

    call  istprm (nph, nat, iphat, rat, iatph, xnatph, novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
    edens, edenvl, idmag, dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
    ixc, rhoint,vint, rs, xf, xmu, xmunew, rnrmav, qtotel, inters, totvol)
    xmunew = xmu
    if (nmix.gt.0) goto 140

    200 continue
    !**************** END OF SCF CYCLE *******************
    !     suspicious exit: run out of iterations (iscmt=nscmt)
    if(master)call wlog (':WARNING Convergence not reached; ran out of iterations.')
    goto 211

    !     clean exit from the loop: self-consistency is achieved
    210 continue
    if (nscmt.gt.0 .and. master) then
      call wlog('Electronic configuration')
      call wlog('  type     l     N_el')

      do iph= 0,nph
        do il = 0,lx
          write (slog,'(2i6, f9.3)') iph,il,xnmues(il,iph)
          call wlog(slog)
        enddo
      enddo
      call wlog('Charge transfer:  type  charge ')
      do iph=0,nph
        write (slog,'(5x,i3, 2f9.3)') iph, -qnrm(iph) + xion(iph)
        call wlog(slog)
      enddo
      write(slog,'(a,i4,a)') 'Convergence reached in ',iscmt,' iterations.'
      call wlog(slog)
    endif
    211 CONTINUE
    
    IF(ramp_scf) THEN
       ! increase rfms1 and restart scf loop
       IF(rfms1.LT.rfms1_end) THEN
         rfms1 = rfms1_start + (rfms1_end - rfms1_start)*DBLE(iramp)/DBLE(nramp)
         iramp = iramp + 1
         GOTO 195
       END IF
    END IF
    if (master) close(28)  ! convergence.scf file


    if (track_dos_convergence) then
      open(49,file="ldos0_history.dat")
      if (iscmt.gt.nscmt) iscmt=nscmt
      do i=1,negx !251
        !         write(49,'(i7,x,100f20.5)') i,edensTrack(i,0:min(12,iscmt))
        write(49,'(101f20.5)') edos(i),dosTrack(i,0:iscmt)
      enddo
      close(49)
    endif


    if (worker) go to 400

    if (nohole.gt.0) then
      !        testing new final state potential
      do 220 j = 1,251
        220    edens(j,0) = edens(j,0) - drho(j)

        !        notice that vclap is actually for the next iteration
        !        in SCMT loop, thus vclap may be wrong if self-consistency has not been reached
        do 230 j = 1,251
          230    vclap(j,0) = vclap(j,0) - dvcoul(j)

          call  istprm (nph, nat, iphat, rat, iatph, xnatph,             &
          &      novr, iphovr, nnovr, rovr, folp, folpx, iafolp,             &
          &      edens, edenvl, idmag,                                       &
          &      dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,            &
          &      ixc, rhoint,vint, rs, xf, xmu, xmunew,                      &
          &      rnrmav, qtotel, inters, totvol)
        endif

        !     correct the excitation energy
        !     emu = emu -vclap(1,0) + vcoul(1,0) done also above
        !     emu = emu+xmu  should be done in principle but leads
        !     to worse estimate of edge position. fix later. ala

        if (ipr1 .ge. 2 .and. master)  call wpot (nph, edens, imt, inrm, rho, vclap, vcoul, vtot, ntitle, title)


        ! Debug: Fer
        ! Changing eorb by hand for now
        !     eorb(1,0) = -14.0050405749
        !     eorb(1,0) = -14.1374268761

        ! Debug: Fer
        !     write(6,fmt='(a,f16.10)') ' eorb in pot: ', eorb(1,0)

        ! Debug: Fer
        ! Write out some of the stuff that goes into pot.bin
        !     print *, 'xmu:    ', xmu
        !     print *, 'vint:   ', vint
        !     print *, 'eorb:   ', eorb(1,0)
        !     print *, 'emu:    ', emu
        !     print *, 'erelax: ', erelax
        !     print *, 'iorb:   ', iorb

        if ( (ChSh_Type == 1) .or. (ChSh_Type == 3) ) then

          !       Set the target center
          iCenter = 0
          sh_iz            = iz(iCenter)
          sh_ihole         = ihole
          sh_rmt           = rmt(iCenter)
          sh_jri           = imt(iCenter)+1
          sh_dx            = dx
          sh_ri            = ri
          sh_p2f           = xmu-vint
          sh_edge          = xmu
          sh_vxc(1:251)    = vtot(1:251,iCenter) - vint
          sh_dgcn(1:251,:) = dgc(:,:,iCenter)
          sh_dpcn(1:251,:) = dpc(:,:,iCenter)
          !       sh_adgc          = adgc(:,:,iCenter)
          !       sh_adpc          = adpc(:,:,iCenter)
          sh_adgc          =  0.0d0
          sh_adpc          =  0.0d0
          sh_eorb          = eorb(:,iCenter)
          sh_kappa         = kappa(:,iCenter)

          call correorb(sh_iz, sh_ihole, sh_rmt, sh_jri, &
          sh_dx,sh_ri, sh_p2f, sh_edge, sh_vxc, &
          sh_dgcn, sh_dpcn, sh_adgc, sh_adpc, &
          sh_eorb, &
          sh_neg, sh_eng, sh_rhoj, sh_kappa, sh_norbp, 0) !KJ iph=0
          e_chsh = -sh_eng(1,1)
          write(6,fmt='(a,f12.8)') 'e_chsh: ', e_chsh
          emu = e_chsh + erelax
        end if

        !     write stuff into pot.bin
        if (master) call wrpot (nph, ntitle, title, rnrmav, xmu, vint, rhoint,        &
        &            emu, s02, erelax, wp, ecv,rs,xf, qtotel,              &
        &            imt, rmt, inrm, rnrm, folp, folpx, xnatph,            &
        &            dgc0, dpc0, dgc, dpc, adgc, adpc,                     &
        &            edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,      &
        &            eorb(1,0), kappa(1,0), iorb, qnrm, xnmues, nhtmp,     &
        &            ihole, inters, totvol, iafolp, xion, iunf, iz, jumprm,ipr1) !KJ added ipr1 3-2012

        !     write misc.dat
        if (ipr1 .ge. 1 .and. master)  then
          open (unit=1, file='misc.dat', status='unknown', iostat=ios)
          call chopen (ios, 'misc.dat', 'potph')
          call wthead(1, ntitle, title)
          close (unit=1)
        endif

        400 call par_barrier

        OPEN(UNIT=1993, FILE='chemical.dat',STATUS='REPLACE', POSITION="APPEND")
        WRITE(1993,*) scf_temperature, scf_temperature/0.000086173423, xmu*hart
        CLOSE(1993)

        !     Deallocate local variables
        deallocate(xnmues,xnvmu,gtr,xrhoce,xrhole,yrhoce,yrhole)

        return
      end
