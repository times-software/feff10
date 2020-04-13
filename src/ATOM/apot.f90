!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: apot.f90,v $:
! $Revision: 1.23 $
! $Author: jorissen $
! $Date: 2013/01/18 02:19:08 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AtomicPotentials
  USE AtomicPotIO
  USE ErrorMod
  USE constants
  USE DimsMod, only: nphx=> nphu, lx
  USE par
  use config,only: IsOccupied, DumpConfig2
  use atoms_inp,only: nat,nph,iphat,iatph,rat
  use potential_inp,only: nohole, ihole, novr, iphovr, nnovr, rovr, xion, iunf, iz, ipr1
  use rixs_inp
  !     Josh Kas - modified pot.f to calculate the atomic potentials separately from
  !                the scf loop. Note also that moveh has been removed.

  !     Cluster code -- multiple shell single scattering version of FEFF
  !     This program (or subroutine) calculates free atom potentials and
  !     saves the information in apot.bin

  IMPLICIT NONE
  INTEGER, parameter :: NRPts = 251

  ! Input:
  ! iz(0:nphx) - atomic number, input
  ! nohole     - control over core-hole treatment
  ! nph        - number of unique pots
  ! ihole      - hole/edge code of absorbing atom
  ! iphat      - given specific atom, which unique pot?
  ! iunf       - freeze f-orbitals or not?
  ! ipr1       - print control
  ! nat        - number of atoms in problem.
  ! xion(0:nphx)  - ionicity, input
  !     rat(3,natx)  -  cartesian coords of specific atom
  !     iatph(0:nphx) - given unique pot, which atom is model?
  !                   (0 if none specified for this unique pot)
  !     novr(0:nphx)  - number of overlap shells for unique pot
  !     iphovr(novrx,0:nphx)  - unique pot for this overlap shell
  !     nnovr(novrx,0:nphx)  -  number of atoms in overlap shell
  !     rovr(novrx,0:nphx)   -  r for overlap shell

  ! Local variables:
  
  ! xionp         - temp variable to hole xion(iph)
  real*8 xionp

  !     overlapped density*4*pi
  real*8, allocatable :: edens(:,:), edenvl(:,:), vclap(:,:)

  !     ATOM output
  !     Note that ATOM output is dimensioned NRPts, all other r grid
  !     data is set to nrptx, currently 250
  !     rho(NRPts,0:nphx+1)     -   density*4*pi
  real*8, allocatable :: rho(:,:)
  !     rnrm - norman radius of each unique pot.
  real*8, allocatable :: rnrm(:), rnrmTmp(:)
  !     vcoul(NRPts,0:nphx+1)   -   coulomb potential
  real*8, allocatable :: vcoul(:,:)
  real*8 dr(NRPts), drho(NRPts), dvcoul(NRPts)

  !     need irregular solution for complex potential. fix later
  real*8 dgc0(NRPts), dpc0(NRPts), dgc00(NRPts), dpc00(NRPts)

  !     additional data needed for relativistic version
  real*8, allocatable, dimension(:,:,:) :: dgc, dpc, adgc, adpc
  real*8, allocatable, dimension(:,:) :: rhoval,dmag, xnvmu, xnval, eorb,xnel !KJ 5-2012 added xnel througout to enable configuration output
  real*8 et, etinit, etfin, efrozn, erelax, emu, s02, hx, x0, Mkkp, gamch  ! JK - added Mkkp matrix elements and gamch for RIXS calcs.
  integer, allocatable :: kappa(:,:), iorb(:,:), norb(:), nqn(:,:)  !KJ 9-2012 added nqn
  
  CHARACTER*512 slog
  !     Josh use nhtmp to save nohole value
  integer nhtmp, nfree, i, ifree, ispinr, iph, itmp, iEdge, ihole0


10 format (4x, a, i5)

  ! Allocate local arrays:
  allocate(edens(NRPts,0:nphx), edenvl(NRPts,0:nphx),vclap(NRPts,0:nphx))
  allocate(rho(NRPts,0:nphx+1))
  allocate(rnrm(0:nphx), rnrmTmp(0:nphx))
  allocate(vcoul(NRPts,0:nphx+1))
  allocate(dgc(NRPts,30,0:nphx+1), dpc(NRPts,30,0:nphx+1))
  allocate(adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1), rhoval(NRPts,0:nphx+1))
  allocate(dmag(NRPts,0:nphx+1), xnvmu(0:lx,0:nphx+1))
  allocate(xnval(30,0:nphx+1), eorb(30,0:nphx+1), xnel(30,0:nphx+1))
  allocate(kappa(30,0:nphx+1), iorb(-4:3,0:nphx+1), norb(0:nphx+1),nqn(30,0:nphx+1))
  CALL rixs_read

  
  ! loop over edges to make edges.dat for rixs calculations.
  OPEN(UNIT = 17, FILE = 'edges.dat', STATUS = 'REPLACE')
  DO iEdge = 1, RixsI%nEdges
     IF(TRIM(ADJUSTL(RixsI%Edges(iEdge))).NE.'VAL') THEN
        CALL setedg(RixsI%Edges(iEdge),ihole)
     ELSE
        ! Write to edges.dat for valence
        WRITE(17,*) 0.d0, 1.d0, RixsI%gam_exp(2)
        EXIT 
     END IF
     ! Initialize variables
     dgc(:,:,:) = 0.d0
     dpc(:,:,:) = 0.d0
     kappa(:,:) = 0
     nqn(:,:) = 0
     edens(:,:) = 0.d0
     edenvl(:,:) = 0.d0
     vclap(:,:) = 0.d0
     rho(:,:) = 0.d0
     rnrm(:) = 0.d0
     rnrmTmp(:) = 0.d0
     adgc(:,:,:) = 0.d0
     adpc(:,:,:) = 0.d0
     rhoval(:,:) = 0.d0
     dmag(:,:) = 0.d0
     xnvmu(:,:) = 0.d0
     xnval(:,:) = 0.d0
     xnel(:,:) = 0.d0
     eorb(:,:) = 0.d0
     iorb(:,:) = 0
     norb(:) = 0
     vcoul(:,:) = 0.d0
     
     !     increase the length of hydrogen bonds for potential only
     call moveh (nat, iphat, iz, rat)
     
     !     Josh - for now if nohole=2 reset to 0 so that regular nohole
     !     potential is used
     nhtmp = nohole
     if (nohole.eq.2) nohole = 0

     nfree = 1
     ! If any of the atoms are ionized, nfree = 2.
     do i=0,nph
        if (abs(xion(i)) .gt. 1.d-3) nfree = 2
     end do
     
     !     Free atom potentials and densities
     !     Final state is (usually) with a core hole, initial state is 
     !     w/o a corehole.
     !     NB wsatom is needed in SUMAX, if changed here, change it there
     
     !     do not save spinors
     !     Call twice if any of xion.neq.0 ( first time with xion=0 to set rnrm)
     
     do ifree = 1, nfree
        
        ispinr = 0
        ! Calculate for absorbing atom (iph = 0)
        iph = 0
        ! Calculate all edges up to ihole - 1
        !      do itmp = 1, ihole - 1
        !         ! Only calculate orbitals that contain at least one electron.
        !         IF(IsOccupied(iz(0),itmp)) THEN
        !            ! If user has specified NOHOLE, run and fill nph+1 elements of arrays
        !            if (nohole.ge.0) then
        !               xionp = xion(iph)               
        !               if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
        !               ! Run self-consistend dirac-fock atomic solver.
        !               call scfdat ( ipr1, nph+1, nph, iz(0), itmp, xionp, iunf,    &
        !                    &     vcoul(1,nph+1), rho(1,nph+1), dmag(1,nph+1), rhoval(1,nph+1),&
        !                    &     ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                    &
        !                    &     s02, efrozn, et, xnvmu(0,nph+1),                             &
        !                    &     xnval(1,nph+1), iorb(-4,nph+1), norb(nph+1),                 &
        !                    &     eorb(1,nph+1), kappa(1,nph+1) )
        
        !            ! Otherwise fill iph elements of arrays.
        !            else
        
        !               xionp = xion(iph)
        !               if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
        
        !               ! Run self-consistent dirac-fock atomic solver.
        !               call scfdat ( ipr1, iph, nph, iz(iph), itmp, xionp, iunf,    &
        !                    &         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),    &
        !                    &         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                &
        !                    &         s02, efrozn, et, xnvmu(0,iph),                           &
        !                    &         xnval(1,iph), iorb(-4,iph), norb(iph),                   &
        !                    &         eorb(1,iph), kappa(1,iph) )
        !            endif
        !            ! etfin is absorbing atom final state total energy, see nohole case below.
        !            etfin(itmp) = et
        !         END IF
        !      end do
        ! Now calculate from 29 to ihole
        !     do itmp = 29, ihole, -1
        !KJ     IF(IsOccupied(iz(0),ihole)) THEN
        IF(IsOccupied(0,ihole)) THEN  !KJ 12-2010 
           if (nohole.ge.0) then
              xionp = xion(0)
              if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
              ! Run self-consistent dirac-fock atomic solver.
              call scfdat ( ipr1, nph+1, nph, iz(0), ihole, xionp, iunf,    &
                   &     vcoul(1,nph+1), rho(1,nph+1), dmag(1,nph+1), rhoval(1,nph+1),&
                   &     ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                    &
                   &     s02, efrozn, et, xnvmu(0,nph+1),                             &                           
                   &     xnval(1,nph+1), xnel(1,nph+1), iorb(-4,nph+1), norb(nph+1),  &
                   &     eorb(1,nph+1), kappa(1,nph+1), nqn(:,nph+1) )
           else
              xionp = xion(iph)
              if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
              ! Run self-consistent dirac-fock atomic solver.
              call scfdat ( ipr1, iph, nph, iz(iph), ihole, xionp, iunf,    &
                   &         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),    &
                   &         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                &
                   &         s02, efrozn, et, xnvmu(0,iph),                           &
                   &         xnval(1,iph), xnel(1,iph), iorb(-4,iph), norb(iph),                   &
                   &         eorb(1,iph), kappa(1,iph), nqn(:,iph) )
           endif
           ! etfin is absorbing atom final state total energy, see nohole case below.
           etfin = et
        ELSE
           WRITE(slog,'(I2)') ihole
           CALL Error('No electrons in initial state specified by ihole = ' // slog)
        END IF
        !     end do
        
        ! Calculate for other potentials.
        do iph = 1, nph
           ! Calculate only if cell has an atom in it. Josh Kas
           IF(iz(iph).gt.0) THEN
              ! Write to log.
              write(slog,10) 'free atom potential and density for atom type', iph
              !call wlog(slog)
              
              itmp = 0
              xionp = xion(iph)
              if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
              ! Run self-consistent dirac-fock atomic solver.
              call scfdat ( ipr1, iph, nph, iz(iph), itmp, xionp, iunf,    &
                   &         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),    &
                   &         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                &
                   &         s02, efrozn, et, xnvmu(0,iph),                           &
                   &         xnval(1,iph), xnel(1,iph), iorb(-4,iph), norb(iph),                   &
                   &         eorb(1,iph), kappa(1,iph), nqn(:,iph) )
           END IF
        end do ! End of loop over iph
        
        ! Now, run the absorbing atom again with no core hole.
        write(slog,10) 'initial state energy'
        !call wlog(slog)
        
        !     Save initial state energy and spinors for core hole orbital,
        !     do not save potentials, except for nohole.
        ispinr = ihole
        itmp = 0
        ! If user specified nohole, run with nohole and fill the iph = 0 element of arrays.
        if (nohole.ge.0) then
           iph = 0
           xionp = xion(iph)
           if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
           call scfdat ( ipr1, iph, nph, iz(iph), itmp, xionp, iunf,      &
                &         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),    &
                &         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                &
                &         s02, efrozn, etinit, xnvmu(0,iph),                       &
                &         xnval(1,iph), xnel(1,iph), iorb(-4,iph), norb(iph),                   &
                &         eorb(1,iph), kappa(1,iph), nqn(:,iph) )
           
           ! Otherwise, run with nohole and fill the iph = nph+1 element of arrays.
        else        
           xionp = xion(0)
           if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
           call scfdat ( ipr1, nph+1, nph, iz(0), itmp, xionp, iunf,      &
                &     vcoul(1,nph+1), rho(1,nph+1), dmag(1,nph+1), rhoval(1,nph+1),&
                &     ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,                    &
                &     s02, efrozn, etinit, xnvmu(0,nph+1),                         &
                &     xnval(1,nph+1), xnel(1,nph+1), iorb(-4,nph+1), norb(nph+1),                 &
                &     eorb(1,nph+1), kappa(1,nph+1), nqn(:,nph+1) )
        endif
        
        !KJ 5-2012.  Now that all calls to scfdat (->inmuat->getorb) are finished, output configuration information.
        !            Moved this here so that it includes all core hole / screening / ionicity contributions.
        call DumpConfig2(112,xnel(:,0:nph+1),xnval(:,0:nph+1),nqn(:,0:nph+1),kappa(:,0:nph+1),norb(0:nph+1),nhtmp)
        
        !     testing new potential for the final state. ala
        hx = 0.05
        x0 = -8.8
        if (nohole.gt.0) then
           do i = 1,NRPts
              dr(i) = exp(x0+hx*(i-1))
           end do
           if (nohole.eq.1) then
              do i = 1,NRPts
                 drho(i) = dgc0(i)**2 + dpc0(i)**2
              end do
           else
              do i = 1,NRPts
                 drho(i)=dr(i)**2 * (rho(i,0)-rhoval(i,0)-rho(i,nph+1)+rhoval(i,nph+1))
              end do
           endif
           call potslw ( dvcoul, drho, dr, hx,NRPts)
           do i=1,NRPts
              !           drho(i) = drho(i)/ dr(i)**2
              !           use 1/2 of core-hole as in transition state
              drho(i) = drho(i)/2.0d0/ dr(i)**2
           end do
        else
           do i=1,NRPts
              drho(i) = 0
              dvcoul(i) = 0
           end do
        endif
        
        ! etinit is absorbing atom initial state (no hole)
        ! efrozn is ionization energy with frozen orbitals (koopman's theorem)
        ! etfin-etinit is ionization energy in adiabatic approximation
        ! Debug: Fer
        !    print *, '-efrozn: ', -efrozn
        !    print *, 'etinit: ', etinit
        !    print *, 'etfin: ', etfin
        erelax = -efrozn - ( etfin - etinit)
        emu = etfin - etinit
        ! Debug: Fer
        !    print *, ' emu 1: ', emu
        ! Josh - added check for low energy edges.
        IF(emu.le.0.d0) emu = -efrozn
        ! Debug: Fer
        !    print *, ' emu 2: ', emu
        ! Find norman radius.
        ! Overlap potentials and densitites
        do iph = 0, nph
           write(slog,10)  'overlapped atomic potential and density for unique potential', iph
           call wlog(slog)
! Modified by FDV
! Split line to fix error with Solaris Studio
           call ovrlp (iph, iphat, rat, iatph, novr, iphovr, nnovr, &
                       rovr, iz, nat, rho, dmag, rhoval, vcoul, edens, &
                       edenvl, vclap, rnrmTmp)
           if (iph.eq.0) emu = emu - vclap(1,0)+vcoul(1,0)
           ! Debug: Fer
           !       print *, ' emu 3: ', emu
        end do
        
        if (ifree.eq.1) then
           ! Set the Norman radii if this is the atomic potential with no ionicity.
           rnrm(0:nph) = rnrmTmp(0:nph)
        endif
        
     end do ! End of loop over ifree
     IF(iEdge.EQ.1) THEN
        CALL WriteAtomicPots(nph, iz(0:nph), ihole, rho, dmag(:,0:nph+1), rhoval, vcoul, dgc0,  &
             & dpc0, dgc(:,:,0:nph+1), dpc(:,:,0:nph+1), adgc(:,:,0:nph+1), adpc(:,:,0:nph+1), &
             & erelax, emu, xnvmu, xnval(:,0:nph+1), norb, eorb, drho, dvcoul, iphat,    &
             & rat, iatph(0:nph), novr(0:nph), iphovr, nnovr, rovr, nat, edens, &
             & edenvl, vclap,  rnrm(0:nph), kappa(:,0:nph+1), iorb(:,0:nph+1), s02)
        dgc00 = dgc0
        dpc00 = dpc0
        ihole0 = ihole
        ! Write to edges.dat
        WRITE(17,*) '# emu, M_kk, gam'
        Mkkp = 1.d0
        CALL setgam(iz(0),ihole,gamch)
    
        WRITE(17,*) emu, Mkkp, gamch/hart
     ELSE
        ! Calculate dipole matrix elements between this edge and the first edge.
        CALL setgam(iz(0),ihole,gamch)
        !CALL DipoleMatrixElements(Mkkp,dgc00,dpc00,ihole0,ihole)
        WRITE(17,*) emu, Mkkp, gamch/hart
     END IF
  END DO

  ! Deallocate local variables
  deallocate(edens, edenvl,vclap)
  deallocate(rho)
  deallocate(rnrm, rnrmTmp)
  deallocate(vcoul)
  deallocate(dgc, dpc)
  deallocate(adgc, adpc, rhoval)
  deallocate(dmag, xnvmu)
  deallocate(xnval, eorb)
  deallocate(kappa, iorb, norb)
  ! Set ihole back to ihole0 - probably not needed, but just in case
  ihole = ihole0
     
END SUBROUTINE AtomicPotentials
