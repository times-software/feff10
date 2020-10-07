!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $Revision: 1.0 $
! $Author: tts $
! $Date: 2020/4/1 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_sommerfeld_scf
  ! This module contains a version of the SCF loop-step that
  ! computes the densities using Sommerfeld expansion.
  !
  ! The formulation here follows the path PhysRevB.52.11502 Equation (7)
  ! 1) Convert finite temperature integral to complex contour
  ! 2) Sommerfeld expand leg (2) -> (3)
  ! 3) Sum of poles are approximate as integrals
  !
  ! T=0K result is consistent with scmtmp.f90 (diff Efermi ~ 0.001 eV)
  ! Speed is similar to scmtmp; could be slower since points are computed
  ! dynamically in scmtmp.
  ! To use, include a TEMPERATURE card in feff.inp with the temperature in eV
  !
  ! The calling procedure is:
  !   call thscf_init(ecv, mu)
  !   call thscf_main(...)
  !   call thscf_deinit
  !
  ! The algorithm proceeds as follows:
  !
  ! I. the energy-dependent density is calculated in the imaginary plane along
  !    a contour connecting the points (ecv,0) -> (ecv,eimmax) -> (??, eimmax)
  !    where ecv is the core-valence separation energy and eimmax is height of
  !    contour.
  !
  ! II. Starting with a guess for the chemical potential mu:
  !     A. The density * fermi distribution is integrated over pre-calculated
  !        energies and correction term are added on
  !     B. mu is adjusted using the secant rule until the total number of
  !        electrons is correct
  !
  ! III. The densities are updated from the previous interaction of the SCF
  !      loop using the broydn algorithm
  !
  ! Paralellization notes:
  !
  ! Step I.  will use one processor per energy point
  ! Step II. will use one processor per energy point for leg (3) -> (4)
  !          but is serial in optimizing mu
  ! Step III. is fast and performed by all processors

  use DimsMod, only: nphx=>nphu, lx, nrptx, istatx
  use constants
  use atoms_inp, only: nat, nph, iatph, iphat, rat
  use potential_inp, only: ixc, xnatph, xion, iunf, iz, ihole, lmaxsc, &
    nohole, nscmt, icoul, ca1, rfms1, lfms1, jumprm, scf_temperature, xntol, nmu

  implicit none
  private
  public sommerfeld_scf_main, sommerfeld_scf_init, sommerfeld_scf_deinit
  ! public energies

  real*8 kT ! temperature
  real*8 deim   ! minimum de along imaginary axis
  real*8 mu0  ! initial guess at chemical potential
  real*8 ecv  ! core valence separation

  double precision ri05(251)
  integer :: iter_step ! Track iscmt
  integer :: ne ! number of energy points for legs 1,2,3
  integer :: np ! number of points in leg (3) -> (4)
  complex*16, allocatable :: energies(:)
  complex*16, allocatable :: rhoe(:,:,:) ! LDOS
  complex*16, allocatable :: rhore(:,:,:) ! spatial density
  DOUBLE PRECISION, ALLOCATABLE :: history_xntot(:), history_xmu(:)
  COMPLEX*16, ALLOCATABLE :: rhoe1(:,:,:), rhore1(:,:,:) ! points for d/dx
  DOUBLE PRECISION:: lower_bound, upper_bound
  COMPLEX*16 :: dstep ! for d/dx
CONTAINS

  subroutine sommerfeld_scf_init(ecv0, mu, iscmt)
    implicit none
    !
    ! Input module parameters: kT, mu0, ecv, ne
    !
    integer, intent(in) :: iscmt
    double precision, intent(in) :: ecv0, mu
    integer :: i

    kT = scf_temperature / hart ! convert to hartree
    mu0  = mu
    ecv  = ecv0
    iter_step = iscmt
    do i = 1,251
      ri05(i) = exp (-8.8+0.05*(i-1)) ! radial grid
    end do

    ! ne = 80
    ne = 220
    call sommerfeld_scf_energies_init()
  end subroutine

  subroutine sommerfeld_scf_energies_init()
    ! initialize integration grid, `energies
    ! Grid connects vertices:
    !    1) (ecv,0)
    !    2) (ecv,eimmax)
    !    3) (emax,eimmax)
    !   [4) (mu,0)] !  leg (3)->(4) will be constructed in the main routine
    !
    ! This routine determines values for e1 and ecv, generates the energy grid
    ! and calculates the
    !   np - number of  points  for leg (3) -> (4)
    !   n1 - number of energy points on vertical leg
    !   n2 - number of points in horizontal leg in case of first and second iteration
    implicit none
    real*8  :: e1, emax, demu, previous_energy, de, eimmax
    integer :: n1, n2, i
    logical :: final_leg = .false.
    allocate(history_xntot(nmu), history_xmu(nmu))

    eimmax = 0.15d0 ! 4.05 eV same as grids.f90
    ! eimmax = 0.30d0
    emax = mu0 + 0.5

    e1 = 2*pi*kT
    np = 1
    if (e1 < eimmax) then
        np = ceiling(eimmax / (2*pi*kT))
        e1 = np * 2 * pi * kT
    endif

    n1 = 30
    ! eimmax = e1*n1**2 ! max height similar to scmtmp
    n2 = ne - n1
    if (final_leg) n2 = ne - n1*2

    allocate(energies(ne+np))
    allocate(rhoe(0:lx,0:nphx,ne+np), rhore(251,0:nphx,ne+np))
    ALLOCATE(rhoe1(0:lx,0:nphx,2), rhore1(251,0:nphx,2))

    ! vertical leg (1) -> (2)
    de = e1/n1**2
    do i = 1,n1
      energies(i) = ecv + coni * de * i**2
      if (final_leg) energies(ne + 1 - i) = emax + de* i**2 * coni
    end do

    ! horizontal leg (2) -> (3) Linear grid
    de = (emax - ecv) / n2
    do i = 1, n2
      energies(n1+i) = ecv + i * de + e1 * coni
    end do

    ! save for leg4
    deim = e1

    ! Track search bounds
    lower_bound = DBLE(energies(1))
    upper_bound = DBLE(energies(ne))
  end subroutine

  subroutine sommerfeld_scf_deinit
    deallocate(energies)
    deallocate(rhoe, rhore)
    deallocate(history_xntot, history_xmu)
    DEALLOCATE(rhoe1, rhore1)
  end subroutine

  subroutine sommerfeld_scf_main(iscmt, ecv, vclap, edens, edenvl, vtot, vvalgs, &
    rmt, rnrm, qnrm, rhoint, vint, xmu, xnferm, xnvmu, xnval, x0, ri, dx, &
    adgc, adpc, dgc, dpc, rhoval, xnmues, ok, rgrd)

    !
    ! This routine implements a temperature dependent version of scmt.f90
    !

    ! Parameters

    use par
    implicit none
    integer, intent(in) :: iscmt
    double precision, intent(in) :: ecv, vclap(251,0:nphx)
    double precision, intent(in) :: vtot (251,0:nphx), vvalgs (251,0:nphx)
    real*8, intent(in) :: rmt(0:nphx),rnrm(0:nphx)
    double precision, intent(inout) :: qnrm(0:nphx), edens(251,0:nphx), edenvl(251,0:nphx)
    double precision, intent(in) :: rhoint, vint, xnferm
    double precision, intent(inout) :: xmu
    real*8, intent(inout) :: xnvmu(0:lx,0:nphx+1)
    double precision, intent(in) :: xnval (30,0:nphx)
    double precision, intent(in) :: x0, ri(nrptx), dx
    double precision, intent(in) :: adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
    double precision, intent(in) :: dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
    double precision, intent(inout) :: rhoval(251,0:nphx+1)
    real*8, intent(inout) :: xnmues(0:lx, 0:nphx)
    logical :: ok
    double precision, intent(in) :: rgrd

    ! local vars
    ! integer nmu
    integer ie, ip, iph, imu, ir, il, iproc, converged
    complex*16 ee
    integer nr05(0:nphx)
    double precision xntot, xntotprev, xmunew, xmuprev
    double precision dq(0:nphx)
    double precision diff

    character*512 slog
    complex,    allocatable :: gtr(:,:)
    complex*16, allocatable, dimension(:,:):: xrhoce, xrhole, ph
    complex*16, allocatable :: yrhole(:,:,:)
    allocate(gtr(0:lx, 0:nphx),xrhoce(0:lx,0:nphx))
    allocate(xrhole(0:lx,0:nphx),yrhole(251,0:lx,0:nphx),ph(lx+1, 0:nphx))

    rhoe(:,:,:) = 0.0d0
    rhore(:,:,:) = 0.0d0
    xrhole(:,:) = 0.0d0
    yrhole(:,:,:) = 0.0d0
    rhoval(:,:) = 0.0d0

    write (slog,10) iscmt
    10  format('THERMAL SCF ITERATION NUMBER',i3)
    call wlog(slog)
    write(29,*) trim(slog)

    !call wlog (' Calculating energy and space dependent l-DOS ....')

    xntotprev = 0
    ok = .false.
    converged = 0 ! use integer instead of logical so we can send with MPI

    ! I. calculate densities along energy contour
    CALL calculate_contour(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, &
                           & vint, xnval, x0, ri, dx, adgc, adpc, dgc, dpc, &
                           & xrhole, yrhole, ph, gtr, nr05, iscmt)

    ! II. find mu which gives correct number of electrons
    ! nmu = 100 ! TO DO make this configurable
    xmunew = xmu
    do imu = 1, nmu
      call par_barrier

      ! A. Calculate vertical leg (3) -> (EF) residues
      !    Note that these are stored at the end of the rhoe,rhore arrays
      do ip = 1, np
        ie = ne + ip
        iproc = mod(ip-1, numprocs)
        if (this_process .eq. iproc) then
          ! Step gets smaller closer to efermi
          ! ee = xmunew + de * ((np-ip+1)**2) *coni
          ee = xmunew + coni * pi * kT * (2*ip - 1)
          energies(ie) = ee
          call rhofmslie(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, vint, &
                        xnval, x0, ri, dx, adgc, adpc, dgc, dpc, ie,       &
                        rhoe(:,:,ie), rhore(:,:,ie), xrhole, yrhole, ph,&
                        gtr, ri05, nr05, ee, iscmt)
          ! send result back to master
          if (worker) then
            call par_send_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), 0, ie)
            call par_send_dc(rhore(1,0,ie), 251 * (nphx+1), 0, ie)
            call par_send_dc(energies(ie), 1, 0, ie)
          end if
        end if

        if (master .and. iproc .ne. 0) then
          call par_recv_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), iproc, ie)
          call par_recv_dc(rhore(1,0,ie), 251 * (nphx+1), iproc, ie)
          call par_recv_dc(energies(ie), 1, iproc, ie)
        end if
      end do
      call par_barrier

      ! Set dstep
      if (master) then
        ! PRINT*, "Number of poles = ", np
          ! step size 0.1% of xmunew
        !   dstep = xmunew*0.01d0
        ! call print_poles()
        dstep = 0.001d0/hart ! step size of 0.02 V
        ! dstep = 1e-6
      endif
      call par_bcast_double(dstep, 1, 0)
      call par_barrier

      !    Compute points needed for derivatives
      do ip = 1, 2
          ie = ip
          iproc = mod(ip-1, numprocs)
          if (this_process .eq. iproc) then
              ! Step gets smaller closer to efermi
              if (ie.EQ.1) then
                !   ee = xmunew + coni*DIMAG(energies(ne)) - dstep
                  ee = xmunew + coni*deim
              else
                  ee = xmunew + dstep + coni*deim
                !   ee = xmunew + coni*DIMAG(energies(ne)) + dstep
              endif
              call rhofmslie(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, vint, &
                          xnval, x0, ri, dx, adgc, adpc, dgc, dpc, ie,       &
                          rhoe1(:,:,ie), rhore1(:,:,ie), xrhole, yrhole, ph,&
                          gtr, ri05, nr05, ee, iscmt)
              ! send result back to master
              if (worker) then
                  call par_send_dc(rhoe1(0,0,ie), (lx+1) * (nphx+1), 0, ie)
                  call par_send_dc(rhore1(1,0,ie), 251 * (nphx+1), 0, ie)
              end if
          end if

        if (master .and. iproc .ne. 0) then
          call par_recv_dc(rhoe1(0,0,ie), (lx+1) * (nphx+1), iproc, ie)
          call par_recv_dc(rhore1(1,0,ie), 251 * (nphx+1), iproc, ie)
        end if
      end do
      call par_barrier
      if (master) then
        ! B. Integrate rhoe,rhore over fermi distribution
        call integrate_rhos_lt(nr05,xmunew,xntot,xnmues,rhoval,.FALSE.)
        history_xntot(imu) = xntot
        history_xmu(imu) = xmunew

        ! C. Check for convergence.
        if (abs(xntot - xnferm) .lt. xntol) then
          converged = 1
        else
          ! We will use secant method to explore the landscape
          ! When secant method failed to converge, we will use
          ! bisection method.
          if (imu.LE.10) then
              call update_mu(imu, xntot-xnferm, xntotprev-xnferm, xmunew, xmuprev, xnferm, iscmt) ! updates xmu and xmuprev
          else
              call bracketing_method(imu,xnferm,xmunew,xmuprev)
          endif
          xntotprev = xntot

          if (ABS(xmunew-xmuprev).LT.1.0e-20) then
              ! Usually not a problem
              converged = 1
              PRINT*, "WARNING: Convergence tolerance not reached but chemical potential stopped changing."
              PRINT*, "         Delta xntot is", xntot-xnferm
          endif
        end if
        ! write(slog,*) imu, xmunew*hart, xntot, xnferm
        ! call wlog(slog)
      end if
      call par_barrier
      call par_bcast_int(converged, 1, 0)
      call par_bcast_int(imu, 1, 0)
      call par_bcast_double(xmunew, 1, 0)

      if (converged.eq.1) then
        ok = .true.
        exit
      else
        ! print warning if convergence in mu isn't reached by end of loop
        if (imu.EQ.nmu .AND. master)then
            WRITE(slog,*) "WARNING: Reached end of cycle but not converged"
            CALL wlog(slog)
            call print_history()
        endif
      end if
    end do

    ! the rest of this routine is quick and is performed by all nodes

    call par_bcast_double(xntot, 1, 0)
    call par_bcast_double(xnmues, (lx+1)*(nphx+1), 0)
    call par_bcast_double(rhoval, 251*(nphx+1), 0)

    ! Print out occupation numbers
    ! TO DO (this is taken from scmt.f90 and could probably be cleaned up)
    write(29,*) '  Electronic configuration'
    write(29,*) '  type     l     N_el'
    310  format (2i6, f9.3)
    do ip= 0,nph
      do il = 0,lx
        write (slog,310) ip,il,xnmues(il,ip)
        write(29,*) trim(slog)
        ! check that occupation numbers are consistent with those set in getorb.f
        diff = abs(xnmues(il,ip) - xnvmu(il,ip))
        if (diff.gt.13.1 .or. (il.eq.2 .and. diff.gt. 9.1) .or.        &
          &   (il.eq.1 .and. diff.gt.5.1) .or. (il.eq.0 .and. diff.gt.1.95)) then
          ! the difference in number of atoms is too large
          call wlog (' Found bad counts.')
          write (slog,311) xnvmu(il,ip)
          311       format('  Occupation number in getorb is ', f9.3)
          call wlog(slog)
          call wlog ('  Will repeat this iteration. ')
          if (iscmt.gt.1) ok = .false.
        endif
      end do
    end do

    ! III. if ok, then update output values using Broyden algorithm. otherwise, scf loop will restart
    if (ok) then
      xmu = xmunew
      ! find rhoval via Broyden algorithm
      call broydn( iscmt, ca1, nph, xnvmu, nr05 , xnatph, rnrm, qnrm, edenvl, rhoval, dq) ! xnatph,qrnm,edenvl,rhoval => rhoval,dq

      ! calculate new vclap - overlap coulomb potential
      call coulom (icoul, nph, nr05 , rhoval, edenvl, edens, nat, rat, iatph, iphat, rnrm, dq, iz, vclap) ! everything => vclap

      ! update array edens
      do iph=0,nph
        do ir=1,nr05 (iph)
          edens(ir,iph)=edens(ir,iph)-edenvl(ir,iph)+rhoval(ir,iph)
        end do
        do ir=nr05 (iph)+1,251
          edens(ir,iph)=0.0d0
          edenvl(ir,iph)=0.0d0
        end do
      end do
    end if

    deallocate(gtr,xrhoce)
    deallocate(xrhole,yrhole,ph)
  end subroutine

  subroutine update_mu(imu, dn, dnprev, xmunew, xmuprev, xnferm, iscmt)
    integer, intent(in) :: imu, iscmt
    double precision, intent(in):: dn, dnprev, xnferm
    double precision, intent(inout):: xmunew, xmuprev

    integer :: inrange, ncycle
    double precision xmutmp, xmutmp2, step_size
    character(512) :: slog

    if (imu .eq. 1) then
      ! first step, pick an arbitrary direction to move in
      ! TO DO make the scale of this depend on dn
      xmuprev = xmunew
      if (dn .lt. 0.0d0) then
        xmunew = xmunew + 0.1;
      else
        xmunew = xmunew - 0.1;
      end if
    else
        xmutmp = xmunew
        inrange = 0
        step_size = 1.d0
        ! will always exit when step_size is close to zero
        do ncycle = 1, 4
            xmutmp2 = xmunew - step_size*dn * (xmunew - xmuprev) / (dn - dnprev)
            if (is_in_grid(xmutmp2)) then
               inrange = 1
               exit
            else
              step_size = 0.5d0*step_size
              write(slog,*) "increase step size ", step_size
              call wlog(slog)
            endif
        enddo
        xmunew = xmutmp2
        xmuprev = xmutmp
        ! if (.NOT.is_in_grid(xmunew)) then
        !     write(slog,*) "xmunew is not in grid"
        !     call wlog(slog)
        !     call par_stop
        ! endif
    end if
  end subroutine

  SUBROUTINE integrate_rhos_lt(nr05,xmu,xntot,xnmues,rhoval,print_flag)
    use par
    ! Integrate densities over energy including Fermi distribution
    ! An implementation of 'integrate_rhos' using interpolation around
    ! the chemical potential.
    ! Integrate densities over energy by finding the fermi level
    ! interpolate the original grid up to fermi level
    ! attach the segment (3)->(4) to the end of array, discard e > mu.
    !
    ! Input:
    !   nr05:   indices of Norman radii for each potential type
    !   xmu:    chemical potential
    ! Output:
    !   xntot:  total number of electrons in cluster
    !   xnmues: number of electrons for each potential type and L value
    !   rhoval: valence electron density vs r for each potential type
    LOGICAL, INTENT(IN) :: print_flag
    INTEGER, INTENT(IN) :: nr05(0:nphx)
    DOUBLE PRECISION, INTENT(IN) :: xmu
    DOUBLE PRECISION, INTENT(OUT) :: xntot, xnmues(0:lx, 0:nphx), rhoval(251,0:nphx+1)


    ! Local variables
    !   bug     - If there is a problem (when  xmu_m3 is outside of array bound)
    !   ne_terp - Number of points to interpolate
    !   nxmu    - Number of points in
    !   xmu_m3  - a point below chemical potential
    !   xmu_p3  - a point above chemical potential
    DOUBLE PRECISION ::  xnmues_old(0:lx, 0:nphx)
    DOUBLE PRECISION :: dxmu
    REAL*8 :: sommerfeld_factor
    COMPLEX*16 :: ee, ep, de, xmu_m3, xmu_p3
    COMPLEX*16, ALLOCATABLE :: rhoe_terp(:,:,:), rhore_terp(:,:,:)
    COMPLEX*16, ALLOCATABLE :: energies_terp(:)
    INTEGER :: ne_terp, ind_m3, ind_p3
    INTEGER :: iep, ip, ir, ie, iph, il
    LOGICAL :: out_of_bound = .FALSE.
    CHARACTER(512) message, filename
    INTEGER, PARAMETER :: nxmu = 2 ! Add a few more points


    sommerfeld_factor = pi*pi*kT*kT/6.d0
    xntot = 0
    xnmues = 0
    rhoval = 0

    ! Find where xmu lies on the contour
    ind_m3 = binarysearch(DBLE(xmu), DBLE(energies(1:ne)), ne)-1 ! E(ind_m3) < xmu
    ind_p3 = binarysearch(DBLE(xmu), DBLE(energies(1:ne)), ne)   ! xmu < E(ind_p3)
    dxmu = (DBLE(xmu) - DBLE(energies(ind_m3)))/DBLE(nxmu)
    xmu_m3 = energies(ind_m3)
    xmu_p3 = xmu + coni*DIMAG(energies(11))

    IF (ind_m3.EQ.0) out_of_bound=.TRUE.

    IF (out_of_bound) THEN
      WRITE(message,'(a,3f10.5)') "ERROR: xmu not in range. ", &
                    & DBLE(xmu*hart), DBLE(energies(1)*hart), DBLE(energies(ne)*hart)
      CALL wlog(message)
      ! If xmu below ecv then, we dont do interpolation, we treat
      ! everything below that as zero
      ! JJK - This should really use the energies of the core-states and
      ! Fermi-occupations of fixed energy core-states.
      ! We could actually change the total density accordingly as well.
      ! That is \rho = \rho_val + \rho_core
      !         \rho_val  = \int_ecv^inf f(mu,E,T)ImG/pi
      !         \rho_core = \sum_i |\phi_i|^2 f(mu,E_i,T)
      ! I think the above indicates a serious problem,
      ! so I'm going to put a stop here
      !
      ! TTS - Quick fix. Decrease gradient size if the new xmu is out of range

      call par_stop
    ELSE
      ! Generate a new grid by removing points above chemical potential

      ! Define number of points for legs (1) -> (3)
      ne_terp = ind_m3 + nxmu
      ALLOCATE(rhoe_terp(0:lx,0:nphx,ne_terp+np))
      ALLOCATE(rhore_terp(251,0:nphx,ne_terp+np))
      ALLOCATE(energies_terp(ne_terp+np))
      energies_terp = 0.d0
      rhoe_terp =  0.d0
      rhore_terp = 0.d0

      ! Define number of points for legs (1) -> mu - epsilon
      DO ie=1, ind_m3
        ! Copy vertical points and all points less than xmu
        energies_terp(ie) = energies(ie)
        rhoe_terp(:,:,ie) = rhoe(:,:,ie)
        rhore_terp(:,:,ie) = rhore(:,:,ie)
      ENDDO

      ! Define number of points for mu - epsilon -> mu
      DO ie=ind_m3+1, ind_m3+nxmu
        ! Interpolate from [ e<xmu, xmu ]
        energies_terp(ie) = energies(ind_m3) + (ie-ind_m3) * dxmu

        DO iph=0, nphx
          DO il = 0, lx
            rhoe_terp(il,iph,ie) = interp1d(DBLE(energies_terp(ie)), DBLE(energies(1:ne)), rhoe(il,iph,1:ne))
          ENDDO
          DO il = 1,251
            rhore_terp(il,iph,ie) = interp1d(DBLE(energies_terp(ie)), DBLE(energies(1:ne)), rhore(il,iph,1:ne))
          ENDDO
        ENDDO
        ! write(message,*) DBLE(energies_terp(ie)), DBLE(rhoe_terp(2,0,ie)), DIMAG(rhoe_terp(2,0,ie))
        ! call wlog(message)
      ENDDO

      ! Define number of points for leg(3) -> leg(4)
      DO ie = 1, np
        ! Attach leg (3) -> (4)
        energies_terp(ne_terp+ie) = energies(ne+ie)
        rhoe_terp(:,:,ne_terp+ie) = rhoe(:,:,ne+ie)
        rhore_terp(:,:,ne_terp+ie) = rhore(:,:,ne+ie)
      ENDDO
    ENDIF

    ! I. Integrate density over energy
    !    The first calculated point has finite imaginary part.
    !    So, include contribution from leg between real axis and
    !    this point by approximating the integrand as constant
    !    along that leg (this is the same as what is done in scmt.f90).
    ep = DBLE(energies(1))
    iep = 1

    ! Smooth integral part
    DO ie = 1, ne_terp
      ee = energies_terp(ie)
      de = ee - ep

      DO iph = 0, nph
        DO il = 0, lx
          if (iunf.EQ.0 .AND. il.GT.2) CYCLE
          xnmues(il,iph) = xnmues(il,iph) + dimag( (rhoe_terp(il,iph,ie) + rhoe_terp(il,iph,iep)) * de )
        END DO

        DO ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( (rhore_terp(ir,iph,ie) + rhore_terp(ir,iph,iep)) * de )
        END DO
      END DO

      ep = ee
      iep = ie
    END DO

    ! Sum of residues
    DO ip = 1, np
      DO iph = 0, nph
        DO il = 0, lx
          if (iunf.EQ.0 .AND. il.GT.2) CYCLE
          xnmues(il, iph) = xnmues(il, iph) + dimag( -4.d0*pi*coni*kT * rhoe_terp(il,iph,ne_terp+ip) )
        END DO
        DO ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( -4.d0*pi*coni*kT * rhore_terp(ir,iph,ne_terp+ip) )
        END DO
      END DO
    END DO


    ! II. Compute the derivative at z = mu + i*eimmax
    !      Evaluate the derivative by linear interpolation
    !      G'(z) = slope between point energies(ind_m3) and energies(ind_p3)
    ! We need to compute the derivatives for rhoe, rhore

    !! de = energies(ind_p3)-energies(ind_m3)
    de = dstep
    DO iph = 0, nph
      DO il = 0, lx
        if (iunf.EQ.0 .AND. il.GT.2) CYCLE
        ! Factor of 2 for degeneracy
        xnmues(il,iph) = xnmues(il,iph) + 2.d0*sommerfeld_factor*DIMAG(&
                            (rhoe1(il,iph,2)-rhoe1(il,iph,1))/de)
      END DO

      DO ir = 1, nr05(iph)
        rhoval(ir, iph) = rhoval(ir, iph) + 2.d0*sommerfeld_factor*DIMAG( &
                            (rhore1(ir,iph,2)-rhore1(ir,iph,1))/de)
      END DO
    END DO

    ! III. Sum over angular momenta and sites to get total electron count
    DO iph = 0, nph
      DO il = 0, lx
        xntot = xntot + xnmues(il, iph) * xnatph(iph)
      END DO
    END DO
    IF (print_flag) PRINT*, "    xntot = ", xntot
    ! call print_rhoe(energies_terp, rhoe_terp, ne_terp)
    ! call print_grid(energies_terp, ne_terp)
    DEALLOCATE(rhoe_terp, rhore_terp)
    DEALLOCATE(energies_terp)
  END SUBROUTINE

  SUBROUTINE calculate_contour(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint,&
                              & vint, xnval, x0, ri, dx, adgc, adpc, dgc, dpc,&
                              & xrhole, yrhole, ph, gtr, nr05, iscmt)
    ! Calculate densities along energy contour.
    ! Note: Original version was hard coded in "thscf_main" subroutine.
    !       It has been moved into a separate subroutine.
    ! Import mpi parameters
    USE PAR
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: edens(251,0:nphx), edenvl(251,0:nphx)
    DOUBLE PRECISION, INTENT(IN) :: vtot (251,0:nphx), vvalgs (251,0:nphx)
    REAL*8, INTENT(IN) :: rmt(0:nphx),rnrm(0:nphx)
    DOUBLE PRECISION, INTENT(IN) :: rhoint, vint, xnval (30,0:nphx)
    DOUBLE PRECISION, INTENT(IN) :: x0, ri(nrptx), dx
    DOUBLE PRECISION, INTENT(IN) :: adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
    DOUBLE PRECISION, INTENT(IN) :: dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
    COMPLEX*16, INTENT(INOUT) :: xrhole(:,:), yrhole(:,:,:)
    COMPLEX*16, INTENT(INOUT) :: ph(:, :)
    COMPLEX, INTENT(INOUT) :: gtr(:, :)
    INTEGER, INTENT(IN) :: iscmt, nr05(0:nphx)

    ! Local variables
    COMPLEX*16 :: ee
    INTEGER :: ie, ip, iph, iproc, il, ir
    CHARACTER*512 :: slog

    DO ie = 1,ne
      ! Spread energy points evenly across processors
      iproc = mod(ie-1, numprocs)
      if (this_process .eq. iproc) then
        ee = energies(ie)
        call rhofmslie(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, vint, &
                      xnval, x0, ri, dx, adgc, adpc, dgc, dpc, ie,          &
                      rhoe(:,:,ie), rhore(:,:,ie), xrhole, yrhole, ph, gtr, &
                      ri05, nr05, ee, iscmt)

        ! send result back to master
        if (worker) then
          call par_send_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), 0, ie)
          call par_send_dc(rhore(1,0,ie), 251 * (nphx+1), 0, ie)
        end if

        ! periodically print out progress update
        if (ie.eq.1 .or. mod(ie,20).eq.0) then
          if (worker) par_type = 3
          write(slog,30) ie, dble(ee)*hart
          30    format('     point # ', i3, '  energy = ', f7.3)
          call wlog(slog)
          if (worker) par_type = 2
        endif
      end if

      if (master .and. iproc .ne. 0) then
        call par_recv_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), iproc, ie)
        call par_recv_dc(rhore(1,0,ie), 251 * (nphx+1), iproc, ie)
      end if
    END DO
    CALL par_barrier
  END SUBROUTINE

  COMPLEX*16 FUNCTION interp1d(x0, x, y)
    ! Linear interpolate (x,y) at point x0
    IMPLICIT none
    INTEGER :: n, idx
    REAL(8), INTENT(IN) :: x0, x(:)
    COMPLEX*16, INTENT(IN) :: y(:)
    COMPLEX*16 :: y1, y2
    REAL(8) :: x1, x2
    n = size(x)
    if ((x(1).LT.x0).and.(x0.LT.x(n))) then
      idx = binarysearch(x0,x,n)
      y2 = y(idx)
      x2 = x(idx)
      y1 = y(idx-1)
      x1 = x(idx-1)
      interp1d = (y2-y1)*(x0-x1)/(x2-x1) + y1
    else
      IF (x(1).EQ.x0) THEN
        interp1d = y(1)
      ELSE IF (ABS(x(n)-x0).LT.1e-10) THEN
        interp1d = y(n)
      ELSE
        PRINT*, "ERROR: Out of range !",x(1)*hart,"<=",x0*hart,"<=",x(n)*hart
        call par_stop
      ENDIF
    endif
    RETURN
  END FUNCTION interp1d

  INTEGER FUNCTION binarysearch(x0, x, n)
    ! Binary Insert: Return index i such that a[i-1] < x0 <= a[i]
    IMPLICIT none
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x0, x(n)

    INTEGER :: L, R, M, P

    IF (n.LE.0) goto 1234

    L = 0
    R = n-1

    ! Put End
    IF (x0.GE.x(n)) THEN
      binarysearch = n+1
      GOTO 222
    ENDIF

    ! Put Front
    IF (x0.LE.x(1)) THEN
      binarysearch = 1
      GOTO 222
    ENDIF

    111 M = (R+L)/2

    IF (M.LE.0) GOTO 1234

    IF(M.EQ.L) THEN
      binarysearch = M + 1
      GOTO 222
    ENDIF

    IF (x0.LT.x(M)) THEN
      GOTO 123 !LEFT
    ELSE IF (x0.GT.x(M)) THEN
      GOTO 321 !RIGHT
    ELSE
      binarysearch = M ! FOUND
      GOTO 222
    ENDIF

    ! LEFT
    123 IF (x0.EQ.x(M)) THEN
          binarysearch = M
          GOTO 222
        ELSE
          R = M
          GOTO 111 ! Divide again
        ENDIF


    ! RIGHT
    321 IF (x0.EQ.x(R)) THEN
          binarysearch = R
          GOTO 222
        ELSE
          L = M
          GOTO 111 ! Divide again
        ENDIF


    1234 binarysearch = -1 ! Not found or too short
    PRINT*, "ERROR: Binary search failed"
    222 return
  END FUNCTION binarysearch

  SUBROUTINE bracketing_method(current_step, xnferm, xmunew, xmuprev)
      ! This subroutine sets the search interval for xmu
      ! based on previous values.
      ! Note: current_step must be > 1
      implicit none
      integer, intent(in) :: current_step
      double precision, intent(in) :: xnferm
      real*8, intent(inout) :: xmunew, xmuprev

      ! Local variable
      integer :: iOrder(current_step), i_left, i_right, imu, ip
      real*8 :: sorted_xmu(current_step), sorted_xntot(current_step)
      character(512) :: slog

      ! quad precision for arithematics
      real*8 :: aa, bb, cc, dd, xnferm_precise

      sorted_xmu(:) = -100.d0
      sorted_xntot(:) = -1.d0
      CALL qsortd(iOrder, current_step, history_xmu(1:current_step))

      do imu = 1, current_step
          ip = iOrder(imu)
          sorted_xmu(imu) = history_xmu(ip)
          sorted_xntot(imu) = history_xntot(ip)
      enddo
      ! open(file="history.dat", status="replace", unit=1818)
      ! write(1818,*) imu, sorted_xmu(imu), sorted_xntot(imu)
      ! close(1818)

      ! Find the search bounds for xmu
      ! The integrated dos is strictly increasing
      ! Picked an off centered search to avoid having xmu on
      ! the bounds.
      xnferm_precise = xnferm
      i_left = binarysearch(xnferm_precise, sorted_xntot, current_step)-1
      i_right = i_left+1
      if (i_right.GT.current_step) then
          write(slog,*) "Bracketing error: i_right out of bound", i_right
          call wlog(slog)
      endif
      lower_bound = sorted_xmu(i_left)
      upper_bound = sorted_xmu(i_right)
      if ((sorted_xntot(i_right)-xnferm)*(sorted_xntot(i_left)-xnferm).GT.0) then
          print*, "Error points are not bounded", i_left, i_right
          call par_stop
      endif
      xmuprev = xmunew

      ! Regula falsi
    !   if (imu.LT.40) then
        !   xmunew = (sorted_xntot(i_right)-xnferm)*lower_bound - (sorted_xntot(i_left)-xnferm)*upper_bound
        !   xmunew = xmunew / (sorted_xntot(i_right) - sorted_xntot(i_left))
    !   else
    !       ! Bisection method
      aa = lower_bound
      bb = upper_bound
      cc = (bb+aa)*0.5d0
      xmunew = cc
  END SUBROUTINE

  SUBROUTINE print_history()
        ! This subroutine sets the search interval for xmu
        ! based on previous values.
        ! Note: current_step must be > 1
        implicit none

        ! Local variable
        integer :: iOrder(nmu), i_left, i_right, imu, ip
        double precision :: sorted_xmu(nmu), sorted_xntot(nmu)
        character(512) :: slog

        sorted_xmu(:) = -100.d0
        sorted_xntot(:) = -1.d0
        CALL qsortd(iOrder, nmu, history_xmu)

        open(file="scf_history_sorted.dat", status="replace", unit=1818)
        open(file="scf_history.dat", status="replace", unit=1919)
        do imu = 1, nmu
            ip = iOrder(imu)
            sorted_xmu(imu) = history_xmu(ip)
            sorted_xntot(imu) = history_xntot(ip)
            write(1818,100) imu, sorted_xmu(imu), sorted_xmu(imu)*hart, sorted_xntot(imu)
            write(1919,100) imu, history_xmu(imu), history_xmu(imu)*hart, history_xntot(imu)
        enddo
        100 format (i4, 1x, f30.22, 1x, f30.22, 1x, f30.22)
        close(1818)
        close(1919)
  END SUBROUTINE

  SUBROUTINE print_grid(energies_terp, ne_terp)
    IMPLICIT none
    INTEGER, INTENT(IN) :: ne_terp
    COMPLEX*16, INTENT(IN) :: energies_terp(ne_terp+np)
    INTEGER :: ie

    ! To generate animation for grid
    OPEN(UNIT=9991, FILE="contour.sommerfeld", STATUS="OLD",action='write',position='append')
    DO ie = 1, ne_terp+np
        WRITE(9991,*) DBLE(energies_terp(ie)), DIMAG(energies_terp(ie))
    ENDDO
    WRITE(9991,*)
    WRITE(9991,*)
    CLOSE(9991)
  END SUBROUTINE

  SUBROUTINE print_poles()
      IMPLICIT none
      INTEGER :: ip, iph, il

      DOUBLE PRECISION :: xn_arr(0:lx, 0:nphx), xntot

      OPEN(UNIT=9991, FILE="poles.dat", STATUS="OLD",action='write',position='append')

      DO ip = 1, np
        xn_arr(:,:) = 0.d0
        xntot = 0.d0
        DO iph = 0, nph
          DO il = 0, lx
            if (iunf.EQ.0 .AND. il.GT.2) CYCLE
            xn_arr(il, iph) = xn_arr(il, iph) + dimag( -4.d0*pi*coni*kT * rhoe(il,iph,ne+ip) )
          END DO
        END DO
        DO iph = 0, nph
          DO il = 0, lx
            xntot = xntot + xn_arr(il, iph) * xnatph(iph)
          END DO
        END DO
        WRITE(9991,100) ip, DBLE(energies(ne+ip)), DIMAG(energies(ne+ip)), xntot
      END DO
      100 format(i4, 1x, f30.25, 1x, f30.25, 1x, f30.25)
      WRITE(9991,*)
      WRITE(9991,*)
      CLOSE(9991)
  END SUBROUTINE

  SUBROUTINE print_rhoe(energies_terp, rhoe_terp, ne_terp)
      IMPLICIT none
      INTEGER, INTENT(IN) :: ne_terp
      COMPLEX*16, INTENT(IN) :: energies_terp(ne_terp+np)
      COMPLEX*16, INTENT(IN) :: rhoe_terp(0:lx,0:nphx,ne_terp+np)
      INTEGER :: ie
      COMPLEX*16 :: aa, bb, cc, ee

      OPEN(UNIT=9991, FILE="rhoe.sommerfeld", STATUS="OLD",action='write',position='append')
      OPEN(UNIT=99911, FILE="rhoe.real.sommerfeld", STATUS="OLD",action='write',position='append')
      OPEN(UNIT=9981, FILE="rhoe.xmu.sommerfeld", STATUS="OLD",action='write',position='append')
      OPEN(UNIT=99811, FILE="rhoe.xmu.real.sommerfeld", STATUS="OLD",action='write',position='append')
      OPEN(UNIT=777, FILE="rhoe.div.sommerfeld", STATUS="OLD",action="write",position='append')
      DO ie = 1, ne
          WRITE(9991,1818) DBLE(energies(ie)), DIMAG(energies(ie)), dimag(rhoe(0,0,ie)), dimag(rhoe(1,0,ie)), dimag(rhoe(2,0,ie))
          WRITE(99911,1818) DBLE(energies(ie)), DIMAG(energies(ie)), DBLE(rhoe(0,0,ie)), DBLE(rhoe(1,0,ie)), DBLE(rhoe(2,0,ie))
      ENDDO

      DO ie = 1, ne_terp
          WRITE(9981,1818) DBLE(energies_terp(ie)), DIMAG(energies_terp(ie)), dimag(rhoe_terp(0,0,ie)), dimag(rhoe_terp(1,0,ie)), dimag(rhoe_terp(2,0,ie))
          WRITE(99811,1818) DBLE(energies_terp(ie)), DIMAG(energies_terp(ie)), DBLE(rhoe_terp(0,0,ie)), DBLE(rhoe_terp(1,0,ie)), DBLE(rhoe_terp(2,0,ie))
      ENDDO
      DO ie = 10, ne-1
          bb = interp1d(DBLE(energies(ie)), &
                & DBLE(energies(10:ne)), rhoe(0,0,10:ne)) &
                & + interp1d(DBLE(energies(ie)), &
                & DBLE(energies(10:ne)), rhoe(1,0,10:ne)) &
                & + interp1d(DBLE(energies(ie)), &
                & DBLE(energies(10:ne)), rhoe(2,0,10:ne))
          aa = interp1d(DBLE(energies(ie+1)), &
                 & DBLE(energies(10:ne)), rhoe(0,0,10:ne)) &
                 & + interp1d(DBLE(energies(ie+1)), &
                 & DBLE(energies(10:ne)), rhoe(1,0,10:ne)) &
                 & + interp1d(DBLE(energies(ie+1)), &
                 & DBLE(energies(10:ne)), rhoe(2,0,10:ne))
          ee = energies(ie+1) - energies(ie)
          cc = (bb-aa)/ee
          WRITE(777,1818) DBLE(energies(ie)), DIMAG(energies(ie)), DBLE(cc), DIMAG(cc)
      ENDDO
      1818 format (f10.5, 1x, 4(1pe13.6,1x))
      WRITE(9991,*)
      WRITE(9991,*)
      WRITE(99911,*)
      WRITE(99911,*)
      WRITE(99811,*)
      WRITE(99811,*)
      WRITE(9981,*)
      WRITE(9981,*)
      WRITE(777,*)
      WRITE(777,*)

      CLOSE(99811)
      CLOSE(99911)
      CLOSE(9991)
      CLOSE(9981)
      CLOSE(777)
  END SUBROUTINE

  LOGICAL FUNCTION is_in_grid(xmu)
      ! Check if xmu is in the grid generated
      implicit none
      real*8, intent(in) :: xmu
      if (xmu.gt.DBLE(energies(1)) .AND. xmu.lt.DBLE(energies(ne))) then
          is_in_grid = .TRUE.
      else
          is_in_grid = .FALSE.
      endif
      return
  END FUNCTION
END MODULE
