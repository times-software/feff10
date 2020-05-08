!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $Revision: 1.0 $
! $Author: tts $
! $Date: 2020/4/1 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_thermal_scf
  ! Note: Thermal_scf with bad count fix
  !     If experience problem, try using more atoms. It will converge.
  !
  ! This module contains a version of the SCF loop-step that uses a finite temperature
  ! Fermi distribution for the valence electron occupation.
  !
  ! To use, include a TEMPERATURE card in feff.inp with the temperature in eV
  !
  ! The calling procedure is:
  !   call thscf_init(ecv, mu)
  !   call thscf_main(...)
  !   call thscf_deinit
  !
  ! The algorithm proceeds as follows:
  !
  ! I. the energy-dependent density is calculated in the imaginary plane along a contour connecting the points (ecv,0) -> (ecv,e1) -> (mu + 10 * kT, e1) where ecv is the core-valence separation energy and e1 is currently the smallest multiple of 2*pi*kT that is larger than 4 eV.
  !
  ! II. Starting with a guess for the chemical potential mu:
  !   A. The density * fermi distribution is integrated over pre-calculated energies and
  !      the contributions from the Matsubara poles are added on
  !   B. mu is adjusted using the secant rule until the total number of electrons is correct
  !
  ! III. The densities are updated from the previous interaction of the SCF loop using the broydn algorithm
  !
  ! Paralellization notes:
  !
  ! Step I. will use one processor per energy point
  ! Step II. will use one processor per Matsubara pole, but is serial in optimizing mu
  ! Step III. is fast and performed by all processors

  use DimsMod, only: nphx=>nphu, lx, nrptx, istatx
  use constants

  use atoms_inp, only: nat, nph, iatph, iphat, rat
  use potential_inp, only: ixc, xnatph, xion, iunf, iz, ihole, lmaxsc, &
    nohole, nscmt, icoul, ca1, rfms1, lfms1, jumprm, scf_temperature

  implicit none

  private

  public thscf_init, thscf_deinit, thscf_main, thscf_debug
  public fermi_distribution
  public energies

  real*8 temp ! temperature
  real*8 mu0  ! initial guess at chemical potential
  real*8 ecv  ! core valence separation

  double precision ri05(251)

  INTEGER :: window_size_terp ! region around chemical potential to interpolate
  integer :: ne ! number of energy points
  integer :: iter_step ! step number
  integer :: np ! number of Matsubara poles within contour
  complex*16, allocatable :: energies(:)
  complex*16, allocatable :: rhoe(:,:,:) ! LDOS
  complex*16, allocatable :: rhore(:,:,:) ! spatial density
  DOUBLE PRECISION :: history(100), historyE(100)
CONTAINS

  subroutine thscf_init(ecv0, mu, iscmt)
    integer i, iter_step
    integer iscmt
    double precision ecv0, mu

    temp = scf_temperature / hart
    mu0  = mu
    ecv  = ecv0
    iter_step = iscmt
    do i = 1,251
      ri05(i) = exp (-8.8+0.05*(i-1))
    end do

    ! this should be dependent on the temperature
    ! added multiplication by 2 then temperatur is less then 1 eV
    ! Was needed for SiC calculations, SCF was failing due to bad counts
    ! for low temperatures, the fermi distribution is sharp, and fine spacing is
    ! required around mu
    ! for high temperatures, the fermi distribution and the DOS (at i*2*pi*kT)
    ! are both quite smooth, and only a smal number of points are needed

    ne = 320 ! Default: 200
    ! if (scf_temperature <= 1.0) then
    !   ne = ne*2
    ! end if
    call thscf_energies_init()
  end subroutine

  subroutine thscf_energies_init()
    ! initialize integration grid, `energies
    !
    ! Input module parameters: temp, mu0, ecv, ne
    !
    ! Grid connects vertices:
    !    1) (ecv,0)
    !    2) (ecv,e1)
    !    3) (emax,e1)
    !   [4) (emax,0)] ! The fermi function vanishes on leg (3)->(4), so it is not included
    !
    ! This routine determines values for e1 and ecv, generates the energy grid and calculates the
    ! number of Matsubara poles (np) lying within the contour
    ! n1 number of energy points on vertical leg
    ! n2 number of points in horizontal leg in case of first and second iteration
    ! n3 number of points around chemical potential.
    real*8  :: e1, emax, de, demu, previous_energy
    integer :: n1, n2, n3, i
    logical :: final_leg = .false.

    ! Determine imaginary part based on temperature (this should be configurable / overridable)
    ! For now, use the smallest multiple of 2*pi*kT that is greater than 0.15 Hartree (~4 eV)
    ! Modified: 1. window_size_terp: Only +/-12kT is important
    !           2. emax: 2.5 Hartree should be large enough for Mu to oscillate

    window_size_terp = 12
    emax = mu0 + 2.5

    e1 = 2 * pi * temp
    np = 1
    if (e1 < 0.15) then
      np = ceiling(0.15 / (2*pi*temp))
      e1 = np * 2*pi*temp
    endif

    n1 = 10
    n2 = ne - n1
    if (final_leg) n2 = ne - n1*2

    allocate( energies(ne+np) )
    allocate( rhoe(0:lx,0:nphx,ne+np), rhore(251,0:nphx,ne+np) )

    ! vertical leg (1) -> (2)
    de = e1 / n1**2
    do i = 1,n1
      energies(i) = ecv + de * i**2 * coni
      if (final_leg) energies(ne + 1 - i) = emax + de * i**2 * coni
    end do

    ! horizontal leg (2) -> (3) Linear grid
    de = (emax - ecv) / n2
    do i = 1, n2
      energies(n1+i) = ecv + i * de + e1 * coni
    end do
  end subroutine

  subroutine thscf_deinit
    deallocate(energies)
    deallocate(rhoe, rhore)
  end subroutine

  subroutine thscf_main(iscmt, ecv, vclap, edens, edenvl, vtot, vvalgs, &
    rmt, rnrm, qnrm, rhoint, vint, xmu, xnferm, xnvmu, xnval, x0, ri, dx, &
    adgc, adpc, dgc, dpc, rhoval, xnmues, ok, rgrd)

    !
    ! This routine implements a temperature dependent version of scmt.f90
    !

    ! Parameters

    use par

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
    integer ie, ip, iph, imu, nmu, ir, il, iproc, converged
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
    nmu = 100 ! TO DO make this configurable
    xmunew = xmu
    do imu = 1,nmu
      call par_barrier

      ! A. Calculate contributions from Matsubara poles
      !    Note that these are stored at the end of the rhoe,rhore arrays
      do ip = 1, np
        ie = ne + ip
        iproc = mod(ie-1, numprocs)

        if (this_process .eq. iproc) then
          ee = xmunew + coni * pi * temp * (2*ip - 1)
          energies(ie) = ee

          call rhofmslie(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, vint, &
                        xnval, x0, ri, dx, adgc, adpc, dgc, dpc, ne+ip,       &
                        rhoe(:,:,ie), rhore(:,:,ie), xrhole, yrhole, ph,&
                        gtr, ri05, nr05, ee, iscmt)
          ! send result back to master
          if (worker) then
            call par_send_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), 0, ie)
            call par_send_dc(rhore(1,0,ie), 251 * (nphx+1), 0, ie)
          end if
        end if

        if (master .and. iproc .ne. 0) then
          call par_recv_dc(rhoe(0,0,ie), (lx+1) * (nphx+1), iproc, ie)
          call par_recv_dc(rhore(1,0,ie), 251 * (nphx+1), iproc, ie)
        end if
      end do
      call par_barrier

      if (master) then
        ! B. Integrate rhoe,rhore over fermi distribution
        call integrate_rhos_lt(nr05,xmunew,xntot,xnmues,rhoval,.FALSE.)
        history(imu) = xntot-xnferm
        historyE(imu) = xmunew*hart

        ! C. Check for convergence. (TODO make tolerance configurable)
        if (abs(xntot - xnferm) .lt. 1e-4) then
          converged = 1
        else
          ! If not converged, adjust mu appropriately and repeat
          call update_mu(imu, xntot-xnferm, xntotprev-xnferm, xmunew, xmuprev, iscmt) ! updates xmu and xmuprev
          xntotprev = xntot
        end if
      end if

      call par_bcast_int(converged, 1, 0)
      call par_bcast_int(imu, 1, 0)
      call par_bcast_double(xmunew, 1, 0)

      if (converged.eq.1) then
        ok = .true.
        exit
      end if
      ! TO DO print warning if convergence in mu isn't reached by end of loop
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

  subroutine update_mu(imu, dn, dnprev, xmunew, xmuprev, iscmt)
    integer, intent(in) :: imu, iscmt
    double precision, intent(in):: dn, dnprev
    double precision, intent(inout):: xmunew, xmuprev

    double precision xmutmp

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
      ! use secant method
      xmutmp = xmunew
      xmunew = xmunew - dn * (xmunew - xmuprev) / (dn - dnprev)
      xmuprev = xmutmp

    end if
  end subroutine

  subroutine integrate_rhos(nr05,xmu,xntot,xnmues,rhoval)
    ! Integrate densities over energy including Fermi distribution

    ! Input:
    !   nr05:   indices of Norman radii for each potential type
    !   xmu:    chemical potential
    ! Output:
    !   xntot:  total number of electrons in cluster
    !   xnmues: number of electrons for each potential type and L value
    !   rhoval: valence electron density vs r for each potential type

    integer, intent(in) :: nr05(0:nphx)
    double precision, intent(in) :: xmu
    double precision, intent(out) :: xntot, xnmues(0:lx, 0:nphx), rhoval(251,0:nphx+1)
    double precision ::  xnmues_old(0:lx, 0:nphx)
    complex*16 ee, ep, de, f, fp
    integer iph, il, ir, ie, iep, ip
    xntot = 0
    xnmues(:,:) = 0
    rhoval(:,:) = 0

    ! I. Integrate density over energy, fermi distribution
    !    The first calculated point has finite imaginary part. So, include contribution from
    !    leg between real axis and this point by approximating the integrand as constant
    !    along that leg (this is the same as what is done in scmt.f90).
    ep = dble(energies(1))
    iep = 1

    open(unit=99,file='rhoe_th.dat')
    open(unit=899,file='fermi_th.dat')
    do ie = 1, ne+np
      f = fermiF(energies(ie), temp, xmu)

      if (ie > ne) f = -2*pi*coni
      ! write(99,'(i4,20e14.8)') ie, energies(ie), rhoe(0,0,ie)*f, rhoe(1,0,ie)*f, rhoe(2,0,ie)*f
      write(899,'(20E20.10E3)') DBLE(energies(ie)), DIMAG(energies(ie)), DBLE(f), DIMAG(f)
      write(99,'(20E20.10E3)') DBLE(energies(ie)), DIMAG(energies(ie)),  DBLE(rhoe(0,0,ie)), DIMAG(rhoe(0,0,ie)), &
                                                    & DBLE(rhoe(1,0,ie)), DIMAG(rhoe(1,0,ie)), DBLE(rhoe(2,0,ie)), DIMAG(rhoe(2,0,ie))
    end do
    write(99,*)
    write(899,*)
    close(99)
    close(899)


    do ie = 1, ne
      ee = energies(ie)
      de = ee - ep

      f = fermiF(ee, temp, xmu)
      fp = fermiF(ep, temp, xmu)

      do iph = 0, nph
        do il = 0, lx
          if (iunf.eq.0 .and. il.gt.2) cycle
          xnmues(il,iph) = xnmues(il,iph) + dimag( (rhoe(il,iph,ie)*f + rhoe(il,iph,iep)*fp) * de )
        end do

        do ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( (rhore(ir,iph,ie)*f + rhore(ir,iph,iep)*fp) * de )
        end do
      end do

      ep = ee
      iep = ie
    end do
    ! xnmues_old = xnmues

    ! II. Add on contribution from Matsubara poles
    do ip = 1, np
      do iph = 0, nph
        do il = 0, lx
          if (iunf.eq.0 .and. il.gt.2) cycle
          xnmues(il, iph) = xnmues(il, iph) + dimag( -4*pi*coni*temp * rhoe(il,iph,ne+ip) )
        end do
        do ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( -4*pi*coni*temp * rhore(ir,iph,ne+ip) )
        end do
      end do
    end do
    ! WRITE(299,'(20F10.5)') xnmues_old
    ! WRITE(299,'(20F20.10)') SUM(xnmues_old(:,0)), SUM(xnmues_old(:,1))
    ! WRITE(299,'(20F20.10)') SUM(xnmues(:,0)-xnmues_old(:,0)), SUM(xnmues(:,1)-xnmues_old(:,1))
    ! Sum over angular momenta and sites to get total electron count
    do iph = 0, nph
      do il = 0, lx
        xntot = xntot + xnmues(il, iph) * xnatph(iph)
      end do
    end do
  end subroutine

  subroutine thscf_debug()
    integer i
    do i = 1,ne
    end do
  end subroutine

  function fermi_distribution(ee, mu)
    complex*16 ee, fermi_distribution
    double precision mu

    fermi_distribution = 1.0d0 / (exp((ee - mu) / temp) + 1.0d0)
  end function

  SUBROUTINE integrate_rhos_lt(nr05,xmu,xntot,xnmues,rhoval,print_flag)
    use par
    ! Integrate densities over energy including Fermi distribution
    ! An implementation of 'integrate_rhos' using interpolation around
    ! the chemical potential.
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
    !   nxmu    - Number of points in window_size_terp
    !   dx      - window_size_terp * kT
    !   xmu_m3  - chemical potential - dx
    !   xmu_p3  - chemical potential + dx
    DOUBLE PRECISION ::  xnmues_old(0:lx, 0:nphx)
    DOUBLE PRECISION :: dxmu, dx
    COMPLEX*16 :: ee, ep, de, f, fp
    COMPLEX*16 :: xmu_m3, xmu_p3
    COMPLEX*16, ALLOCATABLE :: rhoe_terp(:,:,:), rhore_terp(:,:,:), energies_terp(:)
    INTEGER :: ne_terp
    INTEGER :: iep, ip, ir
    INTEGER :: ind_m3, ind_p3, ie, iph, il
    LOGICAL :: out_of_bound = .FALSE.
    CHARACTER(512) message
    INTEGER, PARAMETER :: nxmu = 1000 ! 2000 not working for 6kT


    xntot = 0
    xnmues = 0
    rhoval = 0
    dx = window_size_terp*temp
    dxmu = (2*dx)/nxmu
    xmu_m3 = xmu - dx + coni*DIMAG(energies(11))
    xmu_p3 = xmu + dx + coni*DIMAG(energies(11))

    ! Find where does xmu +/- window_terp lies on the contour
    ind_m3 = binarysearch(DBLE(xmu_m3), DBLE(energies(1:ne)), ne)-1
    ind_p3 = binarysearch(DBLE(xmu_p3), DBLE(energies(1:ne)), ne)
    IF (ind_m3.EQ.0) out_of_bound=.TRUE.

    IF (out_of_bound) THEN
      write(message,'(a,2f10.5)') "ERROR: xmu_m3 less than minimum. ", & 
                    & DBLE(xmu_m3*hart), DBLE(energies(1)*hart)
      CALL wlog(message)
      ! If xmu below ecv then, we dont do interpolation, we treat
      ! everything below that as zero
      ! JJK - This should really use the energies of the core-states and 
      ! Fermi-occupations of fixed energy core-states.
      ! We could actually change the total density accordingly as well. 
      ! That is \rho = \rho_val + \rho_core
      ! \rho_val  = \int_ecv^inf f(mu,E,T)ImG/pi
      ! \rho_core = \sum_i |\phi_i|^2 f(mu,E_i,T)
      ! I think the above indicates a serious problem, so I'm going to put a stop here
      call par_stop
      ne_terp = ne
      ALLOCATE(rhoe_terp(0:lx,0:nphx,ne_terp+np))
      ALLOCATE(rhore_terp(251,0:nphx,ne_terp+np))
      ALLOCATE(energies_terp(ne_terp+np))
      rhoe_terp = rhoe
      rhore_terp = rhore
      energies_terp = energies
    ELSE
      ! Interpolating
      ne_terp = ind_m3 + nxmu + (ne-ind_p3+1)
      ALLOCATE(rhoe_terp(0:lx,0:nphx,ne_terp+np))
      ALLOCATE(rhore_terp(251,0:nphx,ne_terp+np))
      ALLOCATE(energies_terp(ne_terp+np))
      energies_terp = 0.d0
      rhoe_terp =  0.d0
      rhore_terp = 0.d0
      DO ie=1, ind_m3
        ! Copy vertical points and all points less than xmu-3kT
        energies_terp(ie) = energies(ie)
        rhoe_terp(:,:,ie) = rhoe(:,:,ie)
        rhore_terp(:,:,ie) = rhore(:,:,ie)
      ENDDO

      DO ie=ind_m3+1, ind_m3+nxmu
        ! Remove xmu +/- 3kT from original grid and replaced by interpolated points
        energies_terp(ie) = xmu_m3 + (ie-ind_m3-1) * dxmu

        DO iph=0, nphx
          DO il = 0, lx
            rhoe_terp(il,iph,ie) = interp1d(DBLE(energies_terp(ie)), DBLE(energies(1:ne)), rhoe(il,iph,1:ne))
          ENDDO
          DO il = 1,251
            rhore_terp(il,iph,ie) = interp1d(DBLE(energies_terp(ie)), DBLE(energies(1:ne)), rhore(il,iph,1:ne))
          ENDDO
        ENDDO
      ENDDO

      ! Insert xmu_p3 into the thing
      energies_terp(ind_m3+nxmu+1) = DBLE(xmu_p3) + coni*DIMAG(energies_terp(ind_m3+nxmu))

      DO ie=1, (ne-ind_p3)+np
        !  Going back to original grid at ind_p3
        energies_terp(ind_m3+nxmu+1+ie) = energies(ie+ind_p3)
        rhoe_terp(:,:,ind_m3+nxmu+1+ie) = rhoe(:,:,ie+ind_p3)
        rhore_terp(:,:,ind_m3+nxmu+1+ie) = rhore(:,:,ie+ind_p3)
      ENDDO
    ENDIF

    OPEN(unit=99, file='RHOE_TH.dat')
    OPEN(unit=1233, file="rhoe_th.dat")
    OPEN(unit=899, file='FERMI_TH.dat')
    DO ie = 1, ne_terp
      f = fermiF(energies_terp(ie), temp, xmu)
      IF (ie > ne_terp) f = -2*pi*coni
      WRITE(899,'(20E20.10E3)') DBLE(energies_terp(ie)*hart), DBLE(f)
      WRITE(99,'(20E20.10E3)') DBLE(energies_terp(ie)*hart), DIMAG(energies_terp(ie)),  DBLE(rhoe_terp(0,0,ie)*f), DIMAG(rhoe_terp(0,0,ie)*f), &
        & DBLE(rhoe_terp(1,0,ie)*f), DIMAG(rhoe_terp(1,0,ie)*f), DBLE(rhoe_terp(2,0,ie)*f), DIMAG(rhoe_terp(2,0,ie)*f)

    ENDDO
    DO ie=1,ne
      f = fermiF(energies(ie), temp, xmu)
      f = 1.d0
      WRITE(1233,'(20E20.10E3)') DBLE(energies(ie)*hart), DIMAG(energies(ie)),  DBLE(rhoe(0,0,ie)*f), DIMAG(rhoe(0,0,ie)*f), &
            & DBLE(rhoe(1,0,ie)*f), DIMAG(rhoe(1,0,ie)*f), DBLE(rhoe(2,0,ie)*f), DIMAG(rhoe(2,0,ie)*f)
    ENDDO
    WRITE(99,*)
    WRITE(899,*)
    CLOSE(99)
    CLOSE(899)
    CLOSE(1233)


    ! I. Integrate density over energy, fermi distribution
    !    The first calculated point has finite imaginary part. So, include contribution from
    !    leg between real axis and this point by approximating the integrand as constant
    !    along that leg (this is the same as what is done in scmt.f90).
    ep = DBLE(energies(1))
    iep = 1

    DO ie = 1, ne_terp
      ee = energies_terp(ie)
      de = ee - ep

      f = fermiF(ee, temp, xmu)
      fp = fermiF(ep, temp, xmu)

      DO iph = 0, nph
        DO il = 0, lx
          if (iunf.EQ.0 .AND. il.GT.2) CYCLE
          xnmues(il,iph) = xnmues(il,iph) + dimag( (rhoe_terp(il,iph,ie)*f + rhoe_terp(il,iph,iep)*fp) * de )
        END DO

        DO ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( (rhore_terp(ir,iph,ie)*f + rhore_terp(ir,iph,iep)*fp) * de )
        END DO
      END DO

      ep = ee
      iep = ie
    END DO
    xnmues_old = xnmues
    IF (print_flag) PRINT*, "    xnmues (contour) = ", SUM(xnmues)

    ! II. Add on contribution from Matsubara poles
    DO ip = 1, np
      DO iph = 0, nph
        DO il = 0, lx
          if (iunf.EQ.0 .AND. il.GT.2) CYCLE
          xnmues(il, iph) = xnmues(il, iph) + dimag( -4*pi*coni*temp * rhoe_terp(il,iph,ne_terp+ip) )
        END DO
        DO ir = 1, nr05(iph)
          rhoval(ir, iph) = rhoval(ir, iph) + dimag( -4*pi*coni*temp * rhore_terp(ir,iph,ne_terp+ip) )
        END DO
      END DO
    END DO
    IF (print_flag)  PRINT*, "    xnmues (Poles) = ", SUM(xnmues-xnmues_old)

    ! Sum over angular momenta and sites to get total electron count
    DO iph = 0, nph
      DO il = 0, lx
        xntot = xntot + xnmues(il, iph) * xnatph(iph)
      END DO
    END DO
    IF (print_flag) PRINT*, "    xntot = ", xntot

    DEALLOCATE(rhoe_terp, rhore_terp)
    DEALLOCATE(energies_terp)
  END SUBROUTINE

  SUBROUTINE calculate_contour(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, &
                              & vint, xnval, x0, ri, dx, adgc, adpc, dgc, dpc, &
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

  COMPLEX*16 FUNCTION fermiF(x,t,mu)
    ! Fermi function with cutoff to prevent float overflow
    IMPLICIT none
    COMPLEX*16, INTENT(IN) :: x
    REAL(8), INTENT(IN) :: t, mu
    REAL(8) :: b
    COMPLEX*16 :: a
    a = (x-mu)/t
    b = DIMAG((x-mu)/t)
    if (DBLE(x-mu)/t.le.5e2) then
        fermiF =  1/(1+EXP((x-mu)/t))
    else
      fermiF = 0.d0
    endif
    RETURN
  END FUNCTION fermiF
end module
