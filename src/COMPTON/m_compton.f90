module compton_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module handles calculations of compton profiles, J(p_q)
!
! J(p_q) is a projected electron momentum density along the
! direction of momentum transfer (p_q := p . qhat).
!
! (1)  J(p_q) = Integral d^3p rho(p) delta(p_q - p . qhat)
!
! where
!
! (2)  rho(p) = Integral d^3r d^3r' e^(i p . (r - r')) rho(r,r')
!
! and
!
! (3) rho(r,r') = Integral dE (2/pi) Im G(r,r',E) f(E).
!
! Here, G is the Green's function and f is the fermi distribution.
!
! See PRB 85 115135 (2012) for more details.
!
! A slice of the FMS matrix g_nLn'L' with n = 0 is saved by the FMS routines
! into gg_slice.bin. This allows calculation of rho(r,r') for r contained
! within the Voronoi cell of the central atom. If rho for r outside of the
! central-cell is needed, the corresponding portions of g must be saved
! (search for gglu_save_slice in FMS/gglu.f90 for the correct code to modify).
! Note that the full matrix is rather large, which is why we currently only
! save part of it...
!
! If we assume the z-axis lies along qhat, then the p_x and p_y integrals
! in eq. (1) are trivial, giving delta functions that set x'=x and y'=y.
! Thus, after performing the x' and y' integrals, we have
!
! (4)  J(p_q) = Integral dx dy dz dz' rho(r,r') exp(i p_z (z-z'))
!
! where r' is defined as (x,y,z').
! The x-y integral is expressed in radial coordinates, which we call
! (s, phi).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use constants
  use rotation_mod
  use rhorrp_mod
  use compton_inp
  implicit none

  type compton_grid
    integer ns, nphi, nz, nzp
    real*8, pointer, dimension(:) :: s, phi, z, zp
    logical rotate
    real*8, dimension(3,3) :: rotation_matrix
  end type
contains

  subroutine compton_init
    call rhorrp_init
  end subroutine

  subroutine compton_deinit
    call rhorrp_deinit
  end subroutine


  subroutine compton_jzzp(gr, jzzp)
    !
    ! Calculate J(z, zp)
    !
    ! Parameters:
    !   gr:   (compton_grid) parameters for integration grid
    !   jzzp: output matrix
    !
    ! Integrates rho(r,rp) over x-y plane where r=(x,y,z) and rp=(x,y,zp)

    ! TODO: Currently, the trapezoid rule is used for all integration.
    ! This could probably be made more efficient/accurate by switching
    ! to adaptive Gaussian quadrature (e.g., quadpack).
    !

    use par
    !use rotation_mod
    type(compton_grid) :: gr
    real*8, intent(inout) :: jzzp(gr%nz, gr%nzp)

    real*8 s,phi,z,zp, r2, rnrm2
    integer is, iphi, iz, izp
    real*8, dimension(3) :: v, vp

    real*8 ds, dphi
    real*8 :: rho, prev, J, Jprev

    integer row_proc
    character*512 slog

    ! The restrict_r_voronoi flag in rhorrp_mod causes the rhorrp subroutine to
    ! return 0 if the r coordinate is outside the central atom's Voronoi cell.
    restrict_r_voronoi = .true.

    ! Integration loop over z'
    do izp=1,gr%nzp
      ! Decide which process should handle this z' point
      row_proc = mod((izp - 1), numprocs)

      if (row_proc == this_process) then
        write(slog,'(a,i5,a,i5)') ' Calculating j(z,zp) for point ', izp, '/', gr%nzp
        call wlog(slog)

        ! Integration  loop over z
        do iz=1,gr%nz
          jzzp(iz,izp) = 0
          prev = 0

          ! Integration loop over s
          do is=1,gr%ns
            s = gr%s(is)

            J = 0 ! J holds integral over phi at the current z,z',s point
            ! Integration loop over phi
            do iphi=1,gr%nphi
              phi = gr%phi(iphi)

              ! fill in values for current r-vector and r'-vectors
              v(1) = s * cos(phi)
              vp(1) = v(1)
              v(2) = s * sin(phi)
              vp(2) = v(2)
              v(3) = gr%z(iz)
              vp(3) = gr%zp(izp)

              ! at this point, v and vp are in a coordinate system with qhat in the z direction
              ! here, we rotate them to the cluster's axes (if necessary)
              if (gr%rotate) then
                call rotate_in_place(gr%rotation_matrix, v)
                call rotate_in_place(gr%rotation_matrix, vp)
              end if

              ! Calculate the density matrix at this r, r' point
              call rhorrp(v, vp, rho)

              ! multiply by cylindrical coordinate measure 
              rho = rho * s

              ! trapezoid rule for phi integration
              if (iphi > 1) then
                dphi = phi - gr%phi(iphi-1)
                J = J + (prev + rho)/2.0 * dphi / (2*pi) ! (this last factor depends on symmetry
              endif

              prev = rho
            end do

            ! trapezoid rule for s
            ! (jzzp(iz,izp) is integral over s and phi at current z,z' point)
            if (is > 1) then
              ds = s - gr%s(is-1)
              jzzp(iz,izp) = jzzp(iz,izp) + (J + Jprev) / 2.0 * ds
            endif
            Jprev = J
          end do
        end do

        ! send calculation at current point back to master
        if (worker) then
          call par_send_double(jzzp(1,izp), gr%nz, 0, izp)
        end if
      end if ! row_proc == this_proc

      ! receive calculations from workers 
      if (master .and. .not. row_proc .eq. 0) then
        call par_recv_double(jzzp(1,izp), gr%nz, row_proc, izp)
      end if
    end do

    ! turn off this flag.
    restrict_r_voronoi = .false.
    call par_barrier
  end subroutine

  subroutine jpq(gr, jzzp, pq, J)
    !
    ! Calculate Compton profile J(p_q) from J(z,z')
    !
    ! Parameters:
    !   gr: (compton_grid) grid parameters
    !   jzzp: precalculated jzzp matrix (see compton_jzzp subroutine)
    !   pq: value of project momentum
    !   J: output - Compton profile at pq
    !
    ! J(p_q) = \int e^{i p_q (z-z')} J(z,z') dz dz'
    !
    ! For pq = 0, this uses the trapezoid rule
    ! For pq /= 0, a piecewise linear Fourier transform is used instead
    !
    use compton_inp,only: window, window_cutoff
    type(compton_grid), intent(in) :: gr
    real*8, intent(in), dimension(gr%nz,gr%nzp) :: jzzp
    real*8, intent(in) :: pq
    real*8, intent(out) :: J

    ! local variables
    integer iz, izp
    real*8 z, pz, zp, pzp ! pz = 'previous z', pzp = "previous z prime"
    complex*16 :: J1, J2,J2p
    complex*16 :: A, B, D
    real*8 :: cutoff
    real*8 w

    ! Use upper bound of z' integration unless a nonzeo cutoff is specified in the input file
    if (window_cutoff == 0) then
      cutoff = gr%zp(gr%nzp)
    else
      cutoff = window_cutoff
    end if

    ! J1 contains integral over z and z'
    ! J2 is integral over z' for current z point
    ! J2p is integral over z' for previous z point
    J1 = 0
    J2p = 0

    pz = 0
    pzp = 0

    ! Integration loop over z
    do iz=1,gr%nz
      J2 = 0

      z = gr%z(iz)
      pzp = gr%zp(1)

      ! Integration loop over z'
      do izp=2,gr%nzp

        zp = gr%zp(izp)

        ! window fourier transform to reduce aliasing from sharp cutoff
        if (window == 0) then
          w = 1.0
          if (abs(zp) > cutoff) w = 0.0
        else if (window == 1) then
          w = cos(pi * abs(zp) / (2 * cutoff)) ** 2
          if (abs(zp) > cutoff) w = 0
        else
          w = 1.0
        end if

        ! Update integral with current z' point
        if (pq == 0) then
          ! trapezoid
          J2 = J2 + (jzzp(iz,izp) + jzzp(iz,izp-1))/2.0 * (zp - pzp) * w
        else
          ! piecewise linear FT
          D = (jzzp(iz,izp) - jzzp(iz,izp-1)) / (zp - pzp)
          A = jzzp(iz,izp-1) + D * (zp - pzp + coni / pq)
          B = jzzp(iz,izp-1) + coni * D / pq

          J2 = J2 + exp(coni * pq * zp) * A / (coni * pq) * w&
               &  - exp(coni * pq * pzp) * B / (coni * pq) * w
        endif

        pzp = zp
      enddo

      if (iz > 1) then
        ! Update integral with current z point
        if (pq == 0) then
          ! trapezoid
          J1 = J1 + (J2 + J2p) / 2.0 * (z - pz)
        else
          ! piecewise linear fourier transform
          D = (J2 - J2p) / (z - pz)
          A = J2p + D * (z - pz - coni / pq)
          B = J2p - coni * D / pq

          J1 = J1 + exp(-coni * pq * z) * A / (-coni * pq) &
               &  - exp(-coni * pq * pz) * B / (-coni * pq)
        endif
      endif
      J2p = J2
      pz = z
    enddo

    ! return result
    ! XXX check that dimag(J1) is small here...
    J = dble(J1)
  end subroutine

  subroutine compton_build_grid(grid)
    !
    ! set up grid for compton integration
    !

    use compton_inp
    type (compton_grid), intent(out) :: grid

    real*8, dimension(3) :: zhat, axis
    real*8 :: theta

    integer i

    ! Set upper bound on s and z grid to Norman radius if 0 is specified
    if (smax == 0) smax = rnrm(0)
    if (zmax == 0) zmax = rnrm(0)
    grid%ns   = ns
    grid%nphi = nphi
    grid%nz   = nz
    grid%nzp  = nzp

    allocate( grid%s(grid%ns), grid%phi(grid%nphi), grid%z(grid%nz), grid%zp(grid%nzp) )

    ! pre-calculate linear grids for all coordinates
    do i=1,ns
      grid%s(i) = smax/(ns-1) * (i - 1)
    end do

    do i=1,nphi
      grid%phi(i) = phimax/(nphi-1) * (i - 1)
    end do

    do i=1,nz
      grid%z(i) = -zmax + 2*zmax/(nz-1) * (i - 1)
    end do

    do i=1,nzp
      grid%zp(i) = -zpmax + 2*zpmax/(nzp-1) * (i - 1)
    end do

    ! set up rotation matrix
    zhat(:) = 0.0; zhat(3) = 1.0
    call rotation_axis_angle(zhat, qhat, axis, theta)
    if (abs(theta) > 1e-10) then
      grid%rotate = .true.
      call rotation_matrix(axis, theta, grid%rotation_matrix)
    else
      grid%rotate = .false.
      ! set rotation matrix to the identity matrix
      grid%rotation_matrix(:,:) = 0.0
      do i=1,3
        grid%rotation_matrix(i,i) = 1.0
      end do
    end if

  end subroutine
end module
! vim: et
