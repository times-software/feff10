!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ggtf.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggtf( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,  &
     &                 toler1, toler2, lcalc, msord)

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  use stkets

  implicit none
 
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential
  !     TFQMR: Saad, Iterative Methods for Sparse Matrices, p.225 (1996).

  ! Inputs
  integer, intent(in) :: i0 (0:nphx), lipotx(0:nphx)
  integer, intent(in) :: nsp
  integer, intent(in) :: ipi,ipf
  complex, intent(in) :: tmatrx(nspx, istatx)
  real,    intent(in) :: toler1, toler2
  logical, intent(in) :: lcalc(0:lx)
  complex, intent(in) :: g0( istatx, istatx)

  ! Outputs:
  ! Big work matrices
  integer, intent(out) :: msord
  complex, intent(out) :: g0t( istatx, istatx)
  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

  ! Definitions added to satisfy implicit none:
  integer :: i,j,is,is1,is2,ip,nit,nitx  ! Loop indecies
  integer :: m1,m2,l1,ipart,istart,isp1,ist2

  complex, allocatable, dimension(:) :: xvec, uvec, avec, wvec
  complex, allocatable, dimension(:) :: dvec, rvec, vvec

  complex alpha, beta, aa, rho, eta
  real tau, nu, cm, err

  !      notice that in gglu we invert (1-Gt), but here (1-tG).
  !     multiply T and G0 matrices together, construct g0t = 1 - T*G0
  !     notice that the signs below for g0t ARE correct since 1 is the
  !     unit matrix
  !     since t is tri-diagonal, this product can be computed in n^2 time
  !     also fill up some work matrices for use in eigenvalue and
  !     determinant calculations and elsewhere
  !     cycle over dimensions of matrix g0t

  ! Allocate local variables:
  allocate(xvec(istatx), uvec(istatx), avec(istatx), wvec(istatx))
  allocate(dvec(istatx), rvec(istatx), vvec(istatx))

  do j = 1,istatx
     do i = 1,istatx
        g0t(i,j) = 0
     enddo
  enddo

  do j = 1,istate
     do i = 1,istate
        !         T diagonal contribution T(i, i)
        if ( abs( g0(i, j)) .gt. toler2 )                       &
             &    g0t(i,j)=g0t(i,j) - tmatrx(1,i) * g0(i,j) 

        !         T off-diagonal contribution T(ist2, i) in tmatr(2,i)
        !         T off-diagonal contribution T(i, ist2) in tmatr(2,ist2)
        l1   = lrstat(2,i)
        m1   = lrstat(3,i)
        isp1 = lrstat(4,i)
        m2 = m1+isp1
        if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
           !            spin-flip contribution
           ist2 = i + (-1)**isp1
           if ( abs( g0(ist2, j)) .gt. toler2)                     &
                &       g0t(i, j) = g0t(i, j)                          &
                &                   - tmatrx(nsp, ist2) * g0(ist2, j) 
        endif
     enddo

     g0t(j, j) = g0t(j, j) + one
  enddo

  IP_LOOP: do ip=ipi, ipf
     ipart = nsp*(lipotx(ip)+1)**2
     IS_LOOP: do is1 = 1, ipart
        is2 = is1+i0(ip)
        l1   = lrstat(2,is2)
        if (.not.lcalc(l1)) cycle IS_LOOP

        !         start first tier with xvec=0
        istart = -1
        msord = 0
        do is = 1, istate
           rvec(is) = 0
           avec(is) = 0
           xvec(is) = 0
        enddo
        !         RESTART here if necessary
50      continue
        istart = istart+1

        if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
        do is = 1,istate
           uvec(is) = - avec(is)
        enddo
        !         uvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
        uvec(is2) = uvec(is2) + 1
        call matvec( istatx,istate,g0t,uvec,avec,1)
        msord = msord + 1
        do is = 1, istate
           wvec(is) = uvec(is)
           vvec(is) = avec(is)
           dvec(is) = 0
        enddo
        call cdot( istatx, istate, uvec, uvec, aa)
        tau = sqrt(real(aa))
        nu = 0
        eta = 0
        !         choose rvec = uvec /aa so that dot products are about 1
        do is = 1, istate
           rvec(is) = uvec(is) / aa
        enddo
        rho = 1

        !         it seems ran out of precision for nit>150
        nitx = 20
        NIT_LOOP: do nit = 0, nitx
           if (mod(nit,2).eq.0) then
              call cdot( istatx, istate, rvec, vvec, aa)
              alpha = rho / aa
           else
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
           endif

           do  is = 1, istate
              wvec(is) = wvec(is) - alpha * avec(is)
           enddo
           aa = nu**2 * eta / alpha
           do is = 1, istate
              dvec(is) = uvec(is) + aa * dvec(is)
           enddo
           call cdot( istatx, istate, wvec, wvec, aa)
           nu = sqrt(real(aa)) / tau
           cm = 1 / sqrt(1+nu**2)
           tau = tau * nu * cm
           eta = cm**2 * alpha
           do is = 1, istate
              xvec(is) = xvec(is) + eta * dvec(is)
           enddo

           !c          Check convergence criteria: | rvec | < tol
           err = (1.e0 + nit) / istate
           err = tau * sqrt(err) * 10
           if ( abs(err).lt.toler1) goto 700

           if (mod(nit,2) .ne.0) then
              aa = rho
              call cdot( istatx, istate, rvec, wvec, rho)
              beta = rho / aa
              do is = 1, istate
                 uvec(is) = wvec(is) + beta * uvec(is)
                 vvec(is) = beta * ( avec(is) + beta * vvec(is))
              enddo

              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
              do is = 1, istate
                 vvec(is) = avec(is) + vvec(is)
              enddo
           else
              do is = 1, istate
                 uvec(is) = uvec(is) - alpha * vvec(is)
              enddo
           endif
        enddo NIT_LOOP
        ! Restart since ran out of iterations
        goto 50

        ! Exit if tolerance has been achieved
700     continue
        ! End of TFQMR iterations

        ! At this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
        ! pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
        ! for each ipot
        do is2=1,ipart
           gg( is2, is1, ip) = zero
           do is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +                   &
                   &        g0( is2+i0(ip), is) * xvec(is)
           enddo
        enddo
     enddo IS_LOOP
  enddo IP_LOOP

  ! Deallocate local variables:
  deallocate(xvec, uvec, avec, wvec)
  deallocate(dvec, rvec, vvec)

  return
end subroutine ggtf
