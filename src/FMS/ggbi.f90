!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ggbi.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggbi( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,  &
     &                 toler1, toler2, lcalc, msord)

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  use stkets

  implicit none
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential
  !     BiCGStab algorithm: Saad, Iterative methods for ..., p. 220 (1996)

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
  integer :: i,j,m1,m2,l1
  integer :: is,ist2,is1,is2,isp1
  integer :: ipart,ip,istart,ipass,nit,nitx  

  !     Lanczos method variables
  complex, allocatable, dimension(:) :: xvec, yvec, avec, asve
  complex, allocatable, dimension(:) :: rvec, pvec, svec
  complex aa, dd, aw, wa, ww
  complex del, delp, omega, chi, psi

  !      notice that in gglu we invert (1-Gt), but here (1-tG).
  !     multiply T and G0 matrices together, construct g0t = 1 - T*G0
  !     notice that the signs below for g0t ARE correct since 1 is the
  !     unit matrix
  !     since t is tri-diagonal, this product can be computed in n^2 time
  !     also fill up some work matrices for use in eigenvalue and
  !     determinant calculations and elsewhere
  !     cycle over dimensions of matrix g0t

  ! Allocate local variables
  allocate(xvec(istatx), yvec(istatx), avec(istatx), asve(istatx))
  allocate(rvec(istatx), pvec(istatx), svec( istatx))

  ! Set g0t=0
  do j = 1,istatx
     do i = 1,istatx
        g0t(i,j) = 0
     enddo
  enddo

  J_LOOP: do j = 1,istate
     I_LOOP:   do i = 1,istate
        ! T diagonal contribution T(i, i)
        if ( abs( g0(i, j)) .gt. toler2 ) then
           g0t(i,j)=g0t(i,j) - tmatrx(1,i) * g0(i,j) 
        endif
        ! T off-diagonal contribution T(ist2, i) in tmatr(2,i)
        ! T off-diagonal contribution T(i, ist2) in tmatr(2,ist2)
        l1   = lrstat(2,i)
        m1   = lrstat(3,i)
        isp1 = lrstat(4,i)
        m2 = m1+isp1
        if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
           ! spin-flip contribution
           ist2 = i + (-1)**isp1
           if ( abs( g0(ist2, j)) .gt. toler2) then
              g0t(i, j) = g0t(i, j)-tmatrx(nsp, ist2) * g0(ist2, j)
           endif
        endif
     enddo I_LOOP

     g0t(j, j) = g0t(j, j) + one
  enddo J_LOOP

  IP_LOOP: do ip=ipi, ipf
     ipart = nsp*(lipotx(ip)+1)**2
     IS_LOOP: do is1 = 1, ipart
        is2 = is1+i0(ip)
        l1   = lrstat(2,is2)
        if (.not.lcalc(l1)) cycle IS_LOOP

        ! Start first tier with xvec=0
        istart = -1
        msord = 0
        do is = 1, istate
           avec(is) = 0
           xvec(is) = 0
        enddo
        ! RESTART here if necessary
50      continue
        istart = istart+1

        if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)

        do is = 1,istate
           rvec(is) = - avec(is)
        enddo

        ! rvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
        rvec(is2) = rvec(is2) + 1

        ! Check convergence criteria: |r_n+1| < tol
        ipass = 0
        CONVERGE1: do is = 1, istate
           if (   (abs( real(rvec(is))).gt.toler1).OR. &
                & (abs(aimag(rvec(is))).gt.toler1)) then
              ipass=1
              exit CONVERGE1   !We don't need to check the others
           endif
        enddo CONVERGE1

        if (ipass.eq.0) goto 700

        do is = 1, istate
           pvec(is) = rvec(is)
        enddo

        call matvec( istatx,istate,g0t,pvec,avec,1)
        msord = msord + 1

        ! choose yvec that del and delp close to one
        call cdot( istatx, istate, avec, avec, aa)
        call cdot( istatx, istate, rvec, avec, wa)
        aw = real(wa) - coni* aimag(wa)
        call cdot( istatx, istate, rvec, rvec, ww)
        dd = aa*ww - aw*wa
        if (abs(dd/aa/ww) .lt.1.e-8) then
           do is = 1,istate
              yvec(is) = rvec(is) / ww
           enddo
        else
           ww = ( ww - aw ) / dd
           aa = ( wa - aa) / dd
           do is = 1,istate
              yvec(is) = rvec(is) * aa + avec(is) * ww
           enddo
        endif
        call cdot( istatx, istate, yvec, rvec, del)

        ! it seems ran out of precision for nit>150
        nitx = 30
        NIT_LOOP: do nit = 0, nitx
           call cdot( istatx, istate, yvec, avec, delp)
           omega = del / delp

           do  is = 1, istate
              svec(is) = rvec(is) - omega * avec(is)
           enddo
           ! Check convergence criteria: |s_n+1| < tol
           ipass = 0
           CONVERGE2: do is = 1, istate
              if (   (abs(real( svec(is))).gt.toler1).OR.  &
                   & (abs(aimag(svec(is))).gt.toler1)) then
                 ipass=1
                 exit CONVERGE2
              endif
           enddo CONVERGE2

           if (ipass.eq.0)  then
              do is = 1, istate
                 xvec(is) = xvec(is) + omega*pvec(is)
              enddo
              goto 700
           endif

           call matvec( istatx,istate,g0t,svec,asve,1)
           msord = msord + 1
           call cdot( istatx, istate, asve, asve, aa)
           call cdot( istatx, istate, asve, svec, wa)
           chi = wa / aa
           do is = 1, istate
              xvec(is) = xvec(is) + omega*pvec(is) + chi*svec(is)
           enddo

           do is = 1, istate
              rvec(is) = svec(is) - chi* asve(is)
           enddo
           ! Check convergence criteria: |r_n+1| < tol
           ipass = 0
           CONVERGE3: do is = 1, istate
              if   ( (abs(real( rvec(is))).gt.toler1).OR. &
                   & (abs(aimag(rvec(is))).gt.toler1) ) then
                 ipass=1
                 exit CONVERGE3
              endif
           enddo CONVERGE3
           ipass = 0

380        continue
           if (ipass.eq.0) goto 700

           ! Prepare for next iteration
           call cdot( istatx, istate, yvec, rvec, del)
           psi = del / (delp * chi)

           do is = 1, istate
              pvec(is) = rvec(is) + psi * (pvec(is)-chi*avec(is))
           enddo
           call matvec( istatx,istate,g0t,pvec,avec,1)
           msord = msord + 1

        enddo NIT_LOOP
        ! Restart since ran out of iterations
        goto 50

        ! Exit if tolerance has been achieved
700     continue
        ! print*, ' BI iterations:', nit + istart*nitx
        ! end of BI iterations

        ! at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
        ! pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
        ! for each ipot
        do is2=1,ipart
           gg( is2, is1, ip) = zero
           do is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +       &
                   &        g0( is2+i0(ip), is) * xvec(is)
           enddo
        enddo
     enddo IS_LOOP
  enddo IP_LOOP

  ! Deallocate local variables
  deallocate(xvec, yvec, avec, asve)
  deallocate(rvec, pvec, svec)



  return
end subroutine ggbi
