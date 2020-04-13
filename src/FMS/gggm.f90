!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: gggm.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gggm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,  &
     &                 toler1, toler2, lcalc, msord)

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  use stkets

  implicit none
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential
  !     Lanczos algorithm: Graves-Morris,Salam, Num.Algor.21,p.213(1999)

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
  integer :: i,j,ip,is1,is2,is,nit  ! Loop indecies
  integer :: m1,m2,l1,isp1,ist2,ipass,nitx,istart,ipart
  real    :: tol

  !     Lanczos method variables
  complex, allocatable, dimension(:) :: xvec, wvec, x0, x1, avec
  complex, allocatable, dimension(:) :: bvec, r0, r1, t0, t1

  complex aa, dd, aw, wa, ww, e0, e1, alpha, beta, theta, q0, q1

  !      notice that in gglu we invert (1-Gt), but here (1-tG).
  !     multiply T and G0 matrices together, construct g0t = 1 - T*G0
  !     notice that the signs below for g0t ARE correct since 1 is the
  !     unit matrix
  !     since t is tri-diagonal, this product can be computed in n^2 time
  !     also fill up some work matrices for use in eigenvalue and
  !     determinant calculations and elsewhere
  !     cycle over dimensions of matrix g0t

  ! Allocate local variables:
  allocate(xvec(istatx), wvec(istatx), x0(istatx), x1(istatx), avec(istatx))
  allocate(bvec(istatx), r0(istatx),   r1(istatx), t0(istatx), t1(istatx))

  do j = 1,istatx
     do i = 1,istatx
        g0t(i,j) = 0
     enddo
  enddo


  do j = 1,istate
     do i = 1,istate
        !         T diagonal contribution T(i, i)
        if ( abs( g0(i, j)) .gt. toler2 )                       &
             &    g0t(i,j)=g0t(i,j) + tmatrx(1,i) * g0(i,j) 

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
                &                   + tmatrx(nsp, ist2) * g0(ist2, j) 
        endif
     enddo
     !       g0t(j, j) = g0t(j, j) + one
  enddo

  IP_LOOP:  do ip=ipi, ipf
     ipart = nsp*(lipotx(ip)+1)**2
     IS_LOOP: do is1 = 1, ipart
        is2 = is1+i0(ip)
        l1   = lrstat(2,is2)
        if (.not.lcalc(l1)) cycle

        ! Start first tier with xvec=0
        istart = -1
        msord = 0
        do is = 1, istate
           bvec(is) = 0

           xvec(is) = 0
        enddo
        ! rvec = bvec - A*xvec , in our case bvec(is) = delta_{is,is2}
        bvec(is2) = 1

        ! RESTART here if necessary
50      continue
        istart = istart+1

        if (istart.gt.0) then
           do is = 1, istate
              xvec(is) = xvec(is) + x0(is) / q0
           enddo
           call matvec( istatx,istate,g0t,xvec,avec,1)
           do is = 1, istate
              bvec(is) = avec(is) - xvec(is)
           enddo
           bvec(is2) = bvec(is2) + 1
        endif
        do is = 1,istate
           r0(is) = bvec(is)
           x0(is) = 0
           x1(is) = bvec(is) 
        enddo

        call matvec( istatx,istate,g0t,bvec,r1,1)
        msord = msord + 1

        ! Choose wvec that del and delp close to one
        call cdot( istatx, istate, r0, r0, ww)
        call cdot( istatx, istate, r1, r1, aa)
        call cdot( istatx, istate, r0, r1, wa)
        aw = real(wa) - coni* aimag(wa)
        dd = aa*ww - aw*wa
        if (abs(dd/aa/ww) .lt.1.e-8) then
           do is = 1,istate
              wvec(is) = r0(is) / ww
           enddo
        else
           ww = ( ww - aw ) / dd
           aa = ( wa - aa) / dd
           do is = 1,istate
              wvec(is) = r0(is) * aa + r1(is) * ww
           enddo
        endif
        !         update dot products to avoid round off errors
        call cdot( istatx, istate, wvec, r0, e0)
        call cdot( istatx, istate, wvec, r1, e1)
        q0 = 1
        q1 = 1

        !         it seems ran out of precision for nit>150
        nitx = 10
        NIT_LOOP: do nit = 1, nitx
           tol = toler1 * abs(q1) /10
           !c          Check convergence criteria: |r1| < tol / 10
           !c          so mostly code will not exit here
           ipass = 1
           do is = 1, istate
              if ( abs(real(r1(is))).gt.tol) goto 99
              if ( abs(aimag(r1(is))).gt.tol) goto 99
           enddo
           ipass = 0
99         continue
           if (ipass.eq.0) then
              do is = 1, istate
                 xvec(is) = xvec(is) + x1(is) / q1
              enddo
              goto 700
           endif

           alpha = e1 / e0
           do is = 1, istate
              t0(is) = r1(is) - alpha* r0(is)
           enddo
           call matvec( istatx,istate,g0t,t0,t1,1)
           msord = msord + 1

           call cdot( istatx, istate, t0, t1, wa)
           call cdot( istatx, istate, t0, t0, ww)
           call cdot( istatx, istate, t1, t1, aa)
           aw = real(wa) - coni* aimag(wa)
           theta = (wa - aa) / (ww - aw)

           do is = 1, istate
              r0(is) = t1(is) - theta * t0(is)
           enddo
           dd = 1- theta
           do is = 1, istate
              x0(is) = t0(is) + dd * (x1(is) - alpha*x0(is))
           enddo
           q0 = dd * (q1 - alpha*q0)
           tol = toler1 * abs(q0)

           ! Check convergence criteria: |r0| < tol
           ipass = 1
           do is = 1, istate
              if ( abs(real(r0(is))).gt.tol) goto 380
              if ( abs(aimag(r0(is))).gt.tol) goto 380
           enddo
           ipass = 0
380        continue
           if (ipass.eq.0) then 
              do is = 1, istate
                 xvec(is) = xvec(is) + x0(is) / q0
              enddo
              goto 700
           endif

           ! prepare for next iteration
           call cdot( istatx, istate, wvec, r0, e0)
           beta = e0 / e1
           do is = 1, istate
              t0(is) = r0(is) - beta * r1(is)
           enddo
           call matvec( istatx,istate,g0t,t0,avec,1)
           msord = msord + 1
           dd = beta * theta
           do  is = 1, istate
              r1(is) = avec(is) + dd * r1(is)
           enddo
           call cdot( istatx, istate, wvec, r1, e1)

           dd = beta * (1-theta)
           do is = 1, istate
              x1(is) = x0(is) - dd * x1(is) + t0(is)
           enddo
           q1 = q0 - (1-theta) * beta * q1
        enddo NIT_LOOP
        ! Restart since ran out of iterations
        goto 50

        ! Exit if tolerance has been achieved
700     continue
        ! End of GM iterations

        ! At this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
        ! pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
        ! for each ipot
        do is2=1,ipart
           gg( is2, is1, ip) = zero
           do is = 1,istate
              gg(is2,is1,ip) = gg(is2,is1,ip) +  g0(is2+i0(ip),is) * xvec(is)
           enddo
        enddo

     enddo IS_LOOP
  enddo IP_LOOP

  ! Deallocate local variables:
  deallocate(xvec, wvec, x0, x1, avec)
  deallocate(bvec, r0, r1, t0, t1)



  return
end subroutine gggm
