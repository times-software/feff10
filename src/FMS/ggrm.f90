!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ggrm.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggrm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg, toler1, toler2, lcalc, msord)
  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  use stkets
  implicit none
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential
  integer, intent(in) :: i0 (0:nphx), lipotx(0:nphx)
  integer, intent(in) :: nsp
  integer, intent(in) :: ipi,ipf
  complex, intent(in) :: tmatrx(nspx, istatx)
  real,    intent(in) :: toler1, toler2
  logical, intent(in) :: lcalc(0:lx)
  complex, intent(in) :: g0( istatx, istatx)
  integer, intent(out) :: msord
  complex, intent(out) :: g0t( istatx, istatx)
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)
  ! Definitions added to satisfy implicit none:
  integer :: i,j,is,is1,is2,ip,nit,nitx     ! Loop indeces
  integer :: l1,m1,m2,ist2,isp1,ipart,istart,ipass
  real    :: xfnorm
  complex, allocatable :: xvec(:), xket(:), xbra(:), xketp(:), xbrap(:)
  complex, allocatable :: zvec(:), rvec(:), svec(:), tket(:), tbra(:)
  real*8 ::  dum1, dum2
  complex alphac, betac, aa, bb, yy, aac, bbc, gamma
  real alpha, beta

  !      notice that in gglu we invert (1-Gt), but here (1-tG).
  !     multiply T and G0 matrices together, construct g0t = 1 - T*G0
  !     notice that the signs below for g0t ARE correct since 1 is the
  !     unit matrix
  !     since t is tri-diagonal, this product can be computed in n^2 time
  !     also fill up some work matrices for use in eigenvalue and
  !     determinant calculations and elsewhere
  !     cycle over dimensions of matrix g0t

  do j = 1,istatx
     do i = 1,istatx
        g0t(i,j) = 0
     enddo
  enddo

  ! Allocate local variables:
  allocate(xvec(istatx), xket(istatx), xbra(istatx))
  allocate(xketp(istatx), xbrap(istatx))
  allocate(zvec(istatx), rvec(istatx), svec(istatx))
  allocate(tket(istatx), tbra(istatx))

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
        msord=0
        do is = 1, istate
           rvec(is) = 0
           xvec(is) = 0
        enddo
        !         RESTART here if necessary
50      continue
        istart = istart+1

        if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,rvec,1)
        !         rvec = g0t*xvec - bvec, in our case bvec(is) = delta_{is,is2}
        rvec(is2) = rvec(is2) - 1
        do is = 1,istate
           xket(is) = - rvec(is)
        enddo

        call cdot( istatx, istate, xket, xket, bb)

        if (abs(bb).eq.0) goto 700

        xfnorm = 1.e0 / real(dble(bb))
        do is = 1, istate
           xbra(is) = xket(is) * xfnorm
        enddo
        !         |t> = A |n> ; |n> - xket, |n-1> - xketp
        call matvec ( istatx, istate, g0t, xket, tket, 1)
        msord = msord + 1
        call cdot( istatx, istate, xbra, tket, aa)
        aac = real(aa) - coni*aimag(aa)
        bb = 0
        bbc= 0
        betac = aa
        yy = 1
        !         initialize vectors
        do is = 1,istate
           xketp(is) = 0
           xbrap(is) = 0
           zvec(is) = xket(is)
           xvec(is) = xvec(is) + zvec(is)/betac
        enddo


        do is = 1, istate
           svec(is) = tket(is)
           rvec(is) = rvec(is) + svec(is) / betac
        enddo

        !         it seems ran out of precision for nit>150
        nitx = 100
        NIT_LOOP: do nit = 1, nitx
           !           use recursion method to calculate a_n+1, b_n, |n+1>, <n+1|
           do is = 1, istate
              tket(is) = tket(is) - aa*xket(is) - bb*xketp(is)
           enddo
           call matvec ( istatx, istate, g0t, xbra, tbra, 2)
           do is = 1, istate
              tbra(is) = tbra(is) - aac*xbra(is) - bbc*xbrap(is)
           enddo
           call cdot( istatx, istate, tbra, tket, bb)
           if (abs(bb).eq.0) goto 700

           bb = sqrt (bb)
           bbc = real(bb) - coni*aimag(bb)
           do is = 1, istate
              xketp(is) = xket(is)
              xbrap(is) = xbra(is)
           enddo

           do is = 1, istate
              xket(is) = tket(is) / bb
              xbra(is) = tbra(is) / bbc
           enddo

           call matvec ( istatx, istate, g0t, xket, tket, 1)
           msord = msord + 1
           call cdot( istatx, istate, xbra, tket, aa)
           aac = real(aa) - coni*aimag(aa)

           !           update iterative solution xvec, 
           !           and residual rvec = g0t*xvec - |1>
           alphac = bb / betac
           do is = 1, istate
              zvec(is) = xket(is) - alphac * zvec(is)
              svec(is) = tket(is) - alphac * svec(is)
           enddo

           betac = aa - alphac*bb
           yy = - alphac * yy
           gamma = yy / betac
           do is = 1, istate
              xvec(is) = xvec(is) + gamma * zvec(is)
              rvec(is) = rvec(is) + gamma * svec(is)
           enddo
           !c          Check convergence criteria: | rvec | < tol
           !           call vecvec( istatx, istate, rvec, rvec, dum2)
           !           if (dum2.le.tol) goto 700
           !c          Check convergence criteria: | rvec | < tol
           ipass = 1
           do is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 260
              if ( abs(aimag(rvec(is))).gt.toler1) goto 260
           enddo
           ipass = 0
260        continue
           if (ipass.eq.0) goto 700

        enddo NIT_LOOP
        !         restart since ran out of iterations
        goto 50

        !         exit if tolerance has been achieved
700     continue
        !         end of RM iterations

        !         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
        !         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
        !         for each ipot
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
  deallocate(xvec, xket, xbra)
  deallocate(xketp, xbrap)
  deallocate(zvec, rvec, svec)
  deallocate(tket, tbra)

  return
end subroutine ggrm

subroutine cdot ( istatx, istate, abra, aket, cc)
  !     dot product of two vectors
  !     notice that we keep bra  vector as it's complex conjugate
  !     thus need to conjugate abra here
  use constants, only: coni

  implicit none
  integer, intent(in) :: istatx,istate
  complex, intent(in), dimension(istatx) ::  abra, aket
  complex, intent(out) :: cc

  complex ::  aa
  integer :: is

  cc = 0
  do is = 1,istate
     aa = real(abra(is)) - coni*aimag(abra(is))
     cc = cc + aa * aket(is)
  enddo

  return
end subroutine cdot

subroutine vecvec ( istatx, istate, avec, bvec, cc)
  !     dot product of two vectors
  
  implicit none

  integer, intent(in) :: istatx,istate
  complex, intent(in), dimension(istatx) :: avec,bvec
  real*8,  intent(out) :: cc

  real*8 :: aa,bb
  integer :: is

  cc = 0
  do is = 1,istate
     aa = dble(real(avec(is))) * dble(real(bvec(is)))
     bb = dble(aimag(avec(is))) * dble(aimag(bvec(is)))
     cc = cc + aa + bb
  enddo
  return
end subroutine vecvec

subroutine matvec (istatx, istate, amat, bvec, cvec, itrans)
  !     itrans = 1  cvec = amat * bvec
  !     itrans = 2  cvec = amat^+ * bvec
  !     itrans = 3  cvec = amat^T * bvec
  use constants, only: coni

  implicit none

  integer, intent(in) :: istatx,istate, itrans
  complex, intent(in), dimension(istatx) :: bvec
  complex, intent(in), dimension(istatx,istatx) :: amat
  complex, intent(out), dimension(istatx) :: cvec

  complex aa
  integer is,i,j

  !     initialize cvec
  do is = 1,istatx
     cvec(is) = 0
  enddo
  !     cycle over dimensions of amat
  do j = 1,istate
     do i = 1,istate
        if (itrans.eq.1) then
           cvec(i) = cvec(i) + amat(i, j) * bvec(j)
        elseif(itrans.eq.2) then
           aa = real(amat(i, j)) - coni*aimag(amat(i, j))
           cvec(j) = cvec(j) + aa * bvec(i)
        else
           cvec(j) = cvec(j) + amat(i, j) * bvec(i)
        endif
     enddo
  enddo


  return
end subroutine matvec
