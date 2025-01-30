!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: chiklu.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine chiklu( nx, dr, chik, fnew, yvec, maxsize, matsize)
!     fix later : this routine works only with RGRID 0.01 !
!     need to pass grid information and interpolate between
!     0.05 and finer grids
      use dimsmod, only: nrptx,nphasx,lx,nclusx,nspx=>nspu !istatx_dont_use_here=>istatx
	  use constants
      implicit real (a-h,o-z)
      implicit integer (i-n)
!  input:
!     nx - number of points on 0.05 grids
!     dr - now 0.01 radial grid
!     chik - K*chi on coarse 0.05 grid
!     yvec(r,n) = K(r,r')*phi^n_init(r')*phi^n_final(r'),
!        where n - index in  double basis set
!     matsize - maximum n
!  output
!     fnew - (1-K*chi)^-1 r
!     yvec = (1-K*chi)^-1 yvec

!KJ here used to be : include 'xparamtd.h'  - only location in feff where this file needed, so I copy/pasted it.
!KJ then most of them were in xparam.h and hence dimsmod, so only istatx needs to be redefined! (as nclxtd is separate from nclusx)
!KJ Then I noticed that lrstat isn't used and  istate is just an index, so I took those out, too.
!KJ Oh, and since lx is now a variable, so is istatx, and therefore allocate statements are needed ...  Yawn ... You still here?  Want my job?
      integer istatx   ! =(lx+1)**2*nclusx*nspx
      integer, allocatable :: ipiv(:)
!KJ!**** array of state kets at current energy
!KJ      common /stkets/ lrstat(4, istatx), istate   !KJ plus this common doesn't exist anywhere else in FEFF, so what's the point??
      integer istate !KJ

      complex, allocatable :: tmatrx(:,:)
      complex   chik(251,251)
!     big work matrices
      complex, allocatable ::   g0(:,:), g0t(:,:), g0s(:,:)
!     return matrix containing info about each unique potential
      complex,allocatable ::   gg(:,:,:)
      double precision dr(nrptx)
      complex*16 fnew(nrptx), yvec(nrptx, maxsize)

      character*3  cerr
      character*13 trans

 400  format(i4)

      istatx=nx  ! JK - turns out that istatx is never used below, and all of the matrices should be
                                    ! allocated with nx. The normal definition
                                    ! of istatx only works if the number of
                                    ! atoms inside fms happens to be large
                                    ! enough to make nclusx large. This was always the case when
                                    ! nclusx was static.
      !PRINT*, 'istatx, lx, nlusx, nspx'
      !PRINT*, istatx, lx, nclusx, nspx
      !allocate(ipiv(istatx),tmatrx(nspx, istatx),g0( istatx, istatx), &
      ! g0t( istatx, istatx),g0s( istatx, nspx*(lx+1)**2),gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx))
      ! Below g0s really has only one dimension since the second is always set
      ! to 1. Also g0, tmatrx, and gg are never used, so taking them out
      allocate(ipiv(nx), &
	     g0t( istatx, istatx),g0s( istatx, 1))



! -------------------- LU gg
!     multiply T and G0 matrices together, construct g0t = 1 - G0*T
!     notice that the signs below for g0t ARE correct since 1 is the
!     unit matrix
!     since t is tri-diagonal, this product can be computed in n^2 time
!     also fill up some work matrices for use in eigenvalue and
!     determinant calculations and elsewhere
! JK - Note that this is not really T and G0, with indices corresponding to atomic index.
!      Instead, g0t contains 1 - chi0K(r,r'), so that the indices run over radial grids.
      istate = nx
      do 320 icol = 1,istate
        do 310 irow = 1,istate
        g0t(irow, icol) = -chik(irow, icol)
 310    continue
        g0t(icol, icol) = g0t(icol, icol) + one
 320  continue

! --- invert matrix by LU decomposition
!     call cgetrf from lapack.  this performs an LU decomposition on
!     the matrix g0t = 1 - g0*T
      call cgetrf( istate, istate, g0t, istatx, ipiv, info )
      if (info.lt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)abs(info)
          call wlog('        Argument #'//cerr//                        &
     &                ' had an illegal value.')
      elseif (info.gt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)info
          call wlog('        g0t('//cerr// ','//cerr//                  &
     &                ') is exactly 0 -- '//                            &
     &                'this matrix cannot be decomposed.')
      endif

!     now we want g_c = (g0t)^-1 * g0.  Rather than calculating
!     the inverse of g0t from the LU decomposition, we can compute
!     g0t^-1 * g0 directly by backsubstituting the columns of G0.
!     See sect. 2.3 in Numerical Recipes or LAPACK Users' Guide
!     sect. 2.3

!     third arg in number of output columns, istate for full
!     matrix, ipart(ik) for just the parts of the matrix needed
!     to contruct fine structure + DOS functions

      ipart = 1
      do 700 im = 0, matsize
        if (im.eq.0) then
           do 590 is1 = 1, istate
             i = 1+5*(is1-1)
             g0s(is1,ipart) = real(dr(i))
  590      continue
        else
           do 591 is1 = 1, istate
             i = 1+5*(is1-1)
             g0s(is1,ipart) = real(dble(yvec(i,im)))
  591      continue
        endif

        trans = 'NotTransposed'
        call cgetrs(trans, istate, ipart, g0t, istatx,                  &
     &                ipiv, g0s, istatx, info)

        ! JK - g0s now contains [1-K*chi0]^{-1} r
        ! or depending on matsize (which seems to always be 0), so the above.
        ! JK - g0s now contains [1-K*chi0]^{-1} K phi_init*phi_final
        if (info.lt.0) then
            call wlog('    *** Error in cgetrf')
            write(cerr,400) abs(info)
            call wlog('        Argument #'//cerr//                      &
     &              ' had an invalid value.')
        endif

! **** at this point g0s contains the full MS ****

!  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each
!  ipot

        do 600 is=1,istate
          i = 1+5*(is-1)
          if (im.eq.0) then
            fnew(i) = cmplx(dble(real(g0s(is,ipart))),                  &
     &                      dble(aimag(g0s(is,ipart)))) / dr(i)
            if (is.gt.1) then
              do 595 j=1,4
  595         fnew(i-j) = ( j*fnew(i-5) + (5-j)*fnew(i)) / 5
            endif
          else
            yvec(i,im) = cmplx(dble(real(g0s(is,ipart))),               &
     &                      dble(aimag(g0s(is,ipart)))) / dr(i)
            if (is.gt.1) then
              do 596 j=1,4
  596         yvec(i-j,im) = ( j*yvec(i-5,im) + (5-j)*yvec(i,im)) / 5
            endif
          endif
 600    continue
 700  continue

      deallocate(ipiv,g0t,g0s)
      return
      end
