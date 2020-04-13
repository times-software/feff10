!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: gglu.f90,v $:
! $Revision: 1.9 $
! $Author: bmattern $
! $Date: 2012/02/09 18:04:57 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gglu_h(nsp,i0,ipi,ipf,lipotx,g0,tmatrx,tmatrxfull,g0t,gg,ck,ie)

  use DimsMod, only: istatx, nphx=>nphu, nspx=>nspu, lx, nphasx
  use constants
  USE IOMod
  use stkets
  use controls,only : gglu_save_slice

  implicit none

  complex :: ck  !KJ debugging
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential

  ! Inputs
  integer, intent(in) :: i0 (0:nphx), lipotx(0:nphx), ie
  integer, intent(in) :: nsp
  integer, intent(in) :: ipi,ipf
  complex, intent(in) :: tmatrx(nspx, istatx)
  complex, intent(in) :: tmatrxfull(istatx, istatx)
  complex, intent(in) :: g0( istatx, istatx)

  ! Outputs:
  ! Big work matrices
  complex, intent(out) :: g0t( istatx, istatx)
  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

  ! Definitions added to satisfy implicit none:
  integer :: i,j,ip,icase,is1,is2,ist1,ist2  ! Loop indecies
  integer :: isp1,jj1,jj2,iatom1,iatom2,ipart,info,nstates
  integer :: l1,m1,m2, icol,iat1,irow,iat2,l2,ik 

  integer, allocatable :: ipiv(:)

  ! Big work matrices
  complex, allocatable :: g0s(:,:)
  complex, allocatable :: q(:,:) !KJ debug

  character*3  cerr
  character*13 trans
  character*10 fname

  logical,parameter :: makeg=.false.
  logical,parameter :: FullG = .FALSE.

400 format(i4)
  ! Allocate local variables:
  ! JK - added check if allocated. This way we don't have to allocate at each
  ! energy point. 9/2009
  IF(.NOT.ALLOCATED(ipiv)) allocate(ipiv(istatx))
  IF(.NOT.ALLOCATED(g0s)) allocate(g0s( istatx, nspx*(lx+1)**2))
  IF(.NOT.ALLOCATED(q)) allocate(q(nspx*(lx+1)**2*8,nspx*(lx+1)**2*8)) !KJ debug

  ! -------------------- LU gg
  !     multiply T and G0 matrices together, construct g0t = 1 - G0*T
  !     notice that the signs below for g0t ARE correct since 1 is the
  !     unit matrix
  !     since t is tri-diagonal, this product can be computed in n^2 time
  !     also fill up some work matrices for use in eigenvalue and
  !     determinant calculations and elsewhere
      IF(.TRUE.) THEN ! tmatrix is a full istate by istate matrix in g0t
         do icol = 1,istate
           iat1 = lrstat(1, icol)
           l1   = lrstat(2, icol)
!           write(88,*) tmatrx(1,icol) - tmatrxfull(icol,icol)
           do irow = 1,istate
              g0t(irow,icol) = 0.d0
              do ik = 1, istate
                    iat2 = lrstat(1, ik)
                    l2   = lrstat(2, ik)
                    IF((l1.EQ.l2).AND.(iat1.EQ.iat2)) THEN
                       g0t(irow,icol) = g0t(irow,icol) - g0(irow,ik)*  tmatrxfull(ik,icol)
                    END IF
              end do
           end do
              g0t(icol,icol) = g0t(icol,icol)+1
        end do   
      ELSE
          do j = 1,istate
             do i = 1,istate
                !         T diagonal contribution
                g0t(i, j) = - g0(i, j) * tmatrx(1, j)
                !         T off-diagonal contribution
                l1   = lrstat(2,j)
                m1   = lrstat(3,j)
                isp1 = lrstat(4,j)
                m2 = m1+isp1
                if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
                   ist2 = j + (-1)**isp1
                   g0t(i, j) = g0t(i, j) - g0(i, ist2) * tmatrx(nsp, j)
                endif
             enddo
             g0t(j, j) = g0t(j, j) + one
          enddo
      END IF


  ! --- invert matrix by LU decomposition
  !     call cgetrf from lapack.  this performs an LU decomposition on
  !     the matrix g0t = 1 - g0*T
  call cgetrf( istate, istate, g0t, istatx, ipiv, info )
  if (info.lt.0) then
     call wlog('    *** Error in cgetrf when computing G')
     write(cerr,400)abs(info)
     call wlog('        Argument #'//cerr//' had an illegal value.')
  elseif (info.gt.0) then
     call wlog('    *** Error in cgetrf when computing G')
     write(cerr,400)info
     call wlog('        g0t('//cerr// ','//cerr//  ') is exactly 0 -- '//                            &
                          'this matrix cannot be decomposed.')
  endif

  !     now we want g_c = (g0t)^-1 * g0.  Rather than calculating
  !     the inverse of g0t from the LU decomposition, we can compute
  !     g0t^-1 * g0 directly by backsubstituting the columns of G0.
  !     See sect. 2.3 in Numerical Recipes or LAPACK Users' Guide
  !     sect. 2.3

  !     third arg in number of output columns, istate for full
  !     matrix, ipart(ik) for just the parts of the matrix needed
  !     to contruct fine structure + DOS functions

  IP_LOOP: do ip=ipi, ipf
     ipart = nsp*(lipotx(ip)+1)**2
     do is1 = 1, istate
        do is2 = 1, ipart
           g0s(is1,is2) = g0(is1, is2 + i0(ip))
        enddo
     enddo
     trans = 'NotTransposed'
     call cgetrs(trans, istate, ipart, g0t, istatx, ipiv, g0s, istatx, info)

     if (info.lt.0) then
        write(cerr,400) abs(info)
        call wlog('    *** Error in cgetrf')
        call wlog('        Argument #'//cerr// ' had an invalid value.')
     endif

     ! **** at this point g0s contains the full MS ****
     !  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each
     !  ipot
     do is2=1,ipart
        do is1=1,ipart
           gg( is1, is2, ip) = g0s( is1+i0(ip), is2)
        enddo
     enddo
  enddo IP_LOOP


  return
end subroutine gglu_h 
