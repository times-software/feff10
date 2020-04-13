!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: gglu.f90,v $:
! $Revision: 1.10 $
! $Author: bmattern $
! $Date: 2012/07/01 16:22:04 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gglu( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t,gg,ck,ie)

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  USE IOMod
  use stkets
  use fms_mod, only: gg_full

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
  complex, intent(in) :: g0( istatx, istatx)

  ! Outputs:
  ! Big work matrices
  complex, intent(out) :: g0t( istatx, istatx)
  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

  ! Definitions added to satisfy implicit none:
  integer :: i,j,ip,icase,is1,is2,ist1,ist2  ! Loop indecies
  integer :: isp1,jj1,jj2,iatom1,iatom2,ipart,info,nstates
  integer :: l1,m1,m2

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
           g0t(i, j) = g0t(i, j)                          &
                &                   - g0(i, ist2) * tmatrx(nsp, j)
        endif
     enddo

     g0t(j, j) = g0t(j, j) + one

  enddo


  !         open(24,file='tmatrixr.txt')
  !	   write(24,*) tmatrx(1,:)
  !	   close(24)
  !      if(makeg) then
  !         open(24,file='tmatrixr.txt')
  !	   write(24,*) tmatrx(1,:)
  !	   close(24)
  !       matrix inversion
  !         trans = 'NotTransposed'
  !         call cgetrf(istate,istate,g0t,istatx,ipiv,info)
  !	   g0ti=cmplx(0,0)
  !	   do i=1,istatx; g0ti(i,i)=cmplx(1,0); enddo
  !         call cgetrs(trans,istate,istate,g0t,istatx,ipiv,g0ti,istatx,
  !     1       info)
  !         goto 199
  !	endif

  ! --- invert matrix by LU decomposition
  !     call cgetrf from lapack.  this performs an LU decomposition on
  !     the matrix g0t = 1 - g0*T
  call cgetrf( istate, istate, g0t, istatx, ipiv, info )
  if (info.lt.0) then
     call wlog('    *** Error in cgetrf when computing G')
     write(cerr,400)abs(info)
     call wlog('        Argument #'//cerr//  ' had an illegal value.')
  elseif (info.gt.0) then
     call wlog('    *** Error in cgetrf when computing G')
     write(cerr,400)info
     call wlog('        g0t('//cerr// ','//cerr// ') is exactly 0 -- '//'this matrix cannot be decomposed.')
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
     !KJ test      do 620 ip=0,1
     !KJ test      i0(0)=0;i0(1)=9
     ipart = nsp*(lipotx(ip)+1)**2
     do is1 = 1, istate
        do is2 = 1, ipart
           g0s(is1,is2) = g0(is1, is2 + i0(ip))
        enddo
     enddo


     trans = 'NotTransposed'

     call cgetrs(trans, istate, ipart, g0t, istatx,  ipiv, g0s, istatx, info)

     if (info.lt.0) then
        call wlog('    *** Error in cgetrf')
        write(cerr,400) abs(info)
        call wlog('        Argument #'//cerr//  ' had an invalid value.')
     endif

     ! **** at this point g0s contains the full MS ****

     !  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each ipot

     do is2=1,ipart
        do is1=1,ipart
           gg( is1, is2, ip) = g0s( is1+i0(ip), is2)
        enddo
     enddo

  enddo IP_LOOP

  IF(FullG) THEN
     call cgetrs(trans, istate, istate, g0t, istatx, ipiv, g0, istatx, info)
     CALL WriteData('g0.bin',Int1 = ie)
     CALL Write2D('g0.bin',g0(:istate,:istate))
  END IF

  ! if requested, calculate full FMS matrix
  if (allocated(gg_full)) then
     gg_full(:,:) = g0(:,:)
     call cgetrs(trans, istate, istate, g0t, istatx, ipiv, gg_full, istatx, info)
  end if

  if(makeg) then
     ICASE_LOOP: do icase=1,4
        g0s=cmplx(0,0)
        if(nsp.eq.1) nstates=9 !9 diamond !16 Po
        if(nsp.eq.2) nstates=18 !18 diamond !32 Po
        iatom1=1
        if(icase.eq.1) then
           iatom2=7
           open (99,file='g17r175at.txt',position='append')
        elseif(icase.eq.2) then
           iatom2=36
           open (99,file='g136r175at.txt',position='append')
        elseif(icase.eq.3) then
           iatom2=35
           open (99,file='g135r175at.txt',position='append')
        elseif(icase.eq.4) then
           iatom2=58
           open (99,file='g158r175at.txt',position='append')
        endif
        jj1=(iatom1-1)*nstates ; jj2=(iatom2-1)*nstates


        ipart = nsp*(lipotx(1)+1)**2
        do is1 = 1, istate
           do is2 = 1, ipart
              g0s(is1,is2) = g0(is1, is2 + jj2)
           enddo
        enddo

        trans = 'NotTransposed'
        call cgetrs(trans, istate, ipart, g0t, istatx,                 &
             &                ipiv, g0s, istatx, info)
        if (info.lt.0) then
           call wlog('    *** Error in cgetrf')
           write(cerr,400) abs(info)
           call wlog('        Argument #'//cerr//                      &
                &              ' had an invalid value.')
        endif

        ! **** at this point g0s contains the full MS ****
        q=cmplx(0,0)
        do ist1=1,nstates
           do ist2=1,nstates
              q(ist1,ist2)=g0s(jj1+ist1,ist2)
              ! if (cabs(q(ist1,ist2)).lt.0.00001) q(ist1,ist2)=cmplx(0,0)
           enddo
        enddo
        write(99,167) ck,(q(ist1,1:nstates),ist1=1,nstates)
        close(99)

167     format(5000(e12.4,1x,e12.4,3x))
        !KJ
     enddo ICASE_LOOP
  endif


  !	   q=cmplx(0,0)
  !	   do i1=1,8*nstates
  !	   do i2=1,8*nstates
  !	      do i=1,istate
  !	         q(i1,i2)=q(i1,i2)+g0ti(i1,i)*g0(i,i2)
  !	      enddo
  !!	      if (cabs(q(i1,i2)).lt.0.00001) q(i1,i2)=cmplx(0,0)
  !	   enddo
  !	   enddo
  !	   call writematrixfree(q,8*nstates,'q.txt')


  ! Deallocate local variables:
  !deallocate(ipiv, g0s)


  return
end subroutine gglu
