!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: gglufullpot.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gglufullpot( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t,gg,ck)

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx
  use constants
  use stkets

  implicit none
  complex ck  !KJ debugging
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential

  ! Inputs
  integer, intent(in) :: i0(0:nphx), lipotx(0:nphx)
  integer, intent(in) :: nsp
  integer, intent(in) :: ipi,ipf
  complex, intent(in) :: tmatrx(istatx, istatx)
  complex, intent(in) :: g0( istatx, istatx)

  ! Outputs:
  ! Big work matrices
  complex, intent(out) :: g0t( istatx, istatx)
  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

  ! Definitions added to satisfy implicit none:
  integer :: i,ip,icase,is1,is2,ist1,ist2  ! Loop indecies
  integer :: jj1,jj2,iatom1,iatom2,ipart,info,nstates

  integer, allocatable :: ipiv(:)

  ! Big work matrices
  complex, allocatable :: g0s(:,:)
  complex q(nspx*(lx+1)**2*8,nspx*(lx+1)**2*8) !KJ debug
  !     return matrix containing info about each unique potential

  character*3  cerr
  character*13 trans
  logical,parameter :: makeg=.false.
400 format(i4)

  ! Allocate local variables:
  allocate(ipiv(istatx), g0s( istatx, nspx*(lx+1)**2))

  !KJ  calling cgemm is wasteful, because tmatrx is block-diagonal.  Programmer's laziness!
  call cgemm('N','N',istate,istate,istate,-cmplx(1,0),g0,istatx,tmatrx,istatx,&
       &                   cmplx(0,0),g0t,istatx) !computes g0t = - g0 tmatrx
  do i=1,istate
     g0t(i,i)=g0(i,i)+cmplx(1,0)
  enddo
  !     now g0t = 1 - g0 tmatrx


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
     call cgetrs(trans, istate, ipart, g0t, istatx,                  &
          &                ipiv, g0s, istatx, info)
     if (info.lt.0) then
        call wlog('    *** Error in cgetrf')
        write(cerr,400) abs(info)
        call wlog('        Argument #'//cerr//                      &
             &              ' had an invalid value.')
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
              !	 if (cabs(q(ist1,ist2)).lt.0.00001) q(ist1,ist2)=cmplx(0,0)
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
  deallocate(ipiv, g0s)

  return
end subroutine gglufullpot
