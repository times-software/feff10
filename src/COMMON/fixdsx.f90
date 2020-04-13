!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fixdsx.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/12/15 17:13:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fixdsx (iph, dxorg, dxnew, dgc, dpc, dgcn, dpcn)
  ! dgc(xorg) => dgcn(xnew), dpc(xorg) => dpcn(xnew)

  !     This fixes up the dirac spinor components (dgc and dpc) from ATOM
  !     for the xsect and phase codes.
  use DimsMod, only: nphx=>nphu, nrptx
  use constants
  use ifuns

!  implicit double precision (a-h, o-z)
  implicit none

  real*8, dimension(251,30,0:nphx) :: dgc,dpc
  real*8, dimension(nrptx,30) :: dgcn,dpcn
  real*8, dimension(nrptx) :: xorg,xnew

  ! Added to satisfy implicit none
  integer :: iorb,i,imax,j,jmax,jnew,iph
  real*8  :: dxorg,dxnew,rmax

  !     Use linear interpolation in x whether necessary or not.  If
  !     new grid is same as old, it shouldn't make any difference.

  !     relation between x, r, and j.  xx00 = 8.8 for all grids
  !     in this version, change it if more flexibility is necessary.
  !     xx = -xx00 + (j-1)*delta
  !     rr = exp (xx)
  !     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

  !     The dgc and dpc arrays are zero beyond a certain point, usually
  !     inside the muffin tin radius.  Find this distance.

  do j = 1, 251
     xorg(j) = xxx(j,dxorg)
  enddo

  do j = 1, nrptx
     xnew(j) = xxx(j,dxnew)
  enddo

  IORB_LOOP: do iorb = 1, 30
     imax = 0
     I_LOOP: do i = 251, 1, -1
        if ( abs(dgc(i,iorb,iph)) .ge. 1.0d-11 .or. abs(dpc(i,iorb,iph)) .ge. 1.0d-11 )  then
           imax = i
           exit I_LOOP
        endif
     enddo I_LOOP

     if (imax .eq. 0) then
        jnew = 0
        goto 35
     endif
     !        jmax is the first point where both dpc and dgc are zero in
     !        the original grid
     jmax = imax + 1
     if (jmax .gt. 251) jmax = 251

     rmax = rrr(jmax,dxorg)

     !        How far out do we go in the new grid?  To the last new grid
     !        point before jmax.  Everything will be zero beyond jmax.
     jnew = jjj(rmax,dxnew)

     !        interpolate to new grid using x, only inside of rmax
     do j = 1, jnew
        call terp(xorg,dgc(1,iorb,iph),jmax,3, xnew(j),dgcn(j,iorb)) ! makes dgcn(xnew) given dgc(xorg)
        call terp(xorg,dpc(1,iorb,iph),jmax,3, xnew(j),dpcn(j,iorb))
     enddo

     !        and zero the arrays past rmax
35   do j = jnew+1, nrptx
        dgcn(j,iorb) = 0
        dpcn(j,iorb) = 0
     enddo

  enddo IORB_LOOP

  return
end subroutine fixdsx
