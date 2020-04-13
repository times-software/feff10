!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fixdsp.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fixdsp (dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

  !     This fixes up the dirac spinor components (dgc and dpc) from ATOM
  !     for the xsect code.
  use DimsMod, only: nrptx
  use ifuns
  use constants


!  implicit double precision (a-h, o-z)
  implicit none
  
  real*8, dimension(251)   :: dgc0, dpc0
  real*8, dimension(nrptx) :: dgcx, dpcx, xorg, xnew

  ! Added to satisfy implicit none
  integer :: i,imax,j,jmax,jnew
  real*8 :: dxorg,dxnew,rmax


  !     Use linear interpolation in x whether necessary or not.  If
  !     new grid is same as old, it shouldn't make any difference.

  !     relation between x, r, and j.  xx00 = 8.8 for all grids
  !     in this version, change it if more flexibility is necessary.
  !     xx = -xx00 + (j-1)*delta
  !     rr = exp (xx)
  !     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

  !     The dgc and dpc arrays are zero beyond a certain point, usually
  !     inside the muffin tin radius.  Find this distance.
  I_LOOP: do i = 251, 1, -1
     if ( abs(dgc0(i)) .ge. 1.0d-11 .or.                            &
          &        abs(dpc0(i)) .ge. 1.0d-11 )  then
        imax = i
        exit I_LOOP
     endif
  enddo I_LOOP


  !     jmax is the first point where both dpc and dgc are zero in
  !     the original grid
  jmax = imax + 1
  if (jmax.gt.251) jmax = 251

  do j = 1, jmax
     xorg(j) = xxx(j,dxorg)
  enddo

  rmax = rrr(jmax,dxorg)

  !     How far out do we go in the new grid?  To the last new grid
  !     point before jmax.  Everything will be zero beyond jmax.
  jnew = jjj(rmax,dxnew)
  do j = 1, jnew
     xnew(j) = xxx(j,dxnew)
  enddo

  !     interpolate to new grid using x, only inside of rmax
  do j = 1, jnew
     call terp (xorg, dgc0,  jmax, 3, xnew(j), dgcx(j))
     call terp (xorg, dpc0,  jmax, 3, xnew(j), dpcx(j))
  enddo

  !     and zero the arrays past rmax
  do j = jnew+1, nrptx
     dgcx(j) = 0
     dpcx(j) = 0
  enddo

  return
end subroutine fixdsp
