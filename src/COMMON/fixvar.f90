!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fixvar.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/12/15 17:13:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fixvar (rmt, edens, vtot, dmag,                        &
     &             vint, rhoint, dxorg, dxnew, jumprm,            &
     &             vjump, ri, vtotph, rhoph, dmagx)

  ! edens(xorg) => rhoph(xnew), vtot(xorg) => vtotph(xnew), dmag(xorg) => dmagx(xnew)

  use DimsMod, only: nrptx
  use constants
  use ifuns

  implicit none

  real*8, dimension(251)   :: edens,vtot,dmag
  real*8, dimension(nrptx) :: vtotph,rhoph,dmagx,ri,xorg,xnew

  ! Added to satisfy implicit none:
  integer :: j
  integer :: jmtorg,jmtnew
  integer :: jriorg,jrior1,jrinew,jrine1
  integer :: jumprm

  real*8 :: dxnew,xmt,rmt,vjump,vint,rhoint,dxorg,vmt


  !     PHASE needs
  !     vtot = total potential including gs xcorr, no r**2
  !     edens = rho, charge density, no factor of 4*pi, no r**2
  !     From overlapping, vtot = potential only, ok as is
  !                       edens = density*4*pi, so fix this here.
  !     ri = r grid through imt+1

  !     Only values inside the muffin tin are used, except that XCPOT
  !     (in PHASE) uses values at imt+1 and requires these to be the
  !     interstitial values.  So set the last part of the arrays to
  !     interstitial values...

  !     Use linear interpolation in x whether necessary or not.  If
  !     new grid is same as old, it shouldn't make any difference.

  !     relation between x, r, and j.  xx00 = 8.8 for all grids
  !     in this version, change it if more flexibility is necessary.

  !     xx = -xx00 + (j-1)*delta
  !     rr = exp (xx)
  !     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

  jmtorg = jjj(rmt,dxorg)
  jriorg = jmtorg + 1
  jrior1 = jriorg + 1
  do j = 1, jrior1
     xorg(j) = xxx(j,dxorg)
  enddo

  jmtnew = jjj(rmt,dxnew)
  jrinew = jmtnew + 1
  jrine1 = jrinew + 1
  do j = 1, jrine1
     xnew(j) = xxx(j,dxnew)
  enddo

  !     interpolate to new grid using x, only inside of muffintin
  !     jri (first interstitial point) must be set to interstitial value
  do j = 1, jrinew
     call terp (xorg, vtot,  jriorg, 3, xnew(j), vtotph(j))  ! returns vtotph(xnew) given vtot(xorg)
     call terp (xorg, edens, jrior1, 3, xnew(j), rhoph(j))
     call terp (xorg, dmag,  jrior1, 3, xnew(j), dmagx(j))
  enddo

  if (jumprm .eq. 1) then
     xmt = log(rmt)
     call terp (xorg, vtot,  jriorg, 3, xmt, vmt)
     vjump = vint - vmt
  endif
  if (jumprm .gt. 0) then
     do j = 1, jrinew
        vtotph(j) = vtotph(j) + vjump
     enddo
  endif

  do j = 1, nrptx
     ri(j) = rrr(j,dxnew)
  enddo

  do j = 1, jrinew
     rhoph(j) = rhoph(j)/(4*pi)
  enddo
  do j = jrinew+1, nrptx
     vtotph(j) = vint
     rhoph(j) = rhoint/(4*pi)
     ! fix later : need to calculate interstitial dmint
     !      want interpolation beyond mt also
     dmagx(j) = 0.0d0
  enddo

  return
end subroutine fixvar
