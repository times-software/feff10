!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: istval.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine istval (vtot, rholap, rmt, imt, rws, iws, vint, rhoint,&
     &                   ierr)

!     This subroutine calculates interstitial values of v and rho
!     for an overlapped atom.  Inputs are everything except vint and
!     rhoint, which are returned.  vtot includes ground state xc.
!     rhoint is form density*4*pi, same as rholap
!
!     ierr = 0, normal exit
!          =-1, rmt=rws, no calculation possible
        use DimsMod, only: nrptx

      implicit double precision (a-h, o-z)

      parameter (delta = 0.050000000000000)

      dimension vtot (nrptx)
      dimension rholap (nrptx)

!     Integrations are done in x (r = exp(x), see Louck's grid)
!     Trapezoidal rule, end caps use linear interpolation.
!     imt is grid point immediately below rmt, etc.
!     We will integrate over spherical shell and divide by volume of
!     shell, so leave out factor 4pi, vol = r**3/3, not 4pi*r**3/3,
!     similarly leave out 4pi in integration.

!     If rmt and rws are the same, cannot contribute to interstitial
!     stuff, set error flag
      vol = (rws**3 - rmt**3) / 3
      if (vol .le. 0)  then
         ierr = -1
         return
      endif
      ierr = 0

!     Calculation of vint including exchange correlation
!     Trapezoidal rule from imt+1 to iws
      vint = 0
      do 100  i = imt, iws-1
         fr = rr(i+1)**3 * vtot(i+1)
         fl = rr(i)**3   * vtot(i)
         vint = vint + (fr+fl)*delta/2
  100 continue
!     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * vtot(iws+1)
      fl = rr(iws)**3   * vtot(iws)
      vint = vint + (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
!     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * vtot(imt+1)
      fl = rr(imt)**3   * vtot(imt)
      vint = vint - (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
      vint = vint / vol

!     Calculation of rhoint
!     Trapezoidal rule from imt+1 to iws
      rhoint = 0
      do 200  i = imt, iws-1
         fr = rr(i+1)**3 * rholap(i+1)
         fl = rr(i)**3   * rholap(i)
         rhoint = rhoint + (fr+fl)*delta/2
  200 continue
!     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * rholap(iws+1)
      fl = rr(iws)**3   * rholap(iws)
      rhoint = rhoint + (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
!     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * rholap(imt+1)
      fl = rr(imt)**3   * rholap(imt)
      rhoint = rhoint - (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
      rhoint = rhoint / vol

      return
      end
