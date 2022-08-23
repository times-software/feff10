!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdgeom.f90,v $:
! $Revision: 1.10 $
! $Author: hebhop $
! $Date: 2012/01/07 00:55:48 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     Read GEOM.INP and LDOS.INP and SCREEN.INP
!=======================================================================
subroutine rdscreen
!KJ all variables are in geom.inp ie module atoms_inp ; I removed the calling argument ibounc since it's not set here.
!KJ Added Yoshi's 'rdldos' and 'rdscrn' to this file also.
  use DimsMod, only: istatx, lx, nclusx, init_dimensions
  use constants
  use screenw_inp
  implicit none
!KJ rgrd and ixc accessible through potential_inp via ldos_inp ; nph in atoms ; ispin in global via ldos

  ! Initialize needed modules
  call init_dimensions
  call screenw_read  !KJ emin/emax from screen.inp override values set in ldos.inp as in Yoshi's original version
  !ScreenI%maxl=min(ScreenI%maxl,lx+1) !KJ Yoshi didn't have this - not sure why!  His default maxl=4 causes trouble if lx=2 (or smaller). 

  !     transform to code units (bohrs and hartrees - atomic units)
  ScreenI%emin   = ScreenI%emin   / hart
  ScreenI%emax   = ScreenI%emax   / hart
  ScreenI%eimax  = ScreenI%eimax  / hart
  ScreenI%omega_max  = ScreenI%omega_max  / hart
  ScreenI%xmu  = ScreenI%xmu  / hart
  ScreenI%gam1  = ScreenI%gam1  / hart
  ScreenI%gam2  = ScreenI%gam2  / hart

  ScreenI%rfms = ScreenI%rfms/bohr
  ScreenI%ermin  = ScreenI%ermin  / hart


  return

end subroutine rdscreen
