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
subroutine rdgeom
!KJ all variables are in geom.inp ie module atoms_inp ; I removed the calling argument ibounc since it's not set here.
!KJ Added Yoshi's 'rdldos' and 'rdscrn' to this file also.
  use DimsMod, only: istatx, lx, nclusx, init_dimensions
  use constants
  use atoms_inp
  use ldos_inp
  use screen_inp
  use stkets
  use rotx
  use lnlm
  use xstruc
  use t3j
  implicit none
!KJ rgrd and ixc accessible through potential_inp via ldos_inp ; nph in atoms ; ispin in global via ldos

  ! Initialize needed modules
  call init_dimensions
  call init_stkets(istatx)
  call init_rotx(lx,nclusx)
  call init_lnlm(lx,nclusx)
  call init_xstruc(nclusx)
  call init_t3j(lx)

  call atoms_read
  call ldos_read
  call screen_read  !KJ emin/emax from screen.inp override values set in ldos.inp as in Yoshi's original version
  ScreenI%maxl=min(ScreenI%maxl,lx+1) !KJ Yoshi didn't have this - not sure why!  His default maxl=4 causes trouble if lx=2 (or smaller). 

  !     transform to code units (bohrs and hartrees - atomic units)
  rat=rat/bohr

  rfms2  = rfms2  / bohr
  
  rdirec = rdirec / bohr
  ScreenI%emin   = ScreenI%emin   / hart
  ScreenI%emax   = ScreenI%emax   / hart
  ScreenI%eimax  = ScreenI%eimax  / hart

  ScreenI%rfms = ScreenI%rfms/bohr
  ScreenI%ermin  = ScreenI%ermin  / hart


  return

end subroutine rdgeom
