!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_dimsmod.f90,v $:
! $Revision: 1.31 $
! $Author: jorissen $
! $Date: 2013/01/27 21:04:26 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module DimsMod
  ! This module contains dimensions for data arrays

! The file in which dimensions current to the calculation are saved :
  character*20, parameter :: dimFName = '.dimensions.dat'
! Set the following according to max. available memory on your system :
! The hardcoded limit on cluster size that can NEVER be exceeded :
  integer, parameter :: nclusxhardlimit = 3000
! The hardcoded upper limit on l-values that can NEVER be exceeded :
  integer, parameter :: lxhardlimit = 20
! The hardcoded upper limit on the number of potentials that can NEVER be exceeded :
  integer, parameter :: nphxhardlimit = 31
! The hardcoded upper limit on the number of spins that can NEVER be exceeded :
  integer, parameter :: nspxhardlimit = 2  ! duh
  private dimFName, lxhardlimit, nclusxhardlimit, nspxhardlimit

  integer,parameter :: nclxtd = 100     ! Maximum number of atoms for tdlda module.
  integer,parameter :: nspx   = 2      ! Max number of spins: 1 for spin average; 2 for spin-dep
  integer,parameter :: natx   = 4000    ! Max number of atoms in problem for the pathfinder and ffsort
  integer,parameter :: nattx  = 4000    ! Max number of atoms in problem for the rdinp
!  integer,parameter :: nphx   = 14      ! Max number of unique potentials (potph)
  integer,parameter :: ltot   = 24      ! Max number of ang mom (arrays 1:ltot+1)
  integer,parameter :: nrptx  = 1251    ! Loucks r grid used through overlap and in phase work arrays
  integer,parameter :: nex    = 2000     ! Number of energy points genfmt, etc.
  integer,parameter :: lamtot = 15      ! Max number of distinct lambda's for genfmt 15 handles iord 2 and exact ss
  integer,parameter :: mtot   = 4       ! Vary mmax and nmax independently
  integer,parameter :: ntot   = 2 
  integer,parameter :: npatx  = 8       ! Max number of path atoms, used  in path finder, NOT in genfmt
  integer,parameter :: legtot = npatx+1 ! Matches path finder, used in GENFMT
  integer,parameter :: novrx  = 8       ! Max number of overlap shells (OVERLAP card)
  integer,parameter :: nheadx = 20+nphxhardlimit      ! Max number of header lines !KJ 7-09 added term to accomodate large systems in xsect.bin header
  integer,parameter :: MxPole = 1000    ! Max number of poles that can be used to model epsilon^-1 for HL multipole self energy
  integer,parameter :: nwordx = max(100,2+2*nphxhardlimit)     ! An infuriatingly stupid parameter that shows up in a few places. KJ added 7-09.  used to be 20 - must be at least 2*(1+nphx) for feff.bin header.
  integer,parameter :: novp = 50 ! For istprm, movrlp, ovp2mt - an atom list cutoff that should be high enough to include one atom of each potential type.  Added 2-2011 !KJ

! Added March 2014 for path/ff2x:
  integer,parameter :: nx = 5000000 ! size of stack of paths for pathfinder
  integer,parameter :: np1x = nx
  integer,parameter :: npx_path = 5000001  ! max number of retained paths for pathfinder.  I expect it's irrelevant if > nx.  Not used for dimensioning.
  integer,parameter :: npx_ff2x = 2000 ! max number of paths in ff2x.  Note that these are *unique paths after reducing for symmetry.  Hence can usually be much smaller than above limits.
  integer,parameter :: lx_xsph = 6

! NON PARAMETER STATEMENTS
  integer :: nclusx    ! Maximum number of atoms for FMS.
  integer :: lx        ! Max orbital momentum for FMS module.
  integer :: nspu      ! Max number of spin states (1 or 2)
  integer :: nphu      ! Max number of potential types

  ! OLD XPARAM.H MODULE
  integer,parameter :: natxx = natx
  integer,parameter :: nexx = nex
  integer,parameter :: nkmin = 1
  !integer,parameter :: nphasx = nphx
  integer :: nphasx
  integer :: istatx

contains

  subroutine set_dimensions_for_rdinp
    ! We need some non-zero values for the rdinp module
    ! These will be used before we know the true number of atoms etc.
    ! They will be overwritten once we have read the input
    ! by calling write_dimensions
    nclusx = nclusxhardlimit
    lx = lxhardlimit
    nphu = nphxhardlimit
    nspu = nspxhardlimit
    nphasx = nphu
  end subroutine set_dimensions_for_rdinp


  subroutine write_dimensions(nclusxuserlimit,lxuserlimit)
    implicit none
    ! Write dimension data to a file
	integer,intent(in) :: nclusxuserlimit,lxuserlimit
    integer :: ios  ! IO Status

!   3/ Apply hardcoded dimension limits
    if(nclusxuserlimit.gt.0) then
	   nclusx=min(nclusx,nclusxuserlimit)
	else
	   nclusx=min(nclusx,nclusxhardlimit)
	endif
    if(lxuserlimit.ge.0) then
	   lx=min(lx,lxuserlimit)
	else
       lx=min(lx,lxhardlimit)
	endif

    !sanity checks with defaults:
    if (nspu.lt.1 .or. nspu.gt.2) nspu=1
    if (nphu.lt.1) nphu=1
    if (nphu.gt.nphxhardlimit) then
       nphu=nphxhardlimit
       call wlog('Reducing nphx because it was larger than nphxhardlimit')
       call wlog('Bad crash expected.')
    endif

    open(10,FILE=trim(dimFName),STATUS='unknown',FORM='formatted',IOSTAT=ios)
    call chopen(ios,trim(dimFName),'dimsmod')
    if (ios.ne.0) stop "Error writing dimensions.dat.  Quitting."
    
    write(10,*) nclusx,lx,nphu,nspu
    close(10)
  end subroutine write_dimensions



  subroutine dump_dimensions
    implicit none
    ! Write ALL dimensions to diagnostic file for debugging

    open(90,file='dimensions_debug.dat',form='formatted',status='unknown')
    write(90,*) '********************************'
    write(90,*) 'nclusx',nclusx
    write(90,*) 'lx',lx
    write(90,*) 'nphu',nphu
    write(90,*) 'nspu',nspu
    write(90,*) '********************************'
    write(90,*) 'lxhardlimit',lxhardlimit
    write(90,*) 'nphxhardlimit',nphxhardlimit
    write(90,*) 'nclxtd',nclxtd
    write(90,*) 'nspx',nspx
    write(90,*) 'natx',natx
    write(90,*) 'nattx',nattx
!    write(90,*) 'nphx',nphx
    write(90,*) 'ltot',ltot
    write(90,*) 'nrptx',nrptx
    write(90,*) 'nex',nex
    write(90,*) 'lamtot',lamtot
    write(90,*) 'mtot',mtot
    write(90,*) 'ntot',ntot
    write(90,*) 'npatx',npatx
    write(90,*) 'legtot',legtot
    write(90,*) 'novrx',novrx
    write(90,*) 'nheadx',nheadx
    write(90,*) 'MxPole',MxPole
    write(90,*) 'nwordx',nwordx
    write(90,*) 'novp',novp
    write(90,*) 'natxx',natxx
    write(90,*) 'nexx',nexx
    write(90,*) 'nkmin',nkmin
    write(90,*) 'nphasx',nphasx
    write(90,*) 'istatx',istatx
    close(90)
    return
  end subroutine dump_dimensions



  subroutine init_dimensions
    implicit none
    ! Read dimensions from file
    integer :: ios  ! IO Status
    open(10,FILE=trim(dimFName),STATUS='old',FORM='formatted',IOSTAT=ios)
    call chopen(ios,trim(dimFName),'dimsmod')
    read(10,*) nclusx,lx,nphu,nspu
    close(10)
    ! OLD XPARAM DIMENSIONS
    istatx=(lx+1)**2*nclusx*nspx
    nphasx=nphu
  end subroutine init_dimensions


end module DimsMod
