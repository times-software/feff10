!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xsph.f90,v $:
! $Revision: 1.25 $
! $Author: jorissen $
! $Date: 2012/05/21 23:45:17 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program  ffmod2

! cross-section and phase shifts calculations
! coded by a.ankudinov 2000 ; modules added KJ 2009
! INPUT: mod2.inp geom.dat global.inp and pot.bin
! OUTPUT: xsect.bin and xsph.bin

  use DimsMod, only: lx, nclusx, istatx, init_dimensions
  use stkets
  use rotx
  use lnlm
  use xstruc
  use t3j
  use par
  use xsph_inp, only: mphase
  use errorfile
  implicit none
! locals:
  integer ios

  call par_begin
  if(master) call OpenErrorfileAtLaunch('xsph')
  if (worker) go to 400
!  Initialize needed modules !KJ
  call init_dimensions
  call init_stkets(istatx)
  call init_rotx(lx,nclusx)
  call init_lnlm(lx,nclusx)
  call init_xstruc(nclusx)
  call init_t3j(lx)


! open the log file, unit 11.  See subroutine wlog.
  open (unit=11, file='log2.dat', status='unknown', iostat=ios)
  call chopen (ios, 'log2.dat', 'feff')

! read  INPUT data files: geom.dat, global.dat and mod2.inp.
  call rexsph     ! data passed in modules
  if (mphase .eq. 1)  then
     call wlog('Calculating cross-section and phases ...')
     call xsph
     call wlog('Done with module: cross-section and phases (XSPH).'//char(13)//char(10))
  endif

  close (unit=11)
  400 call par_barrier
  call par_end

  if(master)call WipeErrorfileAtFinish
  stop

  end



!KJ ADAPTATIONS FOR FEFFQ :    ! 7-09
! -some extra variables in m_inpmodules and m_nrixs
! -wrxsph generalized a little; slight change in phase.bin
! -if-block in xsph ; calls xsectjas/xsect for nrixs/no nrixs
! -tdlda not possible for nrixs.  probably easy to fix.
! -new type of energy grid generated when xkmax<0.




! OUTPUT: data for the next modules is written in phase.bin
! auxiliary output can be obtained using 'ipr2' (see feff8.2 manual)

! necessary input information from feff.inp file
! see CARDs description in feff8 manual for more details
! CONTROL mphase: 1-run (0-don't run)  the program
! PRINT ipr2: for auxialry output files (default=0)
! ispec: spectroscopy type (EXAFS, XANES, XES, DANES, FPRIME) 
! vixan, xkstep, xkmax: energy grid for chosen spectroscopy
! RDRIG rgrid: radial grid (default=0.05)
! POTENTIAL info
!   nph: number of unique potentials
!   lmaxph: max orbital momentum for xsph calculations
!   potlbl: labels for unique potentials
!  ATOMS
!   nat: number of atoms
!   rat: their coordinates
!   iphat: type of potential for each site
!   iatph: representative atoms indices in atoms list
!  EXCHANGE ixc  vr0  vi0  ixc0 - exchange correlation model
!  RSIGMA (RPHASES) lreal (default=0)
!  FMS  rfms2 lfms2

!  Global data
!    ipol - polarization type (default:0 - polarization average)
!    ispin - spin type (default=0 - spin independent)
!    le2 - include/exclude quad. transitions (default=2 - include)
!    angks - angle between x-ray propagation and spin (default=0)
!    ptz - polarization tenzor (default=0 for ipol=0)
