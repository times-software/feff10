!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdinp.f90,v $:
! $Revision: 1.85 $
! $Author: jorissen $
! $Date: 2013/01/27 23:17:50 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sub-program exchange point
      program rdinp
!     subroutine rdinp (nabs,nss,ceels)

!    reads 'feff.inp' file and writes several files in special format
!    ready for the use by other modules: geom.dat, global.dat,
!    mod1.inp, mod2.inp, mod3.inp mod4.inp mod5.inp mod6.inp ldos.inp .
!    The subroutine output 'nabs' is needed for configurational average
!    The rest of output, passed to wrtall via modules (COMMON/m_inpmodules.f90).

!     coded s. zabinski 1994
!     last modified by a.l.ankudinov march 2001  for new i/o structure
!     introduced k-space and eels ; introduced modules, merged feffq ; improved dynamic allocation - Kevin Jorissen 7-09


!KJ restructured input data 7-09 :
    use par
    use constants
    use geometry_inp
    use global_inp
    use reciprocal_inp
    use potential_inp
    use ldos_inp
    use opcons_inp
    use xsph_inp
    use fms_inp
    use paths_inp
    use genfmt_inp
    use ff2x_inp
    use sfconv_inp
    use eels_inp
    use rixs_inp
    use compton_inp
    use dimsmod, only: nheadx, nrptx, nspx, nspu, lx, nclusx, nphu, &
      nphxhardlimit, nwordx, write_dimensions, set_dimensions_for_rdinp
    use screen_inp,only:screen_inp_parse
    use crpa_inp
    use band_inp, emin_band=>emin,emax_band=>emax,estep_band=>estep,nkp_band=>nkp  !avoid contamination between modules
    use errorfile
    use hubbard_inp
    use fullspectrum_inp

      implicit none
      include '../HEADERS/vers.h'

!     Single scattering path to go with Overlap information
      integer, parameter :: nssx = 16
      integer indss(nssx), iphss(nssx)
      real*8 degss(nssx), rss(nssx)

!     Local stuff
!      integer,parameter :: nwordx = 20
      integer,parameter :: nbr=100  !KJ 12-2011 max number of SCF iterations.  Changed from 30 to 100.
      integer,parameter :: big = 1.0e5
      character*512 slog
      character*150  line
      character*120 words(nwordx)
      character*20 symfil !KJ file that contains symmetry for k-mesh 11-06
      character*12 tmpstr
      character*6 str6 ! dummy
      character*3 str3
      character*2 edgename
      character*10,external :: name_spectroscopy,name_corehole
      character*3000,external :: list_cards,list_features
      integer ltit(nheadx),iss,iatabs,nttl,iabs,nat
      integer icnt,ios,iatom,ifolp,iovrlp,lxnat,nss,mode,jinit,nwords,itok,nph_read
      integer icoord,i,j,k,iph,iovr,ltmp,iatrd,i1,i2,i3,j1,j2,j3,indexabs,iat,icount,iq,iqq
      logical nogeom,userchl
      logical ceels  !KJ for monolithic version 5-6
      integer mpathold !KJ to fix OVERLAP card 7-06
      real*8 sss(3),alatt,xxx(3),distance(nattx),scalelattice,mindist,qvec(3),dummy,dummy2,userChLifetime
      real*8 folpx,rmult,s02h,tmp,rdims,ratmin,ratmax,xinorm,xnat,cosmdff_dum,magnifier,ratomslist
      real*8 shift(3,3)
      integer nshift,lattice_factor
      integer nclusxuserlimit,lxuserlimit !KJ 7-09
      logical cards_set(150) !KJ needs be big enough to have a field for each programmed card 7-09 for consistency checker
      logical cifread !KJ 10-2011 Take crystal structure from .cif file
      character*120 cifname !KJ 10-2011 Name of .cif file
      integer cif_equivalence !KJ 1-2012 for making potential types from .cif file
!     integer nq !KJ now in global_inp
      logical,parameter :: enforce_alexis_exchange_policy = .false. !KJ 11-2010 - I don't know why Aleksi had it in the first place, just keeping this as precaution
      logical,parameter :: debug_cif = .false.
      integer, allocatable :: iatph(:)

! Added by Fer
! Needed for DMDW
      logical            :: Use_DMDW = .FALSE.
      integer            :: DMDW_Order, DMDW_Type, DMDW_Route
      integer            :: iiAtom, jjAtom
      real               :: mxDij2, Dij2
      character(len=256) :: dym_File

!     Functions :
      integer,external :: itoken,istrln
      real*8,external :: dist, getspin
      integer iTmp

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

      call par_begin
      if (worker) go to 400
      call OpenErrorfileAtLaunch('rdinp')

!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log.dat', 'feff')

      tmpstr = vfeff
      call triml (tmpstr)
      call wlog('Launching FEFF version ' // tmpstr)
      ! Josh adding revision number to feff output.
      !call wlog(' ' // revision)

!     Set dimensions temporarily to the largest allowed values
      call set_dimensions_for_rdinp
      ! Later we will figure out "just-big-enough" values
      allocate (iatph(0:nphxhardlimit))

!     initialize all things to be passed
      call iniall
!KJ   this one left over from iniall, doesn't seem to belong anywhere!
      nat = 0

!     initialize local staff
      iatom = 0
      icoord = 3 !KJ default : SPRKKR coordinates
      ifolp = 0
      iovrlp = 0
      lxnat = 0
      folpx = 1.15d0
      nogeom = .false.
      rclabs = big
      rmult = 1.0d0
      s02h = 1.0d0
      nss = 0
      nclusxuserlimit = -1 !KJ -1 means undef
      lxuserlimit = -1 !KJ
      indss(:) = 0
      iphss(:) = 0
      degss(:) = 0
      rss(:) = 0
      iatph(:) = 0
      symfil='                    '  !KJ 11-06
      cards_set(:)=.false. !KJ 7-09
      userchl=.false. !KJ 6-10 get ch lifetime from setgam
      cifread=.false.
      cif_equivalence=1 !KJ 1-12 use default scheme for setting potentials


!     tokens  0 if not a token
!             1 if ATOM (ATOMS)
!             2 if HOLE
!             3 if OVER (OVERLAP)
!             4 if CONT (CONTROL)
!             5 if EXCH (EXCHANGE)
!             6 if ION
!             7 if TITL (TITLE)
!             8 if FOLP
!             9 if RPATH or RMAX
!            10 if DEBY (DEBYE)
!            11 if RMUL (RMULTIPLIER)
!            12 if SS
!            13 if PRIN (PRINT)
!            14 if POTE (POTENTIALS)
!            15 if NLEG
!            16 if CRIT (CRITERIA)
!            17 if NOGEOM
!            18 if IORDER
!            19 if PCRI (PCRITERIA)
!            20 if SIG2
!            21 if XANE (XANES)
!            22 if CORR (CORRECTIONS)
!            23 if AFOL (AFOLP)
!            24 if EXAF (EXAFS)
!            25 if POLA (POLARIZATION)
!            26 if ELLI (ELLIPTICITY)
!            27 if RGRI (RGRID)
!            28 if RPHA (RPHASES), real phase shifts
!            29 if NSTA (NSTAR), n* for co-linear polarization
!            30 if NOHO (NOHOLE), use no hole for potentials
!            31 if SIG3 third and first cumulants for ss paths
!            32 if JUMP (JUMPRM), remove jumps of potential
!            33 if MBCO (MBCONV), do convolution with exitation spectrum
!            34 if SPIN do calculation for spin-up(down) photoelectron
!            35 if EDGE to specify edge by name
!            36 if SCF  do self-consistency loop
!            37 if FMS  use FMS for cluster of the size rfms
!            38 if LDOS print out l-dos for specified energy range
!            39 if INTE how to find interstitial parameters
!            40 if CFAV to do configuration average
!            41 if S02  to specify S_0^2
!            45 if RSIG (RSIGMA), real self-energy
!            46 if XNCD natural dichroism
!            47 if MULT for quadrupolar etc. transitions
!            48 if UNFR unfreeze f-electrons
!            49 if TDLDA use TDLDA background
!            50 if PMBSE use BSE for background
!            51 if PLASMON       - Added by Josh Kas
!                                - PLASMON
!                                - With this card set, ffmod2 will read exc.dat and
!                                - use a multiple pole self energy
!            52 if SFCO (SFCONV) compute spectral function from response function
!                  and convolve output.
!            53 if SELF print on shell self energy as a function of E.
!            54 if SFSE print off shell self energy and spectral function.
!            55 if RCONV print running convolution with spectral function.
!            56 if ELNE calculate ELNES  !KJ 1-06
!            57 if EXEL calculate EXELFS !KJ 1-06
!            58 if MAGI plot magic angle !KJ 1-06
!            59 if ABSO don't normalize spectrum !KJ 3-06
!            60 if SYMM fix value of icase in PATH module !KJ 6-06
!            61 if REAL work in real space  !KJ 8/06
!            62 if RECIPROCAL work in reciprocal space  !KJ 8/06
!            63 if SGROUP
!            74 if EXTPOT use external mt potentials defined in extpot.aip
!            77 if DIMS specify lx and nclusx for dynamical allocation
!            78 if NRIXS
!            79 if LJMAX convergence parameter for NRIXS
!            80 if LDEC output parameter for NRIXS
!            81 if SETE
!            82 if EPS0 specify dielectric constant to correct exc.dat for MPSE
!            83 if OPCONS make feff create loss.dat from database.
!            84 if NUMDENS use with OPCONS card to specify number densities
!            85 if PREP
!            86 if EGAP
!            87 if CHWIDTH set core hole lifetime manually !KJ 6/2010
!            88 if MDFF calculate mdff - hidden secret option !KJ 11/2010
!            89 if RESTART get starting potentials from pot.inp instead of atomic overlap !KJ 12-2010
!            90 if CONFIG use non-default electronic configuration for some atoms !KJ 12-2010
!            101 if Hubbard corrections
!            107 if SCXC (choose XC for the SCF cycle)
!            -1 if END  (end)
!     mode flag  0 ready to read a keyword card
!                1 reading atom positions
!                2 reading overlap instructions for unique pot
!                3 reading unique potential definitions
!                4 reading EELS input  !KJ
!                5 reading lattice vectors
!                6 reading energy grid
!                7 reading density types

!   call to rdline, which will:
!    1. read from feff.inp if found, otherwise will stop and complain
!       (support for reading from standard input would be easy to add)
!    2. handles line processing tasks like
!         = ignoring comment lines and blank lines
!         = tab removal
!    3. allows 'include' files in input file
!    4. for initial call, set jinit = -1, line = input_file_name
!
      mode  = 0
      jinit = -1
      line  = 'feff.inp'
  200 continue
         call rdline(jinit,line)
         if (line .eq. 'read_line_end')    line='END'
         if (line .eq. 'read_line_error')  line='END'
         words=' '
         nwords = nwordx

         call bwords (line, nwords, words)
         itok = itoken (words(1),'feff.inp')
     if (itok.gt.0) cards_set(itok)=.true.  !KJ 7-09

     !write(*,*) nwords,words(1:nwords)

!        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
!              ATOM
!              Following lines are atom postions, one per line
               mode = 1
               iatom  = iatom  +1
            elseif (itok .eq. 2)  then
!              HOLE     1  1.0
!                   holecode s02
               read(words(2),20,err=900)  ihole
               if (nwords.gt.2) read(words(3),30,err=900)  s02h
               mode = 0
            elseif (itok .eq. 3)  then
!              OVERLAP iph
!                  iph  n  r
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               mode = 2
               iovrlp = iovrlp +1
            elseif (itok .eq. 4)  then
!              CONTROL  mphase, mpath, mfeff, mchi
!               0 - do not run modules, 1 - run module
               if (nwords.eq.5) then
!                 feff7 input file
                  read(words(2),20,err=900)  mpot
                  mphase = mpot
                  mfms = mpot
                  read(words(3),20,err=900)  mpath
                  read(words(4),20,err=900)  mfeff
                  read(words(5),20,err=900)  mchi
               else
!                 feff8 input file
                  read(words(2),20,err=900)  mpot
                  read(words(3),20,err=900)  mphase
                  read(words(4),20,err=900)  mfms
                  read(words(5),20,err=900)  mpath
                  read(words(6),20,err=900)  mfeff
                  read(words(7),20,err=900)  mchi
               endif
               mode = 0
            elseif (itok .eq. 5)  then
!              EXCHANGE  ixc  vr0  vi0 (ixc0)
!              ixc=0  Hedin-Lunqvist + const real & imag part
!              ixc=1  Dirac-Hara + const real & imag part
!              ixc=2  ground state + const real & imag part
!              ixc=3  Dirac-Hara + HL imag part + const real & imag part
!              ixc=5  partially nonlocal: Dirac-Fock for core + HL for
!                     valence electrons, + const real & imag part
!              ixc=10 same as ixc=0 with broadened plasmon HL selfenergy
!              ixc=13 same as ixc=3 with broadened plasmon HL selfenergy
!              ixc=15 same as ixc=5 with broadened plasmon HL selfenergy
!              vr0 is const imag part of potential
!              vi0 is const imag part of potential
!              Default is HL. (ixc=0, vr0=0, vi0=0, ixc0 = 2)
               vr0=0.0
               vi0=0.0
               read(words(2),20,err=900)  ixc
         if (ixc.lt.0) ixc=0  !default
!              if (nwords.ge.3) (read(words(3),30,err=900)  vr0
                read(words(3),30,err=900)  vr0
!              if (nwords.ge.4) read(words(4),30,err=900)  vi0
                read(words(4),30,err=900)  vi0
               if (nwords .gt. 4) read(words(5),20,err=900)  ixc0
               !if (ixc0.lt.0) ixc0=2  !default
               mode = 0
            elseif (itok .eq. 6)  then
!              ION  iph xion(iph)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  xion(iph)
               mode = 0
            elseif (itok .eq. 7)  then
!              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. nheadx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call wlog(' Too many title lines, title ignored')
                  call wlog(' ' // line(1:71))
               endif
               mode = 0
            elseif (itok .eq. 8)  then
!              FOLP iph folp (overlap factor, default 1)
               ifolp = 1
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  folp(iph)
               mode = 0
            elseif (itok .eq. 9)  then
!              RPATH rmax (max r for ss and pathfinder)
               read(words(2),30,err=900)  rmax
               !KJ apparently (at least on mac+ifort) this does not crash for an empty "RPATH" card, but just reads rmax=0 and continues on.
               !   rmax is then set to a default (I think 2x nn distance?) somewhere below.
            elseif (itok .eq. 10)  then
!              DEBYE  temp debye-temp ( idwopt )
!                1     2        3          4
!                          + other if idwopt = 4:
!                            dym_File  DMDW_Order DMDW_Type DMDW_Route
!                                5          6         7         8
!                   temps in kelvin
!                   idwopt = 0 use CD model
!                   idwopt = 1 use EM method
!                   idwopt = 2 use RM method
!                   idwopt = 3 use CL method !KJ 7/06
!                   idwopt = 4 use sig2.dat file !JK (via FDV)
!                   idwopt = 5 use Dynamical Matrix method !FDV
!                   idwopt = -1,-2,... don't calculate DW factors
!                   These add to any sig2 from SIG2 card or files.dat
               read(words(2),30,err=900)  tk
               read(words(3),30,err=900)  thetad
               idwopt=0
               if (nwords.gt.3) then
                 read(words(4),20,err=900)  idwopt

! Added by Fer
! Get the options for the Dynamical Matrix Calculation
                 if (idwopt .eq. 5) then
! Activate the DMDW stuff
                   Use_DMDW = .true.
                   dym_File = "feff.dym"
                   if (nwords .gt. 4) then
!   Get the filename in which to find the dynamical matrix
                     read(words(5),10,err=900) dym_File
                   end if
!   Get the Lanczos recursion order
                   DMDW_Order = 2
                   if (nwords .gt. 5) then
                     read(words(6),20,err=900) DMDW_Order
                   end if
!   Get the type of DMDW calculation to do
                   DMDW_Type = 0
                   if (nwords .gt. 6) then
                     read(words(7),20,err=900) DMDW_Type
                   end if
!   Get the route to determine what to calculate
                   DMDW_Route = 0
                   if (nwords .gt. 7) then
                     read(words(8),20,err=900) DMDW_Route
                   end if

                 end if
                 if (idwopt.gt.5) then   !KJ 7/06 changed 2 to 3.
                                         !Josh - Changed 3 to 4
                                         !FDV - Changed from 4 to 5
                    write(slog,'(a,i5,2x,a)')  ' Option idwopt=',idwopt,'is not available.'
                    call wlog(slog)
                    write(slog,'(a)')   '...setting idwopt=2 to use RM.'
                    call wlog(slog)
                    idwopt = 2
                 endif
               endif
               mode = 0
            elseif (itok .eq. 11)  then
!              RMULTIPLIER  rmult
!              Multiples atom coord, rss, overlap and rmax distances by
!              rmult (default 1).  DOES NOT modify sig2g
               read(words(2),30,err=900)  rmult
               mode = 0
            elseif (itok .eq. 12)  then
!              SS index ipot deg rss
               nss = nss + 1
               if (nss .gt. nssx)  then
                  write(slog,'(a,i8)')                                  &
     &               ' Too many ss paths requested, max is ', nssx
                  call wlog(slog)
                  call par_stop('RDINP')
               endif
               read(words(2),20,err=900)  indss(nss)
               read(words(3),20,err=900)  iphss(nss)
               read(words(4),30,err=900)  degss(nss)
               read(words(5),30,err=900)  rss(nss)
               mode = 0
            elseif (itok .eq. 13)  then
!              PRINT  ipr1  ipr2  ipr3  ipr4 ipr5 ipr6
!              print flags for various modules
!              ipr1 potph  0 pot.bin only
!                          1 add misc.dat
!                          2 add pot.dat
!                          5 add atom.dat
!                          6 add central atom dirac stuff
!                          7 stop after doing central atom dirac stuff
!              ipr2 xsph   0 phase.bin only
!                          2 add  phase.dat
!                          3 add  emesh.dat
!              ipr3 fmstot  currently is dummy
!              ipr4 pathfinder  0 paths.dat only
!                               1 add crit.dat
!                               2 keep geom.dat
!                               3 add fbeta files
!                               5 special magic code, crit&geom only
!                                 not paths.dat.  Use for path studies
!              ipr5 genfmt 0 files.dat, feff.dats that pass 2/3 of
!                            curved wave importance ratio
!                          1 keep all feff.dats
!              ipr6 ff2chi 0 chi.dat
!                          1 add sig2.dat with debye waller factors
!                          2 add chipnnnn.dat for each path
!                          3 add feffnnnn.dat for each path, and
!                            do not add chipnnnn.dat for each path
!                          4 add both feffnnnn.dat and chipnnnn.dat
!                            for each path
               if (nwords.eq.5) then
!                 feff7 input file
                  read(words(2),20,err=900)  ipr1
                  ipr2 = ipr1
                  ipr3 = ipr1
                  read(words(3),20,err=900)  ipr4
                  read(words(4),20,err=900)  ipr5
                  read(words(5),20,err=900)  ipr6
               else
!                 feff8 input file
                  read(words(2),20,err=900)  ipr1
                  read(words(3),20,err=900)  ipr2
                  read(words(4),20,err=900)  ipr3
                  read(words(5),20,err=900)  ipr4
                  read(words(6),20,err=900)  ipr5
                  read(words(7),20,err=900)  ipr6
               endif
               mode = 0
            elseif (itok .eq. 14)  then
!              POTENTIALS
!              Following lines are unique potential defs, 1 per line
               mode = 3
            elseif (itok .eq. 15)  then
!              NLEG nlegmax (for pathfinder)
               read(words(2),20,err=900)  nlegxx
               mode = 0
            elseif (itok .eq. 16)  then
!              CRIT critcw critpw
               read(words(2),30,err=900)  critcw
               read(words(3),30,err=900)  critpw
               mode = 0
            elseif (itok .eq. 17)  then
!              NOGEOM (do not write geom.dat) (disabled)
               nogeom = .true.
               mode = 0
            elseif (itok .eq. 18)  then
!              IORDER  iorder (used in genfmt, see setlam for meaning)
               read(words(2),20,err=900)  iorder
               mode = 0
            elseif (itok .eq. 19)  then
!              PCRIT  pcritk pcrith
!                     (keep and heap criteria for pathfinder)
               read(words(2),30,err=900)  pcritk
               read(words(3),30,err=900)  pcrith
               mode = 0
            elseif (itok .eq. 20)  then
!              SIG2  sig2g   global sig2 used by ff2chi, summed with
!              correlated debye model if DEBYE card used, and with
!              sig2 from files.dat if non-zero.
!              Units are Ang**2
               read(words(2),30,err=900)  sig2g
               mode = 0
            elseif (itok .eq. 21)  then
!              XANES ( xkmax  xkstep vixan)
               if (ixc0.lt.0) ixc0 = 2
!              Use extended k range for xanes
               ispec = 1
!              to avoid problems with debye waller factors below the
!              edge, always use complex p for debye waller
!              set the energy grid. xkstep - step in k to use for high
!              energies up to kmax. Near the Fermi level the energy
!              grid is regular in energy with step=vixan
!              the default value is vixan=gamma_ch/2+vi
               if (nwords.gt.1) read(words(2),30,err=900)  xkmax
               if (nwords.gt.2) read(words(3),30,err=900)  xkstep
               if (nwords.gt.3) read(words(4),30,err=900)  vixan

!              sanity checks
               if (xkstep.lt.0.01) xkstep = 0.01d0
               if (xkstep.gt.2.0) xkstep = 0.5d0
               if (abs(xkmax).lt.2) xkmax = 2.d0 !KJ 7-09 added abs for feff8q
               if (abs(xkmax).gt.200) xkmax = 200.d0 !KJ 7-09 added abs for feff8q
               mode = 0
            elseif (itok .eq. 22)  then
!              CORRECTIONS  e0-shift, lambda correction
!              e0 shift is in eV, edge will be edge-e0
!              lambda corr is a const imag energy in eV
!              e0 and lambda corr same as vr0 and vi0 in EXCH card
               read(words(2),30,err=900)  vrcorr
               read(words(3),30,err=900)  vicorr
               mode = 0
            elseif (itok .eq. 23)  then
!              AFOLP use generalized automatic folp
               folpx = 1.15
               if (nwords.ge.2) read(words(2),30,err=900)  folpx
               mode =0
            elseif (itok .eq. 24)  then
!              EXAFS  xkmax for energy grid
               read(words(2),30,err=900)  xkmax
               mode = 0
            elseif (itok .eq. 25)  then
!              POLARIZATION  X Y Z
               ipol = 1
!              run linear polarization code
               read(words(2),30,err=900)  evec(1)
               read(words(3),30,err=900)  evec(2)
               read(words(4),30,err=900)  evec(3)
               mode = 0
            elseif (itok .eq. 26)  then
!              ELLIPTICITY  E incident direction
               read(words(2),30,err=900)  elpty
               read(words(3),30,err=900)  xivec(1)
               read(words(4),30,err=900)  xivec(2)
               read(words(5),30,err=900)  xivec(3)
               mode = 0
            elseif (itok .eq. 27)  then
!              RGRID  rgrd
!              rgrd will be dpas, default is 0.03 in feff7
               read(words(2),30,err=900)  rgrd
               write(slog,'(a,1pe13.5)') ' RGRID, rgrd; ', rgrd
               call wlog(slog)
               i = 1 + int (12.5d0 / rgrd)
               if (mod(i,2) .eq. 0) i = i + 1
               if (i.gt.nrptx) then
                 write(slog,'(a,i6)')                                   &
     &           ' FATAL error in RGRID: increase in m_dimsmod.f90 nrptx to', i
                 call wlog(slog)
                 call par_stop('RDINP')
               endif
               mode = 0
            elseif (itok .eq. 28)  then
!              RPHASES (real phase shifts only)
               call wlog(' Real phase shifts only will be used.  ' //   &
     &                   'FEFF results will be unreliable.')
               lreal = 2
               mode = 0
            elseif (itok .eq. 29)  then
!              NSTAR, write out n* for colinear polarization
               wnstar = .true.
               mode = 0
            elseif (itok .eq. 30)  then
!              NOHOLE
               if (nohole.lt.0) then
                  nohole = 0
                  if (nwords.ge.2) read(words(2),20,err=900)  nohole
               end if
            elseif (itok .eq. 31)  then
!              SIG3 alphat  thetae   first and third cumulants for ss paths
               read(words(2),30,err=900)  alphat
               if (nwords.ge.3) read(words(3),20,err=900)  thetae
               write(slog,'(a,1pe13.5)') ' SIG3, alphat ; ', alphat
               call wlog(slog)
               mode = 0
            elseif (itok .eq. 32)  then
!              JUMPRM remove potential jumps at muffin tin radii
               jumprm = 1
            elseif (itok .eq. 33)  then
!              MBCONV do many body convolution with excitation spectrum
               mbconv = 1
            elseif (itok .eq. 34)  then
!              SPIN  specifies spin direction on central atom
               read(words(2),20,err=900)  ispin
!              set default spin along z axis
               if (ispin.ne.0) spvec(3) = 1.d0
               if (nwords.gt.2) read(words(3),30,err=900)  spvec(1)
               if (nwords.gt.3) read(words(4),30,err=900)  spvec(2)
               if (nwords.gt.4) read(words(5),30,err=900)  spvec(3)
            elseif (itok .eq. 35)  then
!              EDGE     L3
!                   holecode
               call upper(words(2))

               ! For RIXS, read second edge, core or valence
               RixsI%Edges(1) = words(2)
               RixsI%nEdges = 1
               ! Allow for multiple edges for RIXS calculations.
               DO
                  IF(nwords.gt.(RixsI%nEdges+1)) THEN
                     CALL upper(words(RixsI%nEdges+2))
                     IF(TRIM(ADJUSTL(words(RixsI%nEdges+2))).EQ."VAL") THEN
                        RixsI%MBConv = .TRUE.
                        RixsI%nEdges = RixsI%nEdges + 1
                        RixsI%Edges(RixsI%nEdges) = words(RixsI%nEdges+1)
                        EXIT
                     ELSE
                        ! Check that this is an actual edge
                        call setedg (words(RixsI%nEdges+2), ihole)
                     END IF
                  ELSE
                     EXIT
                  END IF
                  RixsI%Edges(RixsI%nEdges+1) = words(RixsI%nEdges+2)
                  RixsI%nEdges = RixsI%nEdges + 1
               END DO
               ! Set the hole to the first edge in case this is not a RIXS calc.
               call setedg (words(2), ihole)
               mode = 0
            elseif (itok .eq. 36)  then
!              SCF    rfms [ lfms nscmt  ca1 nmix  ecv icoul]
!              number of cycles, mode of calculating coulomb potential,
!              convergence accelerator
               nscmt = nbr
               ca1 = 0.2d0

               read(words(2),30,err=900)  rfms1
               if (nwords.gt.2) read(words(3),20,err=900)  lfms1
               if (nwords.gt.3) read(words(4),20,err=900)  nscmt
               if (nwords.gt.4) read(words(5),30,err=900)  ca1
               if (nwords.gt.5) read(words(6),20,err=900)  nmix
               if (nwords.gt.6) read(words(7),30,err=900)  ecv
               if (nwords.gt.7) read(words(8),20,err=900)  icoul
               if (nscmt.le.0 .or. nscmt.gt.nbr) then  !KJ 12-2011 I added the diagnostic message - the user may want to know, y'know?
            call wlog('Invalid number of SCF iterations specified.  Reset to hardwired limit.')
          nscmt = nbr

!            elseif (itok .eq. 107) then
!              SCXC
!               read(words(2),30,err=900) iscfxc !-LC- 11=vBh 12=PZ 21=PDW 22=KSDFT
!write(*,*) "RDINP1", iscfxc
!               !-LC- check if the value of iscfxc is valid
!               if ( (iscfxc .eq. 11) .or. (iscfxc .eq.12) .or. (iscfxc .eq. 21) .or. (iscfxc .eq. 22) ) then
!                  call wlog('Error: iscfxc should take one of the values &
!                11 for vBH, 12 for PZ, 21 for PDW, or 22 for KSDT ... stopping')
!                stop
!               endif

         endif
         if (nwords.gt.5 .and. nmix.gt.30) then !KJ 12-2011 added this so it's done transparently to the user.
            call wlog('Number of Broyden SCF cycles exceeds hardwired maximum of 30; will be reset to 30.')  !KJ This is, I think, vaguely and messily enforced in broydn.f90.  I hate old-school FEFF programming style ...
                  nmix=30
         endif
           if (lfms1.gt.0) lfms1 = 1
!              sanity checks for ca1
               if (ca1.lt.0) ca1 =0
               if (ca1.gt.0.5) then
                 call wlog(' Reduce convergence factors in SCF ')
                 call par_stop(' Cannot run with specified ca1 in SCF card.')
               endif
               if (ecv.ge.0) ecv = -40.0
               if (nmix.le.0) nmix=1
               if (nmix.gt.30) nmix=30
            elseif (itok .eq. 37)  then
!              FMS   rfms2  (lfms2 minv toler1 toler2 rdirec)
!              radius of the cluster to do FMS
               read(words(2),30,err=900)  rfms2
               if (nwords.gt.2) read(words(3),20,err=900)  lfms2
               if (nwords.gt.3) read(words(4),20,err=900)  minv
               if (nwords.gt.4) read(words(5),30,err=900)  toler1
               if (nwords.gt.5) read(words(6),30,err=900)  toler2
               if (nwords.gt.6) read(words(7),30,err=900)  rdirec
               if (rdirec .gt. 2*rfms2 .or. rdirec.lt.0) rdirec=2*rfms2
               if (lfms2.gt.0) lfms2 = 1
         do_fms = 1  !flag that FMS algorithm should run
            elseif (itok .eq. 38)  then
!              LDOS  emin  emax  eimag  neldos
               mldos = 1
               read(words(2),30,err=900)  emin
               read(words(3),30,err=900)  emax
               read(words(4),30,err=900)  eimag
               if (nwords.gt.4) read(words(5),20,err=900)  neldos
               if (neldos.gt.nex) then
                 write (slog, "(a,i4,a)") &
                   "Warning - the number of energy points specified " // &
                   "in the LDOS card is larger than the hardcoded " // &
                   "maximum value (nex = ", nex, "). The maximum value " // &
                   "will be used instead."
                 call wlog(slog)
                 neldos = nex
               end if
            elseif (itok .eq. 39)  then
!              INTERSTITIAL  inters  totvol
!              inters = 1 local V_int (around central atom)
!              inters = 0 extended V_int (average over all atoms)
!              more obscure options described in manual
               read(words(2),20,err=900)  inters
               if (nwords.ge.3) read(words(3),30,err=900)  totvol
            elseif (itok .eq. 40) then
!              CFAV  iphabs nabs rclabs
               read(words(2),20,err=900)  iphabs
               read(words(3),20,err=900)  nabs
               read(words(4),30,err=900)  rclabs
               if (rclabs.lt.0.5) rclabs=big
               mode = 0
            elseif (itok .eq. 41) then
!              S02  s02
               read(words(2),30,err=900)  s02
               mode = 0
            elseif (itok .eq. 42)  then
!              XES ( emin  emax estep)
               if (ixc0.lt.0) ixc0 = 2
!              Use extended k range for xanes
               ispec = 2
!              to avoid problems with debye waller factors below the
!              edge, always use complex p for debye waller
               call wlog('  XES:')
!              keep the same grid variables names as in XANES card
!              with new meaning for ispec=2: xkmax=emin, xkstep=emax
!              and vixan=estep
               xkstep=0.01d0
               if (nwords.gt.1) read(words(2),30,err=900)  xkmax
               if (nwords.gt.2) read(words(3),30,err=900)  xkstep
               if (nwords.gt.3) read(words(4),30,err=900)  vixan
!              sanity checks
               !xkstep = 0.01d0 : JK changed xkstep to max energy.
               if (xkstep.le.xkmax) xkstep=0.01d0
               if (xkmax.ge.0) xkmax = -40.d0
               mode = 0
            elseif (itok .eq. 43)  then
!              DANES ( xkmax  xkstep vixan)
               if (ixc0.lt.0) ixc0 = 2
!              Use extended k range for xanes
               ispec = 3
!              to avoid problems with debye waller factors below the
!              edge, always use complex p for debye waller
               call wlog('  DANES:')
!              set the energy grid. xkstep - step in k to use for high
!              energies up to kmax. Near the Fermi level the energy
!              grid is regular in energy with step=vixan
!              the default value is vixan=gamma_ch/2+vi
               if (nwords.gt.1) read(words(2),30,err=900)  xkmax
               if (nwords.gt.2) read(words(3),30,err=900)  xkstep
               if (nwords.gt.3) read(words(4),30,err=900)  vixan
!              sanity checks
               if (xkstep.lt.0.01) xkstep = 0.01d0
!              if (xkstep.gt.1.0) xkstep = 1.0d0
               if (xkmax.lt.2) xkmax = 2.d0
!              if (xkmax.gt.30) xkmax = 30.d0
               mode = 0
            elseif (itok .eq. 44)  then
!              FPRIME  emin emax estep
               if (ixc0.lt.0) ixc0 = 2
!              Use extended k range for xanes
               ispec = 4
               call wlog(' FPRIME:')
!              set the energy grid.
               read(words(2),30,err=900)  xkmax
               read(words(3),30,err=900)  xkstep
               if (nwords.gt.3) read(words(4),30,err=900)  vixan
!              sanity checks
               if (xkstep.lt.xkmax) xkstep = xkmax
               mode = 0
            elseif (itok .eq. 45)  then
!              RSIGMA  (real self energy only)
               call wlog(' Real self energy only will be used.  FEFF results will be unreliable.')
               if (lreal.lt.1) lreal = 1
               mode = 0
            elseif (itok .eq. 46)  then
!              XNCD or XMCD
               ipol = 2
               mode = 0
            elseif (itok .eq. 47)  then
!              MULTIPOLES le2 (l2lp)
               read(words(2),20,err=900)  le2
               if (nwords.gt.2) read(words(3),20,err=900)  l2lp
               mode = 0
            elseif (itok .eq. 48)  then
!              UNFREEZEF
               iunf = 1
               mode = 0
            elseif (itok .eq. 49)  then
!              TDLDA
               izstd = 1
               if (nwords.gt.1) read(words(2),20,err=900)  ifxc
               mode = 0
            elseif (itok .eq. 50)  then
!              PMBSE
               itdlda = 2
               if (nwords.gt.1) read(words(2),20,err=900)  ipmbse
               if (nwords.gt.2) read(words(3),20,err=900)  nonlocal
               if (nwords.gt.3.and.izstd.eq.0)                          &
     &                          read(words(4),20,err=900)  ifxc
               if (nwords.gt.4) read(words(5),20,err=900)  ibasis
               mode = 0
            elseif (itok .eq. 51)  then ! Added by Josh Kas
!              MPSE [iMP] (alias PLASMON)
!              iMP = 0, use feff default (card does nothing)
!              iMP = 1, use position independent SE, parameterized by the
!                       interstitial density.
!              iMP = 2, use position (density) dependent SE.
               if(nwords.gt.1) then
                  read(words(2),20,err=900) iPlsmn
               else
                  iPlsmn = 1
               end if
               if(nwords.gt.2) read(words(3),20,err=900) NPoles
!              Make old input run like before.
               if(iPlsmn .eq. 4) iPlsmn = 1
            elseif (itok .eq. 52)  then ! Added by Josh Kas
!              SFCONV
               msfconv = 1
            elseif (itok .eq. 53)  then ! Added by Josh Kas
!              SELF (print out on shell self energy Sig(k(E),E) )
               ipse = 1
            elseif (itok .eq. 54)  then ! Added by Josh Kas
!              SFSE k0 (print out self energy Sig(k0,E) )
               ipsk = 1
               read(words(2),30,err=900)  wsigk
            elseif (itok .eq. 55) then ! Added by Josh Kas
!              RCONV (print running convolution with file cfname at energy cen)
!              RCONV cen cname
               read(words(2),30,err=900) cen
               cfname = words(3)(1:12)
            elseif (itok.eq.56) then  !KJ added this card 1-06
!              ELNES
               eels=1   ! switch on ELNES
               absolu=1 !no renormalization in ff2x
!                  now follows the same code as for the XANES card
!              ELNES ( xkmax  xkstep vixan)
               if (ixc0.lt.0) ixc0 = 2
               ispec = 1
!              set the energy grid. xkstep - step in k to use for high
!              energies up to kmax. Near the Fermi level the energy
!              grid is regular in energy with step=vixan
!              the default value is vixan=gamma_ch/2+vi
               if (nwords.gt.1) read(words(2),30,err=900)  xkmax
               if (nwords.gt.2) read(words(3),30,err=900)  xkstep
               if (nwords.gt.3) read(words(4),30,err=900)  vixan
!              sanity checks
               if (xkstep.lt.0.01) xkstep = dble(0.01)
               if (xkstep.gt.2.0) xkstep = dble(0.5)
               if (xkmax.lt.2) xkmax = dble(2)
               if (xkmax.gt.200) xkmax = dble(200)

                 ipol=1   ! override previous entries on POLARIZATION and ELLIPTICITY cards
                 elpty=0
           do i=1,3
                 evec(i)=dble(0)
           enddo
                 mode = 4  ! continue to read the rest of the ELNES card
                 icnt=5  ! number of lines to read
            elseif (itok.eq.57) then  !KJ added this card 1-06
!               EXELFS
               eels=1   ! switch on EXELFS
               absolu=1 !no renormalization in ff2x
!              EXAFS  xkmax for energy grid
               if (nwords.gt.1) read(words(2),30,err=900)  xkmax
                 ipol=1   ! override previous entries on POLARIZATION and ELLIPTICITY cards
                 elpty=0
           do i=1,3
                 evec(i)=dble(0)
           enddo
                 mode = 4  ! continue to read the rest of the EXELFS card
                 icnt=5  ! number of lines to read
            elseif (itok .eq. 58) then !KJ added this card 1-06
!               MAGIC card
                 magic=1
                 read(words(2),30,err=900) emagic
                 icnt=5  ! number of lines to read
            elseif (itok .eq. 59) then !KJ added this card 3-06
!               ABSOLUTE card
                 absolu=1 !KJ end my addition 3-06
            elseif (itok .eq. 60) then !KJ added this card 6-06
!               SYMMETRY card
                 read(words(2),20) ica
           if (ica.lt.1.or.ica.gt.7) ica=-1
           write(slog,'(1x,a,i4,a)') 'SYMMETRY CARD - fixing             &
     &                 icase to ',ica,' in module PATH.'
                 call wlog(slog)
            elseif (itok .eq. 61) then   !KJ 8/06
!               REAL card
                  ispace=1
                  call wlog('Working in real space.')
            elseif (itok .eq. 62) then  !KJ 8/06
!               RECIPROCAL card
                ispace=0
                call wlog('Working in reciprocal space.')
            elseif (itok .eq. 63) then  !KJ 8/06
!               SGROUP card
                  read(words(2),20,err=900) sgroup
            elseif (itok .eq. 64) then  !KJ 8/06
!               LATTICE card
                  mode=5 !read lattice vectors now
            icnt=3 !expecting 3 lines of data
                  read(words(2),10,err=900) latticename
          lattice(1:1)=latticename(1:1)
            if(nwords.gt.2) then
                read(words(3),*,err=900) scalelattice
            else
                scalelattice=dble(1)
            endif
              elseif (itok .eq. 65) then  !KJ 8/06
!               KMESH card
                  ! The number of k-points: could be "KMESH 1000" or "KMESH 10 20 5" or "KMESH 1000 0 0"
                  read(words(2),20,err=900) nkx
          if(nwords.gt.2) then
             read(words(3),20,err=900) nky
           read(words(4),20,err=900) nkz
          endif
          nkp=nkx*nky*nkz
          if(nkp.eq.0) nkp=nkx

          ! Strategy
          ktype=1
          if(nwords.gt.4) read(words(5),20,err=900) ktype
                  ! Apply symmetry to reduce k-mesh?
                  usesym=0
          if(nwords.gt.5) read(words(6),20,err=900) usesym
            elseif (itok .eq. 66) then  !KJ 8/06
!               STRFAC card
                  read(words(2),30,err=900) streta
                  read(words(3),30,err=900) strgmax
                  read(words(4),30,err=900) strrmax
            elseif (itok .eq. 67) then  !KJ 8/06
!               BANDSTRUCTURE card
                  call wlog ('BANDSTRUCTURE card is experimental.')
          mband = 1
          if(nwords.lt.5) then
             call wlog('BANDSTRUCTURE requires at least: emin  emax  estep  ikpath')
             call par_stop(' ')
          endif
          read(words(2),30,err=900) emin_band
          read(words(3),30,err=900) emax_band
          read(words(4),30,err=900) estep_band
          read(words(5),20,err=900) ikpath
          if(nwords.ge.6) read(words(6),20,err=900) nkp_band
          if(nwords.ge.7) read(words(7),*,err=900) freeprop
            elseif (itok .eq. 68) then !JK 8/09
!               COREHOLE hole_treatment    -   default is FSR
                if(nwords.gt.1) then
                   call upper(words(2))
                   if(TRIM(ADJUSTL(words(2))).eq.'NONE') then
                      nohole = 0
                   elseif(TRIM(ADJUSTL(words(2))).eq.'RPA') then
                      nohole = 2
                   elseif((TRIM(ADJUSTL(words(2))).eq.'FSR') .or. (TRIM(ADJUSTL(words(2))).eq.'REGULAR')) then
                      !I'm keeping 'regular' here for compatibility - don't tell John :)
                      nohole = -1
                   else
                      call wlog('Invalid COREHOLE option - choose NONE, RPA, or FSR.')
                      stop
                endif
                ! Set COREHOLE to COREHOLE none if compton calculation is done
                if (do_compton) then
                  nohole = 0
                end if
            end if
            elseif (itok .eq. 71) then  !KJ 8.06
!               TARGET card
                read(words(2),20,err=900) absorber
            elseif (itok .eq. 72) then  !KJ 1.07
!               EGRID card
                if(nwords.gt.1) then
                  ! XXX this case probably doesn't work since phmesh is no longer called...
                    read(words(2),20,err=900) iegrid
                      if(iegrid.eq.2) then
                      read(words(3),10,err=900) egridfile
                         call wlog('Energy grid to be read from file.')
                      elseif(iegrid.eq.3) then
                         read(words(3),20,err=900) egrid3a
                         read(words(4),30,err=900) egrid3b
                         read(words(5),30,err=900) egrid3c
                         call wlog('Energy grid on exponential mesh.')
                      else
                         iegrid=0
                         call wlog('Regular FEFF energy grid.')
                      endif
                else
                   iGrid = 1
                   mode=6
                   open(UNIT=15,FILE='grid.inp',STATUS='UNKNOWN')
                end if
            elseif (itok .eq. 73) then
!              COORDINATES icoord
               read(words(2),20,err=900) icoord
               if(icoord.eq.1) then
                  call wlog('Atom positions are given in Carthesian coordinates.')
            call wlog('The units are Angstrom.')
            call wlog('FEFF-like coordinates.')
               elseif(icoord.eq.2) then
                  call wlog('Atom positions are given in Carthesian coordinates.')
            call wlog('The units are fractions of the resp. lattice vector.')
               elseif(icoord.eq.3) then
                  call wlog('Atom positions are given in Carthesian coordinates.')
            call wlog('The units are fractions of the first lattice vector.')
            call wlog('SPRKKR-like coordinates.')
               elseif(icoord.eq.4) then
                  call wlog('Atom positions are given in lattice coordinates.')
            call wlog('The units are fractions of the resp. lattice vector.')
            call wlog('WIEN2k-like input (beware of funny lattice types).')
               elseif(icoord.eq.5) then
                  call wlog('Atom positions are given in lattice coordinates.')
            call wlog('The units are fractions of the first lattice vector.')
               elseif(icoord.eq.6) then
                  call wlog('Atom positions are given in lattice coordinates.')
            call wlog('The units are Angstrom.')
               else
                  call wlog('Attempt to enter funky lattice coordinates.')
            call wlog('Please stick to one of the formats described in the manual.')
            call wlog('Exiting now.')
            stop
               endif

            elseif (itok .eq. 74) then
!              EXTPOT
               ExternalPot = .TRUE.
            elseif (itok .eq. 75) then
               ! CHBROAD
               read(words(2),20,err=900) iGammaCH
            elseif (itok .eq. 76) then
!              Added by Fer
               read(words(2),*,err=900) ChSh_Type
            elseif (itok .eq. 77) then
               ! DIMS card
               ! First is nclusx, second is lx
               read(words(2),20,err=900) nclusxuserlimit
               read(words(3),20,err=900) lxuserlimit
            elseif (itok .eq. 78)  then
!              NRIXS card    !KJ merged 7/09 from feff8q !KJ 12/2010 merged APS code that treats several q's and MDFF
               read(words(2),20,err=900)  nq  !KJ 12/09 changed format 30 to 20
         if (nq.lt.0) then
            qaverage=.true.
          nq=abs(nq)
         else
            qaverage=.false.
          if(nq.eq.0) nq=1  ! nq=0 was allowed in feff8q and feff8qwithnq
         endif
         call make_qlist(nq)
         !read the first vector from the current line:
         if (qaverage) then
!                    just read one component, assume spherical averaging
                     read(words(3),30,err=900)  qvec(3)
             if(nq.gt.1) then
              read(words(4),30,err=900) dummy
            dummy2=0.d0; if(nwords.ge.5) read(words(5),30,err=900) dummy2
            qw(1)=dcmplx(dummy,dummy2) !weight
           endif
                     qvec(2)=0.0d0
                     qvec(1)=0.0d0
             qn(1)=qvec(3) !norm
                     if (qvec(3).le.0.0d0) then
                        call wlog(' ERROR: momentum transfer negative or zero')
                        call par_stop(' ')
                     end if
         else
                     read(words(3),30,err=900)  qvec(1)
                     read(words(4),30,err=900)  qvec(2)
                     read(words(5),30,err=900)  qvec(3)
             if(nq.gt.1) then
                read(words(6),30,err=900) dummy
                dummy2=0.d0; if(nwords.ge.7) read(words(7),30,err=900) dummy2
                qw(1)=dcmplx(dummy,dummy2)
           endif
             qn(1)=dsqrt(qvec(1)**2+qvec(2)**2+qvec(3)**2)
         end if
         qs(1,:)=qvec
         !read all other vectors, one per line:
         if(nq.gt.1) then
         do i=2,nq
            call rdline(jinit,line)
          if(line.eq.'read_line_end')   line='END'
          if(line.eq.'read_line_error') line='END'
          nwords=nwordx
          call bwords(line,nwords,words)
                  if (qaverage) then
             if(nwords.lt.2) stop 'expecting "q qweight" in feff.inp'
!                    just read one component, assume spherical averaging
                     read(words(1),30,err=900) qvec(3)
           read(words(2),30,err=900) dummy
           dummy2=0.d0; if(nwords.ge.3) read(words(3),30,err=900) dummy2
           qw(i)=dcmplx(dummy,dummy2)
                     qvec(2)=0.0d0
                     qvec(1)=0.0d0
             qn(i)=qvec(3) !norm
                     if (qvec(3).le.0.0d0) then
                        call wlog(' ERROR: momentum transfer negative or zero')
                        call par_stop(' ')
                     end if
                  else
             if(nwords.lt.4) stop 'expecting "qx qy qz qweight" in feff.inp'
                     read(words(1),30,err=900)  qvec(1)
                     read(words(2),30,err=900)  qvec(2)
                     read(words(3),30,err=900)  qvec(3)
           read(words(4),30,err=900) dummy
           dummy2=0.d0; if(nwords.ge.5) read(words(5),30,err=900) dummy2
           qw(i)=dcmplx(dummy,dummy2)
             qn(i)=dsqrt(qvec(1)**2+qvec(2)**2+qvec(3)**2)
                  end if
            qs(i,:)=qvec
         enddo
         qvec=qs(1,:) ! This is a precaution, not sure if it's desirable.  fix later
         endif ! nq>1
         do_nrixs=1
               mode = 0
      elseif (itok .eq. 79) then  !KJ 7-09 merged from feff8q
!              LJMAX lj   abs(LJMAX) gives the number of terms in the expansion
!               of e^{iqr} in terms of spherical Bessel functions.
!               Traditionally it was negative. Do not know if that is needed in this version.
               read(words(2),20,err=900)  lj
      elseif (itok .eq. 80) then  !KJ 7-09 merged from feff8q
!              LDEC ldecmx  : Calculate contributions from different l-final states.
               read(words(2),20,err=900)  ldecmx
            elseif (itok .eq. 81) then
            !SETEDGE - set excitation energies based on elam/mcmasters table
               lopt=.true.
            elseif (itok .eq. 82) then
!           EPS0 - set dielectric constant for MPSE calculation.
               read(words(2),30,err=900) Eps0
            elseif (itok .eq. 83) then
!           OPCONS - create loss.dat file for MPSE from internal database.
               run_opcons = .TRUE.
            elseif (itok .eq. 84) then
!           NUMDENS - set the number densities for creating loss.dat from database.
               read(words(2),20,err=900) iph
               IF(iph.gt.nphxhardlimit) then
                  call wlog("iph > nphxhardlimit in feff.inp")
                  call wlog(TRIM(ADJUSTL(line)))
                  call par_stop
               END IF
               read(words(3),30) NumDens(iph)
            elseif (itok .eq. 85) then
!           PREPS - Print out epsilon from database.
               print_eps = .TRUE.
            elseif (itok .eq. 86) then
!           EGAP - Set gap energy for self-energy calculation.
               read(words(2),*) EGap
            elseif (itok .eq. 87) then
!           CHWIDTH - Set corehole lifetime manually instead of using the tables in COMMON/setgam.f90
               read(words(2),*) userChLifetime
               userchl=.true.
      elseif (itok .eq. 88) then
!           MDFF - Calculate the mixed dynamic form factor  (was called "ADAM" in first version, lol)
         imdff=1
         if(nwords.ge.2) read(words(2),20,err=900) imdff
         if(imdff.eq.2) then
            if(nwords.eq.2) then  !use q and q' from the NRIXS list
             qqmdff=-1.d0
             cosmdff_dum=0.d0
          elseif(nwords.eq.4) then  ! use only q from the NRIXS list; generate q' using:
               read(words(3),30,err=900) qqmdff
               read(words(4),30,err=900) cosmdff_dum
          else  !invalid syntax
             stop "fatal error in feff.inp - expecting:   MDFF 2  q'  angle              or     MDFF 2"
          endif
         endif
         if(imdff.le.0) then
            call wlog('MDFF calculation disabled.')
         elseif(imdff.eq.3) then
            call wlog("EELS type MDFF calculation selected - summed over all q,q' pairs")
         elseif(imdff.eq.2) then
            call wlog("NRIXS type MDFF calculation selected - for a single q,q' pair only.")
               elseif(imdff.eq.1) then
            call wlog("NRIXS type MDFF calculation selected - summed over all q,q' pairs.")
         else
            call wlog('Invalid MDFF option selected.')
          call par_stop('RDINP-2')
         endif
      elseif (itok .eq. 89) then
!           RESTART - get the initial potentials for SCF from a pot.bin file
         StartFromFile = .true.
            elseif (itok .eq. 90) then
!           CONFIG - use non-standard electron configuration for some atoms

         if(words(2) .eq. 'file') then
            configtype=2
         elseif(words(2) .eq. 'feff7') then
            configtype=7
               elseif(words(2) .eq. 'card') then
            configtype=2
          !simply dump whatever's in the card to a file 'config.inp'.
          !if there are mistakes, the user will find out later.
          read(words(3),*) j
          open(62,file='config.inp',form='formatted',status='unknown') !,access='append')
          do i=1,j
             call rdline(jinit,line)
                     write(62,'(a)') line
          enddo
          close(62)
         else
            call wlog('Dubious use of the CONFIG card; calculation will proceed with defaults.')
         endif
      elseif (itok .eq. 91) then
!           SCREEN - pass on some options to the facultative screen.inp file
               if (nwords.lt.3) then
            stop 'SCREEN card must be followed by precisely two arguments, e.g. "SCREEN rfms 5.5"'
         else
            str3=words(2)  !takes first 3 letters
          read(words(3),*) dummy
          call screen_inp_parse(str3,dummy) !KJ 1-2012 used to be "call screen_inp_parse_and_write(str3,dummy)"
          call wlog(":INFO  User provides options for screen.inp")
         endif
      elseif (itok .eq. 92) then
!      CIF - read crystal structure from .cif file
         if (nwords.lt.2) stop 'Error - CIF card must be followed by filename e.g. file.cif'
         read(words(2),'(a)') cifname
         cifread=.true.
      elseif (itok .eq. 93) then
!           EQUIVALENCE - governs choosing of potential types from crystallographic information
         if (nwords.lt.2) call wlog('No equivalence type specified in EQUIVALENCE card - using default.')
               read(words(2),*) cif_equivalence
      elseif (itok .eq. 94) then
!           COMPTON - calculates Compton profile
               ispec = 5
               do_compton = .true.
               ! Set COREHOLE to none for compton calculation
               nohole = 0
               save_gg_slice = .true.
         ltmp=0  !KJ 10-2012 bugfix for Win "ltmp used without being defined"
               if (nwords.gt.1) read(words(2),30,err=900) pqmax
               if (nwords.gt.2) read(words(3),20,err=900) npq
               if (nwords.gt.3) read(words(4),20,err=900) ltmp
               if (ltmp.gt.0) force_jzzp = .true.
      elseif (itok .eq. 95) then
!           RHOZZP - calculate rho(z,z') along a slice
               do_rhozzp = .true.
      elseif (itok .eq. 96) then
!           CGRID - grid parameters for COMPTON,RHOZZP
              if (nwords.gt.1) read(words(2),30,err=900) zpmax
              if (nwords.gt.2) read(words(3),20,err=900) ns
              if (nwords.gt.3) read(words(4),20,err=900) nphi
              if (nwords.gt.4) read(words(5),20,err=900) nz
              if (nwords.gt.5) read(words(6),20,err=900) nzp
      elseif (itok .eq. 97) then
!      CORVAL - set minimum energy for core-valence separation energy search
!             This card is a temporary fix until we fix the core-valence problem properly (FEFF9.7?)
        if (nwords.gt.1) then
           read(words(2),30,err=900) corval_emin ! eV
        else
           call wlog('Ignoring CORVAL card without parameter corval_emin')
        endif
      elseif (itok .eq. 98) then
!           SIGGK - multiply fine structure chi(k) by energy-dependent but otherwise global Debye-Waller factor
!                 i.e. chi(k) * exp{-(sig_gk k)^2}
       if (nwords.gt.1) then
          read(words(2),30,err=900) sig_gk  ! in Angstrom (not Angstrom^2)
        else
           call wlog('Ignoring SIGGK card without parameter sig_gk')
        endif

      elseif (itok .eq. 99) then
!           TEMPERATURE - valence electron temperature for SCF loop
              if (nwords.gt.1) read(words(2),30,err=900) scf_temperature
              ! if (nwords.gt.2) read(words(3),20,err=900) scf_thermal_vxc
              if (nwords.gt.2) read(words(3),20,err=900) iscfxc
              electronic_temperature = scf_temperature
              ! JK - Old version of FEFF used vbh when temperature=0, which
              ! corresponds to iscfxc=11. Let's keep that default
              ! if (scf_temperature.GT.0.d0) THEN
              !    if (iscfxc.ne.22) iscfxc=21 !-LC- Default PDW when using electronic temperature
              ! else
              !    if (iscfxc.ne.12) iscfxc=11 ! JK - Default vbh if temperature = 0
              ! end if
	            ! iscfxc=22
      elseif (itok .eq. 100) then
!           DENSITY - valence electron density
              ispec = 5
              save_gg_slice = .true.
              mode = 7
              open(16, file='density.inp', status='unknown')
      elseif (itok .eq. 101) then
!           RIXS - set options for RIXS calculations.
               RixsI%m_run = 1
               ! experimental broadening in excitation energy
               IF(nwords.GT.1) read(words(2), *, err=900) RixsI%gam_exp(1)
               ! experimental broadening in energy loss
               IF(nwords.GT.2) read(words(3), *, err=900) RixsI%gam_exp(2)
               ! set Fermi level
               IF(nwords.GT.3) read(words(4), *, err=900) RixsI%xmu
      elseif (itok .eq. 102) then
               ! RLPRINT - Make XSPH print wavefunctions.
                PrintRl = .TRUE.
      elseif(itok .eq. 103) then
!           ICORE - set core state to use for matrix elements
               read(words(2), *, err=900) iCoreState
      elseif(itok .eq. 104) then
!            HUBBARD
               i_hubbard = 2
               mldos_hubb = 2
         !     l_hubbard=3
               read(words(2),30,err=900)  U_hubbard
               read(words(3),30,err=900)  J_hubbard
               read(words(4),30,err=900)  fermi_shift
               read(words(5),20,err=900)  l_hubbard
      elseif (itok .eq. 105) then
!            CRPA
               CRPAI%do_CRPA=1
               read(words(2),20,err=900) CRPAI%l_crpa
               read(words(3),30,err=900)  CRPAI%rcut
      elseif (itok .eq. 106) then
            ! FULLSPECTRUM
               mFullSpectrum = 1

      elseif (itok .eq. 107) then
!              SCXC
               read(words(2),20,err=900) iscfxc
!               read(words(2),30,err=900) iscfxc !-LC- 11=vBh 12=PZ 21=PDW 22=KSDFT
!               !-LC- check if the value of iscfxc is valid
               if ( (iscfxc .ne. 11) .and. (iscfxc .ne.12) .and. (iscfxc .ne. 21) .and. (iscfxc .ne. 22) ) then
                  call wlog('Error: iscfxc should take one of the values &
                11 for vBH, 12 for PZ, 21 for PDW, or 22 for KSDT ... stopping')
                stop
               endif
      elseif (itok .eq. 108) then ! JK
              ! HIGHZ: Use finite nucleus to calculate atomic wavefunctions
              ! No real change for low z elements, but high z needs this for
              ! atomic energies and wavefunctions.
              FiniteNucleus = .TRUE.
              CALL wlog('Using finite nucleus.')
      elseif (itok .eq. 109) then
              emaxscf = 5 !eV
              negrid = 400
              xntol = 1e-4
              nmu = 100
              read(words(2),20,err=900) iscfth
              if (nwords.gt.2) read(words(3),30,err=900) emaxscf
              if (nwords.gt.3) read(words(4),20,err=900) negrid
              if (nwords.gt.4) read(words(5),20,err=900) nmu
              if (nwords.gt.5) read(words(6),30,err=900) xntol
      elseif (itok .eq. -1)  then
!              END
               goto 220
      else
               write(slog,'(1x,a)') line(1:70)
               call wlog(slog)
               write(slog,'(1x,a)') words(1)
               call wlog(slog)
               write(slog,'(a,i8)') ' Token ', itok
               call wlog(slog)
               call wlog(' Keyword unrecognized.')
               call wlog(' See FEFF document -- some old features are no longer available.')
               call par_stop('RDINP-2')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
!              We're done reading atoms.
!              Change mode and process current card.
               mode = 0
               goto 210
            endif
            natt = natt+1
            if (natt.gt. nattx)  then
               write(slog,'(a,i8)') 'Too many atoms, maximum is ', nattx
               call wlog(slog)
               call par_stop('RDINP-3')
            endif
            read(words(1),*,err=900)  ratx(1,natt)
            read(words(2),*,err=900)  ratx(2,natt)
            read(words(3),*,err=900)  ratx(3,natt)
            read(words(4),*,err=900)  iphatx(natt)
            if (iatph(iphatx(natt)) .le. 0) iatph(iphatx(natt)) = natt
      !KJ I wonder if the above choice for iatph can be problematic if the atoms list is very long,
      ! and the first entries don't make it into the final cluster (FMS),
      ! or end up being very far from the absorber ...
      !Answer: ffsort redetermines iatph appropriately before writing it to geom.dat

         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
!              We're done reading these overlap instructions.
!              Change mode and process current card.
               mode = 0
               goto 210
            endif
            novr(iph) = novr(iph)+1
            iovr = novr(iph)
            if (iovr .gt. novrx)  then
               write(slog,'(a,i8)') 'Too many overlap shells, maxiomum is ',  novrx
               call wlog(slog)
               call par_stop('RDINP-5')
            endif
            read(words(1),20,err=900) iphovr(iovr,iph)
            read(words(2),20,err=900) nnovr(iovr,iph)
            read(words(3),30,err=900) rovr(iovr,iph)
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
!              We're done reading unique potential definitions
!              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphxhardlimit)  then
               write(slog,'(a,i8)')  'Unique potentials must be between 0 and ',nphxhardlimit
               call wlog(slog)
               write(slog,'(i8,a)') iph, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('RDINP')
            endif
            read(words(2),20,err=900)  iz(iph)
            if (iz(iph).lt. 6) then
               lmaxsc(iph) = 1
            elseif (iz(iph).lt.55) then
               lmaxsc(iph) = 2
            else
               lmaxsc(iph) = 3
            endif
!           No potential label if user didn't give us one
!           Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
            if (nwords .ge. 4)  then
              read(words(4),20,err=900) ltmp
!KJ lx now dynamic             if (ltmp.ge.1 .and. ltmp.le.lx) lmaxsc(iph) = ltmp
              !if (ltmp.ge.1) lmaxsc(iph) = ltmp
              if (ltmp.ge.0) lmaxsc(iph) = ltmp
              if (ltmp.eq.0) call wlog('WARNING - You requested lmax=0 in ' // &
                'POTENTIALS.  Please make sure that this is intentional, ' // &
                'as it may be unstable or give very poor results otherwise.')
            endif
            lmaxph(iph) = 3
            if (iz(iph).lt.6) lmaxph(iph) = 2
            if (nwords .ge. 5)  then
              read(words(5),20,err=900) ltmp
!KJ lx now dynamic              if (ltmp.ge.1 .and. ltmp.le.lx) lmaxph(iph) = ltmp
              !if (ltmp.ge.1) lmaxph(iph) = ltmp
              if (ltmp.ge.0) lmaxph(iph) = ltmp
              if (ltmp.eq.0) call wlog('WARNING - You requested lmax=0 in ' // &
                'POTENTIALS.  Please make sure that this is intentional, ' // &
                'as it may be unstable or give very poor results otherwise.')
            endif
            if (nwords .ge. 6) then
              read(words(6),30,err=900) xnatph(iph)
              lxnat = 1
            endif
            if (nwords .ge. 7) then
              read(words(7),30,err=900) spinph(iph)
            endif
      nph_read=iph
           elseif (mode.eq.4) then  !KJ 1-06 this mode added to read ELNES card
             if(icnt.eq.5) then
                 call fixlinenow(words,nwords)
                 read(words(1),*,err=900) ebeam   ! read beam energy in keV
                 ebeam=ebeam * dble(1000)  ! convert to eV
                 if (nwords.ge.2) read(words(2),20,err=900) aver ! average over sample to beam orientation?
               if (aver.eq.1) icnt=icnt-1 !skip the line for beam orientation
                 if (nwords.ge.3) read(words(3),20,err=900) cross ! calculate cross terms?
                 if (nwords.ge.4) read(words(4),20,err=900) relat ! use relativistic q-vector?
               if (nwords.ge.5) read(words(5),20,err=900) iinput ! read xmu.dat or opconsKK.dat or ... ?   !KJ 5/6
                 if (nwords.ge.6) read(words(6),20,err=900) spcol !column that has spectrum
              elseif(icnt.eq.4) then
                 read(words(1),*,err=900) xivec(1)  ! read direction of incoming beam
                 read(words(2),*,err=900) xivec(2)  ! in arbitrary units
                 read(words(3),*,err=900) xivec(3)
                 xinorm=dsqrt(xivec(1)**2+xivec(2)**2+xivec(3)**2)
               if (xinorm.gt.0.0) then
                  do i=1,3
                       xivec(i)=xivec(i)/xinorm    ! normalize this vector.
                  enddo
               elseif(.not.(aver.eq.1)) then
                  call wlog('WARNING : beam direction unspecified in orientation sensitive EELS calculation.      &
     &                  Please correct before running EELS module.')
               endif
             elseif(icnt.eq.3) then
                 read(words(1),*,err=900) acoll  ! collection semiangle in mrad
                 read(words(2),*,err=900) aconv  ! convergence semiangle in mrad
                 acoll=acoll/dble(1000);aconv=aconv/dble(1000) ! convert from mrad to rad
             elseif(icnt.eq.2) then
                 read(words(1),*,err=900) nqr    ! specify q-mesh, radial parameter
                 read(words(2),*,err=900) nqf    ! specify q-mesh, angular parameter
               if(nqr*nqf.eq.0) then
                  call wlog('WARNING : zero q-mesh points specified for EELS calculation.  Please correct before       &
     &               running EELS module.')
                 endif
             elseif(icnt.eq.1) then
                 read(words(1),*,err=900) thetax ! detector position in plane perpendicular to beam ; angle in mrad
                 read(words(2),*,err=900) thetay ! detector position in plane perpendicular to beam ; angle in mrad
                 thetax=thetax/dble(1000);thetay=thetay/dble(1000) !mrad to rad
                 mode=0  ! finished reading ELNES card
             endif
             icnt=icnt-1    ! now read the next line
         !KJ end my changes
         elseif (mode.eq.5) then  !KJ 04/2007 read lattice vectors
              if(icnt.eq.3) then
                  read(words(1),*,err=900) a1(1)
                  read(words(2),*,err=900) a1(2)
                  read(words(3),*,err=900) a1(3)
            a1=a1*scalelattice
              elseif(icnt.eq.2) then
                  read(words(1),*,err=900) a2(1)
                  read(words(2),*,err=900) a2(2)
                  read(words(3),*,err=900) a2(3)
            a2=a2*scalelattice
              elseif(icnt.eq.1) then
                  read(words(1),*,err=900) a3(1)
                  read(words(2),*,err=900) a3(2)
                  read(words(3),*,err=900) a3(3)
            a3=a3*scalelattice
              endif
              icnt=icnt-1
              if(icnt.eq.0) mode=0  !finished reading LATTICE card
         elseif (mode.eq.6) then ! Josh Kas 11/2009 egrid input
              ! Just write this to grid.inp
              if(itok.eq.0) then
                 write(15,*) (TRIM(ADJUSTL(words(i))) // ' ', i=1, nwords)
              else
                 ! Done reading egrid input.
                 mode=0
                 close(15)
                 goto 210
              end if
         elseif (mode.eq.7) then ! density input
              if (itok.ne.0) then
                 mode=0
                 close(16)
                 goto 210
               end if
               write(16,'(a)') line(1:istrln(line))
         else
            write(slog,'(a,i8)') 'Mode unrecognized, mode ', mode
            call wlog(slog)
            call par_stop('RDINP-6')
         endif
      goto 200
  220 continue
! DONE READING INPUT FILE,
!#{mn
! call rdline with jinit=0 to clean up all input files
       jinit = 0
       call rdline(jinit,line)
!#mn}

! ##########################################################################################


!     Fix up defaults, error check limits, figure out free atoms, etc.


!      if (ispace.eq.1 .and. cifread)  call wlog('CIF input option ignored for real space calculation.')

    if (cifread .and. cards_set(1))  call wlog('CIF and ATOMS cards used: ATOMS card will be ignored.')


! !KJ 8/06 copy data from ATOMS card to dedicated arrays :
      if (ispace.eq.0 .and. (.not.cifread)) then

     call wlog('Taking crystal structure from feff.inp.  Note: .cif input is now recommended.')
     ! Convert the allowed coordinate system input to the one used internally:
     nats=natt !KJ bugfix 2-2012
         if(icoord.eq.1) then
            alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
            ratx=ratx/alatt
         elseif(icoord.eq.2) then
            alatt=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
            ratx(2,:)=ratx(2,:)*alatt
            alatt=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
            ratx(3,:)=ratx(3,:)*alatt
            alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
            ratx(2:3,:)=ratx(2:3,:)/alatt
         elseif(icoord.eq.3) then
! no action required ; this is default
         elseif(icoord.eq.4) then
            alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
            do iatrd=1,nats
               xxx(1)=a1(1)*ratx(1,iatrd)+a2(1)*ratx(2,iatrd)+a3(1)*ratx(3,iatrd)
               xxx(2)=a1(2)*ratx(1,iatrd)+a2(2)*ratx(2,iatrd)+a3(2)*ratx(3,iatrd)
               xxx(3)=a1(3)*ratx(1,iatrd)+a2(3)*ratx(2,iatrd)+a3(3)*ratx(3,iatrd)
               ratx(:,iatrd)=xxx/alatt
            enddo
         elseif(icoord.eq.5) then
            alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
            ratx(2:3,:)=ratx(2:3,:)*alatt
            alatt=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
            ratx(2,:)=ratx(2,:)/alatt
            alatt=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
            ratx(3,:)=ratx(3,:)/alatt
            do iatrd=1,nats
               xxx(1)=a1(1)*ratx(1,iatrd)+a2(1)*ratx(2,iatrd)+a3(1)*ratx(3,iatrd)
               xxx(2)=a1(2)*ratx(1,iatrd)+a2(2)*ratx(2,iatrd)+a3(2)*ratx(3,iatrd)
               xxx(3)=a1(3)*ratx(1,iatrd)+a2(3)*ratx(2,iatrd)+a3(3)*ratx(3,iatrd)
               ratx(:,iatrd)=xxx/alatt
            enddo
         elseif(icoord.eq.6) then
            alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
            ratx(1,:)=ratx(1,:)/alatt
            alatt=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
            ratx(2,:)=ratx(2,:)/alatt
            alatt=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
            ratx(3,:)=ratx(3,:)/alatt
            do iatrd=1,nats
               xxx(1)=a1(1)*ratx(1,iatrd)+a2(1)*ratx(2,iatrd)+a3(1)*ratx(3,iatrd)
               xxx(2)=a1(2)*ratx(1,iatrd)+a2(2)*ratx(2,iatrd)+a3(2)*ratx(3,iatrd)
               xxx(3)=a1(3)*ratx(1,iatrd)+a2(3)*ratx(2,iatrd)+a3(3)*ratx(3,iatrd)
               ratx(:,iatrd)=xxx/alatt
            enddo
         endif

           call init_struct(natt)
           nats=natt
           ppos(1:3,1:nats)=ratx(1:3,1:nats)
           ppot(1:nats)=iphatx(1:nats)
!!! Disable the next section
!!    It is expected that coordinates are given in FRACTIONAL COORDINATES!!
!               do i=1,3
!                 do iatrd=1,nats
!!     Reduce atom positions to first unit cell [0,1]^3
!                   if(dabs(ppos(i,iatrd)).gt.dble(1))                   &
!     &            ppos(i,iatrd)=ppos(i,iatrd)-int(ppos(i,iatrd))
!                   if(ppos(i,iatrd).lt.dble(0))                         &
!     &              ppos(i,iatrd)=ppos(i,iatrd)+dble(1) !KJ fix later
!                 enddo
!               enddo

    elseif ((ispace.eq.0 .and. cifread) .or. cifread) then
       call wlog('Taking crystal structure from .cif file.')
       call importcif(cifname,cif_equivalence)
           !importcif creates ppos,ppot,a1,a2,a3 and more.  init_struct is called inside.

       !Now deal with the settings of the POTENTIALS card:
           if (cards_set(14)) then  !POTENTIALS card in feff.inp
          ! First check that POTENTIALS and cif file are compatible:
!              call wlog(':WARNING  You are using CIF import and POTENTIALS card.  Make sure the two correspond perfectly!')
        if(nph_read.ne.nphstr) then
        call wlog(':WARNING  POTENTIALS card contains different number of potentials than cif file.  Ignoring POTENTIALS card.')
        cards_set(14)=.false.
        endif
        nph=nphstr
        do iph=0,nph
           if(iz(iph).ne.izatom(iph)) then
          call wlog(':WARNING   POTENTIALS card contains different atomic number than cif file.  Ignoring POTENTIALS card.')
          cards_set(14)=.false.
         endif
        enddo
        if (cards_set(14)) then
           !All compatibility checks passed!
         !Now the data from the POTENTIALS card can just be kept as is :
         ! iz,potlbl,lmaxsc,lmaxph,xnatph,spinph
        endif
       endif
       if (.not.cards_set(14)) then  !no POTENTIALS card in feff.inp
          nph=nphstr
        iz(:)=-1 !careful: this must not be initialized to "0"!
        iz(0:nph)=izatom(0:nph)
        do iph=0,nph  ! carefully convert 2-string to 6-string
           potlbl(iph)='      '
         str6(1:2)=label(iph)
         str6(3:6)='    '
         potlbl(iph)=str6
        enddo
        do iph=0,nph
           ! lmaxsc and lmaxph will be set to defaults:
               if (iz(iph).lt. 6) then
                   lmaxsc(iph) = 1
           lmaxph(iph) = 2
                 elseif (iz(iph).lt.55) then
                   lmaxsc(iph) = 2
           lmaxph(iph) = 3
                 else
                   lmaxsc(iph) = 3
           lmaxph(iph) = 3
                 endif
        enddo
        ! set stoichiometry in xnatph:
        lxnat=1
        xnatph(0)=0.01d0
        xnatph(1:nph)=natom(1:nph)
        ! spinph CANNOT be set in the current implementation - the user must have a POTENTIALS card for that to work.
        spinph(:)=0.d0
       endif

     endif

     if ((ispace.eq.0) .or. (cifread .and. ispace.eq.1)) then

!  The FMS routines need only the information above.  However, the path expansion still needs a real space list of coordinates.
!  This list must now be generated, and placed in the arrays corresponding to the ATOMS card.
!  After that, initialization can continue.
         if ((absorber.lt.1.or.absorber.gt.nats)) then
               call wlog ('No absorber specified - assigning it to first position.')
               absorber=1
         endif

!  Now we replace atoms-list from feff.inp by one generated from the atoms in pos by periodic repetition.

         !List atoms up to at least rmax + 33%
         magnifier=1.33d0
         ratomslist=max(20.d0,magnifier*rmax)
         i1=int(ratomslist/dsqrt(a1(1)**2+a1(2)**2+a1(3)**2))+1
         i2=int(ratomslist/dsqrt(a2(1)**2+a2(2)**2+a2(3)**2))+1
         i3=int(ratomslist/dsqrt(a3(1)**2+a3(2)**2+a3(3)**2))+1
     if(lattice.eq.'P'.or.lattice.eq.'H') then
        lattice_factor=1
     elseif(lattice.eq.'F') then
        lattice_factor=4
     else
        lattice_factor=2
     endif
4245     continue  !Come back here to try again if the first list was too long
         j=nats*(2*i1+1)*(2*i2+1)*(2*i3+1)*lattice_factor
         alatt=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
!     write(*,*) 'i1,i2,i3',i1,i2,i3
!     write(*,*) 'ratomslist,lattice_factor',ratomslist,lattice_factor
!     write(*,*) 'j,nattx',j,nattx
!     write(*,*) 'a1,a2,a3',dsqrt(a1(1)**2+a1(2)**2+a1(3)**2),dsqrt(a2(1)**2+a2(2)**2+a2(3)**2),dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
         if(j.gt.nattx) then
          if(i1+i2+i3.eq.0) then
                 call par_stop ('WARNING - current value of nattx does not allow to calculate up to rmax as specified in feff.inp')
        else
           i1=i1-1
         i2=i2-1
         i3=i3-1
         goto 4245
        endif
         endif
          nshift=0
          shift(:,:)=0.d0
              ! Careful here -- bugfix KJ May 2014
              ! lattice contains the 1st letter of the lattice name, it is the "lattice type"
              ! e.g. F, R, C
              ! While latticename contains the full thing
              ! e.g. F, R, CXY
        if(lattice.eq.'F') then
           nshift=3
         shift(1,1)=0.5d0 ; shift(2,1)=0.5d0
         shift(1,2)=0.5d0 ; shift(3,2)=0.5d0
         shift(3,3)=0.5d0 ; shift(2,3)=0.5d0
        elseif(lattice.eq.'R') then
           nshift=1
         !stop 'eek'
                 nshift=0
        elseif(latticename.eq.'CXY') then
           nshift=1
         shift(1,1)=0.5d0 ; shift(2,1)=0.5d0
        elseif(latticename.eq.'CXZ') then
           nshift=1
         shift(1,1)=0.5d0 ; shift(3,1)=0.5d0
        elseif(latticename.eq.'CYZ') then
           nshift=1
         shift(3,1)=0.5d0 ; shift(2,1)=0.5d0
        elseif(lattice.eq.'I'.or.lattice.eq.'B') then
                 nshift=1
                 shift(1,1)=0.5d0 ; shift(2,1)=0.5d0 ; shift(3,1)=0.5d0
        endif  !If P or H, no need to make extra atoms

         natt=0
         if(debug_cif) then
            write(*,*) 'Now making real-space atoms list from cif input.'
            write(*,*) 'nats= ',nats,' lattice=',lattice,' nshift= ',nshift
            write(*,*) 'using ',2*i1+1,'x',2*i2+1,'x',2*i3+1, &
              ' unit cells or ',nats*(2*i1+1)*(2*i2+1)*(2*i3+1)*lattice_factor,' atoms'
            write(*,*) 'original list:'
            do i=1,nats
               write(*,'(a,i4,3(f12.5,x),a2,3x,f12.5)') '   ',i,ppos(:,i)
            enddo
         endif
         do j1=-(i1),i1
         do j2=-(i2),i2
         do j3=-(i3),i3
         do i=1,nats
            natt=natt+1  ! create one more atom
            ratx(:,natt)=ppos(:,i)*alatt+dble(j1)*a1+dble(j2)*a2+dble(j3)*a3
            iphatx(natt)=ppot(i)
              if((j1.eq.0.and.j2.eq.0.and.j3.eq.0).and.i.eq.absorber) then
        iphatx(natt)=0  ! absorber
            indexabs=natt
              endif
              if (iatph(iphatx(natt)).le.0) iatph(iphatx(natt))=natt
        if(nshift.gt.0) then
           do j=1,nshift
              natt=natt+1
            ratx(:,natt)=ratx(:,natt-j)+shift(1,j)*a1+shift(2,j)*a2+shift(3,j)*a3
            iphatx(natt)=ppot(i)
!          write(*,*) j1,j2,j3,i,j,ratx(:,natt),ratx(:,natt-1)
         enddo
        endif  !If P or H, no need to make extra atoms
         enddo
         enddo
         enddo
         enddo
! ratx is now an array in Carthesian coordinates and in Angstrom units.
         iz(0)=iz(ppot(absorber)) !KJ 11-2009 avoids triggering overlap part below


         ! The next section reorders the atom list "ratx" according to distance from absorber.
     ! It also corrects arrays distance(:), iphatx(:), and iatph(:).
         distance=dble(0)
         do i=1,natt
             distance(i)=dsqrt((ratx(1,i)-ppos(1,absorber)*alatt)**2+(ratx(2,i)-ppos(2,absorber)*alatt)**2&
     &       +(ratx(3,i)-ppos(3,absorber)*alatt)**2)
         enddo
         do i=1,natt
            k=i
            mindist=distance(i)
            do j=i,natt
               if(distance(j).lt.mindist) then
                  k=j
                  mindist=distance(j)
               endif
            enddo
            sss=ratx(:,i)
            ratx(:,i)=ratx(:,k)
            ratx(:,k)=sss
            distance(k)=distance(i)
            distance(i)=mindist
            j=iphatx(i)
            iphatx(i)=iphatx(k)
            iphatx(k)=j
            if(i.eq.indexabs) indexabs=k
            if(k.eq.indexabs) indexabs=i
            if(iatph(iphatx(i)).eq.i) iatph(iphatx(k))=k
            if(iatph(iphatx(k)).eq.k) iatph(iphatx(k))=i
         enddo
     iatph(0)=1  !KJ bugfix 6-10-2013

     ! Now TRUNCATE the list so that it only shows "full shells".  I.e. above we have a
     ! "boxy" list created by stacking unit cells.  Now cut out a sphere such that
     ! for each distance to the absorber in that list, all atoms at that distance are given.
     ! FEFF doesn't give a damn (since we have made sure we have all atoms up to rmax by design),
     ! but it confuses users to give incomplete shells at large distance.
     ! The radius of that sphere is the minimum of the 3 supercell lattice parameters
     mindist = min(dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)*i1, dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)*i1, dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)*i1)
     do i=1,natt
         if (distance(i).gt.mindist) then
            j=i
            exit
         endif
     enddo
     ! Now throw away all the atoms > j.
     ratx(:,j:natt)=0.d0
     distance(j:natt)=0.d0
     iphatx(j:natt)=-1
     natt=j-1
     !There really should be way more checks here ...

      endif
! !KJ


!KJ empty titles sometimes cause problems later on:
      if (ntitle.gt.0) then
       if((title(1).eq.'') .or. title(1).eq.' ') title(1)="no title"
    endif


!KJ check on treatment of core hole in case of k-space calculations :
      if(ispace.eq.0.and.nohole.ne.0.and.nohole.ne.2) then
         call wlog('COREHOLE RPA  is recommended for k-space.')
      endif

!KJ prepare true nohole calculation in k-space (NOHOLE 2 counts as a core hole calculation as far as kspace is concerned)
      if(ispace.eq.0.and.nohole.eq.0) then
       icorehole = 0   ! this will set the variable 'corehole' in FMS through reciprocal.inp
                     ! for NOHOLE 2, reapot will choose to ignore the core hole - see reapot.f90
    endif


! !KJ added this check 1-06
      if(magic.eq.1.and.(eels.ne.1)) then
          call wlog('To use MAGIC card you must have ELNES card.  Ignoring MAGIC card.')
          magic=0
        endif
! !KJ

!  !KJ another check for eels 1-06
      if((eels.eq.1).and.(aver.eq.1).and.(cross.eq.1)) then
          call wlog('WARNING : you have asked to calculate an orientation' // &
            ' averaged spectrum, but you have also asked         &
     &   to calculate cross-terms.  Averaging kills the cross terms. ' // &
       'Hence the program ignores your request and does not calculate cross terms.')
      endif
!  !KJ

!  !KJ  set up a variable needed for elnes 1-06
        if(eels.eq.1) then
          if(aver.eq.1) then
             ipstep=1
             ipmin=10
             ipmax=10
          else
            ipmin=1
            ipmax=9
            if(cross.eq.1) then
               ipstep=1
            else
               ipstep=4
            endif
          endif
        endif
!  !KJ

!KJ:   For MDFF calculations:  12-2010 and 03-2011
      !SANITY CHECKS
      if(do_nrixs.eq.1 .and. (imdff.eq.1 .or. imdff.eq.2)) then
       mixdff=.true.
    elseif(imdff.eq.1 .or. imdff.eq.2) then
       call wlog('ERROR - the selected MDFF option is only available with the NRIXS card.')
     call par_stop('RDINP')
    else
       mixdff=.false.
    endif

    if((imdff.eq.2).and.(nq.ne.2)) then
       call wlog("Current version of this type of MDFF calculation requires that you set nq=2 in the NRIXS card.")
     call par_stop(" ")
      endif
    if((imdff.eq.2).and.((dabs(cosmdff_dum)-dble(1)).lt.0.01d0)) then
       call wlog("Just letting you know - you're calculating a DFF using' // &
         ' MDFF technology.  Should work fine, but quite unnecessary.")
      endif
    if((imdff.eq.3) .and. .not.(cards_set(56).or.cards_set(57))) then
       call wlog("The selected MDFF option requires the ELNES or EXELFS card.  Aborting now.")
    endif

      !Initialize q' for "type 2" MDFF calculations
    if(imdff.eq.2) then
       if(qqmdff.ge.0.d0) then
        !make a new q' based on the parameters given in the MDFF card
      !q' == q, scaled to length qqmdff, rotated around x-axis by angle cosmdff_dum
      !I've freely chosen the rotation axis since user input only fixes 2 degrees of freedom for q'
      !If user wants to choose all 3, they should use the NRIXS list (i.e. don't specify qqmdff in feff.inp)
      dummy=qqmdff/qn(1)  ! q' / q
      qs(2,1)=qs(1,1)*dummy
      qs(2,2)=dummy*(  qs(1,2)*dcos(cosmdff_dum *pi/180.d0) + qs(1,3)*dsin(cosmdff_dum *pi/180.d0) )
      qs(2,3)=dummy*( -qs(1,2)*dsin(cosmdff_dum *pi/180.d0) + qs(1,3)*dcos(cosmdff_dum *pi/180.d0) )
     else
        !just use the q' as specified in the NRIXS card
     endif
    endif

    !initialize q,q' angles for MDFF+NRIXS calculation
    if(mixdff) then
     do iq=1,nq
     do iqq=1,nq
      cosmdff(iq,iqq)=dcos(pi/180.d0 * ((qs(iq,1)*qs(iqq,1)+qs(iq,2)*qs(iqq,2)+qs(iq,3)*qs(iqq,3))/(qn(iq)*qn(iqq))) )   ! cos<q,q'> = cos( q.q' /q /q')
     enddo
     enddo
    endif
!:KJ

!KJ :
!     NRIXS sanity checks and initialization:
      if (do_nrixs.eq.1) then
         xivec=qvec   !overwriting variables is OK since nrixs can never be combined with pola,elli,...
     if(enforce_alexis_exchange_policy) then
        call wlog('EXCHANGE ignored for NRIXS calculation.')
        vr0=0.d0
        vi0=0.d0
     endif
     le2=lj
     l2lp = 30 ! JK - 15 is the number used in the old NRIXS code. This may need to be variable in future.
     if(qaverage) then
        elpty=-nq ! JK - passing nrixs info on number of q point in elpty for now
     else
        elpty=nq
     end if
     if(xkmax.lt.0.d0) call wlog('Uniform energy mesh selected.')
         call init_feffq   !calculate xivnorm
     if (xivnorm.lt.0.01) call wlog('Warning - NRIXS calculation with very small q-vector.  Results may be bad.')
     if(ica.gt.0 .and. ica.lt.8) call wlog('Warning - SYMMETRY card ignored because of polarized NRIXS calculation.')
     if(.not. qaverage) then
        ica=7  !disables all symmetry in the Path Expansion.
     else
        ica=5  !disables most symmetry in the Path Expansion ; any rotation around z allowed.
     endif
       call make_qtrig
!OK     * set vr0=vri=0
!OK     * disable ELLIP, POLA, NSTAR, SPIN, CFAV, XNCD, RPHASES, TDLDA, XES, PMBSE
!OK     * warn about constant step energy grids for negative kmax - first find out how it works :)
!OK     * get lj into le2 (what Aleksi used - ripped from MULT card - and written to global.inp by wrtall)
!OK     * get nq and qvec into whatever Aleksi "stole" from the ?? card
!OK     * give warnings if momentum transfer small or large.  Aleksi put these in mkptz, but this gets tricky
!       since in the general version of feff, the vectors are used for different purposes, and I want to avoid
!       mess in low-level routines if possible (ie, better to have the "if nrixs" loop here than in mkptz).
!       Aleksi warns if smaller than 0.01 - random value.
!OK     * NOTE TO SELF : fix mkptz evec block
      endif

!-KJ


!     need smaller rgrid for nonlocal exchange
      if (ixc0.lt.0) ixc0 = 0 !for EXAFS, EXELFS (maybe NRIXS also?)
      if (mod(ixc, 10).ge.5 .and. rgrd.gt.0.03) rgrd=0.03d0
      if (mod(ixc0,10).ge.5 .and. rgrd.gt.0.03) rgrd=0.03d0
!     must use linear polarization to use nstar
      if (wnstar)  then
         if (ipol.ne.1)  then
            call wlog(' Must have linear polarization to use NSTAR.  NSTAR will be turned off.')
            wnstar = .false.
         endif
      endif

!     Do not use ihole .le. 0
      if (ihole .le. 0)  then
         call wlog(' Use NOHOLE to calculate without core hole.  Only ihole greater than zero are allowed.')
         call par_stop('RDINP')
      endif

!     Find out how many unique potentials we have
!     in POTENTIAL card
      nph = 0
      do 300  iph = nphxhardlimit, 0, -1
         if (iz(iph) .gt. -1)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue


!     cannot use OVERLAP and ATOMS cards together
      if (iatom .gt. 0 .and. iovrlp .gt. 0)  then
        call wlog(' Cannot use ATOMS and OVERLAP in the same feff.inp.')
        call par_stop('RDINP')
      endif

!     cannot use OVERLAP and CFAVERAGE   cards together
      if (novr(0) .gt. 0) then
!        OVERLAP is used, cannot do configuration average
         iphabs = 0
         nabs = 1
         rclabs = big
      endif

!     Must do LDOS-calculcation if HUBBARD card is chosen
      if (i_hubbard .eq. 2 .and. mldos .eq. 0) then
!KJ 7-2014 disabling this -- what if I already ran that stuff and just want to re-run "eels" or "ff2x"?
!        Could have a slightly smarter check based on CONTROL card.
!         mldos=1
!         emin=-20
!         emax=20
!         eimag=0.05
!         neldos=nex
!         l_hubbard=2
      endif

!     Must have nspx=2 if doing HUBBARD
      if (i_hubbard .eq. 2 .and. nspx .lt. 2) then
         call wlog('ERROR -- HUBBARD calculation requires compiling code with nspx=2, but you have nspx < 2.')
         call wlog('(nspx is the number of spins)')
         call wlog('Quitting now.  Please recompile the code appropriately (ask the authors for help), or remove the HUBBARD card.')
         stop
      endif

!     Set nspu to 1 (spin-averaged calculation) or 2 (spin-polarized calculation)
      if (i_hubbard .eq. 2) then
         nspu=2
      elseif (ispin.ne.0) then  ! A more sophisticated criterium may be needed here
         nspu=MIN(2,nspx)
      else
         nspu=1
      endif

!     Must have central atom
      if (iz(0) .le. 0)  then
         if (iphabs .gt. 0) then
!           central atom is of the iphabs type
            iz(0) = iz(iphabs)
            potlbl(0) = potlbl(iphabs)
            lmaxsc(0) = lmaxsc(iphabs)
            lmaxph(0) = lmaxph(iphabs)
            xion(0) = xion(iphabs)
         else
            call wlog(' No absorbing atom (unique pot 0) was defined.')
            call par_stop('RDINP')
         endif
      endif

!     No gaps allowed in unique pots.  Make sure we have enough
!     to overlap all unique pots 0 to nph.
      if (iphabs.gt.0 .and. iatph(0).le.0)   iatph(0) = iatph(iphabs)
      do 340  iph = 0, nph
         if (iatph(iph) .le. 0  .and.  novr(iph) .le. 0)  then
!           No model atom, no overlap cards, can't do this unique pot
            write(slog,'(a,i8)') ' No atoms or overlap cards for unique pot ', iph
            call wlog(slog)
            call wlog(' Cannot calculate potentials, etc.')
            call par_stop('RDINP-')
         endif
!        by default freeze f-electrons and reset lmaxsc=2
         if (iunf.eq.0 .and. lmaxsc(iph).gt.2) then
        write(slog,'(a,i4,a)') 'Resetting lmaxsc to 2 for iph = ',iph,'.  Use  UNFREEZE to prevent this.'
        call wlog(slog)
        lmaxsc(iph)=2
     endif
  340 continue

!     Need number of atoms of each unique pot, count them.  If none,
!     set to one. Do statistics for all atoms in feff.inp.
      do iph = 0, nph
        if (lxnat.eq.0) then
          xnatph(iph) = 0
          do iat = 1, natt
              if (iphatx(iat) .eq. iph)  xnatph(iph) = xnatph(iph)+1
          enddo
          if (iph.gt.0 .and. iph.eq.iphabs) xnatph(iph) = xnatph(iph)-1
        else
          if (xnatph(iph).le. 0.01) then
            if (iph.eq.0) then
              xnatph(iph) = 0.01d0
            else
              write (slog,'(a,i4)') ' Inconsistency in POTENTIAL card is detected for unique pot ', iph
              call wlog (slog)
              call wlog (' Results might be meaningless.')
            endif
          endif
        endif
        if (xnatph(iph) .le. 0)  xnatph(iph) = 1
      enddo
      if (lxnat.ne.0) then
!        normalize statistics to have one absorber
         do 351 iph = 1, nph
  351    xnatph(iph) = xnatph(iph) /xnatph(0)
         xnatph(0) = 1
      endif
      xnat = 0
      do 352 iph = 0,nph
  352 xnat = xnat + xnatph(iph)

!     Find distance to nearest and most distant atom (use overlap card
!     if no atoms specified.)
      if (natt .lt. 2)  then
         ratmin = rovr(1,0)
         ratmax = rovr(novr(0),0)
      else
         ratmax = 0
         ratmin = 1.0e10
         iatabs = iatph(0)
         icount = 0
         if (iatabs.le.0) iatabs = iatph( iphabs)
         if (iatabs.le.0) call par_stop('RDINP fatal error: iatabs=NaN')

         do 412  iat = 1, natt
           if (iphatx(iat) .eq. iphabs .or. iphatx(iat).eq.0)  icount = icount +1
           if (iat.ne.iatabs) then
!           skip absorbing atom
            tmp = dist (ratx(1,iat), ratx(1,iatabs))
            if (tmp .gt. ratmax)  ratmax = tmp
            if (tmp .lt. ratmin)  ratmin = tmp
           endif
  412    continue
         if (nabs.le.0) nabs = icount
      endif

!     Set total volume
      if (totvol.gt.0) totvol = totvol * ratmin**3 * xnat

!     Set rfms if they are too small
      if (rfms1 .lt. ratmin) rfms1 = -1.e0
      if (rfms2 .lt. ratmin) rfms2 = -1.e0
      if (rfms2 .lt. ratmin .and. ispec.lt.2) ispec = - ispec
      if (rfms2 .lt. ratmin .and. ispec.eq.3) ispec = - ispec
!     if ispec.le.0 MS expansion will be used, else - FMS method.


!     Set rmax if necessary
      if (rmax.le.0 .and. nss.le.0 .and. ispec.le.0)  then
!        set to min (2+ times ratmin, ratmax) (magic numbers to
!        avoid roundoff, note that rmax is single precision).
         rmax = min (2.2 * ratmin, 1.01 * ratmax)
      endif

!     Set core hole lifetime (central atom quantity) and s02
!     KJ added 'if' construction and userChLifetime
      iph = 0
      if (userchl) then
         if (userChLifetime.gt.dble(0)) then
            gamach = userChLifetime
         else
            call setgam(iz(iph),ihole,gamach)
            gamach=min(gamach,abs(userChLifetime))
         endif
      else
         call setgam (iz(iph), ihole, gamach)
      endif
      if (s02.eq.1.d0) s02=s02h
      write(slog,'(a,f7.3,a)') 'Core hole lifetime is ',gamach,' eV.'
      call wlog(slog)

!KJ NOW DEAL WITH DIMENSIONS FOR DYNAMICAL ALLOCATION
! These will be written to dimensions.dat using call to WriteDimensions
! Set the appropriate values here :
!   1/ code figures out dimensions to corresponding input options
!   2/ code truncates if a/ exceeds user limits specified in DIMS card b/ hardcoded limits in COMMON/m_dimsmod.f90
!   (the DIMS card overrides the values in COMMON/m_dimsmod.f90)

!   1/
      lx=0
      do iph=0,nph
       lx=max(lx,lmaxph(iph))
       lx=max(lx,lmaxsc(iph))
      enddo

!   1/ :
      rdims=max(rfms1,rfms2,rmax)
      nclusx=0
      do iat = 1, natt
         if (dist(ratx(:,iat),ratx(:,iatabs)).le.rdims) then
            nclusx=nclusx+1
         endif
      enddo
      nphu=nph  !number of unique potential types
      if (abs(ispin).eq.1 .or. abs(ispin).eq.2) then  !SPIN card set; do spin-polarized calculation:
         nspu=2
         ! JK - set spinph here if it was not set in POTENTIALS card.
         DO iph = 0, nph
            IF(spinph(iph).LT.-1.d8) THEN
               IF(iph.EQ.0) THEN
                  PRINT*, 1
                  spinph(iph) = getspin(iz(iph),ihole,xion(iph),iph)
               ELSE
                  spinph(iph) = getspin(iz(iph),0,xion(iph),iph)
               END IF
               call wlog('No spin set in POTENTIALS card. Using default spins:')
               call wlog('iph   spinph')
               WRITE(slog,'(I3,X,f3.1)') iph, spinph(iph)
               call wlog(slog)
            END IF
         END DO

      elseif(i_hubbard .ge. 2) then  !HUBBARD-U calculation:
         nspu=2
         spinph = 0.d0
      else    ! regular spin-averaged calculation:
         nspu=1
         spinph = 0.d0
      endif

!    2/ : happens inside write_dimensions
      call write_dimensions(nclusxuserlimit,lxuserlimit)

!    Now fix lmax values to final lx-value:
      do iph=0,nph
       lmaxph(iph)=min(lx,lmaxph(iph))
       lmaxsc(iph)=min(lx,lmaxsc(iph))
      enddo

!    Rmax can't really be fixed here since nclusx is a number, not a distance.
!    We must rely on subsequent programs to use nclusx as input and cut off.
!    I.e., lmaxph/sc can be trusted 'blindly' from here on ; but rfms1/2 cannot.

!KJ done with dimensions for dynamical allocation


!KJ   CHECK THAT NO INVALID COMBINATION OF CARDS IS USED :
      call consistency_checker(cards_set)
    if((cards_set(9) .and. rmax.gt.2.5d0) .and. (cards_set(37) .and. rfms2.gt.2.5d0)) &
       call wlog("WARNING  You are using RPATH and FMS.  This is syntactically permitted, but it's almost always a bad idea.")


!     Convert everything to code units, and use rmult factor
!     rmax is for pathfinder, so leave it in Ang.
      rmax = rmax * rmult
      rfms1 = rfms1 * rmult
      rfms2 = rfms2 * rmult
      totvol = totvol * rmult**3
!     Use rmult factor.  Leave distances in Ang.
      do 430  iat = 1, natt
         do 420  i = 1, 3
            ratx(i,iat) = ratx(i,iat) * rmult
  420    continue
  430 continue
      do 460  iph = 0, nph
         do 450  iovr = 1, novr(iph)
            rovr(iovr,iph) = rovr(iovr,iph) * rmult
  450    continue
  460 continue
      do 462  iss = 1, nss
!        rss used only to make paths.dat, so leave it in Angstroms.
         rss(iss) = rss(iss) * rmult
  462 continue

!     Clean up control flags
      if (mpot .ne. 0)  mpot = 1
      if (mphase .ne. 0)  mphase = 1
      if (mfms .ne. 0)  mfms = 1
      if (mpath  .ne. 0)  mpath = 1
      if (mfeff  .ne. 0)  mfeff = 1
      if (mchi   .ne. 0)  mchi = 1
      if (nss    .le. 0)  ms = 1
      if (ifolp  .ne. 0)  iafolp = -1
      if (natt.le.0) then
!       Overalp geometry
        mfms = 0
        mpathold=mpath !KJ 7-06 for writing paths.dat
    ! mpath : will path module actually be run?  !KJ
    ! mpathold : do we want to create a new list of paths?  !KJ
        mpath = 0
        ms = 0
!       no SCF loop
        nscmt = 0
        do 464 iph = 0, nph
          if (novr(iph).le.0) call par_stop('Bad OVERLAP cards.')
  464   continue
      endif

     !!KJ Dec 2013 I don't quite understand some of the control flag logic above.  It looks like FMS is turned on for many EXAFS calculations.
     !!KJ It doesn't *really* matter since rfms<0 but it looks confusing ...
     !if (rfms2 .lt. 0.01d0) mfms=0 !KJ
     !KJ Jan 2014: I take it back.  Turning fms off means fms.bin does not get written.  In some cases this leads to problems.
     !KJ E.g. the XES BN example on Windows crashes in this case (as "ne" is not initialized properly in ff2xmu).
     !KJ I would still like a more elegant solution than relying on a file full of 0s to get written.
     !KJ But until I have the time to test it properly -- let's just turn the old option back on.

      if (iafolp .ge. 0) folp(0:nphxhardlimit)=folpx

      if (ntitle .le. 0)  then
         ntitle = 1
         title(1) = 'Once upon a time ...' !KJ ;)
      endif
      do i = 1, ntitle
         ltit(i) = istrln (title(i))
      enddo
      nttl = ntitle
      ! Check RIXS values for nEdges.
      IF(RixsI%m_run.EQ.0) RixsI%nEdges = 1
!     write atoms.dat, global.inp, modN.inp and ldos.inp
      call wrtall

! Write the dmdw.inp file
      open(unit=65,file='dmdw.inp',status='unknown',iostat=ios)
      if ( Use_DMDW ) then
        write(65,fmt='(i4)') 1
        write(65,fmt='(i4)') DMDW_Order
! Modifying this to match the new more flexible format
!       write(65,fmt='(i4,2f11.3)') 1, tk, tk
        write(65,fmt='(i4,f11.3)') 1, tk
        write(65,fmt='(i4)') DMDW_Type
        write(65,fmt='(a)') trim(dym_File)
! Now we write the path selectors for the standalone run of the dmdw module
! We will be adding possible choices in the future, maybe make it as
! flexible as the input of dmdw itself.
! Routes:
!     0              Don't do anything
!     1              All SS paths from absorber (assumed to be atom 1, for now)
!     2              Same as 1 + all DS paths
!     3              Same as 2 + all TS paths
!    11              All SS paths
!    12              Same as 1 + all DS paths
!    13              Same as 2 + all TS paths
! NOTE: This is not "pretty" code, will fix later
! Calculate the maximum distance within the input cluster to detemine
! "safe" path cutoffs.
        mxDij2 = 0.0
        do iiAtom = 1,natt-1
          do jjAtom = iiAtom+1,natt
            Dij2 = sum((ratx(:,iiAtom)-ratx(:,jjAtom))**2)
            if ( Dij2 > mxDij2 ) then
              mxDij2 = Dij2
            end if
          end do
        end do
        if ( DMDW_Route == 0 ) then
          write(65,fmt='(i4)')  0
        end if
        if ( DMDW_Route == 1 ) then
          write(65,fmt='(i4)')  1
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 1, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
        end if
        if ( DMDW_Route == 2 ) then
          write(65,fmt='(i4)')  2
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 1, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
          write(65,fmt='(4i4,4x,f7.2)') &
                3, 1, 0, 0,    1.1*sqrt(mxDij2)*2.0*1.8897
        end if
        if ( DMDW_Route == 3 ) then
          write(65,fmt='(i4)')  3
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 1, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
          write(65,fmt='(4i4,4x,f7.2)') &
                3, 1, 0, 0,    1.1*sqrt(mxDij2)*2.0*1.8897
          write(65,fmt='(5i4,   f7.2)') &
                4, 1, 0, 0, 0, 1.1*sqrt(mxDij2)*3.0*1.8897
        end if
        if ( DMDW_Route == 11 ) then
          write(65,fmt='(i4)')  1
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 0, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
        end if
        if ( DMDW_Route == 12 ) then
          write(65,fmt='(i4)')  2
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 0, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
          write(65,fmt='(4i4,4x,f7.2)') &
                3, 0, 0, 0,    1.1*sqrt(mxDij2)*2.0*1.8897
        end if
        if ( DMDW_Route == 13 ) then
          write(65,fmt='(i4)')  3
          write(65,fmt='(3i4,8x,f7.2)') &
                2, 0, 0,       1.1*sqrt(mxDij2)*1.0*1.8897
          write(65,fmt='(4i4,4x,f7.2)') &
                3, 0, 0, 0,    1.1*sqrt(mxDij2)*2.0*1.8897
          write(65,fmt='(5i4,   f7.2)') &
                4, 0, 0, 0, 0, 1.1*sqrt(mxDij2)*3.0*1.8897
        end if
      else
        write(65,fmt='(i4)') -999
      end if
      close(65)

!     In case of OVERLAP and SS calculateions write 'paths.dat'
!     without invoking the pathfinder. Single scattering paths only.
      if (nss .gt. 0  .and.  mpathold .eq. 1)  then !KJ 7-06 : fix bug
         open (unit=1, file='paths.dat', status='unknown', iostat=ios)
         call chopen (ios, 'paths.dat', 'rdinp')
         do 750  i = 1, ntitle
            write(1,748)  title(i)(1:ltit(i))
  748       format (1x, a)
  750    continue
         write(1,751)
  751    format (' Single scattering paths from ss lines cards in feff input')
         write(1,706)
  706    format (1x, 71('-'))
         do 760  iss = 1, nss
            if (rmax.le.0  .or.  rss(iss).le.rmax)  then
!              NB, rmax and rss are in angstroms
               write(1,752) indss(iss), 2, degss(iss), rss(iss)
  752          format ( 2i4, f8.3,'  index,nleg,degeneracy,r=', f8.4)
               write(1,766)
  766          format (' single scattering')
               write(1,754) rss(iss), zero, zero, iphss(iss), potlbl(iphss(iss))
               write(1,753) zero, zero, zero, 0, potlbl(0)
  753          format (3f12.6, i4,  1x, '''', a6, '''', '  x,y,z,ipot')
  754          format (3f12.6, i4,  1x, '''', a6, '''')
            endif
  760    continue
         close (unit=1)
      endif

      call wlog('Your calculation:')
      do i = 1, ntitle
         call wlog(' ' // title(i)(1:ltit(i)))
      enddo
      edgename='  '
      call setedg(edgename,ihole)
      write(slog,'(a,x,a,x,a,x,a,x,a,x,a,x,a,x,a)') trim(potlbl(0)),trim(edgename),'edge',trim(name_spectroscopy(cards_set)),'using',trim(name_corehole(nohole)),'corehole.'
      call wlog(slog)
      call wlog('Using:  '//list_features(cards_set))
      call wlog('Using cards:  '//list_cards(cards_set))
      call wlog(' ')


!     if user doesn't want geom.dat, don't do it
      if (nogeom)  then
!        don't delete geom.dat when done with it either...
         if (ipr4 .lt. 2)  ipr4 = 2
         if (nabs.gt.1) call par_stop('NOGEOM and CFAVERAGE are incompatible')
      else
        iabs = 1
!    !KJ 1-06 : If the user does EELS and doesn't calculate cross terms for an
!       orientation sensitive calculation, FEFF mustn't change the
!       coordinate system, as this would lead to the appearance of
!       cross terms after all.  Therefore, I added an argument to the
!       calling sequence of ffsort.
!       To be precise, giving '.false.' disables the call of ffsort to mkptz.
!       Giving '.true.' makes ffsort work exactly as it always has.
        if((eels.eq.1)) then
           call ffsort(iabs,nss,.false.) !KJ 7-06 added nss
        else
           call ffsort(iabs,nss,.true.) !KJ 7-06 added nss
        endif   !KJ end my changes
       endif
!KJ fix later 8/06 : if ispace=0, ffsort will make path expansion use different coordinates from fms !!

       ceels=(eels.eq.1) !KJ 5-6 for monolithic version

  400 call par_barrier
      call par_end

!     sub-program exchange
      if(master)call WipeErrorfileAtFinish
      stop
!     return

!     normal end of rdinp

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('RDINP fatal error.')

      end program rdinp

      subroutine phstop (iph,line)
    use dimsmod, only: nphxhardlimit
      implicit double precision (a-h, o-z)
      character*(*) line
      character*512 slog
      if (iph .lt. 0  .or.  iph .gt. nphxhardlimit)  then
         write(slog,10) iph, nphxhardlimit, line
         call wlog(slog)
   10    format (' Unique potential index', i5, ' out of range.',       &
     &           ' Must be between 0 and', i5, '.  Input line:',        &
     &           1x, a)
         call par_stop('RDINP - PHSTOP')
      endif
      return
      end
