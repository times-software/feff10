!     Common blocks with all input data
!     the common
!c    atoms.dat
      integer  natt
      integer iphatx(nattx)
      double precision  ratx(3,nattx)
      common /geom/ ratx, iphatx, natt
!c    geom.dat
!       integer  nat
!       integer iatph(0:nphx)
!       integer iphat(natx)
!       double precision  rat(3,natx)
!       common /geom/ ratx, iphatx, natt
!c    global.inp
!       configuration average
      integer iphabs
!     global polarization data
      integer  ipol, ispin, le2
      double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
      complex*16 ptz(-1:1, -1:1)
        integer ldecmx !KJ 7-09 for feff8q
  
	common /global/ ptz, evec, xivec, spvec, elpty, angks, rclabs,    &
     &     ipol, ispin, le2, iphabs,ldecmx !KJ 7-09 added ldecmx for feff8q
!     c    mod1.inp
      character*80 title(nheadx)
!     integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec,
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, iunf,             &
     &     nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx)
      integer lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), spinph(0:nphx)
      double precision  xion(0:nphx)
      logical ExternalPot
!     for OVERLAP option
      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
      double precision  rovr(novrx,0:nphx)
      common /mod1/ title, xion, xnatph, spinph, folp, gamach, rgrd,    &
     &     ca1, ecv, totvol, rovr, rfms1, iz, lmaxsc, mpot, nph, ntitle,&
     &     ihole, ipr1, iafolp, nmix,nohole,jumprm, inters,             &
     &     nscmt, icoul, lfms1, novr, iphovr, nnovr, iunf, ExternalPot
!     c    ldos.inp
      integer mldos, lfms2
      double precision emin, emax, eimag, rfms2
      common /mod7/ emin, emax, eimag, rfms2, mldos, lfms2
!c    mod2.inp
!     integer mphase, ipr2, ixc, ixc0, vr0, vi0, ispec, lreal, lfms2
      integer mphase, ipr2, ixc, ixc0, ispec, lreal, l2lp, iPlsmn
      integer lmaxph(0:nphx), iGammaCH, iGrid
      character*6  potlbl(0:nphx)
!     double precision rgrd, rfms2, gamach, xkstep, xkmax, vixan
      double precision xkstep, xkmax, vixan, vr0, vi0
      common /mod2/ xkstep, xkmax, vixan, vr0, vi0,                     &
     &     lmaxph, mphase, ipr2, ixc, ixc0, ispec, lreal, l2lp,         &
     &     izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis, iPlsmn,       &
     &     iGammaCH, iGrid, potlbl
!     c    mod3.inp
      integer mfms, idwopt, minv
!     integer lmaxph(0:nphx)
!     real rfms2, rprec, rdirec, toler1, toler2
      real rprec, rdirec, toler1, toler2
      double precision   tk, thetad, sig2g
      common /mod3/ tk, thetad, sig2g, rprec, rdirec, toler1,           &
     &       toler2,  mfms, idwopt, minv
!     c    mod4.inp
      integer  mpath, ms, nncrit, nlegxx, ipr4, ica  !KJ added ica 6-06
!     real critpw, pcritk, pcrith,  rmax, rfms2
      real critpw, pcritk, pcrith,  rmax
      common /mod4/ critpw, pcritk, pcrith,  rmax,                      &
     &       mpath, ms, nncrit, nlegxx, ipr4, ica !KJ added ica 6-06
!     c    mod5.inp
      integer  mfeff, ipr5, iorder
      logical  wnstar
      double precision critcw
      common /mod5/ critcw, mfeff, ipr5, iorder, wnstar
!     c    mod6.inp
!     integer  mchi, ispec, idwopt, ipr6, mbconv
!     double precision  vrcorr, vicorr, s02, alphat, sig2g
      integer  mchi, ipr6, mbconv, absolu !KJ added absolu 3-06
      double precision  vrcorr, vicorr, s02, alphat, thetae
      common /mod6/ vrcorr, vicorr, s02, alphat, thetae,                &
     &     mchi, ipr6, mbconv, absolu   !KJ added absolu 3-06
!     c    so2.inp  
      integer  msfconv, ipse, ipsk
      double precision wsigk, cen
      character(12) cfname
      common /sfconvinp/ wsigk, cen, cfname, msfconv, ipse, ipsk
      
!     c    eels.inp
!     EELS variables  !KJ 1-06 this section added for ELNES, EXELFS, MAGIC cards
      real*8 ebeam, aconv, acoll, thetax, thetay, emagic
      integer eels, relat, aver, cross, iinput,spcol
      integer nqr,nqf,magic
      integer ipmin,ipmax,ipstep
      common /eelsva/ ebeam,aconv,acoll,thetax,thetay,emagic,magic,     &
     &     nqr, nqf, aver, cross, relat, iinput, spcol,ipmin, ipmax,    &
     &     ipstep, eels
!     !KJ end
!KJ for the energy grid card EGRID :
      integer iegrid,egrid3a
        real*8 egrid3b,egrid3c
        character*100 egridfile
          common /egridvars/ egrid3b,egrid3c,iegrid,egrid3a,egridfile
!KJ

! Added by Fer
! Used to correct the excitation energy for chemical shifts
      integer  ChSh_Type
      common /Chem_Shft/ ChSh_Type

