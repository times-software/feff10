!      subroutine rdpotp_fs (path, ntitle, title, rnrmav, xmu, vint, rhoint,&
!     &                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
!     &                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
!     &                  dgc0, dpc0, dgc, dpc, adgc, adpc,               &
!     &                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
!     &                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
!     &                  inters, totvol, iafolp, xion, iunf, iz, jumprm, &
!     &                  nph)

      subroutine rdpotp_fs (path, ntitle, title, xnatph, rnrm, iz, nph)

!  opens pot.bin file at pathname path and reads following information
!  General:
!     ntitle - number of title lines
!     title  - title itself
!     emu    - edge position (x-ray energy for final state at Fermi level)
!  Muffin-tin geometry
!     rmt    - muffin-tin radii
!     imt    - index of radial grid just below muffin-tin radii
!     rnrm   - Norman radii
!     inrm   - index of radial grid just below Norman radii
!     rnrmav - average Norman radius
!     folp   - overlap parameter for rmt
!     folpx  - maximum value for folp
!     xnatph - number of atoms of each potential type
!  Atomic orbitals info (Dirac spinors)
!     dgc0   - upper component for initial orbital
!     dpc0   - lower component for initial orbital
!     dgc    - upper components for all atomic orbitals
!     dpc    - lower components for all atomic orbitals
!     adgc   - development coefficient for upper components
!     adpc   - development coefficient for lower components
!     xnval  - number of valence electrons for each atomic orbital
!              used for core-valence separation and non-local exchange
!     eorb  - atomic enrgy of each orbital for the absorber
!  Electron density information 
!     rhoint - interstitial density
!     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
!     xf     - estimate of momentum at Fermi level from rhoint
!     edens  - total electron density
!     edenvl - density from valence electrons
!     dmag   - density for spin-up minus density for spin-down
!     qtotel - total charge of a cluster
!     qnrm   - charge accumulated inside Norman sphere as result of SCF
!     xnmues - occupation numbers of valence orbitals from SCF procedure
!  Potential information
!     xmu    - Fermi level position
!     ecv    - core-valence separation energy
!     vint   - muffin-tin zero energy (interstitial potential)
!     vclap  - Coulomb potential
!     vtot   - vclap + xc potential from edens
!     vvalgs - vclap + xc potential from edenvl (EXCHANGE 5 model)
!  Specific data for convolution with excitation spectrum (see mbconv)
!     s02    - many body reduction factor S_0^2 
!     erelax - estimate of relaxation energy = efrozen - emu, where
!              efrozen is edge position estimate from Koopmans theorem
!     wp     - estimate of plasmon frequency from rhoint
      use dimsmod, only: nheadx,nex,nphx=>nphu

      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 41, 0:nphx), dpc(251, 41, 0:nphx)
      dimension adgc(10, 41, 0:nphx), adpc(10, 41, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(41,0:nphx), qnrm(0:nphx) !, xnmues(0:lx,0:nphx)
      dimension eorb(41)
      dimension iorb(-5:4,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)
      double precision kappa(41) !KJ bugfix previously implicitly declared as integer

      character*80 title(nheadx)
      character*512 path, slog

      dimension dum(13)
      logical v84 !true if this pot.bin written in feff84 format

  10  format(a)
   20 format (bn, i15)

      !initialize values that may not be saved in pot.bin
      do i=1,41
        eorb(i)=0.0
        kappa(i)=0.0
      enddo

!     open (unit=3, file='pot.bin', status='old')
      open (unit=3, file=path, status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,     &
     &            jumprm, iunf
  30  format(9(1x,i4))
!     nph and npadx are not passed to calling subroutine
!     not sure why this was done, but i'm going to pass nph to the
!     calling routine for use there. mpp 11/2005
      do 133  i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
  133 continue
!     Misc double precision stuff from pot.bin
      call rdpadd(3, npadx, dum(1), 13)
      rnrmav = dum(1)
      xmu    = dum(2)
      vint   = dum(3)
      rhoint = dum(4)
      emu    = dum(5)
      s02    = dum(6)
      erelax = dum(7)
      wp     = dum(8)
      ecv    = dum(9)
      rs     = dum(10)
      xf     = dum(11)
      qtotel = dum(12)
      totvol = dum(13)

!     read imt
      read (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      call rdpadd(3, npadx, rmt(0), nph+1)
!     read inrm
      read (3, 40) (inrm(i),i=0,nph)
      read (3, 40) (iz(i),i=0,nph)

      !read in the next line to see if it is pad ascii data (as in
      !feff82 format pot.bin) or a string of integers (as in feff84
      !format pot.bin)
      read (3,10) slog
      !backspace so we can reread the last line to extract the data
      backspace 3
      !determine if we should read the kappa array or not
      call triml(slog)
      if (slog(1:1).ne.'!') then
        !we have a feff84 format pot.bin and must read kappa
        v84=.true.
        read (3, 40) (kappa(i),i=1,41)
      else
        v84=.false.
      endif
!     write(slog,fmt="('v84: ',i3)") v84
!     call wlog (slog)
        

!     call wlog('rnrm')
      call rdpadd(3, npadx, rnrm(0), nph+1)
!     call wlog('folp')
      call rdpadd(3, npadx, folp(0), nph+1)
!     call wlog('folpx')
      call rdpadd(3, npadx, folpx(0), nph+1)
!     call wlog('xnatph')
      call rdpadd(3, npadx, xnatph(0), nph+1)
!     call wlog('xion')
      call rdpadd(3, npadx, xion(0), nph+1)
!     call wlog('dgc0')
      call rdpadd(3, npadx, dgc0(1), 251)
!     call wlog('dpc0')
      call rdpadd(3, npadx, dpc0(1), 251)
!     call wlog('dgc')
      call rdpadd(3, npadx, dgc(1,1,0), 251*41*(nph+1) )
!     call wlog('dpc')
      call rdpadd(3, npadx, dpc(1,1,0), 251*41*(nph+1) )
!     call wlog('adgc')
      call rdpadd(3, npadx, adgc(1,1,0), 10*41*(nph+1) )
!     call wlog('adpc')
      call rdpadd(3, npadx, adpc(1,1,0), 10*41*(nph+1) )
!     call wlog('edens')
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
!     call wlog('vclap')
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
!     call wlog('vtot')
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
!     call wlog('edenvl')
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
!     call wlog('vvalgs')
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
!     call wlog('dmag')
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
!     call wlog('xnval')
      call rdpadd(3, npadx, xnval(1,0), 41*(nph+1) )
!     call wlog ('this is where i should read eorb')
      if (v84) call rdpadd(3, npadx, eorb(1), 41)
!     call wlog ('just did! (or didnt'))
      do 50 iph=0,nph
 50   read (3, 60) (iorb(i,iph),i=-5,4)
 60   format(8(1x,i2))
      call rdpadd(3, npadx, qnrm(0), nph+1 )
!      nn = (lx+1)*(nph+1)
!      call rdpadd(3, npadx, xnmues(0,0), nn )
      close (unit=3)
      
      return
      end
