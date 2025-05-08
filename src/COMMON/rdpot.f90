!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdpot.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2010/11/30 19:41:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,      &
                       emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
                       imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
                       dgc0, dpc0, dgc, dpc, adgc, adpc,               &
                       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
                       eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
                       inters, totvol, iafolp, xion, iunf, iz, jumprm)
!  opens pot.bin file and reads following information
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

!     KJ 11-2010
!     Added option to read feff8(q) potentials.
!     Only difference is that kappa and eorb are not present by default.
!     I've changed feff8q to output these two to a file called forf9.dat .
!     If this file is present, rdpot will use it.
!     There is currently no runtime option for the user to influence this behavior.

      use dimsmod, only: nphx=>nphu, nheadx, lx
      implicit double precision (a-h, o-z)
!Changed the dimensions to 41 to account for superheavy elements. Pavlo Baranov 07/2016

      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 41, 0:nphx), dpc(251, 41, 0:nphx)
      dimension adgc(10, 41, 0:nphx), adpc(10, 41, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(41,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(41), kappa(41)
      dimension iorb(-5:4,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)

      character*80 title(nheadx)

      logical :: feff9format=.true. !KJ
      dimension dum(13)

  10  format(a)
   20 format (bn, i15)

!KJ:
      open(155,file='forf9.dat',form='formatted',status='old',err=132)
	  read(155,*,end=132,err=132)
	  read(155,*,end=132,err=132) kappa
	  read(155,*,end=132,err=132)
	  read(155,*,end=132,err=132) eorb
	  close(155)
	  feff9format=.false. !if anything went wrong, we'll try to treat the calculation as "feff9"
 132  continue
      if(.not.feff9format) then
	     call wlog(':INFO FEFF9 is using a pot.bin and forf9.dat file from a FEFF8 source in routine rdpot.')
      end if
       
!:KJ

      open (unit=3, file='pot.bin', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp, jumprm, iunf

!     nph and npadx are not passed to calling subroutine
      do i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
      enddo
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

      read (3, 40) (imt(i),i=0,nph)
      call rdpadd(3, npadx, rmt(0), nph+1)
      read (3, 40) (inrm(i),i=0,nph)
      read (3, 40) (iz(i),i=0,nph)
      if(feff9format) read (3, 40) (kappa(i),i=1,41)  !KJ added "if"
      call rdpadd(3, npadx, rnrm(0), nph+1)
      call rdpadd(3, npadx, folp(0), nph+1)
      call rdpadd(3, npadx, folpx(0), nph+1)
      call rdpadd(3, npadx, xnatph(0), nph+1)
      call rdpadd(3, npadx, xion(0), nph+1)
      call rdpadd(3, npadx, dgc0(1), 251)
      call rdpadd(3, npadx, dpc0(1), 251)
      call rdpadd(3, npadx, dgc(1,1,0), 251*41*(nph+1) )
      call rdpadd(3, npadx, dpc(1,1,0), 251*41*(nph+1) )
      call rdpadd(3, npadx, adgc(1,1,0), 10*41*(nph+1) )
      call rdpadd(3, npadx, adpc(1,1,0), 10*41*(nph+1) )
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, xnval(1,0), 41*(nph+1) )
      if(feff9format) call rdpadd(3, npadx, eorb(1), 41) !KJ added "if"
      do iph=0,nph
         read (3, 60) (iorb(i,iph),i=-5,4)
      enddo
      call rdpadd(3, npadx, qnrm(0), nph+1 )
      nn = (lx+1)*(nph+1)
      call rdpadd(3, npadx, xnmues(0,0), nn )

      close (unit=3)
	  
 30   format(9(1x,i4))
 40   format(20(1x,i4))
 60   format(8(1x,i2))

      return
      end



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdpot.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2010/11/30 19:41:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   module pot_bin
!
!      use dimsmod, only: nphx=>nphu, nheadx, lx
!      implicit none
!
!      integer imt(0:nphx), inrm(0:nphx)
!      real*8 folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
!      real*8 dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
!      real*8 adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
!      real*8 edens(251, 0:nphx), vclap(251, 0:nphx)
!      real*8 vtot(251, 0:nphx), edenvl(251, 0:nphx)
!      real*8 vvalgs(251, 0:nphx), dmag(251, 0:nphx)
!      real*8 xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
!      real*8 eorb(30), xion(0:nphx), rmt(0:nphx),  rnrm(0:nphx)
!      integer kappa(30)
!      integer iorb(-4:3,0:nphx), iz(0:nphx)
!      real*8 xnatph(0:nphx)
!
!      character*80 title(nheadx)
!      integer ntitle, nohole, ihole, inters, iafolp, iunf, jumprm
!      real*8 rnrmav, xmu, vint, rhoint, emu, s02, erelax, wp, ecv, rs, xf, qtotel, totvol, xion
!
!
!      contains
!
!      subroutine rdpot_mod
!
!!  opens pot.bin file and reads following information
!!  General:
!!     ntitle - number of title lines
!!     title  - title itself
!!     emu    - edge position (x-ray energy for final state at Fermi level)
!!  Muffin-tin geometry
!!     rmt    - muffin-tin radii
!!     imt    - index of radial grid just below muffin-tin radii
!!     rnrm   - Norman radii
!!     inrm   - index of radial grid just below Norman radii
!!     rnrmav - average Norman radius
!!     folp   - overlap parameter for rmt
!!     folpx  - maximum value for folp
!!     xnatph - number of atoms of each potential type
!!  Atomic orbitals info (Dirac spinors)
!!     dgc0   - upper component for initial orbital
!!     dpc0   - lower component for initial orbital
!!     dgc    - upper components for all atomic orbitals
!!     dpc    - lower components for all atomic orbitals
!!     adgc   - development coefficient for upper components
!!     adpc   - development coefficient for lower components
!!     xnval  - number of valence electrons for each atomic orbital
!!              used for core-valence separation and non-local exchange
!!     eorb  - atomic enrgy of each orbital for the absorber
!!  Electron density information 
!!     rhoint - interstitial density
!!     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
!!     xf     - estimate of momentum at Fermi level from rhoint
!!     edens  - total electron density
!!     edenvl - density from valence electrons
!!     dmag   - density for spin-up minus density for spin-down
!!     qtotel - total charge of a cluster
!!     qnrm   - charge accumulated inside Norman sphere as result of SCF
!!     xnmues - occupation numbers of valence orbitals from SCF procedure
!!  Potential information
!!     xmu    - Fermi level position
!!     ecv    - core-valence separation energy
!!     vint   - muffin-tin zero energy (interstitial potential)
!!     vclap  - Coulomb potential
!!     vtot   - vclap + xc potential from edens
!!     vvalgs - vclap + xc potential from edenvl (EXCHANGE 5 model)
!!  Specific data for convolution with excitation spectrum (see mbconv)
!!     s02    - many body reduction factor S_0^2 
!!     erelax - estimate of relaxation energy = efrozen - emu, where
!!              efrozen is edge position estimate from Koopmans theorem
!!     wp     - estimate of plasmon frequency from rhoint
!
!!     KJ 11-2010
!!     Added option to read feff8(q) potentials.
!!     Only difference is that kappa and eorb are not present by default.
!!     I've changed feff8q to output these two to a file called forf9.dat .
!!     If this file is present, rdpot will use it.
!!     There is currently no runtime option for the user to influence this behavior.
!
!
!      logical :: feff9format=.true. !KJ
!      dimension dum(13)
!
!  10  format(a)
!   20 format (bn, i15)
!
!!KJ:
!      open(155,file='forf9.dat',form='formatted',status='old',err=132)
!	  read(155,*,end=132,err=132)
!	  read(155,*,end=132,err=132) kappa
!	  read(155,*,end=132,err=132)
!	  read(155,*,end=132,err=132) eorb
!	  close(155)
!	  feff9format=.false. !if anything went wrong, we'll try to treat the calculation as "feff9"
! 132  continue
!      if(.not.feff9format) then
!	     call wlog(':INFO FEFF9 is using a pot.bin and forf9.dat file from a FEFF8 source in routine rdpot.')
!	  endif
!!:KJ
!
!      open (unit=3, file='pot.bin', status='old')
!      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp, jumprm, iunf
!
!!     nph and npadx are not passed to calling subroutine
!      do i  = 1, ntitle
!         read(3,10) title(i)
!         call triml(title(i))
!      enddo
!!     Misc double precision stuff from pot.bin
!      call rdpadd(3, npadx, dum(1), 13)
!      rnrmav = dum(1)
!      xmu    = dum(2)
!      vint   = dum(3)
!      rhoint = dum(4)
!      emu    = dum(5)
!      s02    = dum(6)
!      erelax = dum(7)
!      wp     = dum(8)
!      ecv    = dum(9)
!      rs     = dum(10)
!      xf     = dum(11)
!      qtotel = dum(12)
!      totvol = dum(13)
!
!      read (3, 40) (imt(i),i=0,nph)
!      call rdpadd(3, npadx, rmt(0), nph+1)
!      read (3, 40) (inrm(i),i=0,nph)
!      read (3, 40) (iz(i),i=0,nph)
!      if(feff9format) read (3, 40) (kappa(i),i=1,30)  !KJ added "if"
!      call rdpadd(3, npadx, rnrm(0), nph+1)
!      call rdpadd(3, npadx, folp(0), nph+1)
!      call rdpadd(3, npadx, folpx(0), nph+1)
!      call rdpadd(3, npadx, xnatph(0), nph+1)
!      call rdpadd(3, npadx, xion(0), nph+1)
!      call rdpadd(3, npadx, dgc0(1), 251)
!      call rdpadd(3, npadx, dpc0(1), 251)
!      call rdpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
!      call rdpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
!      call rdpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
!      call rdpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
!      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
!      call rdpadd(3, npadx, xnval(1,0), 30*(nph+1) )
!      if(feff9format) call rdpadd(3, npadx, eorb(1), 30) !KJ added "if"
!      do iph=0,nph
!         read (3, 60) (iorb(i,iph),i=-4,3)
!      enddo
!      call rdpadd(3, npadx, qnrm(0), nph+1 )
!      nn = (lx+1)*(nph+1)
!      call rdpadd(3, npadx, xnmues(0,0), nn )
!
!      close (unit=3)
!	  
! 30   format(9(1x,i4))
! 40   format(20(1x,i4))
! 60   format(8(1x,i2))
!
!      return
!      end subroutine rdpot_mod
!
