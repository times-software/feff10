!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wrpot.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/03/27 01:29:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrpot ( nph, ntitle, title, rnrmav, xmu, vint, rhoint, &
     &                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
     &                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
     &                  dgc0, dpc0, dgc, dpc, adgc, adpc,               &
     &                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
     &                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
     &                  inters, totvol, iafolp, xion, iunf, iz, jumprm,iprint) !KJ added iprint 3-2012
!  opens pot.bin file and writes following information
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
!  Atomic orbitals info (Dirac spinors)
!     dgc0   - upper component for initial orbital
!     dpc0   - lower component for initial orbital
!     dgc    - upper components for all atomic orbitals
!     dpc    - lower components for all atomic orbitals
!     adgc   - development coefficient for upper components
!     adpc   - development coefficient for lower components
!     xnval  - number of valence electrons for each atomic orbital
!     eorb   - atomic energy of occupied orbitals
!     kappa  - quantum number kappa of occupied orbitals
!     iorb   - index of last occupied orbital for each kappa
!              used for core-valence separation and non-local exchange
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

        use DimsMod, only: nphx=>nphu, nheadx, lx

      implicit double precision (a-h, o-z)

      parameter (npadx=8)
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx) 
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)
      character*80 title(nheadx)

      real*8, intent(inout) :: xnmues(0:lx,0:nphx)
	  integer,intent(in) :: iprint

      dimension dum(13)
      character*75 wfmt

  10  format(a)
      write (wfmt, 35) nph + 1
  35  format( '(', i3, '(1x,i4))' )

      open (unit=3, file='pot.bin', status='unknown', iostat=ios)
      call chopen (ios, 'pot.bin', 'pot')
      write(3,20) ntitle, nph, npadx, nohole, ihole, inters, iafolp,    &
     &            jumprm, iunf
  20  format (9(1x,i4))
      do 133  i  = 1, ntitle
         ll = istrln(title(i))
         write(3,10) title(i)(1:ll)
  133 continue
!     Misc stuff from pot.bin
      dum(1)  = rnrmav
      dum(2)  = xmu
      dum(3)  = vint
      dum(4)  = rhoint
      dum(5)  = emu
      dum(6)  = s02
      dum(7)  = erelax
      dum(8)  = wp
      dum(9)  = ecv
      dum(10)  = rs
      dum(11)  = xf
      dum(12)  = qtotel
      dum(13)  = totvol
      call wrpadd(3, npadx, dum(1), 13)

      write (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      
      call wrpadd(3, npadx, rmt(0), nph+1)

      write (3, 40) (inrm(i),i=0,nph)
      write (3, 40) (iz(i),i=0,nph)
      write (3, 40) (kappa(i),i=1,30)

      call wrpadd(3, npadx, rnrm(0), nph+1)
      call wrpadd(3, npadx, folp(0), nph+1)
      call wrpadd(3, npadx, folpx(0), nph+1)
      call wrpadd(3, npadx, xnatph(0), nph+1)
      call wrpadd(3, npadx, xion(0), nph+1)
      call wrpadd(3, npadx, dgc0(1), 251)
      call wrpadd(3, npadx, dpc0(1), 251)
      call wrpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call wrpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call wrpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call wrpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call wrpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call wrpadd(3, npadx, eorb(1), 30)

      do 50 iph=0,nph
        write (3, 60) (iorb(i,iph),i=-4,3)
 50   continue
 60   format(8(1x,i2))
      call wrpadd(3, npadx, qnrm(0), nph+1 )
      call wrpadd(3, npadx, xnmues(0,0), (lx+1)*(nph+1) )
      close (unit=3)


      if(iprint.lt.4) return
!KJ debugging
! write the same info to a formatted file for inspection

      open (unit=3, file='pot.dat', form='formatted', status='unknown', iostat=ios)
      call chopen (ios, 'pot.dat', 'pot')
      write(3,*)  'ntitle, nph, npadx, nohole, ihole, inters, iafolp, jumprm, iunf'
      write(3,20) ntitle, nph, npadx, nohole, ihole, inters, iafolp, jumprm, iunf
!  20  format (9(1x,i4))
      do   i  = 1, ntitle
         ll = istrln(title(i))
         write(3,10) title(i)(1:ll)
      enddo
!     Misc stuff from pot.bin
	  write(3,*) 'rnrmav',rnrmav
	  write(3,*) 'xmu',xmu
	  write(3,*) 'vint',vint
	  write(3,*) 'rhoint',rhoint
      write(3,*) 'emu',emu
	  write(3,*) 's02',s02
	  write(3,*) 'erelax',erelax
	  write(3,*) 'wp',wp
	  write(3,*) 'ecv',ecv
	  write(3,*) 'rs',rs
	  write(3,*) 'xf',xf
	  write(3,*) 'qtotel',qtotel
	  write(3,*) 'totvol',totvol
      write(3,*) 'imt'
      write (3, 40) (imt(i),i=0,nph)
!  40  format(20(1x,i4))
      write(3,*) 'rmt'
	  write(3,*) rmt(0:nph)
      write(3,*) 'inrm'
      write (3, 40) inrm(0:nph)
	  write(3,*) 'iz'
      write (3, 40) iz(0:nph)
	  write(3,*) 'kappa'
      write (3, 40) kappa(1:30)
      write(3,*) 'rnrm'
      write(3,*) rnrm(0:nph)
      write(3,*) 'folp'
	  write(3,*)  folp(0:nph)
      write(3,*) 'folpx'	  
      write(3,*)  folpx(0:nph)
	  write(3,*) 'xnatph'
      write(3,*)  xnatph(0:nph)
	  write(3,*) 'xion'
      write(3,*)  xion(0:nph)
	  write(3,*) 'dgc0'
      write(3,*)  dgc0(1:251)
	  write(3,*) 'dpc0'
      write(3,*)  dpc0(1:251)
	  write(3,*) 'dgc'
      write(3,*)  dgc(1:251,1:30,0:nph)
	  write(3,*) 'dpc'
      write(3,*)  dpc(1:251,1:30,0:nph)
	  write(3,*) 'adgc'
      write(3,*)  adgc(1:10,1:30,0:nph)
	  write(3,*) 'adpc'
      write(3,*)  adpc(1:10,1:30,0:nph)
	  write(3,*) 'edens'
      write(3,*)  edens(1:251,0:nph)
	  write(3,*) 'vclap'
      write(3,*)  vclap(1:251,0:nph)
	  write(3,*) 'vtot'
      write(3,*)  vtot(1:251,0:nph)
	  write(3,*) 'edenvl'
      write(3,*)  edenvl(1:251,0:nph)
	  write(3,*) 'vvalgs'
      write(3,*)  vvalgs(1:251,0:nph)
	  write(3,*) 'dmag'
      write(3,*)  dmag(1:251,0:nph)
	  write(3,*) 'xnval'
      write(3,*)  xnval(1:30,0:nph)
	  write(3,*) 'eorb'
      write(3,*)  eorb(1:30)

      write(3,*) 'iorb'
      do iph=0,nph
        write (3, 60) iorb(-4:3,iph)
      enddo
 !60   format(8(1x,i2))
	  write(3,*) 'qnrm'
      write(3,*)  qnrm(0:nph)
	  write(3,*) 'xnmues'
      write(3,*)  xnmues(0:lx,0:nph)
      close (unit=3)


      return
      end
