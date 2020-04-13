!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdpotp.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdpotp ( vtot) 
!  opens pot.ch file and reads following information
!  Muffin-tin geometry
!     imt    - index of radial grid just below muffin-tin radii
!  Electron density information 
!     rhoint - interstitial density
!     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
!     xf     - estimate of momentum at Fermi level from rhoint
!     edens  - total electron density
!     edenvl - density from valence electrons
!     qnrm   - charge accumulated inside Norman sphere as result of SCF
!     xnmues - occupation numbers of valence orbitals from SCF procedure
!  Potential information
!     xmu    - Fermi level position
!     vint   - muffin-tin zero energy (interstitial potential)
!     vclap  - Coulomb potential
!     vtot   - vclap + xc potential from edens
      use dimsmod, only: nphx=>nphu, nheadx, lx
      implicit double precision (a-h, o-z)
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)

      character*80 title(nheadx)

      dimension dum(13)

  10  format(a)
   20 format (bn, i15)

      open (unit=3, file='pot.ch', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,     &
     &            jumprm, iunf
  30  format(9(1x,i4))
!     nph and npadx are not passed to calling subroutine
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
      read (3, 40) (kappa(i),i=1,30)
      call rdpadd(3, npadx, rnrm(0), nph+1)
      call rdpadd(3, npadx, folp(0), nph+1)
      call rdpadd(3, npadx, folpx(0), nph+1)
      call rdpadd(3, npadx, xnatph(0), nph+1)
      call rdpadd(3, npadx, xion(0), nph+1)
      call rdpadd(3, npadx, dgc0(1), 251)
      call rdpadd(3, npadx, dpc0(1), 251)
      call rdpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call rdpadd(3, npadx, eorb(1), 30)
      do 50 iph=0,nph
 50   read (3, 60) (iorb(i,iph),i=-4,3)
 60   format(8(1x,i2))
      call rdpadd(3, npadx, qnrm(0), nph+1 )
      nn = (lx+1)*(nph+1)
      call rdpadd(3, npadx, xnmues(0,0), nn )
      close (unit=3)
      imt0 = imt(0)
      rmt0 = rmt(0)
      inrm0 = inrm(0)
      rnrm0 = rnrm(0)

      return
      end
