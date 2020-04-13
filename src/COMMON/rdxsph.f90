!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdxsph.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/02/04 04:55:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdxsph ( ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,    &
     &               ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      use dimsmod, only: nex, ltot, nphx=>nphu, nspx=>nspu
      implicit double precision (a-h, o-z)
!     reads file 'phase.bin' 
!  Energy grid information
!     em   - complex energy grid
!     eref - V_int + i*gamach/2 + self-energy correction
!     ne   - total number of points in complex energy grid
!     ne1  - number of points on main horizontal axis
!     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
!     ne3  - number of points on auxiliary horizontal axis (need for f')
!     xmu  - Fermi energy
!     edge - x-ray frequency for final state at Fermi level
!     ik0  - grid point index at Fermi level
!  Potential type information
!     nph - number of potential types
!     iz  - charge of nuclei (atomic number)
!     potlbl - label for each potential type
!     lmax - max orb momentum for each potential type
!     ihole - index of core-hole orbital for absorber (iph=0)
!     rnrmav - average Norman radius (used in headers only)
!  Main output of xsect and phases module (except that in xsect.bin)
!     ph  - complex scattering phase shifts
!     rkk - complex multipole matrix elements

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
!     kinit, linit, ilinit,  - initial state kappa and ang. mom.
!     lmaxp1  -largest lmax in problem + 1

!     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

!     Local staff
!     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
!     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp 
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp 
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

!     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
!        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.            &
     &          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end
