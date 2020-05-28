!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdxsph.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/02/04 04:55:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdxsph_h ( ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,    &
                    ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1,    &
                     Trans, InvTrans, UseTFrm, aph)
      use dimsmod, only: nex, ltot, nspx=>nspu, nphx=>nphu, lx
      use hubbard_inp
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
!      Hubbard Definitions
      complex  Trans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx),InvTrans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
      complex*16  aph(nex,lx+1,(lx+1)**2,2,0:nphx)
      logical     UseTFrm(0:lx,0:nphx)
      character*30  fname
!     Local staff
!     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')

      ! read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
      read(1,*) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx! TS added
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

! Begin of Hubbard Input
! read Hubbard phase shifts from XSPH:
open(27,file='aphase_hubbard.bin',form='unformatted',status='old')
read(27) aph
close(27)
!        do iph=0,nph
!         do iss=1,2
!         if(iss.eq.1) then
! 2011      format('aphase_up', i2.2, '.dat')
!           write(fname,2011) iph
!           open (unit=27, file=fname, status='old', iostat=ios)
!           call chopen (ios,'aphase.dat','rdxsph')
!             do  ie=1,ne
!              do  lll=0,lx
!               do  mmm=(lll**2+1),(lll+1)**2
!                 read(27,*) aa,aa,aa,aa,aph(ie,lll+1,mmm,iss,iph)
!               enddo
!              enddo
!             enddo
!           close(27)
!         elseif(iss.eq.2) then
! 2012      format('aphase_down', i2.2, '.dat')
!           write(fname,2012) iph
!           open (unit=28, file=fname, status='old', iostat=ios)
!           call chopen (ios,'aphase.dat','rdxsph')
!             do  ie=1,ne
!              do  lll=0,lx
!               do  mmm=(lll**2+1),(lll+1)**2
!                 read(28,*) aa,aa,aa,aa,aph(ie,lll+1,mmm,iss,iph)
!               enddo
!              enddo
!             enddo
!           close(28)
!         endif
!         enddo
!        enddo

! read transformation matrices to diagonalize Hubbard matrices:
open(63, file='transformation_hubbard.bin',form='unformatted',status='old')
read(63) Trans
read(63) InvTrans
close(63)
!      open(file='transf.dat',unit=63,status='old',iostat=ios)
!      open(file='Invtransf.dat',unit=64,status='old',iostat=ios)
!      write(*,*) 'Reading Transformation Matrices'
!       do is=1,2
!         do ip=0,nphx
!           do il=0,lx
!             do im1=1,2*l_hubbard+1
!              read(63,*)(Trans(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
!              read(64,*)(InvTrans(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
!             enddo
!           enddo
!         enddo
!        enddo
!       close(63)
!       close(64)

! Initialize a few control variables:
       UseTFrm(:,:)= .False.
       UseTFrm(l_hubbard,1)=.True.

      return
      end
