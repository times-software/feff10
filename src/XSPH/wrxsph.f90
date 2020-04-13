!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wrxsph.f90,v $:
! $Revision: 1.17 $
! $Author: jorissen $
! $Date: 2012/02/03 00:45:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrxsph (nsp, ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,&
     &                   ik0, em, eref, lmax, iz, potlbl, ph, rkk,kfinmax,indmax,iprint)

!KJ NOTE : updated to include feffq, 7-09.  kfinmax set to 8 by calling routine if no nrixs => identical to old routine
!          except for 2 extra fields at end first line.  Also size of temp array before packing different ; not sure if that'll
!          make a difference.

      use DimsMod, only: nphx=>nphu, ltot, nspx=>nspu, nex
      use global_inp,only: nq
      implicit double precision (a-h, o-z)
!     writes down file 'phase.bin' to be read by rphbin
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

      integer,intent(in) :: kfinmax,indmax  !KJ new feffq variables
      integer,intent(in) :: iprint !KJ 7-09 verbosity = ipr2 from xsph.inp or PRINT card
      character*6 potlbl(0:nphx)
      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex, nq, kfinmax, nspx) !KJ 12/10 added nq
      dimension lmax(0:nphx)
      dimension iz(0:nphx)

!     Local staff
!     npadx control padlib precision (see padlib package)
      parameter (npadx=8)
!     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      complex*16 temp2(nex*kfinmax*nsp) !KJ added nq ; nex*kfinmax*nqloc in withnqwithmdff
      dimension dum(3)
      character*25 innerform
      character*150 phasedatheader

      open (unit=1, file='phase.bin', status='unknown', iostat=ios)
      call chopen (ios, 'phase.bin', 'wrxsph')

      write(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0,npadx,kfinmax,indmax,nq  !KJ last 3 added for feffq
  10  format (11(1x,i4))

      dum(1) = rnrmav
      dum(2) = xmu 
      dum(3) = edge
      call wrpadd(1, npadx, dum(1), 3)

      call wrpadx(1, npadx, em(1), ne)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        temp(ii) = eref (ie, isp)
  60  continue
      call wrpadx (1, npadx, temp(1), ii)

      do 80  iph = 0, nph
         write(1, 20) lmax(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)
         do 75  isp = 1, nsp
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax(iph), lmax(iph)
               ii = ii+ 1
               temp(ii) = ph(ie, ll, isp, iph)
   70       continue
            call wrpadx (1, npadx, temp(1), ii )
   75    continue
   80 continue

!KJ old code :
!        ii = 0
!        do isp = 1, nsp
!        do kdif = 1, 8
!        do ie=1, ne
!            ii = ii + 1
!            temp(ii) = rkk (ie, kdif, isp)
!        enddo
!        enddo
!        enddo
!        call wrpadx (1, npadx, temp(1), ii)
!KJ new code from feffq :
        do iq=1, nq 
        ii = 0
        do isp = 1, nsp
        do kdif = 1, indmax
        do ie=1, ne
            ii = ii + 1
            temp2(ii) = rkk (ie, iq, kdif, isp)
        enddo
        enddo
        enddo
        call wrpadx (1, npadx, temp2(1), ii)
        enddo
        close (unit=1)

!KJ 7-09 next section writes formatted version of phase.bin in phase.dat. 
!KJ  would like to replace this with neater, cleaner iomod routine calls ... FIX LATER 
!KJ  order of variables etc. is DIFFERENT from phase.bin b/c rather than passing on arrays, I'm thinking of plotting stuff here ...
      if (iprint .gt. 0) then
          open (unit=1, file='phase.dat', status='unknown', iostat=ios)
          call chopen (ios, 'phase.dat', 'wrxsph')

          write(1,*) '# nsp, ne, ne1, ne3, nph, ihole, ik0,npadx,kfinmax,indmax'  !KJ last 2 added for feffq
          write(1,'(a,100(1x,i4))') '#',nsp, ne, ne1, ne3, nph, ihole, ik0,npadx,kfinmax,indmax  !KJ last 2 added for feffq
          write(1,'(a,3(1x,e14.7))') '# rnrmav, xmu, edge ',rnrmav,xmu,edge

          write(1,*) '# iph, lmax(iph), iz(iph), potlbl(iph)'
          do iph=0,nph
             write(1,*) '#',lmax(iph), iz(iph), potlbl(iph)
          enddo
!                                   1         2         3         4         5         6         7         8         9         0         1         2         3         4         5      
          phasedatheader= '# em(ie),(eref(ie,isp),isp=1,n),' // &
            '((rkk(ie,iq=1,kdif,isp),kdif=1,in),isp=1,n),' // &
            '(((ph(ie,ll,isp,iph),ll=-lmax(iph),lmax(iph)),isp=1,n),iph=0,nph)         '
          write(phasedatheader(30:30),'(i1)') nsp
          write(phasedatheader(64:65),'(i2)') indmax
          write(phasedatheader(74:74),'(i1)') nsp
          write(phasedatheader(129:129),'(i1)') nsp
          write(phasedatheader(138:140),'(i3)') nph
          write(1,'(a150)') phasedatheader
          do ie=1,ne
             write(1,'(100(1x,e14.7))') em(ie),(eref(ie,isp),isp=1,nsp),((rkk(ie,1,kdif,isp),kdif=1,indmax),isp=1,nsp), &
                        (((ph(ie,ll,isp,iph),ll=-lmax(iph),lmax(iph)),isp=1,nsp),iph=0,nph)
          enddo

          close (unit=1)

      endif


      return
      end
