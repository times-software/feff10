!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2afs.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ff2afs (ipr4, idwopt, critcw, s02, sig2g,              &
     &                   tk, thetad, mbconv, absolu,  &!KJ added absolu 4-06
     &                   vrcorr, vicorr, alphat, thetae, iabs, nabs,    &
     &            elnes,ipmin,ipmax,ipstep)   !KJ added this line  4-06
!     calculate anomalous scattering amplitude for a given edge
!     Writes down main output: chi.dat and xmu.dat
      use dimsmod, only: npx=>npx_ff2x, nheadx, nex, nphx=>nphu, legtot
	  use constants
      implicit double precision (a-h, o-z)

      parameter (eps4 = 1.0e-4)
      integer ipmin,ipmax,ipstep,elnes !KJ my variables  4-06    
      integer absolu !KJ 4-06  


!     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)
      complex*16 gtr(nex),gtrful(ipmin:ipmax,nex) !KJ added gtrful 4-06

!KJ      parameter (npx=15000)
      !parameter (npx = 1200) !now in dimsmod
!     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      complex*16 cchi(nex), ckp
!     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nex), xkp(nex)
      complex*16 chia(nex)

      logical dwcorr
      character*512 slog
      character*2 coment
      parameter (coment='# ')

!     Stuff from feff.bin, note that floating point numbers are
!     single precision.  Be careful throughout this routine, especially
!     when passing things to subroutines or intrinsic functions.
      real rnrmav, xmu, edge
      character*80 title(nheadx), titfms
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
!     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx)
      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)

!     stuff from xsect.bin
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex), fpp(nex)
      
! !KJ locals  1-06
      integer iip,nip
      logical cross 
      character*9 f1,f2
      character*10 f0,f3
      complex*16 gtrtemp(nex*(1+ipmax-ipmin))
      complex*16 kxsec(nex)
      real*8 del,delta_relative
! !KJ end my variables




!     get gtr - result of FMS
      do 112 ie =1,nex
      gtr(ie)= 0
      do 112 iip =ipmin,ipmax  !KJ I added this variable ip 1-06
  112 gtrful(iip,ie) = 0
      ntfms = 0
      nip=ipmax-ipmin+1 !KJ 1-06
      
      open (unit=1, file='fms.bin', status='old', iostat=ios)
      if (ios.le.0) then
         ntfms = 1
         read(1, 113) titfms
  113    format(a)
         read(1, 115) ne, ne1, ne3, nph, npadx
  115    format(5(1x,i3))
         call rdpadx(1, npadx, gtrtemp, ne*nip)  !KJ I added *nip, changed gtr to gtrtemp  1-06
      endif
      close (unit=1)
!KJ Next lines my addition to read several spectra at once. 1-06
        i=1
        do iip=ipmin,ipmax
        do j=1,ne
        gtrful(iip,j)=gtrtemp(i+j-1)
        enddo
          i=i+ne
        enddo
!KJ Now we don't need gtrtemp anymore.  End my changes



! !KJ loop over iip added to process several spectra at once  4-06
! !KJ reading of feff.bin and list.dat moved inside the loop (used to be before reading
! !KJ xsect.bin      
      do iip=ipmin,ipmax,ipstep
        cross=(.not.(iip.eq.1.or.iip.eq.10.or.iip.eq.5.or.iip.eq.9))
      
! !KJ choose different filename for each spectrum.
        if(iip.eq.1) then
          f1(1:9)='chi.dat  '
          f2(1:9)='xmu.dat  '
          f0(1:10)='feff.bin  '
          f3(1:10)='list.dat  ' 	  
        elseif(iip.eq.10) then
          f1(1:9)='chi10.dat'
          f2(1:9)='xmu10.dat'
          f0(1:10)='feff10.bin'
          f3(1:10)='list10.dat'	  
        elseif(iip.gt.1.and.iip.lt.10) then
          f1(1:4)='chi0'
          f1(5:5)= char(48+iip)
          f1(6:9)='.dat'
          f2(1:4)='xmu0'
          f2(5:5)= char(48+iip)
          f2(6:9)='.dat'
          f0(1:5)='feff0'
          f0(6:6)= char(48+iip)
          f0(7:10)='.bin'	
          f3(1:5)='list0'
          f3(6:6)= char(48+iip)
          f3(7:10)='.dat'
        else
          stop 'crazy iip in ff2xmu'
        endif
        do i=1,nex
        gtr(i)=gtrful(iip,i)
        enddo
!KJ


!     open list.dat and read list of paths we want
      open (unit=1, file= f3, status='old', iostat=ios)!KJ changed 'list.dat' to f3 1-06
      ntotal = 0
      if (ios.le.0) then
        call chopen (ios, f3, 'ff2afs')  !KJ id.
        nhead = nheadx
        call rdhead (1, nhead, head, lhead)
!       skip a label line
        read(1,*)
!       ip is index of path, sig2u is debye-waller from user
        do 100  i = 1, npx
           read(1,*,end=110)  ip(i), sig2u(i)
           ntotal = i
  100   continue
  110   continue
      endif
      close (unit=1)


!     lines below allow to skip FMS module for DANES 
!     after XANES calculations
      open (unit=1, file='phase.bin', status='old', iostat=ios)
      if (ios.le.0) then
        read(1,*) ne3, ne3, ne3, ne3
      endif
      close (unit=1)

       call rdfbin (f0, nphx, nex, npx, legtot,   &!KJ changed 'feff.bin' to f0  4-06
     &     nptot, ne, npot, ihole, iorder, ilinit,                      &
     &     rnrmav, xmu, edge, potlbl, iz, phc, ck, xk,                  &
     &     index, nleg, deg, reff,                                      &
     &     crit, ipot, rat, beta, eta, ri, achi, phchi)
!     read xsect.dat file
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,     &
     &  emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)
     !c !!KJ these comments copied from ff2chi and not necessarily relevant  4-06.
! !KJ I have put rdxbin inside the loop since omega is 'recycled' below,
! !KJ which is a problem if the loop executes more than once.
! !KJ Simply reading the file again and again is the lazy solution,
! !KJ but it avoids confusing changes to the code (eg., new variables).



!     make combined title
      if (ntfms.eq.1) then
        ntitle = ntitle + 1
        title(ntitle) = titfms
      endif
      do 120 ihead = 1, nhead
  120 title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

!     write feffnnnn.dat
      if (ipr4.ge.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne,npot,              &
     &        ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,         &
     &        iz,phc,ck,xk,index,                                       &
     &        nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

!     If there is a vicorr, will need a mean free path factor xlam0.
!     Use it as  chi(ie) * exp (2 * reff * xlam0)
!     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 180  ie = 1, ne
            ckp = sqrt (ck(ie)**2 + 2*coni*vicorr)
            xlam0 = aimag(ck(ie)) - dimag(ckp)
            do 170  ipath = 1, nptot
               achi(ie,ipath) = achi(ie,ipath) *                        &
     &               exp (2 * reff(ipath) * xlam0)
  170       continue
  180    continue
      endif

!     k'**2 = k**2 + vr. If there is no real correction
!     (vrcorr = 0), these two grids will be the same.
!           k' is value for output,  k is  value used for
!           interpolations with original grid.

!     vrcorr shifts the edge and the k grid
      if (abs(vrcorr) .gt. eps4)  then
         edge = edge - vrcorr
      endif

!     ik0 is index at fermi level
      do 250  i = 1, ne
         temp = xk(i)*abs(xk(i)) + 2*vrcorr
         if (temp.ge. 0) then
           xkp(i) = sqrt(temp)
         else
           xkp(i) = - sqrt(-temp)
         endif
  250 continue
     

      dwcorr = .false.
      if (tk .gt. 1.0e-3)  dwcorr = .true.

!     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file=f1, status='unknown', iostat=ios) !KJ changed chi.dat to f1
         call chopen (ios, f1, 'ff2afs') !KJ id.  4-06
         open (unit=8, file=f2, status='unknown', iostat=ios) !KJ changed xmu.dat to f2
         call chopen (ios, f2, 'ff2afs') !KJ id.

!        write miscellaneous staff into headers  !KJ corrected typo
         call wrhead (8, ntitle, title, dwcorr, s02,                    &
     &     tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

!        also write information on the screen
         if (alphat .gt. zero)  then
            write(slog,322) alphat
  322       format ('    1st and 3rd cumulants, alphat = ', 1pe20.4)
            call wlog(slog)
         endif
         if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
            write(slog,343) vrcorr*hart, vicorr*hart
  343       format ('    Energy zero shift, vr, vi ', 1p, 2e14.5)
            call wlog(slog)
         endif

         write(slog,370) critcw
         call wlog(slog)
  370    format ('    Use all paths with cw amplitude ratio', f7.2, '%')
         if (dwcorr)  then
            write(slog,380) s02, tk, thetad, sig2g
            call wlog(slog)
         else
            write(slog,381) s02, sig2g
            call wlog(slog)
         endif
  380    format('    S02', f7.3, '  Temp', f8.2, '  Debye temp', f8.2,  &
     &           '  Global sig2', f9.5)
  381    format('    S02', f7.3, '  Global sig2', f9.5)
      endif


!     make chi and sum it
      do 400  i = 1, nex
         cchi(i) = 0
  400 continue
      do 402  ik = 1, ne
         cchi(ik)= s02 * gtr(ik)
  402 continue
!     add correction due to vicorr
      if (vicorr.gt.eps4) then
         call conv(omega,cchi,ne1,vicorr)
!        call conv(omega,xsec,ne1,vicorr)
      endif


!     add Debye-Waller factors
      ispec = 3
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,&
     &  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,&
     &  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, ne, xk, xkp, &
     &  xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,                     &
     &  title, vrcorr, vicorr,  nused)

!     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nex
            chia(ie) = 0
  635    continue
      else
         open (unit=1, file='chia.bin', status='old',                   &
     &   access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,ne
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

      if(iabs.eq.1) then
!        compare grids in xsect.dat and feff.bin
         do i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           delta_relative = (xk(i)**2 - xkxs(i)**2) / (xk(i)**2 + xkxs(i)**2)
           !KJ Jan 2014 added delta_relative test
           !KJ for DANES or FPRIME calculations in which we go to really high energies, e.g. k=200.0
           ! At this point, errors larger than 0.0001 occur - but they are maybe in the 7th or 8th significant digit.
           ! Thus a better test is needed.
           ! Both xk and xkxs originate from the same energy array "em" in module xsph, but they get propagated differently.
           ! Small errors occur on the one hand due to some variables retardedly being stored in single precision for nonexistent performance benefits,
           ! as well as limiting the number of digits in e.g. xsect.dat.
           ! The main source of discrepancy is likely the conversion from energy=k^2 to k=sqrt(energy) then back to energy=k^2.
           if (abs(del) .gt.  10*eps4 .and. delta_relative .gt. 1.d-6 )  then
             call wlog(' Emesh in feff.bin and xsect.dat different.')
!            stop
           endif
         enddo
      endif

!     add contribution from an absorber iabs 
!     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, ne
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (iabs.lt.nabs) then
!        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',               &
     &   access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,ne
  760    write(1) chia(ie)
         close(unit=1)
      endif

      if (iabs.eq.nabs) then
!        The loop over absorbers is finished. Write out the results.
         write(8,600) coment, nused, ntotal
  600    format (a2,1x, i4, '/', i4, ' paths used')
  610    format (a2,1x, 71('-'))

         do 702 ik = 1, ne
            rchtot(ik) = dimag (chia(ik))
  702    continue
!        prepare the output grid omega
         efermi = edge + omega(1) - dble(emxs(1))

!        do convolution with excitation spectrum
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv                                                &
     &      (omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv                                                &
     &      (omega, ne1, efermi, s02p, erelax, wp, rchtot)
         endif

!        normalize to xsec at 50 eV above edge
!        and prepare the output energy grid omega
         edg50 = efermi + 50 / hart
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
         if (absolu.eq.1) xsedge=dble(1)  !KJ 4-06 don't normalize	 
         write(8,660)  coment, xsedge 
  660    format (a2,' xsedge+ 50, used to normalize mu ', 1pe20.4)
         write(8,610) coment
         write(8,665) coment
  665    format (a2,' omega    e    k    mu    mu0     chi     @#')

!        transform from cross section in Angstrom**2 to f"/m*c**2
         do 670 ie = 1,ne
            energy = dble(emxs(ie)) + efermi 
            prefac = 4 * pi * alpinv / energy * bohr**2
!           add alpha**2 to convert to units for f'
            xsec(ie) = xsec(ie) / prefac * alpinv**2
            xsnorm(ie) = xsnorm(ie) / prefac * alpinv**2
  670    continue

         do i=1,nex
         if (.not.cross) then !KJ I added this block 4-06
           kxsec(i)=xsec(i)
         else
           kxsec(i)=dcmplx(0,0)
         endif !KJ end my code
         enddo


!        do correction using brouder method
         ne2 = ne - ne1 - ne3
         call fprime(efermi, emxs, ne1, ne3, ne, ik0, kxsec,xsnorm,chia, &!KJ changed xsec to kxsec 4-06
     &       vrcorr, vicorr, cchi)
         do 850 ie=1,ne1
           fpp(ie)=xsnorm(ie) + dimag(xsnorm(ie)*chia(ie))
           rchtot(ie)=dble(xsnorm(ie)*chia(ie)+cchi(ie))
  850    continue
         do 855 ie=1,ne
           chia(ie) = 0
  855    continue
         call fprime(efermi, emxs, ne1, ne3, ne, ik0, kxsec,xsnorm,chia, &!KJ id.
     &       vrcorr, vicorr, cchi)
         do 860 ie = 1, ne1
            em0 = dble(emxs(ie))
            xsec0 = dble( cchi(ie))
            chi0  = (rchtot(ie) - xsec0)
            if (ne2.gt.0) then
!             DANES
!             - signs to comply with Cromer-Liberman notation for f', f"
              write(8,700)  omega(ie)*hart, em0*hart, xkp(ie)/bohr,     &
     &             -rchtot(ie), -xsec0, -chi0
  700         format (1x, 2f11.3, f8.3, 1p, 3e13.5)
            else
!             FPRIME
              write(8,710)  omega(ie)*hart, em0*hart,                   &
     &              -rchtot(ie), -xsec0, fpp(ie), xsnorm(ie)
  710         format (1x, 2f11.3, 4e13.5)
            endif
  860    continue

         close (unit=8)
         close (unit=3, status='delete')
      endif
!     for if (iabs=abs); or the last absorber

      enddo !KJ of my iip=ipmin,ipmax,ipstep loop  1-06


      return
      end
