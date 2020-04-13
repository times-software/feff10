!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2chi.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ff2chi (ispec, ipr4, idwopt, critcw, s02, sig2g,       &
     &                   tk, thetad, mbconv, absolu,  &!KJ added absolu 3-06
     &                   vrcorr, vicorr, alphat, thetae, iabs, nabs,    &
     &            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
!     adds the contributions from each path and absorber, including
!     Debye-Waller factors. Writes down main output: chi.dat and xmu.dat

      use dimsmod, only: npx=>npx_ff2x, nheadx, nex, nphx=>nphu, legtot 
	  use constants
      implicit double precision (a-h, o-z)

      parameter (eps4 = 1.0e-4)
      integer ipmin,ipmax,ipstep,elnes !KJ my variables  1-06    
      integer absolu !KJ 3-06  

!     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)

!KJ      parameter (npx=15000)
      !parameter (npx = 1200) !now in dimsmod
!     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

!KJ   2014      Since we don't have the same memory / speed concerns as in the 90s,
!     and since nex is now big enough to form a "fine grid",
!     let's see if this simple solution works:
!     parameter (nfinex = 601)  !old statement
      parameter (nfinex = nex)
      complex*16 cchi(nfinex), ckck(nfinex), ccc, ckp
!     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nfinex)
      complex*16 chia(nfinex)
      dimension xkp(nfinex), xk0(nfinex)

      logical dwcorr
      character*512 slog
      character*2 coment
      parameter (coment='# ')

!     Stuff from feff.bin, note that floating point numbers are
!     single precision.  Be careful throughout this routine, especially
!     when passing things to subroutines or intrinsic functions.
      real rnrmav, xmu, edge
      character*80 title(nheadx)
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
      dimension omegax(nfinex)
!#mn
      external getxk

! !KJ locals  1-06
      integer iip,nip
      logical cross 
      character*9 f1,f2
      character*10 f0,f3
      complex*16 kxsec(nex)     
!KJ end my variables      
      

!     lines below allow to skip FMS module for DANES
!     after XANES calculations
      open (unit=1, file='phase.bin', status='old', iostat=ios)
      if (ios.le.0 .and. abs(ispec).eq.3) then
        read(1,*) ne3, ne3, ne3, ne3
      endif
      close (unit=1)

! !KJ loop over iip added to process several spectra at once  1-06
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

!     open list.dat and read list of paths we want
      open (unit=1, file = f3, status='old', iostat=ios) !KJ changed 'list.dat' to f3 1-06
      call chopen (ios, f3, 'ff2chi') !KJ id.
      nhead = nheadx
      call rdhead (1, nhead, head, lhead)
!     skip a label line
      read(1,*)
      ntotal = 0
!     ip is index of path, sig2u is debye-waller from user
      do 100  i = 1, npx
         read(1,*,end=110)  ip(i), sig2u(i)
         ntotal = i
  100 continue
  110 continue
      close (unit=1)
      
      
!     read 'xsect.bin'
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,     &
     &  emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)

! !KJ I have put rdxbin inside the loop since omega is 'recycled' below,
! !KJ which is a problem if the loop executes more than once.
! !KJ Simply reading the file again and again is the lazy solution,
! !KJ but it avoids confusing changes to the code (eg., new variables).

       call rdfbin (f0, nphx, nex, npx, legtot,  &!KJ changed 'feff.bin' to f0  1-06
     &      nptot, ne, npot, ihole, iorder, ilinit,                     &
     &      rnrmav, xmu, edge, potlbl, iz, phc, ck, xk, index,          &
     &      nleg, deg, reff, crit, ipot,                                &
     &      rat, beta, eta, ri, achi, phchi)

!     make combined title
      do 120 ihead = 1, nhead
  120 title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

!     write feffnnnn.dat
      if (ipr4.ge.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne1,npot,             &
     &        ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,         &
     &        iz,phc,ck,xk,index,                                       &
     &        nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

      if (iabs.eq.1) then
!        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           if (abs(ispec).ne.3 .and. abs(del) .gt. 10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
             call wlog                                                  &
     &       (' Results may be meaningless, check input files.')
             call wlog                                                  &
     &       (' Either use XANES card or remove xsect.bin file.')
             write(slog,670)  i, xk(i)/bohr, xkxs(i)/bohr, del
             call wlog(slog)
  670        format(i7, 1p, 3e13.5)
             call par_stop('FF2CHI-1') 
           endif
  680    continue
      endif

!     If there is a vicorr, will need a mean free path factor xlam0.
!     Use it as  chi(ie) * exp (2 * reff * xlam0)
!     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 170  ipath = 1, nptot
            do 180  ie = 1, ne
               ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
               xlam0 = aimag(ck(ie)) - dimag(ckp)
               achi(ie,ipath) = achi(ie,ipath) *                        &
     &              exp (2 * reff(ipath) * xlam0)
 180        continue
 170     continue
      endif

!     Decide on fine grid.  We need two, k' evenly spaced by 
!     delk (0.05 invA) and k0 being the place in the original k 
!     grid corresponding to each k'.  k0 will not in general be on 
!     an original grid point.  Define k' by k'**2 = k**2 + vr.
!     If there is no real correction (vrcorr = 0), these two grids
!     will be the same.
!           k' is value for output, k0 is k value used for
!           interpolations with original grid.

!     vrcorr shifts the edge and the k grid
      if (abs(vrcorr) .gt. eps4)  then
         edge = edge - vrcorr
      endif

!     Find xkmin, beginning of k' grid
      delk = 0.05 * bohr
      tmp = sign (real(one), xk(1))
      e = tmp * xk(1)**2 / 2 + vrcorr
      xkpmin = getxk (e)
      n = xkpmin / delk
!     need 1st int ABOVE xkpmin/delk
      if (xkpmin .gt. 0)  n = n + 1
!     First k grid point moved by vrcorr
      xkmin = n * delk

!     Make xkp (k') and xk0 (k0) fine grids
!     ik0 is index at fermi level
      if (abs(ispec).ne.3.AND.(ispec.ne.-1)) ik0 = 1
      ik0p = 1
      do 250  i = 1, nfinex
         xkp(i) = xkmin + delk * (i - 1)
         tmp = sign (one, xkp(i))
         e = tmp * xkp(i)**2 /2 - vrcorr
         xk0(i) = getxk(e)
         if (xk0(i).lt.eps4)  ik0p = i
         if (xk0(i) .gt. xk(ne1)+eps4)  goto 260
         nkx = i
  250 continue
  260 continue
  !KJ I think at this point nkx is the size of the fine grid, and nkx <= nfinex
  !   note that it is possible that nfinex < nex ...  (i.e. the fine grid may actually be coarse ...)

  !write(*,*) 'Fine grid, nkx =',nkx, ' and ne= ',ne,' and ne1=',ne1

      dwcorr = .false.
      if (tk .gt. 1.0e-3)  dwcorr = .true.

!     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file=f1, status='unknown', iostat=ios) !KJ changed chi.dat to f1 1-06
         call chopen (ios, f1, 'ff2chi') !KJ id.
         open (unit=8, file=f2, status='unknown', iostat=ios) !KJ changed xmu.dat to f2  1-06
         call chopen (ios, f2, 'ff2chi') !KJ id.

!        write miscellaneous staff into headers  !KJ corrected typo
         call wrhead (3, ntitle, title, dwcorr, s02,                    &
     &        tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

         call wrhead (8, ntitle, title, dwcorr, s02,                    &
     &        tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

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
  380       format('    S02', f7.3, '  Temp', f8.2, '  Debye temp',f8.2,&
     &           '  Global sig2', f9.5)
            call wlog(slog)
         else
            write(slog,381) s02, sig2g
  381       format('    S02', f7.3, '  Global sig2', f9.5)
            call wlog(slog)
         endif
      endif


!     make chi and sum it
      do 400  i = 1, nfinex
         cchi(i) = 0
  400 continue
!     add Debye-Waller factors
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,&
     &  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot,rat, tk,thetad,&
     &  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, nkx, xk, xk0,&
     &  xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,                     &
     &  title, vrcorr, vicorr, nused)
!     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nfinex
            chia(ie) = 0
  635    continue
      else
         open (unit=1, file='chia.bin', status='old',                   &
     &   access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,nkx
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

!     add contribution from an absorber iabs 
!     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, nkx
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (iabs.lt.nabs) then
!        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',               &
     &   access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,nkx
  760    write(1) chia(ie)
         close(unit=1)
      endif

      if (iabs.eq.nabs) then
!        the loop over absorbers finished, ready to report results

!        Write it out
         write(3,600)  coment, nused, ntotal
         write(8,600)  coment, nused, ntotal
  600    format (a2, 1x, i4, '/', i4, ' paths used')
         write(3,610) coment
  610    format (a2, 1x, 71('-'))
         write(3,620) coment
  620    format(a2,                                                     &
     &         '      k          chi          mag           phase @#')

         do 702 ik = 1, nkx
           if (abs(ispec).ne.3) then
            rchtot(ik) = dimag (chia(ik))
           else
            rchtot(ik) = dble (chia(ik))
           endif
  702    continue
!        prepare the output grid omegax
         efermi = edge + omega(1) - dble(emxs(1))
         do 590  ik = 1, nkx
            if (xkp(ik) .lt. 0.0) then
               omegax(ik) = - xkp(ik) * xkp(ik) / 2  + efermi
            else
               omegax(ik) = xkp(ik) * xkp(ik) / 2  + efermi
            endif
  590    continue

!        do convolution with excitation spectrum
!        it is currently screwed up since xsnorm is rewritten
!        fix later
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv                                                &
     &      (omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv                                                &
     &      (omegax, nkx, efermi, s02p, erelax, wp, rchtot)
         endif


!        write to 'chi.dat'
         do 660 ik = 1, nkx
            ccc = chia(ik)
            phase = 0
            if (abs(ccc) .gt. 0)  then
               phase = atan2 (dimag(ccc), dble(ccc))
            endif
            if (ik .gt. 1)  call pijump (phase, phase0)
            phase0 = phase
            if (ipr4.ne.4) then
              write(3,630)  xkp(ik)/bohr, rchtot(ik), abs(ccc), phase0
  630         format (1x, f10.4, 3x, 3(1pe13.6,1x))
            else
!             need to report ck into chi.dat for Conradson's program
!             complex*16 should be used in terpc
!KJ    2014 BUGFIX replacing ne by ne1 below -- I guess in the old code we didn't have a branch
!         along the complex axis yet??  In any case, it "traps" terpc, therefore we need to exclude
!         those energy points here.
!  Note that I'm not sure that the "ckp" output is now necessarily correct ...
!  In fact, I think it's not.  What's needed is a proper version of terpc, in which not only ckck but also
!  xkxs can be complex.
              do i=1, ne1
                 ckck(i) = dble(real(ck(i))) + coni*dble(aimag(ck(i)))
              enddo
              call terpc (xkxs, ckck, ne1 , 3, xk0(ik), ckp)
              write(3,650)  xkp(ik)/bohr, rchtot(ik), abs(ccc), phase0, &
     &        dble(ckp)/bohr, dimag(ckp)/bohr
  650         format (1x, f10.4, 3x, 5(1pe13.6,1x))
            endif
  660    continue
         close (unit=3)
   
!        write to 'xmu.dat'
!        normalize to xsec at 50 eV above edge
!        and prepare the output energy grid omegax
         edg50 = efermi + 50 / hart
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
         if (absolu.eq.1) xsedge=dble(1) !KJ 1-06 don't normalize
         write(8,690)  coment, xsedge 
  690    format (a2, ' xsedge+50, used to normalize mu ', 1pe20.4)
         write(8,610) coment
         write(8,695) coment
  695    format (a2,' omega    e    k    mu    mu0     chi     @#')

         do i=1,nex
         if (.not.cross) then !KJ I added this block 1-06
           kxsec(i)=xsec(i)
         else
           kxsec(i)=dcmplx(0,0)
         endif !KJ end my code
         enddo
         

!        do edge correction and write down results to xmu.dat, chi.dat
         do 710 ie = 1, ne
  710    chia(ie) = 0
         if (abs(ispec).eq.3) then
!          transform from cross section in Angstrom**2 to f"/m*c**2
           do 697 ie = 1,ne
             energy = dble(emxs(ie)) + efermi
             prefac = 4 * pi * alpinv / energy * bohr**2
!            add alpha**2 to convert to units for f'
             kxsec(ie) = kxsec(ie) / prefac * alpinv**2   !KJ changed xsec to kxsec  1-06
             xsnorm(ie) = xsnorm(ie) / prefac * alpinv**2
  697      continue
           ne2 = ne - ne1 - ne3
           call fprime(efermi, emxs, ne1, ne3,ne,ik0, kxsec,xsnorm,chia,&
     &       vrcorr, vicorr, cchi)  !KJ changed xsec to kxsec 1-06
           do 850 ie=1,ne1
             omega(ie) = dble(cchi(ie))
  850      continue
         else
           call xscorr (ispec, emxs, ne1, ne, ik0, kxsec, xsnorm, chia, &
     &       vrcorr, vicorr, cchi) !KJ xsec to kxsec 7/06
!          omega is not used as energy array, but as xsec array below
           do 711 ie = 1, ne1
  711      omega(ie) = dimag(kxsec(ie)+cchi(ie))  !KJ xsec to kxsec 7/06
         endif

         do 750  ik = 1, nkx
            em0 = omegax(ik) - efermi + edge
            call terp (xkxs, omega,  ne1, 1, xk0(ik), xsec0)
            call terp (xkxs, xsnorm,  ne1, 1, xk0(ik), xsnor0)
            if (omegax(ik).ge.efermi) then
              chi0 = xsnor0 * rchtot(ik)
            else
              chi0 = xsnor0 * rchtot(ik0p)
            endif
            if (abs(ispec).ne.3) then
              write(8,700)  omegax(ik)*hart, em0*hart, xkp(ik)/bohr,    &
     &              ( chi0 + dble(xsec0) )/xsedge,                      &
     &              xsec0 /xsedge, rchtot(ik)
            else
!             signs to comply with Cromer-Liberman notation for f', f"
              write(8,700)  omegax(ik)*hart, em0*hart, xkp(ik)/bohr,    &
     &             -(xsec0+chi0), -xsec0, -chi0
            endif
  700       format (1x, 2f11.3, f8.3, 1p, 3e13.5)
  750    continue
         close (unit=8)
      endif
!     for if (iabs=abs); or the last absorber


      enddo !KJ of my iip=ipmin,ipmax,ipstep loop   1-06

      return
      end
