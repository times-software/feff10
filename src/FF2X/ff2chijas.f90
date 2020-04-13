!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2chijas.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2011/12/11 01:11:14 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ff2chijas(ispec, ipr4, idwopt, critcw, s02, sig2g,tk, thetad, mbconv, &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs,kfinmax,xivec,ldecmx,ldecs)
!     adds the contributions from each path and absorber, including
!     Debye-Waller factors. Writes down main output: chi.dat and xmu.dat
      use dimsmod, only: npx=>npx_ff2x, nheadx, nex, nphx=>nphu, legtot 
	  use constants
      implicit double precision (a-h, o-z)

      parameter (eps4 = 1.0e-4)
      parameter (npadx=8)

!     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)

      complex*16 gtrl(0:ldecs,0:ldecs,nex)
      !parameter (npx = 1200) !now in dimsmod
!     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

!KJ      parameter (nfinex = 601) !KJ 12-2011 see remarks in subr dwadd and ff2afsjas
	  parameter (nfinex = nex)
      complex*16 cchi(nfinex), ckck(nfinex), ccc, ckp
      integer ldecmx,ldecs
      complex*16 lcchi(nfinex,0:ldecs,0:ldecs) 
!     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nfinex),rchtotl(nfinex,0:ldecs,0:ldecs)
      double precision rch0
      complex*16 chia(nfinex)
      complex*16 lchia(nfinex,0:ldecs,0:ldecs)
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
      real eorg(nex) 
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
	  real,allocatable :: rat(:,:,:), beta(:,:), eta(:,:), ri(:,:), achi(:,:), phchi(:,:), lgachi(:,:,:,:), lgphchi(:,:,:,:)  !KJ see comments in ff2xmujas
      double precision om2x0
 


      double precision lgom2x0(0:ldecs,0:ldecs)
      double precision lgxsec0(0:ldecs,0:ldecs)
      double precision lgomega(nex,0:ldecs,0:ldecs)
      double precision lgom2x(nex,0:ldecs,0:ldecs)
      double precision lgchi0(0:ldecs,0:ldecs) 
!
!
!
!
      integer ii,ilm1,ilm2,id1,id2
      complex*16 xsecl(kfinmax,nex),lgxsec(nex,0:ldecs)
      integer  kiind(kfinmax), lgind(kfinmax),ljind(kfinmax)
      integer lind(kfinmax)

!
!     Regular deltaE
!
      double precision de,emin,emax,e0,ep

!     stuff from xsect.bin
      integer lg1,lg2
      complex*16 emxs(nex), xsec(nex)
      complex*16 xseczero(nex)
      dimension omega(nex),om2x(nex), xkxs(nex), xsnorm(nex)
      dimension omegax(nfinex)
!#mn
      external getxk
!
!
!
      double precision xivec(3)
!
!     Debugging
!
      
       complex, allocatable :: ss(:,:),ss2(:,:) !KJ
       character*25 innerform
       integer ltotch

!KJ:
      allocate( rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx), &
                ri(legtot,npx), achi(nex,npx), phchi(nex,npx), &
                lgachi(nex,0:ldecs,0:ldecs,npx), lgphchi(nex,0:ldecs,0:ldecs,npx) )
	  allocate( ss(nex,npx),ss2(nex,npx))



!      write(6,*) "Here we are ",kfinmax
!     open list.dat and read list of paths we want
      open (unit=1, file='list.dat', status='old', iostat=ios)
      call chopen (ios, 'list.dat', 'ff2chi')
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

       call rdfbinl ('feff.bin', nphx, nex, npx, legtot, &
            nptot, ne, npot, ihole, iorder, ilinit,  &
            rnrmav, xmu, edge, potlbl, iz, phc, ck, xk, index, &
            nleg, deg, reff, crit, ipot, &
            rat, beta, eta, ri, achi, phchi,ldecmx,ldecs, &
            lgachi,lgphchi)


!
!
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0, &
        emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)
!
!     read 'xsecl.bin'
!

      if (ldecmx.ge.0) then 
         open (unit=1,file='xsecl.bin',status='unknown',iostat = ios)
         call chopen(ios,'xsecl.bin','ff2xmuq')
         read(1,*) id1,id2
!     write(6,*) id1,id2
         do ii=1,id2
            read(1,'(4I5)')  kiind(ii),lgind(ii),ljind(ii),lind(ii)
         end do
         do ie=1,nex
            call rdpadx (1, npadx, xsecl(1,ie),kfinmax )
         end do
         close(1)
         do lg1=0,ldecmx
            do ie=1,nex
               lgxsec(ie,lg1)=0.0d0
            end do
            do ii=1,id2
               if (lgind(ii).eq.lg1) then
                  do ie=1,nex
                     lgxsec(ie,lg1)=lgxsec(ie,lg1)+xsecl(ii,ie)
                  end do
               end if
            end do
         end do
      end if
      
!      open(unit=99,file='chisect1.dat',form='formatted',
!     1     status='unknown')(
!      rewind(99)
!      do ie=1,ne
!         write(99,*) dble(emxs(ie)),dble(xsec(ie)),dimag(xsec(ie))
!      end do
!      close(99)

!     make combined title
      do 120 ihead = 1, nhead
  120 title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

!     write feffnnnn.dat
      if (ipr4.ge.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne1,npot, &
              ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl, &
              iz,phc,ck,xk,index, &
              nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

      if (iabs.eq.1) then
!        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
! JAS commenting but leaving it be
!            if (i.gt.1) then 
!              write(6,'(I5,4f12.6)') i,xk(i)-xk(i-1),
!     1             xkxs(i)-xkxs(i-1),
!     1             xk(i)*xk(i)-xk(i-1)*xk(i-1),    
!     1             xkxs(i)*xkxs(i)-xkxs(i-1)*xkxs(i-1)    
!           end if
           if (abs(del) .gt. 10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
             call wlog(' Results may be meaningless, check input files.')
             call wlog(' Either use XANES card or remove xsect.bin file.')
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
!      write(6,*) "Doing vicorr"
!         write(6,*)
         do 170  ipath = 1, nptot
            do 180  ie = 1, ne
               ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
               xlam0 = aimag(ck(ie)) - dimag(ckp)
               achi(ie,ipath) = achi(ie,ipath) * exp (2 * reff(ipath) * xlam0)
 180        continue
            if (ldecmx.ge.0) then 
               do ilm1=0,ldecmx
                  do ilm2=0,ldecmx
                     do ie=1,ne
                        ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
                        xlam0 = aimag(ck(ie)) - dimag(ckp)
                        lgachi(ie,ilm2,ilm1,ipath) = lgachi(ilm2,ilm1,ie,ipath) * exp (2 * reff(ipath) * xlam0)
                     end do
                  end do
               end do
            end if

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
!     JAS modified, instead of constant delta k use
!     constant delta E if this was used earlier, 
!     and go all the way up to xk(ne1)*xk(ne1)/2.0d0
!
      if (abs(xk(ne1)+xk(ne1-2)-2.0d0*xk(ne1-1)).lt.eps4) then
!
!     Here do constant delk grid as before
!
!     Find xkmin, beginning of k' grid
         write(6,*) "   Constant delta k grid"
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
!
         
         
         ik0 = 1
         ik0p = 1
         do 250  i = 1, nfinex
            xkp(i) = xkmin + delk * (i - 1)
            tmp = sign (one, xkp(i))
            e = tmp * xkp(i)**2 /2 - vrcorr
            xk0(i) = getxk(e)
            if (xk0(i).lt.eps4)  ik0p = i
            if (xk0(i) .gt. xk(ne1)+eps4)  goto 260
            nkx = i
 250     continue
 260     continue
      else
!         
!     Do constant delE grid 
!
         write(6,*) "   Constant delta E grid"
         do i=1,ne1
             if (xk(i).gt.0.0d0) then 
               eorg(i)=xk(i)*xk(i)/2.0d0
            else
               eorg(i)=-xk(i)*xk(i)/2.0d0
            end if
         end do
         de=(eorg(ne1)-eorg(1))/dble(nfinex)
         emin=eorg(1)+de/4.0d0
         do i=1,nfinex
            e0=emin+de*dble(i-1)
            xk0(i)=getxk(e0)
            ep=e0+vrcorr
            xkp(i)=getxk(ep)
         end do
         nkx=nfinex
      end if
      dwcorr = .false.
!      write(6,*) "chi this chi that"
      if (tk .gt. 1.0e-3)  dwcorr = .true.

!     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file='chi.dat', status='unknown', iostat=ios)
         call chopen (ios, 'chi.dat', 'ff2chi')
         open (unit=8, file='xmu.dat', status='unknown', iostat=ios)
         call chopen (ios, 'xmu.dat', 'ff2chi')

!        write miscellanious staff into headers
         call wrhead (3, ntitle, title, dwcorr, s02, tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw ) !KJ lazy !,xivec)

         call wrhead (8, ntitle, title, dwcorr, s02, tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw ) !KJ lazy !,xivec)

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
  380       format('    S02', f7.3, '  Temp', f8.2, '  Debye temp',f8.2, '  Global sig2', f9.5)
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
      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do 401  i = 1, nfinex
                  lcchi(i,lg2,lg1)=0   
 401           continue
               
            end do
         end do
      end if


      
!      do i=1,npx
!         do j=1,nex
!            ss2(j,i)=cmplx(0.0d0,0.0d0)
!            ss(j,i)=achi(j,i)*exp(coni*phchi(j,i))
!            do lg1=0,2
!               do lg2=0,2
!            ss2(j,i)=ss2(j,i)+lgachi(j,lg2,lg1,i)
!     1                 *exp(coni*lgphchi(j,lg2,lg1,i))
!               end do
!            end do
!         end do
!      end do
!       do i=1,npx
!          do j=1,nex
!             write(68,'(2I5,4e15.6)') i,j,aimag(ss(j,i)),aimag(ss2(j,i))
!     1            ,real(ss(j,i)),real(ss2(j,i))
!             if (abs(ss(j,i)).gt.0.0d0) then
!             write(6,*) i,j,(abs(ss2(j,i))-abs(ss(j,i)))/abs(ss(j,i))
!     1            ,abs(ss(j,i))
!          end if
!          end do
!       end do

!      stop
!     add Debye-Waller factors
      call dwaddl (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g, &
        sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot,rat, tk,thetad, &
        alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, ldecmx,ldecs,lgachi,lgphchi,  &
           nkx, xk, xk0,xkp, cchi,lcchi, iabs, nabs, ispec, ipr4, ntitle,title, vrcorr, vicorr, nused)

!
!     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nfinex
            chia(ie) = 0
  635    continue
         if (ldecmx.ge.0) then 
            do ie=1,nfinex
               do lg1=0,ldecmx
                  do lg2=0,ldecmx
                     lchia(ie,lg2,lg1) = 0
                  end do
               end do
            end do
         end if
       else
         open (unit=1, file='chia.bin', status='old', access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,nkx
  640    read(1) chia(ie)
!         write(6,*) "from read",ie,chia(ie) 
         close (unit=1, status='delete')
      endif

!     add contribution from an absorber iabs 
!     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, nkx
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ik=1,nkx
                  lchia(ik,lg2,lg1)   = lchia(ik,lg2,lg1) + lcchi(ik,lg2,lg1)/ nabs
               end do
            end do
         end do
      end if
      if (iabs.lt.nabs) then
!        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown', access='sequential', form='unformatted', iostat=ios)
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
  620    format(a2,'      k          chi          mag           phase @#')

         do 702 ik = 1, nkx
            rchtot(ik) = dimag (chia(ik))
!            write(6,*) "Here",ik,rchtot(ik)
  702    continue
         if (ldecmx.ge.0) then 
            do lg1=0,ldecmx
               do lg2=0,ldecmx
                  do ik=1,nkx
                     rchtotl(ik,lg2,lg1) = dimag (lchia(ik,lg2,lg1))
!     write(6,*) "rchtot", ik,lg2,lg1, rchtotl(ik,lg2,lg1)
                  end do
               end do
            end do
         end if
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
            call  exconv(omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv(omegax, nkx, efermi, s02p, erelax, wp, rchtot)
            if (ldecmx.ge.0) then 
               do lg1=0,ldecmx
                  do lg2=0,ldecmx
                     call  exconv(omegax, nkx, efermi, s02p, erelax, wp,  &
                          rchtotl(1,lg2,lg1))
                  end do
               end do
            end if
         endif


!        write to 'chi.dat', but not for the decomposition
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
!KJ 2014 ne => ne1 ; see comments in ff2chi.f90
              do 625 i=1,ne1
  625         ckck(i) = dble(real(ck(i))) +coni*dble(aimag(ck(i)))
              call terpc (xkxs, ckck, ne1, 3, xk0(ik), ckp)
              write(3,650)  xkp(ik)/bohr, rchtot(ik), abs(ccc), phase0, &
              dble(ckp)/bohr, dimag(ckp)/bohr
  650         format (1x, f10.4, 3x, 5(1pe13.6,1x))
            endif
  660    continue
         close (unit=3)
   
!        write to 'xmu.dat'
!        normalize to xsec at 50 eV above edge
!        and prepare the output energy grid omegax
         edg50 = efermi + 50 / hart
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
!         write(8,690)  coment, xsedge 
!  690    format (a2, ' xsedge+50, used to normalize mu ', 1pe20.4)
         write(8,690)  coment 
  690    format (a2, ' Contribution to S(q,w) from a single electron')
         write(8,610) coment
         write(8,695) coment
!  695    format (a2,' omega    e    k    mu    mu0     chi     @#')
 695     format (a2,' omega    e    k   S(qw)  S^0(qw)  chi_q*S^0(qw) & 
        &      @#')

!        do edge correction and write down results to xmu.dat, chi.dat
         do 710 ie = 1, ne
  710    chia(ie) = 0 
         do ie=1,nex
            xseczero(ie)=0
         end do

         call xscorr(ispec, emxs, ne1, ne, ik0, xsec, xsnorm, chia, vrcorr, vicorr, cchi)
!        omega is not used as energy array, but as xsec array below
         do 711 ie = 1, ne1
  711    omega(ie) = dimag(xsec(ie)+cchi(ie))
         do 741 ie = 1, ne1
  741    om2x(ie) = dimag(xsec(ie))
         
         do 750  ik = 1, nkx
            em0 = omegax(ik) - efermi + edge
            call terp (xkxs, omega,  ne1, 1, xk0(ik), xsec0)
            call terp (xkxs, om2x,  ne1, 1, xk0(ik), om2x0)
            call terp (xkxs, xsnorm,  ne1, 1, xk0(ik), xsnor0)
            if (omegax(ik).ge.efermi) then
              chi0 = xsnor0 * rchtot(ik)
            else
              chi0 = xsnor0 * rchtot(ik0p)
            endif 

            write(8,700)  omegax(ik)*hart, em0*hart, xkp(ik)/bohr, ( chi0 + dble(xsec0) )/hart, &
                    xsec0 /hart, chi0/hart, om2x0/hart
!     1              ( chi0 + dble(xsec0) )/xsedge,
!     1              xsec0 /xsedge, rchtot(ik)

  700       format (1x, 2f11.3, f8.3, 1p, 4e13.5)
  750    continue
         close (unit=8)
!
!     Now do it for the decompostions 
!

         if (ldecmx.ge.0) then 
            open (unit=8, file='xmul.dat', status='unknown', iostat=ios)
            call chopen (ios, 'xmul.dat', 'ff2chi')
            write(8,890)  coment 
 890        format (a2, ' Decomposition of S(q,w) for a single electron')
            write(8,895) coment
!     695    format (a2,' omega    e    k    mu    mu0     chi     @#')
 895        format (a2,' omega    k   S^0(qw)  S_{l=0,...,ldecmx}^0(qw)  &
     &           chi^q_{l=0,..ldecmx,l^*=0,...,ldecmx} ')
            write(8,898) coment,ldecmx
 898     format (a2,'and ldecmx= ', i5)        
            
            do lg1=0,ldecmx
               do lg2=0,ldecmx
                  do ie = 1, ne
                     lchia(ie,lg2,lg1) = 0
                  end do 
                  if (lg1.eq.lg2) then 
                  call xscorr(ispec, emxs, ne1, ne, ik0, lgxsec(1,lg2), xsnorm, lchia(1,lg2,lg1), &
                      vrcorr, vicorr, lcchi(1,lg2,lg1))
                  
                  do ie = 1, ne1
!     write(6,*) "Sampsa",  lcchi(ie,lg2,lg1), 
!     1                    lchia(1,lg2,lg1)
                     lgomega(ie,lg2,lg1) = dimag(lgxsec(ie,lg2)+lcchi(ie,lg2,lg1))
!     write(6,*) "Sampsa", ie,lg2,lg1,lgomega(ie,lg2,lg1)
                  end do
                  do ie = 1, ne1
                     lgom2x(ie,lg2,lg1) = dimag(lgxsec(ie,lg2))
                  end do               
                  
               else 
                  call xscorr(ispec, emxs, ne1, ne, ik0, xseczero, xsnorm, lchia(1,lg2,lg1), &
                      vrcorr, vicorr, lcchi(1,lg2,lg1))
                  do ie = 1, ne1
                     lgomega(ie,lg2,lg1) = dimag(lcchi(ie,lg2,lg1))
!     write(6,*) "Sampsa2", ie,lg2,lg1,
!     1                    lgomega(ie,lg2,lg1)
                  end  do
                  do ie = 1, ne1
                     lgom2x(ie,lg2,lg1) = dimag(xseczero(ie))
                  end do               
               end if
               
               
               
            end do              ! lg2
            
         end do                 ! lg1
         
!
!     Now loop over energy and write it out
!
         
          do ik = 1, nkx
             em0 = omegax(ik) - efermi + edge
             call terp (xkxs, xsnorm,  ne1, 1, xk0(ik), xsnor0)
             call terp (xkxs, omega,  ne1, 1, xk0(ik), xsec0)
             call terp (xkxs, om2x,  ne1, 1, xk0(ik), om2x0)
             call terp (xkxs, xsnorm,  ne1, 1, xk0(ik), xsnor0)
            if (omegax(ik).ge.efermi) then
              chi0 = xsnor0 * rchtot(ik)
            else
              chi0 = xsnor0 * rchtot(ik0p)
            endif 

             do lg1=0,ldecmx
                do lg2=0,ldecmx
                   call terp (xkxs, lgomega(1, lg2,lg1),  ne1, 1, xk0(ik),lgxsec0(lg2,lg1))
!                   write(6,*) "Jaj", lgxsec0(lg2,lg1)
                   call terp (xkxs, lgom2x(1,lg2,lg1),  ne1, 1, xk0(ik),lgom2x0(lg2,lg1))
!                   write(6,*) "Juj", lgom2x0(lg2,lg1)
                   if (omegax(ik).ge.efermi) then
                      lgchi0(lg2,lg1) = xsnor0 * rchtotl(ik,lg2,lg1)
                   else
                      lgchi0(lg2,lg1) = xsnor0 * rchtotl(ik0p,lg2,lg1)
                   endif 
                end do ! lg2
             end do !lg1 
!
!     Now we have something to print, so print it out. 
!
!             c
!     total number of channels 
!     1. atomic background: one column
!     2. l-channel contr. to atomic background ldecmx+1 channels
!     3. fine structure (ldecmx+1)^2
!
         
         ltotch=ldecmx+2+(ldecmx+1)*(ldecmx+1)
         if (ltotch.le. 9) then 
            write(innerform,'(''(1x,2f11.3,'',I1,''e11.3)'')') ltotch
         elseif (ltotch.le. 99) then
            write(innerform,'(''(1x,2f11.3,'',I2,''e11.3)'')') ltotch
         else
            call par_stop('Too many shannels in ff2xmu') 
         end if 
                   
!             write(8,'(1x,2f11.3,9e11.3)') omegax(ik)*hart,xkp(ik)/bohr,
!     1            ((( 
!     1            lgchi0(lg2,lg1) 
!     1            )/hart,
!     1            + dble(lgxsec0(lg2,lg1)) )/hart,
!     1            lg2=0,2),lg1=0,2)
             write(8,innerform) omegax(ik)*hart, &
                  xkp(ik)/bohr, &
                  xsec0 /hart, &
                  (dble(lgxsec0(lg2,lg2))/hart,lg2=0,ldecmx), &
                  ((( lgchi0(lg2,lg1) )/xsec0,lg2=0,ldecmx),lg1=0,ldecmx)
          end do
         close(8)
      end if
      endif
     

!     for if (iabs=abs); or the last absorber

      return
      end
