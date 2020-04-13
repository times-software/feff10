!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2xmujas.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2012/02/03 07:17:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ff2xmujas (ispec, ipr4, idwopt, critcw, s02, sig2g, tk, thetad, mbconv, &
                         vrcorr, vicorr, alphat, thetae, iabs, nabs, kfinmax, xivec, ldecmx,ldecs)
!     adds the contributions from each path and absorber, including
!     sDebye-Waller factors. Writes down main output: chi.dat and xmu.dat
      use constants
	  use dimsmod, only: npx=>npx_ff2x, nheadx, nex, nphx=>nphu, legtot
      implicit double precision (a-h, o-z)

      parameter (eps4 = 1.0e-4)
      !parameter (npx = 1200) !now in dimsmod
	  
!     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)
      complex*16 gtr(nex)
      integer kfinmax
      integer ii,ilm1,ilm2,id1,id2,lg1,lg2
      integer ldecmx,ldecs
!
!
      complex*16 lcchi(nex,0:ldecs,0:ldecs) 
      dimension rchtotl(nex,0:ldecs,0:ldecs)
      complex*16 lchia(nex,0:ldecs,0:ldecs)
      double precision lgom2x0(0:ldecs,0:ldecs)
      double precision lgxsec0(0:ldecs,0:ldecs)
      double precision lgomega(nex,0:ldecs,0:ldecs)
      double precision lgom2x(nex,0:ldecs,0:ldecs)
      double precision lgchi0(0:ldecs,0:ldecs) 
      integer  kiind(kfinmax), lgind(kfinmax),ljind(kfinmax)
      integer lind(kfinmax)!



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
	  
      double precision xivec(3)

!     stuff from xsect.bin
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex)
      character*25 innerform
      integer ltotch

	  !KJ 12-2011 Some arrays here are really big.  They require increased stack size, which is annoying for users.
	  ! Using allocate statements they can go on the heap instead.
      real,allocatable :: rat(:,:,:), beta(:,:), eta(:,:)
      real,allocatable :: ri(:,:)
      real,allocatable :: achi(:,:), phchi(:,:)
      real,allocatable :: lgachi(:,:,:,:), lgphchi(:,:,:,:)
      complex*16,allocatable :: gtrl(:,:,:), xsecl(:,:), lgxsec(:,:), xseczero(:)

      
!KJ:	  
	  allocate( lgachi(nex,0:ldecs,0:ldecs,npx), lgphchi(nex,0:ldecs,0:ldecs,npx), achi(nex,npx), phchi(nex,npx) )
	  allocate( rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx), ri(legtot,npx) )
      allocate( gtrl(0:ldecs,0:ldecs,nex), xsecl(kfinmax,nex),lgxsec(nex,0:ldecs), xseczero(nex) )
!	  write(*,*) 'dimensions'
!	  write(*,*) 'nex= ',nex
!	  write(*,*) 'npx= ',npx
!	  write(*,*) 'ldecs= ',ldecs
!	  write(*,*) 'kfinmax= ',kfinmax

!     open list.dat and read list of paths we want
      open (unit=1, file='list.dat', status='old', iostat=ios)
      ntotal = 0
      if (ios.le.0) then
        call chopen (ios, 'list.dat', 'ff2chi')
!       read title line for paths and genfmt.
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


!     get gtr - result of FMS
      do 112 ie =1,nex
  112 gtr(ie) = 0
      ntfms = 0
      open (unit=1, file='fms.bin', status='old', iostat=ios)
      if (ios.le.0) then
         ntfms = 1
         read(1, 113) titfms
  113    format(a)
         read(1, 115) ne, ne1, ne3, nph, npadx
  115    format(5(1x,i3))
         call rdpadx(1, npadx, gtr, ne)
      endif
      close (unit=1)
!
!     get the decompostion of gtr
      if (ldecmx.ge.0) then 
         open (unit=1, file='fmsl.bin', status='old', iostat=ios)
         if (ios.le.0) then
            do ie=1,ne
               call rdpadx(1,npadx,gtrl(0,0,ie),(ldecmx+1)*(ldecmx+1))
            end do
         end if
         close(1)
      end if
      
       call rdfbinl ('feff.bin', nphx, nex, npx, legtot, &
          nptot, ne, npot, ihole, iorder, ilinit, &
          rnrmav, xmu, edge, potlbl, iz, phc, ck, xk, &
          index, nleg, deg, reff, &
          crit, ipot, rat, beta, eta, ri, achi, phchi,&
          ldecmx,ldecs,lgachi, lgphchi)
!     read xsect.bin file
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,&
       emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)
!
!     read the xsectatom
      if (ldecmx.ge.0) then 
         open (unit=1,file='xsecl.bin',status='unknown',iostat = ios)
         call chopen(ios,'xsecl.bin','ff2xmuq')
         read(1,*) id1,id2
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

!
!     make combined title
      if (ntfms.eq.1) then
        ntitle = ntitle + 1
        title(ntitle) = titfms
      endif
      do 120 ihead = 1, nhead
 120  title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

!     write feffnnnn.dat
      if (ipr4.ge.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne,npot,&
             ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,&
             iz,phc,ck,xk,index,&
             nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

!     If there is a vicorr, will need a mean free path factor xlam0.
!     Use it as  chi(ie) * exp (2 * reff * xlam0)
!     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 180  ie = 1, ne
            ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
            xlam0 = aimag(ck(ie)) - dimag(ckp)
            do 170  ipath = 1, nptot
               achi(ie,ipath) = achi(ie,ipath) * &
                    exp (2 * reff(ipath) * xlam0)
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
         open (unit=3, file='chi.dat', status='unknown', iostat=ios)
         call chopen (ios, 'chi.dat', 'ff2chi')
         open (unit=8, file='xmu.dat', status='unknown', iostat=ios)
         call chopen (ios, 'xmu.dat', 'ff2chi')

!        write miscellaneous staff into headers
         call wrhead (8, ntitle, title, dwcorr, s02, tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw ) !,xivec) !KJ lazy

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
  380    format('    S02', f7.3, '  Temp', f8.2, '  Debye temp', f8.2,&
                '  Global sig2', f9.5)
  381    format('    S02', f7.3, '  Global sig2', f9.5)
      endif


!     make chi and sum it
      do 400  i = 1, nex
         cchi(i) = 0
  400 continue
      do 402  ik = 1, ne
         cchi(ik)= s02 * gtr(ik)
  402 continue

!     make decomposition
!
      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ik=1,nex
                  lcchi(ik,lg2,lg1) = 0
               end do
            end do
         end do
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ik=1,ne
                  lcchi(ik,lg2,lg1) = s02*gtrl(lg2,lg1,ik)
               end do
            end do
         end do
      end if

!     add Debye-Waller factors
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,&
       sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,&
       alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, ne, xk, xkp,&
       xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,&
       title, vrcorr, vicorr,  nused)
!     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nex
            chia(ie) = 0
  635    continue
      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ie=1,nex
                  lchia(ie,lg2,lg1) = 0
               end do
            end do
         end do
      end if
      else
         open (unit=1, file='chia.bin', status='old',&
        access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,ne
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

      if(iabs.eq.1) then
!        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           if (abs(del) .gt.  10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
             write(6,*)  xk(i),xkxs(i)
             call par_stop('FF2XMU-1') 
           endif
  680    continue
      endif

!     add contribution from an absorber iabs 
!     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, ne
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ik=1,ne
                  lchia(ik,lg2,lg1)   = lchia(ik,lg2,lg1)   &
                      + lcchi(ik,lg2,lg1)/ nabs
               end do
            end do
         end do
      end if
      if (iabs.lt.nabs) then
!        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',&
        access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,ne
  760    write(1) chia(ie)
         close(unit=1)
      endif
      if (iabs.eq.nabs) then
!        The loop over absorbers is finished. Write out the results.
         write(8,600)  coment, nused, ntotal
  600    format ( a2, 1x, i4, '/', i4, ' paths used')
  610    format ( a2, 1x, 71('-'))

         do 702 ik = 1, ne
            rchtot(ik) = dimag (chia(ik))
  702    continue

      if (ldecmx.ge.0) then 
         do lg1=0,ldecmx
            do lg2=0,ldecmx
               do ik=1,ne
                  rchtotl(ik,lg2,lg1) = dimag (lchia(ik,lg2,lg1))
               end do
            end do
         end do
      end if

!        prepare the output grid omega
         efermi = edge + omega(1) - dble(emxs(1))

!        do convolution with excitation spectrum
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv(omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv(omega, ne1, efermi, s02p, erelax, wp, rchtot)
            if (ldecmx.ge.0) then 
               do lg1=0,ldecmx
                  do lg2=0,ldecmx
                     call  exconv(omega, ne1, efermi, s02p, erelax, wp, rchtotl(1,lg2,lg1)) !KJ it said "omegax" instead of omega; but that doesn't exist - typo?? 2-2012
                  end do
               end do
               end if
         endif

!        normalize to xsec at 50 eV above edge
!        and prepare the output energy grid omega
         edg50 = efermi + 50 / hart
         if (ispec.eq.2) edg50 = efermi
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
         write(8,660)  coment
  660    format (a2, ' Contribution to S(q,w) from a single electron')
         write(8,610) coment 
         write(8,670) coment
  670    format (a2,' omega    e    k   S(qw)  S^0(qw)  chi_q*S^0(qw)       @#')
!  670    format (a2,' omega    e    k    mu    mu0     chi     @#')
!        do correction using brouder method
         vi0 = 0
         call xscorr(ispec,emxs, ne1, ne, ik0, xsec,xsnorm,chia,vrcorr, vi0, cchi)

         do 850 ie=1,ne1
           rchtot(ie)=dimag( xsec(ie)+xsnorm(ie)*chia(ie)+cchi(ie))
  850    continue

         if (ldecmx.ge.0) then 
            do lg1=0,ldecmx
               do lg2=0,ldecmx
                  if (lg1.eq.lg2) then 
                     call xscorr(ispec, emxs, ne1, ne, ik0, &
                         lgxsec(1,lg2), xsnorm, lchia(1,lg2,lg1), &
                         vrcorr, vicorr, lcchi(1,lg2,lg1))
                  else 
                     call xscorr(ispec, emxs, ne1, ne, ik0, &
                         xseczero, xsnorm, lchia(1,lg2,lg1), &
                         vrcorr, vicorr, lcchi(1,lg2,lg1))
                  end if
                  
               end do
            end do


         
            do ie=1,ne1
               do lg1=0,ldecmx
                  do lg2=0,ldecmx
                     if (lg2.eq.lg1) then 
                        rchtotl(ie,lg2,lg1)=dimag( lgxsec(ie,lg1) &
                            +xsnorm(ie)*lchia(ie,lg2,lg1) &
                            +lcchi(ie,lg2,lg1))
                     else
                        rchtotl(ie,lg2,lg1)=dimag( &
                            xsnorm(ie)*lchia(ie,lg2,lg1) &
                            +lcchi(ie,lg2,lg1))
                     end if
                  end do
               end do
            end do
         end if

         do ie=1,nex
            xseczero(ie)=0
         end do



         do 855 ie=1,ne
           chia(ie) = 0
  855    continue
         call xscorr(ispec, emxs, ne1, ne, ik0, xsec,xsnorm,chia, &
            vrcorr, vi0, cchi)
         do 856 ie = 1, ne1
 856     cchi(ie) = dimag(xsec(ie)+cchi(ie)) * coni+rchtot(ie)

!
!     
!

         if (ldecmx.ge.0) then
            do lg1=0,ldecmx
               do lg2=0,ldecmx
                  do ie=1,ne
                     lchia(ie,lg2,lg1) = 0
                  end do
                  if (lg1.eq.lg2) then 
                     call xscorr(ispec, emxs, ne1, ne, ik0, &
                         lgxsec(1,lg2), xsnorm, lchia(1,lg2,lg1), &
                         vrcorr, vicorr, lcchi(1,lg2,lg1))
                     do ie=1,ne1
                        lcchi(ie,lg2,lg1) = dimag(lgxsec(ie,lg1)+ &
                            lcchi(ie,lg2,lg1))  &
                            * coni+rchtotl(ie,lg2,lg1)
                     end do
                     
                  else 
                     call xscorr(ispec, emxs, ne1, ne, ik0, &
                         xseczero, xsnorm, lchia(1,lg2,lg1), &
                         vrcorr, vicorr, lcchi(1,lg2,lg1))
                     do ie=1,ne1
                        lcchi(ie,lg2,lg1) = dimag( &
                            lcchi(ie,lg2,lg1))  &
                            * coni+rchtotl(ie,lg2,lg1)
                     end do
                     
                  end if
               end do
            end do
         end if


         if (vicorr.gt.eps4 .and. ntotal.eq.0) then
!           add correction due to vicorr
            call conv(omega,cchi,ne1,vicorr)

            if (ldecmx.ge.0) then 
               do lg1=0,ldecmx
                  do lg2=0,ldecmx          
                     call conv(omega,lcchi(1,lg2,lg1),ne1,vicorr)
                  end do
               end do
            end if
         endif

         do 860 ie = 1, ne1
            em0 = dble(emxs(ie))
            xsec0 = dimag(cchi(ie))
            rchtot(ie) = dble (cchi(ie))
            chi0  = (rchtot(ie) - xsec0)
            write(8,700)  omega(ie)*hart, em0*hart, xkp(ie)/bohr, &
!     1              rchtot(ie)/xsedge, xsec0/xsedge, chi0
                   rchtot(ie)/hart, xsec0/hart, chi0/hart

!           if you want f'' at the output in el. units use next line
!    1          rchtot(ie)*omega(ie)*prefac, xsec0*omega(ie)*prefac, chi0
!   with        prefac = alpinv / 4 / pi /bohr**2

  700       format (1x, 2f11.3, f8.3, 1p, 3e13.5)
  860    continue

         close (unit=8)
!
!     Now do the same but for the decomposition
!
!     
         if (ldecmx.gt.0) then 
         open (unit=8, file='xmul.dat', status='unknown', iostat=ios)
         call chopen (ios, 'xmul.dat', 'ff2xmu')
         write(8,890) coment
  890    format (a2, ' Decomposition of S(q,w) for a single electron')
         write(8,895) coment
  895     format (a2,' omega    k   S^0(qw)  S_{l=0,...,ldecmx}^0(qw)  &
        &     chi^q_{l=0,..ldecmx,l^*=0,...,ldecmx} ')
         write(8,898) coment,ldecmx
 898     format (a2,'and ldecmx= ', i5)
          
!
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
          do ie = 1, ne1
             em0 = dble(emxs(ie))
             xsec0 = dimag(cchi(ie))
             rchtot(ie) = dble (cchi(ie))
             chi0  = (rchtot(ie) - xsec0)
             do lg1=0,ldecmx
                do lg2=0,ldecmx
                   lgxsec0(lg2,lg1) = dimag(lcchi(ie,lg2,lg1))
                   rchtotl(ie,lg2,lg1) = dble (lcchi(ie,lg2,lg1))
                   lgchi0(lg2,lg1)  = (rchtotl(ie,lg2,lg1) &
                       - lgxsec0(lg2,lg1))
                end do
             end do
              write(8,innerform) omega(ie)*hart, &
                 xkp(ie)/bohr, &
                 xsec0/hart, &
                 (lgxsec0(lg2,lg2)/hart,lg2=0,ldecmx), &
                 (((lgchi0(lg2,lg1) )/xsec0, &
                 lg2=0,ldecmx),lg1=0,ldecmx)
          end do ! ie
          close(unit=8)
          end if
         close (unit=3, status='delete')
      endif
!     for if (iabs=abs); or the last absorber

      return
      end
