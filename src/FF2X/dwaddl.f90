!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dwaddl.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/12/10 23:17:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dwaddl (ntotal,nptot,idwopt,ip,index,crit,critcw,sig2g, &
        sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad, &
        alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, &
        ldecmx,ldecs,lgachi,lgphchi, &
        nkx, xk, xk0, &
        xkp, cchi,lcchi, & 
        iabs, nabs, ispec, ipr4, nhead, &
        head, vrcorr, vicorr, nused)
	  use constants
	  use dimsmod, only: npx=>npx_ff2x, nheadx, legtot, nex, nphx=>nphu
      implicit double precision (a-h, o-z)


      parameter (eps4 = 1.0e-4)
      character*80  head(nheadx)
      !parameter (npx = 1200) !now in dimsmod
!     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

!KJ      parameter (nfinex = 601)  !KJ 12-2011 see remarks in subroutine dwadd and ff2afsjas
	  parameter (nfinex = nex)
      complex*16 cchi(nfinex), ccpath(nfinex), ccc, ckp
!     to keep Im part of cchi 11.18.97 ala
      complex*16 dw, dw1, dw3
      dimension xkp(nfinex), xk0(nfinex)

      logical dwcorr
      dimension rattmp(3,0:legtot)
      dimension iztmp(0:legtot)
      character*512 slog
      character*12 fname
      real rnrmav
      dimension iz(0:nphx)
!     central atom phase shift at l0
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
      integer ldecs,ldecmx
      real lgachi(nex,0:ldecs,0:ldecs,npx),lgphchi(nex,0:ldecs,0:ldecs,npx)
      complex*16 lcchi(nfinex,0:ldecs,0:ldecs),lccpath(nfinex,0:ldecs,0:ldecs)
      integer ilm1,ilm2

      dimension sig2x(0:nphx, 0:nphx)
      character*2 coment
      parameter (coment='# ')
      
      integer istrln
      external istrln 
!
!     Debugging 
!
      complex   ss(nex,npx),ss2(nex,npx)

!
!     DENDS
!




!     Keep stats on paths used to make chi



      nused = 0
      xkref = dble(ck(1)**2) - xk(1)*abs(xk(1)) 
!     open the files for sigrm and sigem
      if (idwopt.eq.1) then
         iem = 111
         open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_em.dat', 'sigem')
      elseif (idwopt.eq.2) then
         irm1 =111
         open(unit=irm1,file='s2_rm2.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm2.dat', 'sigrm')
         irm2 = 112
         open(unit=irm2,file='s2_rm1.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm1.dat', 'sigrm')
      endif
      if (alphat .gt. 0) then
        icum = 113
        open(unit=icum, file='cum.dat', status='unknown', iostat=ios)
        call chopen (ios, 'cum.dat', 'sig3')
        Write(icum, 363)
  363  format('# first and third icumulant for single scattering paths')
        write(icum,364) thetae, alphat
        write(icum,365)
  364   format ('# Einstein-Temp. =',f9.2 ,'   ', 'alpha=',f9.5)
  365   format ('#       file   sig1    sig2    sig3 ')
      endif

      if (idwopt.ge.1) then
!        initialize statistics for max DW for sigrm
         sig2mx=0
         do 400 iph1=0,nphx
         do 400 iph2=0,nphx
  400    sig2x(iph1, iph2) = 0
      endif
     

!     cycle over all paths in the list
      do 560  ilist = 1, ntotal
!        find index of path
         do 410  j = 1, nptot
            if (ip(ilist) .eq. index(j))  then
               ipath = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, ip(ilist)
  420    format (' did not find path i, ip(i) ', 2i10)
         call wlog(slog)
  430    continue
!        Path ipath is the path from feff.bin that corresponds to
!        the path ilist in list.dat.  The index of the path is
!        ip(ilist) and index(ipath).

!        Use this path if it passes critcw filter
         if (crit(ipath) .lt. critcw)  goto 550
!         write(6,*) "Now here"
!        do debye-waller factors, get sig2d from correlated debye 
!        model if required
!        A note about units:  sig2g, sig2u() and sig2d are all in
!        Angs**2.  Convert to code units after we've summed them.
         sig2 = sig2g + sig2u(ilist)
         if (dwcorr .and. idwopt.ge.0)  then
!            write(6,*) "JH"
!           note that stuff from feff.bin is single precision and
!           mostly in multidim. arrays.  sigms is double precision
!           and its arrays are dimensioned for a single path, so
!           use tmp variables to call it.  tk, thetad and sig2d are 
!           all dp, and therefore OK.  Also note that sigms takes
!           inputs in angstroms, except for rs which is in bohr.
            rs = rnrmav
            do 460  ileg = 1, nleg(ipath)
!               write(6,*) "ileg",ileg,nleg(ipath)
               
               iztmp(ileg) = iz(ipot(ileg,ipath))
               do 450  j = 1, 3
                  rattmp(j,ileg) = rat(j,ileg,ipath) * bohr
  450          continue
  460       continue
            iztmp(0) = iztmp(nleg(ipath))
            do 470  j = 1,3
               rattmp(j,0) = rattmp(j,nleg(ipath))
  470       continue
            if (idwopt.eq.0) then 
!             use CD model
              call sigms (tk, thetad, rs, legtot, nleg(ipath), rattmp, iztmp, sig2d)
            elseif (idwopt.eq.1) then 
!             use EM method
              call sigem(sig2mx,sig2x,iem,tk,ipath,nleg(ipath),rattmp,sig2d)
            else 
!             use RM
              call sigrm(sig2mx,sig2x,irm1,irm2,tk,ipath,nleg(ipath),rattmp,sig2d)
            endif
            sig2 = sig2 + sig2d
!            write(6,*) "JH done"
         endif
         sig2 = sig2 / (bohr**2)
         
!        Do first and third cumulants
         sig1 = 0
         sig3 = 0
!         write(6,*) "Check"
         if (alphat .gt. zero  .and. nleg(ipath) .eq. 2)  then
           if (thetae.le.0.d0) then
!            call sig3  to get sig1 and sig3 for single scattering paths
!           use reff(ipath) for r, note that reff is single precision
             iz1 = iztmp(nleg(ipath))
             iz2 = iztmp(1)
             call sigte3(iz1, iz2, sig2, alphat, thetad, reff(ipath), sig1, sig3)
           else
!            this gets sig1 and sig3 for single scattering paths
!            using Morse potential
             call sigm3(sig1, sig2, sig3, tk, alphat, thetae)
           endif
           write(icum,475) index(ipath),  sig1 * bohr, sig2*(bohr**2), sig3*(bohr**3)
  475      format( i10,f9.5,f9.5,' ',f9.7)
         endif

!        put the debye-waller factor and other cumulants into 
!        achi and phchi
         if (mbconv .gt. 0) s02 = 1.0
         do 480  i = 1, ne1
!            write(6,*) i
            dw = exp(-2 * sig2 * ck(i)**2)
            dw1 = exp (2 * coni * ck(i) * sig1)
            dw3 = exp ((-4 * coni * ck(i)**3 * sig3) / 3)
            dw = dw * dw1 * dw3
            phdw = 0.0
            if (abs(dw).gt.0) phdw = atan2 (dimag(dw), dble(dw))
            achi(i,ipath) = achi(i,ipath) * abs(dw) * s02 * deg(ipath)
            phchi(i,ipath) = phchi(i,ipath) + phdw
!
!     make same for the other stuff
!
            if (ldecmx.ge.0) then 
               
               do ilm1=0,ldecmx
                  do ilm2=0,ldecmx
!     write(6,*) ilm2,ilm1,i,ipath
                     lgachi(i,ilm2,ilm1,ipath) = lgachi(i,ilm2,ilm1,ipath) * abs(dw) * s02 * deg(ipath)
                     lgphchi(i,ilm2,ilm1,ipath) = lgphchi(i,ilm2,ilm1,ipath) + phdw
                  end do              
               end do
            end if

  480    continue

!      do i=1,npx
!         do j=1,nex
!            ss2(j,ipath)=cmplx(0.0d0,0.0d0)
!            ss(j,ipath)=achi(j,ipath)*exp(coni*phchi(j,ipath))
!            do lg1=0,2
!               do lg2=0,2
!            ss2(j,ipath)=ss2(j,ipath)+lgachi(j,lg2,lg1,ipath)
!     1                 *exp(coni*lgphchi(j,lg2,lg1,ipath))
!               end do
!            end do
!         end do
!      end do
!       do i=1,npx
!          do j=1,nex
!             write(67,'(2I5,4e15.6)') i,j,aimag(ss(j,ipath)),
!c     1            aimag(ss2(j,ipath))
!     1            ,real(ss(j,ipath)),real(ss2(j,ipath))
!             if (abs(ss(j,i)).gt.0.0d0) then
!             write(6,*) i,j,(abs(ss2(j,i))-abs(ss(j,i)))/abs(ss(j,i))
!     1            ,abs(ss(j,i))
!          end if
!          end do
!       end do

!         write(6,*) "done here"
!        make sure no 2pi jumps in phase
         do 490  i = 2, ne1
!           phchi is single precision, so use tmp variables
            curr = phchi (i, ipath)
            old = phchi (i-1, ipath)
            call pijump (curr, old)
            phchi (i, ipath) = curr
  490    continue
!         write(6,*) "H1"
         if (ldecmx.ge.0) then 
            do ilm1=0,ldecmx
               do ilm2=0,ldecmx
                  do i = 2, ne1
!     phchi is single precision, so use tmp variables
                     curr = lgphchi (i,ilm2,ilm1, ipath)
                     old = lgphchi (i-1,ilm2,ilm1,ipath)
                     call pijump (curr, old)
                     lgphchi (i,ilm2,ilm1, ipath) = curr
                  end do
               end do
            end do
         end if
!         write(6,*) "H2"
         
         do 500  ik = 1, nkx
            call terp1 (xk, achi(1,ipath),  ne1, xk0(ik), achi0)
            call terp1 (xk, phchi(1,ipath), ne1, xk0(ik), phchi0)
            ccpath(ik) = achi0 * exp (coni * (2 * xk0(ik) * reff(ipath) + phchi0))
!           note that this already includes s02, deg, sig2, etc.
!           sum total complex chi
            cchi(ik) = cchi(ik) + ccpath(ik)
  500    continue
         nused = nused + 1
!
!     
!
         if (ldecmx.ge.0) then 
            do ilm1=0,ldecmx
               do ilm2=0,ldecmx
                  do ik = 1, nkx
                     call terp1 (xk, lgachi(1,ilm2,ilm1,ipath),ne1, xk0(ik), achi0)
                     call terp1 (xk, lgphchi(1,ilm2,ilm1,ipath),ne1, xk0(ik), phchi0)
                  
!                  write(6,*) "Here",ik,ilm2,ilm1,achi0
                     lccpath(ik,ilm2,ilm1) = achi0 * exp (coni * (2 * xk0(ik) * reff(ipath) + phchi0))
!     note that this already includes s02, deg, sig2, etc.
!     sum total complex chi
                     lcchi(ik,ilm2,ilm1) = lcchi(ik,ilm2,ilm1) + lccpath(ik,ilm2,ilm1)
!     write(6,*) "Here",ik,ilm2,ilm1, lcchi(ik,ilm2,ilm1)
                     
                  end do
               end do
            end do
         end if
!         write(6,*) "H3"

         if (iabs.eq.nabs) then
!           Put path into chi.dat, xmu.dat as required
            if (abs(sig2u(ilist)) .gt. 0.000001)  then
              write(3,515)  coment, index(ipath), sig2*(bohr**2), &
                crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, sig2u(ilist)
              write(8,515)  coment, index(ipath), sig2*(bohr**2), &
                crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, sig2u(ilist)
            else
              write(3,515) coment, index(ipath), sig2*(bohr**2), &
                crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
              write(8,515) coment, index(ipath), sig2*(bohr**2), &
                crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
            endif
  515       format(a2, 1x, i10, 5x, f9.5, 2f10.2, i6, f9.4, f9.5)
         endif

!        write out a chinnnn.dat for this path, if necessary.
         if (ipr4 .ge. 2 .and. iabs.eq.nabs .and. ispec.eq.0)  then
!           make filename chipnnnn.dat
            write(fname,520)  index(ipath)
  520       format('chip', i4.4, '.dat')
            open (unit=9, file=fname, status='unknown',iostat=ios)
            call chopen (ios, fname, 'ff2chi')
            do 530  ihead = 1, nhead
               lhead = istrln(head(ihead))
               if (lhead .gt. 0)  then
                  write(9,320) head(ihead)(1:lhead)
  320             format (a)
               endif
  530       continue
            if (dwcorr)  then
               write(9,340)  s02, tk, thetad, sig2g
  340          format (' S02', f7.3, '  Temp', f8.2,'  Debye temp',f8.2,'  Global sig2', f9.5)
            else
               write(9,341)  s02, sig2g
  341          format (' S02', f7.3,'                                        Global sig2', f9.5)
            endif
            if (alphat .gt. zero)  then
               write(9,321)  alphat
  321          format (' 1st and 3rd cumulants, alphat = ', 1pe20.4)
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,342)  vrcorr, vicorr
  342          format (' Energy zero shift, vr, vi ', 1p, 2e14.5)
            endif
            write(9,*) 'Debye-waller factor ', sig2, sig3

            write(9,610)
  610       format (1x, 71('-'))
            write(9,535)
  535       format ('       k         chi           mag          ','phase        phase-2kr  @#')
            do 540  i = 1, nkx
               ckp = sqrt (xkp(i)*abs(xkp(i)) + xkref)
!              it would be better to use interpolation for ckp
!              fix later if complaints about chipnnn.dat files, ala
               xlam0 =  - dimag(ckp)
               ccc = ccpath(i) * exp(2 * reff(ipath) * xlam0)
               phase = 0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,630)  xkp(i)/bohr, dimag(ccc), abs(ccc), phase,phase-2*xk0(i)*reff(ipath)
  630          format (1x, f10.4, 3x, 4(1pe13.6,1x))
  540       continue
            close (unit=9)
         endif

  550    continue
  560 continue

!     close files opened for sigem and sigrem
      if (idwopt.eq.1) then
        close (unit=iem)
      elseif (idwopt.eq.2) then
        close (unit=irm1)
        close (unit=irm2)
      endif
      if (alphat .gt. 0) then
        close (unit=icum)
      endif


      return
      end
