!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: genfmtjas.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2011/12/11 02:39:17 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine genfmtjas
      use dimsmod, only: legtot, lamtot, ltot, nex, nspx=>nspu, mtot, nphx=>nphu
	  use constants
	  use nrixs_inp,only: indmax,jinit,jmax,kfinmax,ljmax,ldecs,lind
	  use global_inp,only: ipol,ispin,le2,l2lp,ldecmx,angks,elpty,evec,xivec,ptz,nq,qaverage,qweights=>qw,qtrig
	  use genfmt_inp,only: ipr3=>ipr5,critcw,iorder,wnstar
      use nlm
      use lambda
      use clmz
      use fmatrx
      use rotmat
      use pdata
      use str 
      implicit double precision (a-h, o-z)

!  altered by matt newville (jan 1999): 
!  format of feff.bin changed to packed-ascii, and all writes changed.
!  altered by alex ankudinov(feb 2000); disabled use of paths in LDOS.
      include '../HEADERS/vers.h'

!     Aleksi
      double precision elpty0
!

      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cfac
      complex*16  cchi(nex)

      complex*16  eref2(nex,nspx), ph4, bmati
      dimension   ph4(nex,-ltot:ltot, nspx, 0:nphx)
      dimension   bmati(-mtot:mtot, 8, -mtot:mtot, 8)
      complex*16  hbmatr(-jinit:jinit,nq,-mtot:mtot, kfinmax)
      complex*16  hbmatl(-jinit:jinit,nq,-mtot:mtot, kfinmax)

      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex), ckp
      dimension   ffmag(nex)
      dimension   eps1(3), eps2(3), vec1(3), vec2(3)
      integer is2
      character*128 string
      character*512 slog
      logical done

!
!     Aleksi
!
      integer mj
      double precision hbmat(0:1,kfinmax,-jinit:jinit)
      dimension kiind(kfinmax), lgind(kfinmax),ljind(kfinmax)
      !integer lind(kfinmax) JK - lind is passed through module nrixs_inp now
      complex*16 rkk(nex,nq,kfinmax)
      complex*16 rkk2(nex,nq,kfinmax,nspx)
	  real*8 qbeta(nq)
	  complex pha(nq)

!      padlib staff
       double precision phff(nex), amff(nex), xkr(nex)
       complex*16 lgcchi
       integer  mpadx
       parameter (mpadx = 8)
       character*75 wfmt, atsym*2
       external atsym, cwig3j, istrln
       complex*16  fmatl(-jinit:jinit,nq,lamtot),fmatr(-jinit:jinit,nq,lamtot)
       

!
!     Stuff for the decomposition
!
       complex*16  lgfmatl(-jinit:jinit,nq,0:ldecs,lamtot), lgfmatr(-jinit:jinit,nq,0:ldecs,lamtot)
       complex*16  pgtrl(0:ldecs,0:ldecs,nex)
       integer ilm1,ilm2
       integer ix1,ix2,mj1,mj2,ig1,ig2,lg1,lg2
       integer ims1,ims2
	   real*8 qcst(nq),qsnt(nq),qcsf(nq),qsnf(nq)
!
!
!     These are for the case when spherical averaging is done
!
      double precision hbmats(0:1,0:1,kfinmax,-jinit:jinit)
      complex*16  fmats(-jinit:jinit,0:1,lamtot,lamtot)
      complex*16  lgfmats(-jinit:jinit,0:1,0:ldecs,lamtot,lamtot)
      complex*16  hbmatrs(-jinit:jinit,0:1,-mtot:mtot,-mtot:mtot, kfinmax)

	  complex*16 fmatsum


       integer idum,indmaxt

       

!     Input flags:
!     iorder, order of approx in f-matrix expansion (see setlam)
!             (normal use, 2.  Do ss exactly regardless of iorder)

!     used for divide-by-zero and trig tests
       parameter (eps = 1.0e-16)
       external getxk, xstar
       elpty0=0.0d0
!
!     Have to read indmax here to save some room.
!     kfinmax is too large
!

!write(*,*) 'dimension for NRIXS:'
!write(*,*) 'nex',nex
!write(*,*) 'nphx',nphx
!write(*,*) 'ltot',ltot
!write(*,*) 'mtot',mtot
!write(*,*) 'jinit',jinit
!write(*,*) 'nq',nq
!write(*,*) 'kfinmax',kfinmax
!write(*,*) 'lamtot',lamtot
!write(*,*) 'ldecs',ldecs


       open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')
     
      read(1,11) idum, idum, idum, idum, idum, idum, idum, idum, idum,indmaxt
      close(1)
 11   format (10(1x,i4))
!     Read phase calculation input, data returned via commons
      call rdxsphjas (ne, ne1, ne3, npot, ihole, rnrmav, xmu, edge, ik0, em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1)

       call setkap (ihole, kinit, linit)
!
!     zero gtrl
!

!     read in the angles of the rotations for the different q-vectors
!     KJ : this used to be in a subroutine "angread" ; I've substituted the code below.  3-2011
      if(.not. qaverage) then   !KJ if (nq.gt.0) then
	     qcst=qtrig(:,1)
		 qsnt=qtrig(:,2)
		 qcsf=qtrig(:,3)
		 qsnf=qtrig(:,4)
!         call angread(nqloc,pha,beta,qweights)
         do iq=1,nq
		    pha(iq)=cmplx(qcsf(iq),qsnf(iq))
			pha(iq)=conjg(pha(iq))
			beta(iq)=atan2(qsnt(iq),qcst(iq))
!			if(beta(iq).lt.0.0d0 .or. beta(iq).gt.pi) stop 'something wrong with angles in getgtrjas'
		 enddo
      else
         pha(1)=cmplx(1.0d0,0.0d0)
         beta(1)=0.0d0
         qweights(1)=1.0d0
      end if

!write(*,*) 'indmax',indmax
!write(*,*) 'indmaxt',indmaxt

       ilinit = linit + 1
       is = 1
       if (ispin.eq.1) is = nspx
!      need to sum over spin-up and -down for |ispin|=1 (fix later)
!      as now |ispin|= 1 and 2 should give same answer with path
!      expansion, but not with FMS
!      for ispin=2 the variables already written into is=1 positions
       do 10 ie = 1, ne
  10   eref(ie) = eref2(ie,is)
       do 20 iph = 0, npot
       do 20 ie = 1, ne
       do 20 il = -lmax(ie, iph), lmax(ie, iph)
  20   ph(ie,il, iph) = ph4(ie, il, is, iph)
       do ie=1,ne
          do kdif=1,indmax
		  do iq=1,nq 
             rkk(ie,iq,kdif) = rkk2(ie,iq,kdif,is)
		  enddo
          enddo
       enddo
!     Open path input file (unit in) and read text .  Use unit 1.
       ntext  = 5
       open (unit=1, file='paths.dat', status='old', iostat=ios)
       call chopen (ios, 'paths.dat', 'genfmt')
       call rdhead (1, ntext , text, ltext)
       if (ntext  .le. 0)  then
          text (1) = ' '
       endif
!     Save indices of paths for use by ff2chi
       open (unit=2, file='list.dat', status='unknown', iostat=ios)
       call chopen (ios, 'list.dat', 'genfmt')
!     Put phase header on top of list.dat
       call wthead (2, ntext , text )
       write(2, 125)
 125   format (1x, 71('-'))
       write(2, 135)
 135   format ('  pathindex     sig2   amp ratio    ','deg    nlegs  r effective')
       
!     Open nstar.dat if necessary
       if (wnstar)  then
          open (unit=4,file='nstar.dat', status='unknown', iostat=ios)
          call chopen (ios, 'nstar.dat', 'genfmt')
          write(4,'(1x,a,f8.4)' ) ' polarization', evec
          write(4,'(1x,a)' ) ' npath     n*'
       endif
!     Set crit0 for keeping feff.dat's
       if (ipr3 .le. 0)  crit0 = 2*critcw/3
!     Make a header for the running messages.
       write(slog, 155) critcw
 155   format ('    Curved wave chi amplitude ratio', f7.2, '%')
       call wlog(slog)
       if (ipr3 .le. 0)  then
         write(slog,165) crit0
         call wlog(slog)
       endif
 165   format ('    Discard feff.dat for paths with cw ratio <', f7.2, '%')
       write(slog,195)
 195   format ('    path  cw ratio     deg    nleg  reff')
       call wlog(slog)

!     open feff.bin for storing path info
!     for now, use double precision.  After it's working, try
!     single precision.
!     Use single precision for all fp numbers in feff.bin
      open (unit=3, file='feff.bin', status='unknown', iostat=ios)
      call chopen (ios, 'feff.bin', 'genfmt')
!
!      open file for decomposition
!
      if (ldecmx.ge.0) then 
         open (unit=7,file='feffl.bin',status='unknown',iostat = ios)
         call chopen(ios,'feffl.bin','genfmt')
      end if
!     put label line in feff.bin so other programs know it really
!     is a feff.bin file
       string = '#_feff.bin v03: ' // vfeff
       jstr   = istrln(string)
       write(3, '(a)')  string(1:jstr)

!     save stuff that is the same for all paths
!     header, ck, central atom phase shifts
       write(3, '(a2,6(1x,i4))') '#_', npot, ne, mpadx

!     Misc stuff from phase.bin and genfmt call
 345   format(a2,3(1x,i7), 3(1x,g14.7))
       write(3, 345) '#=', ihole, iorder, ilinit, rnrmav, xmu, edge
       do 380 i = 0, npot
          if (potlbl(i).eq.' ') potlbl(i)  = atsym(iz(i))
          if (potlbl(i).eq.' ') potlbl(i)  = 'null'
 380   continue 
 395   format('(',i3,'(1x,a6),',i3,'(1x,i3))')
       write(wfmt, 395) npot+1, npot+1
       write(string,wfmt) (potlbl(i),i=0,npot) , (iz(i),i=0,npot)
       jstr = istrln(string)
       write(3, '(a2,a)') '#@',string(:jstr)

!     Central atom phase shifts
      ll = linit+1
      if (kinit.lt.0) ll = -ll
      call wrpadx(3,mpadx, ph(1,ll, 0),ne)

!     Set nlm factors in common /nlm/ for use later
      call snlm (ltot+1, mtot+1)

!     Make xk and ck array for later use
       do 850  ie = 1, ne
!        real momentum (k)
         xk(ie) = getxk (dble(em(ie)) - edge)
!        complex momentum (p)
         ck(ie) = sqrt (2*(em(ie) - eref(ie)))
         ckmag(ie) = abs(ck(ie))
         xkr(ie) = real(xk(ie))
 850   continue
       call wrpadx(3,mpadx, ck,ne)
       call wrpadd(3,mpadx, xkr,ne)
!     While not done, read path, find feff.
       npath  = 0
       ntotal = 0
       nused  = 0
       xportx = -1
	   xport = 0.d0
 1000  continue

!        Read current path
         call rdpath (1, done, ipol)
         icalc = iorder
         if (.not.done)  then
            npath = npath + 1
            ntotal = ntotal + 1
            if (wnstar)  then
!              should be ipol=1
               do 1150 ic =1,3
                  vec1(ic) = rat(ic,1) - rat(ic,0)
                  vec2(ic) = rat(ic,nleg-1) - rat(ic,0)
                  eps1(ic) = evec(ic)
 1150          continue
               if (elpty0.ne.0.0) then
                  eps2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
                  eps2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
                  eps2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
               endif
               ndeg = nint (deg)
               xxstar = xstar (eps1, eps2, vec1, vec2, ndeg, elpty0)
               write(4,'(1x,i6,f10.3)')  npath, xxstar
            endif
            
!        Need reff
         reff = 0
         do 1200  i = 1, nleg
            reff = reff + ri(i)
 1200    continue 
         reff = reff/2

!        Set lambda for low k
         call setlam (icalc, 1)

!        Calculate and store rotation matrix elements
         do 1300  isc = 1, nleg
            call rot3i (lmaxp1, mmaxp1, isc)
 1300    continue

!
!     Zero pgtrl if needed
!
         if (ldecmx.ge.0) then 
            do ie=1,nex
               do ilm1=0,ldecmx
                  do ilm2=0,ldecmx
                     pgtrl(ilm2,ilm1,ie)=dcmplx(0.0d0,0.0d0)
                  end do
               end do   
            end do  
         end if

         if (ipol.gt.0)  then
!           one more rotation in polarization case
!           NEED MORE rot3j FOR CENTRAL ATOM ( l \pm 1 )
            call rot3i (lmaxp1, mmaxp1, nleg+1) ! JK - I don't think there is any reason that we should be limited to ilinit+1 for NRIXS.
			!KJ 3-2011  I'm unsure of this, but Aleksi Soininen changed the rot3i call as follows in "withnqwithmdff" :
			!call rot3i(ilinit+1,ilinit+1,nleg+1)
			!Possibly we need more in order to rotate the Green's matrix as in mkgtr??
         endif 
         if (elpty.ge.0) then 
            call  mmtrjas( hbmatr,hbmatl, ipol, ispin, ptz, nspx,pha,qbeta) !KJ 3-2011 last 2 added for withnqwithmdff
!KJ            call  mmtrjas( hbmatr,hbmatl, ipol, ispin, le2, angks, ptz, &
!KJ                 nspx,jinit,jmax,kfinmax,ljmax,indmaxt,lind)
         else
            call  mmtrjas0( hbmatrs,ipol, ispin, ptz, nspx, ltot)
!KJ            call  mmtrjas0( hbmatrs,ipol, ispin, le2, angks, ptz, &
!KJ                 nspx,jinit,jmax,kfinmax,ljmax,ltot,indmaxt,lind)
         end if
         if (indmaxt.ne.indmax) then 
            write(6,*) "indmax mismatch in genfmtjas"
            stop
         end if

!
!        Big energy loop
         do 5000  ie = 1, ne
!           complex rho
            do 2010  ileg = 1, nleg
               rho(ileg) = ck(ie) * ri(ileg)
 2010       continue
!           if ck is zero, xafs is undefined.  Make it zero and jump
!           to end of calc part of loop.
            if (abs(ck(ie)) .le. eps)  then
               cchi(ie) = 0
               write(slog,2055)  ie, ck(ie)
 2055          format (' genfmt: ck=0.  ie, ck(ie)', i5, 1p, 2e14.5)
               call wlog(slog)
               goto 4990
            endif
!           Calculate and store spherical wave factors c_l^(m)z^m/m!
!           in a matrix clmi(il,im,ileg), ileg=1...nleg.
!           Result is that common /clmz/ is updated for use by fmtrxi.
!
!           zero clmi arrays
            do 2100  ileg = 1, legtot
               do 2100 im = 1, mtot+ntot+1
                  do 2100  il = 1, ltot+1
                     clmi(il,im,ileg) = 0
 2100       continue
            mnmxp1 = mmaxp1 + nmax
            do 2150  ileg = 1, nleg
               isc0 = ileg-1
               if (isc0.eq.0) isc0=nleg
               isc1 = ileg
               lxp1 = max (lmax(ie,ipot(isc0))+1, lmax(ie,ipot(isc1))+1)
               mnp1 = min (lxp1, mnmxp1)
               call sclmz (rho, lxp1, mnp1, ileg)
 2150       continue

!           Calculate and store scattering matrices fmati.
!           First matrix
            call fmtrxi (lamx, laml0x, ie, 2, 1)
!           Last matrix if needed
            if (nleg .gt. 2)  then
               call fmtrxi (laml0x, lamx, ie, nleg, nleg-1)
            endif
!           Intermediate scattering matrices
            do 2200  ilegp = 2, nsc-1
               ileg = ilegp + 1
               call fmtrxi (lamx, lamx, ie, ileg, ilegp)
 2200       continue

!           Big matrix multiplication loops.
!           Calculates trace of matrix product
!           M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
!           We will calculate the trace over lambda_N, working from
!           right to left.
!           Use only 2 pmati arrays, alternating indp (index p)
!           1 and 2.

!           to start f(2,1) -> pmat(1)
            indp = 1
            do 2250 lmp = 1, laml0x
            do 2250 lm = 1, lamx
               pmati(lm,lmp,indp)= fmati(lm,lmp,1)
 2250       continue
!           f(N,N-1) * ... * f(3,2) * [f(2,1)]
!           Term in [] is pmat(1)
            do 2900 isc = 2, nleg-1
!              indp is current p matrix, indp0 is previous p matrix
               indp = 2 - mod(isc,2)
               indp0 = 1 + mod(indp,2)
               do 2850  lmp = 1, laml0x
               do 2850  lm = 1, lamx
                  pllp=dcmplx(0.0d0,0.0d0)
                  do 2800 lmi = 1, lamx
                     pllp = pllp + fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
 2800             continue
 2850             pmati(lm,lmp,indp) = pllp
 2900       continue

!           srho=sum pr(i), prho = prod pr(i)
            srho=0
            prho=1
            do 3200  ileg = 1, nleg
               srho = srho + rho(ileg)
               prho = prho * rho(ileg)
 3200       continue

!           Termination matrix, fmati(...,nleg)
!           Polarization enters only this matrix
!           this will fill fmati(...,nleg) in common /fmtrxi/
            if (elpty.ge.0.0d0) then 

               call mmtrxijas(kfinmax,indmax,l2lp, &
                    rkk,laml0x,jinit,hbmatl, &
                    hbmatr, &
                    ie,1,nleg,lind,fmatl,fmatr,ldecmx,ldecs, &
                    lgfmatr,lgfmatl) ! JK - switched lgfmatr and lgfmatl to be consistent with notation in subroutine.
                                     !      Doesn't actually matter though.
!     Final trace over matrix, first for the total
               ptrac=0
               do 4400  lm = 1, laml0x
                  do 4400  lmp = 1, laml0x
                     
                     fmatsum=dcmplx(0.0d0,0.0d0)
					 do iq=1,nq 
                     do mj=-jinit,jinit,2
                        fmatsum=fmatsum+fmatl(mj,iq,lm)*fmatr(mj,iq,lmp)
                     end do
					 end do
                     ptrac = ptrac+fmatsum*pmati(lmp,lm,indp)
              
 4400          continue
            else
               call mmtrxijas0(kfinmax,indmax,l2lp,rkk,laml0x,jinit,hbmatrs,ie,1,nleg,lind,fmats,ldecmx,ldecs,lgfmats)
               
!     Final trace over matrix, first for the total
               ptrac=0
               do 4500  lm = 1, laml0x
                  do 4500  lmp = 1, laml0x
                     fmatsum=dcmplx(0.0d0,0.0d0)
                     do is2=0,1
                        do mj=-jinit,jinit,2
                           fmatsum=fmatsum+fmats(mj,is2,lmp,lm)
                        end do
                     end do
                   ptrac = ptrac+fmatsum*pmati(lmp,lm,indp)
 4500          continue 
            end if
! 
!           Calculate xafs
!           Complex chi (without 2kr term)
!           ipot(nleg) is central atom
!           cdel1 = exp(2*coni*ph(ie,ilinit+1,0))
!           central atom phase shift are included in normalized
!           reduced matrix elements rkk(....)
            cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho

!           now factor 1/(2*l0+1) is inside termination matrix
!           cchi(ie) = ptrac * cfac/(2*l0+1)
            cchi(ie) = ptrac * cfac

!           When ck(ie)=0, xafs is set to zero.  Calc above undefined.
!           Jump to here from ck(ie)=0 test above.

!
!  
!     
            if (ldecmx.ge.0) then 
               if (elpty.ge.0d0) then 
                  lgcchi=dcmplx(0.0d0,0.0d0)
                  do lg1=0,ldecmx
                     do lg2=0,ldecmx
                        ptrac=0
                        
                        do 4401  lm = 1, laml0x
                           do 4401  lmp = 1, laml0x
                              fmatsum=dcmplx(0.0d0,0.0d0)
							  do iq=1,nq 
                              do mj=-jinit,jinit,2
                                 fmatsum=fmatsum+lgfmatr(mj,iq,lg1,lm) *lgfmatl(mj,iq,lg2,lmp)
                              end do
							  end do
                              ptrac = ptrac+fmatsum*pmati(lmp,lm,indp)
 4401                   continue
                         cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho
                         pgtrl(lg2,lg1,ie) = ptrac * cfac
                         lgcchi=lgcchi+ pgtrl(lg2,lg1,ie)
                      end do    !lg2
                   end do       ! lg1
                else
                   lgcchi=dcmplx(0.0d0,0.0d0)
                   do lg1=0,ldecmx
                      ptrac=0
                      do 4501  lm = 1, laml0x
                         do 4501  lmp = 1, laml0x
                            fmatsum=dcmplx(0.0d0,0.0d0)
                            do is2=0,1
                               do mj=-jinit,jinit,2
                                  fmatsum=fmatsum+lgfmats(mj,is2,lg1,lmp,lm)
                               end do
                            end do
                            ptrac = ptrac+fmatsum*pmati(lmp,lm,indp)
                            
 4501                 continue
                      cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho
                      pgtrl(lg1,lg1,ie) = ptrac * cfac
                      lgcchi=lgcchi+ pgtrl(lg1,lg1,ie)
                      
                   end do       ! lg1
                end if
             end if
 4990      continue

 5000    continue
!        end of energy loop

!        Make importance factor, deg*(integral (|chi|*d|p|))
!        make ffmag (|chi|)
!        xport   importance factor
         do 6810  ie = 1, ne1
            ckp = ck(ie)
            xlam0 = dimag(ck(ie)) - dimag(ckp)
            ffmag(ie) = abs( cchi(ie) * exp(2*reff*xlam0) )
 6810    continue

!        integrate from edge (ik0) to ne
         nemax = ne1 - ik0 + 1
         call trap (ckmag(ik0), ffmag(ik0), nemax, xport)
         xport = abs(deg*xport)

         if (xportx.le.0)  xportx = xport
         crit = 100 * xport / xportx
!        Write path data to feff.bin if we need it.
         if (ipr3 .ge. 1  .or.  crit .ge. crit0)  then
!           write path info
 7225       format('(i6,1x,i3,1x,f7.3,1x,f11.7,1x,f9.4,',i3,'(1x,i2))')
            write(wfmt, 7225) nleg
            write(string,wfmt) ipath, nleg, deg, reff*bohr,crit, (ipot(i),i=1, nleg)
            jstr = istrln(string)
            write(3,'(a2,a)') '##',string(:jstr)
            call wrpadd(3,mpadx, rat(1,1),3*nleg)
            call wrpadd(3,mpadx, beta,nleg)
            call wrpadd(3,mpadx, eta,nleg)
            call wrpadd(3,mpadx, ri,nleg)
            phffo = 0
            do 7700  ie = 1, ne
               phff(ie) = 0
               if (abs(cchi(ie)) .ge. eps) then
                  phff(ie) = atan2 (dimag(cchi(ie)),dble(cchi(ie)))
               end if

!  remove 2 pi jumps in phase
               if (ie.gt.1) call pijump (phff(ie), phffo)
               phffo    = phff(ie)
               amff(ie) = dble(abs(cchi(ie)))
 7700       continue
            call wrpadd(3,mpadx, amff,ne)
            call wrpadd(3,mpadx, phff,ne)
!
!     assuming that everything from now on is LINEAR 
!     we output the same infromation for all the 9*9 
!     different parts of pgtrl
!
            if (ldecmx.ge.0) then 
               do ilm1=0,ldecmx
                  do ilm2=0,ldecmx
                     phffo = 0
                     do 7701  ie = 1, ne
                        phff(ie) = 0
                        
                        if (abs(pgtrl(ilm2,ilm1,ie)) .ge. eps) then
                           phff(ie) = atan2 (dimag(pgtrl(ilm2,ilm1,ie)),dble(pgtrl(ilm2,ilm1,ie)))
                        end if
                        
!     remove 2 pi jumps in phase
                        if (ie.gt.1) call pijump (phff(ie), phffo)
                        phffo    = phff(ie)
                        amff(ie) = dble(abs(pgtrl(ilm2,ilm1,ie)))
 7701                continue
                     call wrpadd(7,mpadx, amff,ne)
                     call wrpadd(7,mpadx, phff,ne)
                     
                  end do
               end do
            end if



!           Put feff.dat and stuff into list.dat
!           zero is debye-waller factor column
            write(2,8215) ipath, zero, crit, deg, nleg, reff*bohr
 8215       format(1x, i8, f12.5, 2f10.3, i6, f9.4)

!           Tell user about the path we just did
            write(slog, 8225) ipath, crit, deg, nleg, reff*bohr
 8225       format (3x, i4, 2f12.5, i6, f9.4)
            call wlog(slog)
            nused = nused+1
         else
!           path unimportant, tell user
            write(slog, 8235) ipath, crit, deg, nleg, reff*bohr
 8235       format (3x, i4, 2f12.5, i6, f9.4, ' neglected')
            call wlog(slog)
         endif
!  goto next path
         goto 1000
!  done with loop over paths
       end if
!     close paths.dat, list.dat, feff.bin, nstar.dat
       close (unit=1)
       close (unit=2)
       close (unit=3)
       if (ldecmx.ge.0) then 
          close (unit=7)
       end if
       if (wnstar) close (unit=4)
       write(slog,'(1x,i4,a,i4,a)') nused,' paths kept, ',ntotal,' examined.'
       call wlog(slog)
       return
       end
