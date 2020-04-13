      subroutine ldos_h ( nph, edens_sp, edenvl_sp, dmag, vtot_sp, vvalgs_sp, &
                  rmt,rmt_sp,  rnrm, ixc, rhoint, vint,vint_sp, xmu, jumprm,  &
                  x0, dx, rgrd, xion, iunf, iz, xnval, adgc, adpc, dgc,dpc,  &
                  ihole, qnrm, xnmues, emin, emax, eimag, rfms2, lfms2, lmaxph,  &
                  nat, iphat, rat, minv, rdirec, toler1, toler2)

      use controls,only: ispace
      use DimsMod, only: nphx=>nphu, lx, nrptx, nex, ltot, natx, nspx=>nspu
      use constants
      use par
      use ldos_inp,only: neldos
      use hubbard_inp,only: i_hubbard

      implicit none
      real*8,intent(in) :: vint_sp(2)
      real*8,intent(in) :: dmag(251,0:nphx)
      real*8,intent(in) :: vtot_sp(251,0:nphx,2), vvalgs_sp(251,0:nphx,2)
      real*8,intent(in) :: rnrm(0:nphx)
      real*8,intent(out) :: rmt(0:nphx)
      real*8,intent(in) :: rmt_sp(0:nphx,2) !, rnrm_sp(0:nphx,2)
      real*8,intent(in) :: xnval(30,0:nphx)
      integer,intent(in) :: iz(0:nphx)
      real*8,intent(in) :: xnmues(0:lx,0:nphx), qnrm(0:nphx)
      integer,intent(in) :: lmaxph(0:nphx), iphat(natx)
      real*8,intent(in) :: rat(3,natx), xion(0:nphx)
      real,intent(in) ::   rfms2, rdirec, toler1, toler2
      real*8,intent(in) :: edens_sp(251,0:nphx,2), edenvl_sp(251,0:nphx,2)
      integer,intent(in) :: lfms2,nph,nat,iunf,ihole,ixc,minv
      real*8,intent(in) :: emin,emax,rhoint,dx,x0,xmu,rgrd
      real*8,intent(out) :: eimag,vint
      integer,intent(inout) :: jumprm

!     work space
      real*8 ri(nrptx)
      real*8 Vnlm(0:lx,(lx+1)**2,2,0:nphx)
      real*8 gap_up(0:nphx,0:lx,(lx+1)**2), gap_down(0:nphx,0:lx,(lx+1)**2)
      real*8 vtot(251,0:nphx), vvalgs(251,0:nphx)
      real*8 vtot_tmp(251,0:nphx), vvalgs_tmp(251,0:nphx)
      integer lmax(0:nphx)
      real*8 edens(251,0:nphx), edenvl(251,0:nphx),dmagx(nrptx)
      real*8 dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      real*8 vtotph_sp(nrptx,2),vvalph_sp(nrptx,2)
      real*8,intent(in) :: dgc(251,30,0:nphx), dpc(251,30,0:nphx)
      real*8,intent(in) :: adgc(10,30,0:nphx), adpc(10,30,0:nphx)
      real*8 dgcn(nrptx,30), dpcn(nrptx,30)
      real*8 xrhoce(0:lx,nex,0:nphx)
      real*8 xrhoce_sp(0:lx,2,nex,0:nphx)
      real*8 xmrhoce(0:lx,nspx*(lx+1)**2,nex,0:nphx)
      real*8 xmrhoce_sp(0:lx,nspx*(lx+1)**2,2,nex,0:nphx)
      complex*16 xrhole(0:lx,nex,0:nphx)
      complex*16 xrhole_sp(0:lx,2,nex,0:nphx)
      complex*16 xmrhole(0:lx,nspx*(lx+1)**2,nex,0:nphx)
      complex*16 xmrhole_sp(0:lx,nspx*(lx+1)**2,2,nex,0:nphx)
      complex*16 ph(nex,ltot+1,0:nphx)
      complex*16 ph_sp(nex,ltot+1,2,0:nphx)
      complex*16 aph(nex,lx+1,(lx+1)**2,0:nphx)
      complex*16 aph_sp(nex,lx+1,(lx+1)**2,2,0:nphx)
      integer inclus(0:nphx)
      integer iph
      complex*16 em(nex)
      complex*16 eref(nex)
      complex*16 eref_sp(nex,2)
      character*30  fname
      character*512 slog
      ! added for implicit none
      integer jri,jri1,itmp,i,ie,im,il,ne3,iph0,idwopt,iss,lll,ios,ne1,mmm,is,i_opt,msapp,ik0,ne,lmaxsc
      real*8 vjump,tk,thetad,critcw,sig2g,aa,edge,enext,de

!     msapp=2 - G_c + FMS only (no paths added)
      msapp=2
      vint=0.0d0

  if (emin.lt.emax) then
     !c    plot DOS between emin and emax with eimag above real axis
     !if(master)call wlog('              LDOS calculation for specified grid')
     lmaxsc = lx

     ne = neldos
	 if (ne .gt. nex) then
	   ne = nex
	   write(slog,*) "Warning: neldos = ", neldos, "is larger than nex = ", nex, "."
	   call wlog(slog)
	 endif

     !KJ debugging        ne = min (401, nex)
     de = (emax-emin)/(ne-1)
     if (eimag.lt.0) eimag=3*de
     enext=emin
     do i=1,ne
        em(i) = enext + coni*eimag
        enext= enext + de
     enddo

     ! ik0 is a starting point for path filters in energy
     ik0 = ne-45
     edge = xmu

     if(ispace.eq.0) call kprep(em,ne,.true.)  !KJ for rec space code


      Vnlm(:,:,:,:)=0.0d0

      do i_opt=1,i_hubbard
        !if(master)call wlog ('Calculating energy and space dependent l-DOS.')
        do iph = 0, nph
           write(slog,30) iph
30         format('     potential type ', i2)
           !if(master)call wlog(slog)
           lmax(iph) = lx

            do is = 1, 2

               vint = vint_sp(is)
               vtot(:,iph)= 0.0d0
               vvalgs(:,iph)=0.0d0
               edens(:,iph)=0.0d0
               edenvl(:,iph)=0.0d0

               rmt(iph)=rmt_sp(iph,is)
               vtot(:,iph)=vtot_sp(:,iph,is)
               vvalgs(:,iph)=vvalgs_sp(:,iph,is)
               edens(:,iph)=edens_sp(:,iph,is)
               edenvl(:,iph)=edenvl_sp(:,iph,is)
               call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph), vint, rhoint, dx, rgrd, jumprm, vjump, ri, vtotph, dum, dmagx)
               if (mod(ixc,10) .ge.5) then
                  if (jumprm .gt. 0) jumprm = 2
                  call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), dmag(1,iph), vint, rhoint, dx, rgrd , jumprm, vjump, ri, vvalph, dum, dmagx)
                  if (jumprm .gt. 0) jumprm = 1
               endif
               call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)

               jri = (log(rmt(iph)) + x0) / rgrd + 2
               jri1 = jri+1
               eref(1) = vtotph(jri1)

               do ie =1, ne
                  eref_sp(ie,is) = eref(1)
                  eref(ie) = eref(1)
               enddo
               do i = 1, jri1
                   vtotph(i) = vtotph(i) - eref(1)
               enddo
               if (ixc.ge.5) then
                   do i = 1, jri1
                      vvalph(i) = vvalph(i) - eref(1)
                  enddo
               else
                  do i = 1, jri1
                     vvalph(i) = vtotph(i)
                  enddo
               endif

               do i=1,251
                  vtotph_sp(i,is)=vtotph(i)
                  vvalph_sp(i,is)=vvalph(i)
               enddo
               itmp = 0
               call rhol_h(iph,rgrd, x0, ri, ne, em, ixc, rmt(iph), rnrm(iph),vtotph, vvalph, dgcn, dpcn, eref(1), adgc(1,1,iph), adpc(1,1,iph), xrhole(0,1,iph), xmrhole(0,1,1,iph), xrhoce(0,1,iph), xmrhoce(0,1,1,iph),ph(1,1,iph),aph(1,1,1,iph), iz(iph), xion(iph), iunf, itmp, lmaxsc, xnval(1,iph), Vnlm(0,1,1,iph), i_opt,is, xmu)

               ne1 = ne
               ne3 = 0
               do ie=1,ne
               do il=0,lx
                  xrhole_sp(il,is,ie,iph)=xrhole(il,ie,iph)
                  xrhoce_sp(il,is,ie,iph)=xrhoce(il,ie,iph)
                  ph_sp(ie,il+1,is,iph)= ph(ie,il+1,iph)
                  do im=(il**2+1),(il+1)**2
                     xmrhole_sp(il,im,is,ie,iph)=xmrhole(il,im,ie,iph)
                     xmrhoce_sp(il,im,is,ie,iph)=xmrhoce(il,im,ie,iph)
                     aph_sp(ie,il+1,im,is,iph)=aph(ie,il+1,im,iph)
                  enddo
               enddo
               enddo
 
            enddo ! is=1,2

         enddo !iph=0,nph

         !call fms for a cluster around central atom
         if (lfms2.ne.0) then
            iph0 = 0
            call fmsdos_h(2, rfms2, lfms2, iph0, idwopt, tk,thetad,sig2g, lmaxph, nat, iphat, rat, inclus(0), ne, ne1, ne3, nph, em, eref, eref_sp, iz, ph_sp, aph_sp, minv, rdirec, toler1, toler2,Vnlm,i_opt,xmu)
            inclus(1:nph) = inclus(0)
         else
            do iph0 = 0, nph
               call fmsdos_h(2, rfms2, lfms2,iph0,idwopt,tk,thetad,sig2g, lmaxph, nat, iphat, rat, inclus(iph0), ne, ne1, ne3, nph, em,eref, eref_sp, iz, ph_sp, aph_sp, minv, rdirec, toler1, toler2,Vnlm,i_opt,xmu)
            enddo
         endif


         do iph = 0,nph
            write (slog,'(a, i5)') ' Calculating chi and rho...', iph
            call wlog(slog)
            !KJ added if master to call ff2rho_h because I don't think it does anything for workers
            if (master) call ff2rho_h (critcw,xrhoce_sp(0,1,1,iph), xmrhoce_sp(0,1,1,1,iph),  xrhole_sp(0,1,1,iph), xmrhole_sp(0,1,1,1,iph),iph, msapp, em,lfms2, qnrm, xnmues, xmu, inclus,i_opt)

            if(master .and. i_opt.eq.1) then
               write(fname,"('hubbard', i2.2, '.dat')") iph
               open (unit=24, file=fname, status='old', iostat=ios)
               call chopen (ios,'hubbard.dat','phase')
               read(24,*)
               do iss=1,2
               do lll=0,lx
               do mmm=(lll**2+1),(lll+1)**2
                  read(24,*) aa,aa,aa,aa,Vnlm(lll,mmm,iss,iph)
               enddo
               enddo
               enddo
               close(24)
            endif

            do lll=0,lx
            do mmm=(lll**2+1),(lll+1)**2
                gap_up(iph,lll,mmm)= abs(Vnlm(lll,mmm,1,iph))
                gap_down(iph,lll,mmm)= abs(Vnlm(lll,mmm,2,iph))
            enddo
            enddo

         enddo !iph=0,nph


      enddo !i_opt=1,i_hubbard

	  call par_barrier

    endif

      return
      end


