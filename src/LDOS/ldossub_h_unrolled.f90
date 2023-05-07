
      subroutine ldos_h_unrolled ( nph, edens_sp, edenvl_sp, dmag, vtot_sp, vvalgs_sp, &
                  rmt,rmt_sp,  rnrm, ixc, rhoint, vint,vint_sp, xmu, jumprm,  &
                  x0, dx, rgrd, xion, iunf, iz, xnval, adgc, adpc, dgc,dpc,  &
                  ihole, qnrm, xnmues, emin, emax, eimag, rfms2, lfms2, lmaxph,  &
                  nat, iphat, rat, minv, rdirec, toler1, toler2)

      use controls,only: ispace
      use DimsMod, only: nphx=>nphu, lx, nrptx, nex, ltot, natx, nspx=>nspu
      use constants
      use par
      use ldos_inp,only: neldos
      use hubbard_inp,only: i_hubbard, l_hubbard

      implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      real*8,intent(in) :: vint_sp(2)
      real*8,intent(in) :: dmag(251,0:nphx)
      real*8,intent(in) :: vtot_sp(251,0:nphx,2), vvalgs_sp(251,0:nphx,2)
      real*8,intent(in) :: rnrm(0:nphx)
      real*8,intent(out) :: rmt(0:nphx)
      real*8,intent(in) :: rmt_sp(0:nphx,2) !, rnrm_sp(0:nphx,2)
      real*8,intent(in) :: xnval(41,0:nphx)
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
      real*8 gap_up(0:nphx,0:lx,(lx+1)**2), gap_down(0:nphx,0:lx,(lx+1)**2)
      real*8 vtot(251,0:nphx), vvalgs(251,0:nphx)
      real*8 vtot_tmp(251,0:nphx), vvalgs_tmp(251,0:nphx)
      real*8 edens(251,0:nphx), edenvl(251,0:nphx),dmagx(nrptx)
      real*8 dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      real*8 vtotph_sp(nrptx,2),vvalph_sp(nrptx,2)
      real*8,intent(in) :: dgc(251,41,0:nphx), dpc(251,41,0:nphx)
      real*8,intent(in) :: adgc(10,41,0:nphx), adpc(10,41,0:nphx)
      real*8 dgcn(nrptx,41), dpcn(nrptx,41)
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
      integer jri,jri1,i,ie,im,il,ne3,iph0,idwopt,iss,lll,ios,ne1,mmm,is,ik0,ne,j 
      real*8 vjump,tk,thetad,critcw,sig2g,aa,edge,enext,de
      ! added KJ to avoid communicating through files:
      real*8   Vnlm(0:lx,(lx+1)**2,2,0:nphx)
      complex TFrm(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx), TFrmInv(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
      complex, allocatable, dimension(:,:,:,:) :: gtr
      complex, allocatable, dimension(:,:,:,:,:) :: gtr_m
      complex, allocatable, dimension(:,:,:,:,:,:) :: gtr_off


     if (emin.ge.emax) return
write(*,*)'in ldossub_h_unrolled, I am ',master
     ne = neldos
	 if (ne .gt. nex) then
	   ne = nex
	   write(slog,*) "Warning: neldos = ", neldos, "is larger than nex = ", nex, "."
	   call wlog(slog)
	 endif
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

     Vnlm(:,:,:,:)=0.0d0

     allocate(gtr(0:lx, 2, 0:nphx, nex))
     allocate(gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex))
     allocate(gtr_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,0:nphx,nex))

 ! ********** FIRST PASS

        do iph = 0, nph
           do is = 1, 2

               vint = vint_sp(is)
               rmt(iph)=rmt_sp(iph,is)
               vtot(:,iph)=vtot_sp(:,iph,is)
               vvalgs(:,iph)=vvalgs_sp(:,iph,is)
               edens(:,iph)=edens_sp(:,iph,is)
               edenvl(:,iph)=edenvl_sp(:,iph,is)
               call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph), &
                            vint, rhoint, dx, rgrd, jumprm, vjump, ri, &
                            vtotph, dum, dmagx)
               ! if (mod(ixc,10) .ge.5) then
               if ((mod(ixc,10).ge.5).and.(ixc.ne.6).and.(ixc.ne.7)) then ! TTS
                  if (jumprm .gt. 0) jumprm = 2
                  call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), &
                               dmag(1,iph), vint, rhoint, dx, rgrd , &
                               jumprm, vjump, ri, vvalph, dum, dmagx)
                  if (jumprm .gt. 0) jumprm = 1
               endif
               call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)

               jri = (log(rmt(iph)) + x0) / rgrd + 2
               jri1 = jri+1
               eref(1) = vtotph(jri1)

               eref_sp(1:ne,is) = eref(1)
               eref(1:ne) = eref(1)
               vtotph(1:jri1) = vtotph(1:jri1) - eref(1)
               ! if (ixc.ge.5) then
               if ((ixc.ge.5).and.(ixc.ne.6).and.(ixc.ne.7)) then ! TTS
                  vvalph(1:jri1) = vvalph(1:jri1) - eref(1)
               else
                  vvalph(1:jri1) = vtotph(1:jri1)
               endif
               vtotph_sp(1:251,is)=vtotph(1:251)
               vvalph_sp(1:251,is)=vvalph(1:251)
write(*,*)'calling rhol_h_step1, I am ',master

               ! makes xrhole, xrhoce, ph :
               call rhol_h_step1(iph,rgrd, x0, ri, ne, em, ixc, rmt(iph), &
                 rnrm(iph),vtotph, vvalph, dgcn, dpcn, eref(1), &
                 adgc(1,1,iph), adpc(1,1,iph), xrhole(0,1,iph), &
                 xrhoce(0,1,iph), ph(1,1,iph), iz(iph), xion(iph), iunf, 0, &
                 lx, xnval(1,iph), is, xmu)

               do ie=1,ne
               do il=0,lx
                  xrhole_sp(il,is,ie,iph)=xrhole(il,ie,iph)
                  xrhoce_sp(il,is,ie,iph)=xrhoce(il,ie,iph)
                  ph_sp(ie,il+1,is,iph)= ph(ie,il+1,iph)
               enddo
               enddo
 
            enddo ! is=1,2
         enddo !iph=0,nph

         !call fms for a cluster around central atom
         !calculates gtr, gtr_m, and gtr_off
         if (lfms2.ne.0) then
            iph0 = 0
            call fmsdos_h_step1(rfms2, lfms2, iph0, lmaxph, nat, iphat, rat, &
              inclus(0), ne, nph, em, eref_sp, iz, ph_sp, rdirec, toler1, &
              toler2, gtr, gtr_m, gtr_off)
            inclus(1:nph) = inclus(0)
         else
            do iph0 = 0, nph
               call fmsdos_h_step1(rfms2, lfms2,iph0, lmaxph, nat, iphat, &
                 rat, inclus(iph0), ne, nph, em, eref_sp, iz, ph_sp, rdirec, &
                 toler1, toler2, gtr, gtr_m, gtr_off)
            enddo
         endif


         do iph = 0,nph
            write (slog,'(a, i5)') ' Calculating chi and rho...', iph
            call wlog(slog)
            ! calculates xmrhoce_sp, Vnlm, TFrm, and TFrmInv
            ! writes ldos.dat, lmdos.dat, rhocm.dat, transf.dat, invtransf.dat, hubbard.dat
            if (master) call ff2rho_h_step1 (xrhoce_sp(:,:,:,iph), &
              xmrhoce_sp(:,:,:,:,iph),  xrhole_sp(:,:,:,iph), ne, iph, &
              em, qnrm, xnmues, xmu, inclus, Vnlm, TFrm, TFrmInv, gtr, &
              gtr_m, gtr_off)
         enddo

 
 ! ********** SECOND PASS

        do iph = 0, nph
           do is = 1, 2

               vint = vint_sp(is)
               rmt(iph)=rmt_sp(iph,is)
               vtot(:,iph)=vtot_sp(:,iph,is)
               vvalgs(:,iph)=vvalgs_sp(:,iph,is)
               edens(:,iph)=edens_sp(:,iph,is)
               edenvl(:,iph)=edenvl_sp(:,iph,is)
               call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph), &
                 vint, rhoint, dx, rgrd, jumprm, vjump, ri, vtotph, dum, dmagx)
               if (mod(ixc,10) .ge.5) then
                  if (jumprm .gt. 0) jumprm = 2
                  call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), &
                    dmag(1,iph), vint, rhoint, dx, rgrd , jumprm, vjump, ri, &
                    vvalph, dum, dmagx)
                  if (jumprm .gt. 0) jumprm = 1
               endif
               call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)

               jri = (log(rmt(iph)) + x0) / rgrd + 2
               jri1 = jri+1
               eref(1) = vtotph(jri1)

               eref_sp(1:ne,is) = eref(1)
               eref(1:ne) = eref(1)
               vtotph(1:jri1) = vtotph(1:jri1) - eref(1)
               if (ixc.ge.5) then
                  vvalph(1:jri1) = vvalph(1:jri1) - eref(1)
               else
                  vvalph(1:jri1) = vtotph(1:jri1)
               endif
               vtotph_sp(1:251,is)=vtotph(1:251)
               vvalph_sp(1:251,is)=vvalph(1:251)

               ! makes  xmrhole, xmrhoce, aph :
               call rhol_h_step2(iph,rgrd, x0, ri, ne, em, ixc, rmt(iph), &
                 rnrm(iph),vtotph, vvalph, dgcn, dpcn, eref(1), adgc(:,:,iph), &
                 adpc(:,:,iph), xmrhole(:,:,:,iph), xmrhoce(:,:,:,iph), &
                 aph(:,:,:,iph), iz(iph), xion(iph), iunf, 0, lx, xnval(:,iph), &
                 Vnlm(:,:,:,iph), is)

               do ie=1,ne
               do il=0,lx
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
         !calculates gtr_m
         if (lfms2.ne.0) then
            iph0 = 0
            call fmsdos_h_step2(rfms2, lfms2, iph0, lmaxph, nat, iphat, rat, &
              inclus(0), ne, nph, em, eref_sp, aph_sp, rdirec, toler1, toler2, &
              TFrm, TFrmInv, gtr_m)
            inclus(1:nph) = inclus(0)
         else
            do iph0 = 0, nph
               call fmsdos_h_step2(rfms2, lfms2, iph0, lmaxph, nat, iphat, rat, &
                 inclus(iph0), ne, nph, em, eref_sp, aph_sp, rdirec, toler1, &
                 toler2, TFrm, TFrmInv, gtr_m)
            enddo
         endif

         if(master) then
             do iph = 0,nph
                write (slog,'(a, i5)') ' Calculating chi and rho...', iph
                call wlog(slog)
                ! Writes rhocm.dat and lmdos.dat :
                call ff2rho_h_step2 (xmrhoce_sp(:,:,:,:,iph), &
                  xmrhole_sp(:,:,:,:,iph), ne, iph, em, qnrm, xnmues, xmu, &
                  inclus, gtr_m)
             enddo

             ! pass Vnlm to XSPH:
             open(21, file='v_hubbard.bin',form='unformatted')
             write(21) Vnlm
             close(21)
             ! pass TFrm to FMS:
             open(63, file='transformation_hubbard.bin',form='unformatted')
             write(63) TFrm
             write(63) TFrmInv
             close(63)
         endif
         
	  call par_barrier
      deallocate(gtr, gtr_off, gtr_m)


      return
      end



