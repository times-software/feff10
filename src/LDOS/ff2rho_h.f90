      subroutine ff2rho_h(critcw, xrhoce, xmrhoce, xrhole, xmrhole, &
                        iph, msapp, em, lfms, qnrm, xnmues, xmu, inclus,i_opt,   &
                        Vnlm, TFrm, TFrmInv )

      use par
      use constants
      use DimsMod, only: lx, nex, nspx=>nspu, nphx=>nphu
      use hubbard_inp, only: l_hubbard, u_hubbard, J_hubbard, fermi_shift

      !implicit double precision (a-h, o-z)
      implicit none
!     the  output is l-dos in xrhoce
      integer, intent(in) :: i_opt, iph, msapp, lfms
      real*8, intent(in) :: critcw, xmu
      real*8, intent(in) ::   xrhoce(0:lx,2, nex)
      real*8, intent(out) ::   xmrhoce(0:lx,nspx*(lx+1)**2,2,nex)
      complex*16, intent(in) ::  xrhole(0:lx, 2, nex)
      complex*16, intent(in) ::  xmrhole(0:lx, nspx*(lx+1)**2, 2, nex)
      complex*16, intent(in) :: em(nex)
      real*8, intent(in) :: qnrm (0:nphx), xnmues (0:lx,0:nphx)

      real*8   xrhoce_sc(0:lx,2, nex)
      real*8   xmrhoce_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,nex)
      real*8   occup(0:nphx,0:lx,(lx+1)**2,2)
      real*8   Vnlm(0:lx,(lx+1)**2,2,0:nphx)
      real*8   rho_lm(nex), enm(nex),onlm(600),en(600)
      real*8   occ(0:lx,nspx*(lx+1)**2,2)
      real*8   occ_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2)
      complex gtr(0:lx, 2, 0:nphx,nex)
      complex gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex)
      complex gtr_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,0:nphx,nex)
      integer inclus(0:nphx)
      real*8 xnnmues (0:lx,(lx+1)**2,0:nphx)
      complex*16 cchi(0:lx,2,nex)
      complex*16 mchi(0:lx,nspx*(lx+1)**2,2,nex)
      complex*16 mchi_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,nex)
      integer nph, ne 
      !double precision max_gap_up, max_gap_down
      character*30 fname
      character*2, parameter :: comment='# '
      complex TFrm(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx), TFrmInv(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
      real*8 J_stoner 
      real a_spin1((2*l_hubbard+1),(2*l_hubbard+1))
      real a_spin2((2*l_hubbard+1),(2*l_hubbard+1))
      real evi_spin1(2*l_hubbard+1)
      real evi_spin2(2*l_hubbard+1)
      real evr_spin1(2*l_hubbard+1) 
      real evr_spin2(2*l_hubbard+1) 
      integer indic(2*l_hubbard+1) 
      real t1, t2,re1,re2,dim1,dim2 
      real veci_spin1((2*l_hubbard+1),(2*l_hubbard+1))
      real veci_spin2((2*l_hubbard+1),(2*l_hubbard+1))
      real vecr_spin1((2*l_hubbard+1),(2*l_hubbard+1))
      real vecr_spin2((2*l_hubbard+1),(2*l_hubbard+1))
      character JOBZ, UPLO
      integer i, j, k1, k2, N, LDA, LWORK, INFO
      real WORK(21)
      complex ctmp
      integer ios, ipp, il, ie, im, ne1, ne3, im1, im2, iph0, is, ip, ik, imm, ifms
      real*8 xmuO, xsum2, xsum1, U_hubb, tmp, su, x1, y1, cout, del3, xsum, aver_n0

      n=0
      
      !KJ I don't think anything useful happens here for slave nodes
      ! Sure, we read some stuff from file, but they're all local arrays, so ...
      if (.not.master) return

      N = 2*l_hubbard+1
      LWORK = 3*N
      LDA = 2*l_hubbard+1
      JOBZ = 'V'
      UPLO = 'U'
      xmuO=xmu-fermi_shift/hart
      t1 = 23.0E+00
      t2 = 23.0E+00

!     initialize TFrm and TFrmInv by identity matrix
       do 600 ip = 0,nphx
         do 600 il = 0,lx
           do 600 is = 1,2
             do 600 im1=1,2*l_hubbard+1
               do 600 im2=1,2*l_hubbard+1
                 if(im1.eq.im2) then
                  TFrm(im1,im2,is,il,ip)=cmplx(1.0,0.0) 
                  TFrmInv(im1,im2,is,il,ip)=cmplx(1.0,0.0) 
                 else
                  TFrm(im1,im2,is,il,ip)=cmplx(0.0,0.0)  
                  TFrmInv(im1,im2,is,il,ip)=cmplx(0.0,0.0)  
                 endif
600     continue

      if (msapp.eq.0 .or. msapp.eq.2) then  ! it is always =2
!        read gtrNN.bin

         iph0 = 0
         if (lfms.le.0) iph0 = iph
         write(fname,705)  iph0
705      format('gtr', i2.2, '.bin')
         open (unit=3, file=fname, status='old', access='sequential', form='unformatted', iostat=ios)
         read (3) ne, ne1, ne3, nph, ifms
         read (3) ((((gtr(il,is,ipp,ie), il=0,lx), ipp=0,nph), ie=1,ne),is=1,2)

         if (lfms.le.0) iph0 = iph
         write(fname,715)  iph0
715      format('gtr_m', i2.2, '.bin')

         open (unit=13, file=fname, status='old', access='sequential', form='unformatted', iostat=ios)
         read (13) ne, ne1, ne3, nph, ifms
         read (13) (((((gtr_m(il,im,is,ipp,ie), im=(il**2)+1,(il+1)**2),il=0,lx),ipp=0,nph), ie=1,ne),is=1,2)

         if (lfms.le.0) iph0 = iph
         write(fname,725)  iph0
725      format('gtr_off', i2.2, '.bin')


         open (unit=31, file=fname, status='old', access='sequential', form='unformatted', iostat=ios)
         read (31) ne, ne1, ne3, nph, ifms
         read(31)((((((gtr_off(il,im1,im2,is,ipp,ie)  ,im1=1,(l_hubbard+1)**2),im2=1,(l_hubbard+1)**2),&
                       ipp=0,nph),ie=1,ne),is=1,2),il=0,lx)

         if (iph0.eq.0) then
           close (unit=3)
           close (unit=13)
           close (unit=31)
         else
           close (unit=3)
           close (unit=13)
           close (unit=31)
         endif

      else  ! never happens
!        gtr.bin was not written. set gtr to zero
         gtr(:,:,:,:)=cmplx(0)
         gtr_m(:,:,:,:,:)=cmplx(0)
         gtr_off(:,:,:,:,:,:)=cmplx(0)
      endif

!     Open ldos.dat,lmdos.dat and rhoc.dat an rhocm.dat(output) and start header
      if (master) then
           if(i_opt.eq.1) then
               905   format('ldos', i2.2, '.dat')
               write(fname,905) iph
               open (unit=3, file=fname, status='unknown', iostat=ios)
               call chopen (ios, 'ldos.dat', 'ff2rho')
               write (3,915) comment, ' Fermi level (eV): ', xmu*hart
               915   format( 2a, f7.3)
               write (3,915) comment, ' Charge transfer : ', qnrm(iph)
               write (3,915)  comment, '   Electron counts for each orbital momentum:'
               935   format(a,6x,i1,3x,f8.3)
               do i=0,lx
                  write (3,935) comment, i, xnmues(i,iph)
               enddo
               write (3,975) comment, inclus(iph)
               975   format (a, ' Number of atoms in cluster: ', i3)
           endif

           976   format('lmdos', i2.2, '.dat')
           write(fname,976) iph
           open (unit=15, file=fname, status='unknown', iostat=ios)
           call chopen (ios, 'lmdos.dat', 'ff2rho')
           write (15,977) comment, ' Fermi level (eV): ', xmu*hart
           977   format( 2a, f7.3)
           write (15,977) comment, ' Charge transfer : ', qnrm(iph)
           write (15,977)  comment, 'Electron counts for each magnetic-orbital momentum:'
           978   format(a,6x,i1,3x,i3,3x,f8.3)
           do i=0,lx
           do im=(i**2)+1,(i+1)**2
              write (15,978) comment, i,im,xnmues(i,iph)
!KJ            The original instruction writes xnnmues rather than xnmues.  However, this is an uninitialized variable.
!              I am guessing it's a typo for "xnmues".  It's only a header anyway.
!              write (15,978) comment, i,im,xnnmues(i,im,iph)
           enddo
           enddo
           write (15,980) comment, inclus(iph)
           980   format (a, ' Number of atoms in cluster: ', i3)
           if(i_opt.eq.1) then
               write(fname,985) iph
               985   format('rhoc', i2.2, '.dat')
               open (unit=4, file=fname, status='unknown', iostat=ios)
               call chopen (ios, 'rhoc.dat', 'ff2rho')
           endif
           write(fname,986) iph
           986   format('rhocm', i2.2, '.dat')
           open (unit=14, file=fname,status='unknown', iostat=ios)
           call chopen (ios, 'rhocm.dat', 'ff2rho')

           ! chi from fms is contained in gtr
           do is=1,2
           do i = 1, ne
           do j = 0, lx
              cchi(j,is,i) =  dble( real( gtr(j,is,iph,i) ))  +coni* dble(aimag( gtr(j,is,iph,i) ))
           enddo
           enddo
           enddo
           do is=1,2
           do i = 1, ne
           do j = 0, lx
           do im = (j**2)+1,(j+1)**2
              mchi(j,im,is,i) =  dble( real( gtr_m(j,im,is,iph,i) ))  + coni * dble(aimag( gtr_m(j,im,is,iph,i) ))
           enddo
           enddo
           enddo
           enddo
           do is=1,2
           do ie = 1, ne
           do j = 0,lx
           do im1 = 1, (lx+1)**2
           do im2 = 1, (lx+1)**2
            mchi_off(j, im1,im2,is,ie) = dble(real(gtr_off(j,im1,im2,is,iph,ie)))  + coni * dble(aimag(gtr_off(j,im1,im2,is,iph,ie)))
           enddo
           enddo
           enddo
           enddo
           enddo

           ! Write it out
           if(i_opt.eq.1) then
               write(3,1805) comment, dimag(em(1))*hart
               1805    format (a,' Lorentzian broadening with HWHH ', f10.4, ' eV')
               write(3,1815) comment
               1815  format (a, 71('-'))
               write(3,1825) comment
               1825 format(a,'     e        sDOS(up)   pDOS(up)      dDOS(up) ', &
                       '   fDOS(up)   sDOS(down)    pDOS(down)   dDOS9(down)',&
                       '   fDOS(down)    @#')
           endif
           write(15,1835) comment, dimag(em(1))*hart
           1835    format (a,' Lorentzian broadening with HWHH ', f10.4, ' eV')
           write(15,1845) comment
           1845  format (a, 71('-'))
           write(15,1855) comment
           1855 format(a,'     e      s(0)DOS-up   p(-1)DOS-up   p(0)DOS-up ', &
                '  p(+1)DOS-up  d(-2)DOS-up  d(-1)DOS-up   d(0)DOS-up ', &
                '  d(+1)DOS-up  d(+2)DOS-up  f(-3)DOS-up   f(-2)DOS-up ',&
                '  f(-1)DOS-up  f(0)DOS-up   f(+1)DOS-up   f(+2)DOS-up ',&
                '  f(+3)DOS-up  s(0)DOS-dn   p(-1)DOS-dn   p(0)DOS-dn ', &
                '  p(+1)DOS-dn  d(-2)DOS-dn  d(-1)DOS-dn   d(0)DOS-dn ', &
                '  d(+1)DOS-dn  d(+2)DOS-dn  f(-3)DOS-dn   f(-2)DOS-dn ',&
                '  f(-1)DOS-dn  f(0)DOS-dn   f(+1)DOS-dn   f(+2)DOS-dn ',&
                '  f(+3)DOS-dn   @#')
           ! write  l-dos to 'ldosNN.dat'
           cout=0
           do ik = 1, ne
              cout=cout+1
              if(i_opt.eq.1) write(4,1935)  dble(em(ik))*hart, ((xrhoce(il,is,ik), il=0, lx),is=1,2)
              if(i_opt.eq.2) write(14,1936) dble(em(ik))*hart, (((xmrhoce(il,im,is,ik), im=(il**2)+1,(il+1)**2),il=0,lx),is=1,2)
              if (msapp.ne.1) then  ! always true
                 do is = 1,2
                 do il = 0,lx
                    do im = (il**2)+1, (il+1)**2
                       if(i_opt.eq.1)  xmrhoce(il,im,is,ik)=xrhoce(il,is,ik)/(2*il+1)+ &
                         dimag(mchi(il,im,is,ik)*xrhole(il,is,ik))
                       if(i_opt.eq.2)   xmrhoce(il,im,is,ik)=xmrhoce(il,im,is,ik)/(2*il+1)+ &
                         dimag(mchi(il,im,is,ik)*xmrhole(il,im,is,ik))
                    enddo
                    if(i_opt.eq.1) then
                       if(il.eq.l_hubbard) then
                           xmrhoce_off(:,:,:,is,ik)=0.0d0
                           do im1=1,(lx+1)**2
                              xmrhoce_off(il,im1,im1,is,ik)=xrhoce(il,is,ik)/(2*il+1)
                              do im2=1,(lx+1)**2
                                 xmrhoce_off(il,im1,im2,is,ik)=xmrhoce_off(il,im1,im2,is,ik) &
                                   + dimag(mchi_off(il,im1,im2,is,ik)*xrhole(il,is,ik))
                              enddo
                           enddo
                       endif
                    endif
                    xrhoce_sc(il,is,ik)=xrhoce(il,is,ik)+ dimag(cchi(il,is,ik)*xrhole(il,is,ik))
                 enddo
                 enddo
              endif
              if(i_opt.eq.1) then
                 write(3,1935)  dble(em(ik))*hart, ((xrhoce_sc(il,is,ik),il=0, lx),is=1,2)
                 1935    format (1x, f10.4, 1x, 9(1pe13.6,1x))
              endif
              write(15,1936) dble(em(ik))*hart,(((xmrhoce(il,im,is,ik), &
                  im=(il**2)+1,(il+1)**2),il=0,lx),is=1,2)
              1936    format (1x, f10.4, 1x, 33(1pe13.6,1x))
           enddo !ik

           if(i_opt.eq.1) close (unit=3)
           close (unit=4)
           close (unit=15)
           close (unit=14)
           ! Josh - set xmu0 by hand here to test things.
           ! xmu0 = -15.d0/hart
           del3=dble(xmuO-em(1))/600.0d0
           do il = 0,lx
            do is = 1,2
               do im = (il**2)+1, (il+1)**2
                  do ik=1,ne
                     rho_lm(ik)=xmrhoce(il,im,is,ik)
                     enm(ik)=dble(em(ik))
                  enddo
                  en(1)=dble(em(1))
                  do j=1,600
                     en(j)=en(1)+(j-1)*del3
                     x1=en(j)
                     call terp(enm,rho_lm,ne,3,x1,y1)
                     onlm(j)=y1
                  enddo
                  call  trap(en,onlm,600,su)
                  occ(il,im,is)=su*hart/2.0d0
               enddo
               if(il.eq.l_hubbard.and.iph.eq.1) then
                  do im1 = 1,(lx+1)**2
                  do im2 = 1,(lx+1)**2
                     do ik=1,ne
                        rho_lm(ik)=xmrhoce_off(il,im1,im2,is,ik)
                        enm(ik)=dble(em(ik))
                        IF((im1.EQ.5).AND.(im2.EQ.6)) write(97,*) enm(ik),rho_lm(ik)
                     enddo
                     en(1)=dble(em(1))
                     do j=1,600
                        en(j)=en(1)+(j-1)*del3
                        x1=en(j)
                        call terp(enm,rho_lm,ne,3,x1,y1)
                        onlm(j)=y1
                        if((iph.EQ.1).AND.(il.EQ.2).AND.(im1.GE.5).AND.(im2.GE.5))   write(98,*) en(j), onlm(j)
                     enddo
                     if((iph.EQ.1).AND.(il.EQ.2).AND.(im1.GE.5).AND.(im2.GE.5))   write(98,*)
                     call  trap(en,onlm,600,su)
                     occ_off(il,im1,im2,is)=su*hart/2.0d0
                  enddo !im2
                  enddo !im1
               endif
            enddo !is
           enddo !il

           !Make occ_off symmetric. It should be anyway, but isn't
           !possibly due to numerical issues. - J.K.
           do im1 = 1, (2*l_hubbard+1)
              do im2 = im1+1, (2*l_hubbard+1)
                 occ_off(l_hubbard,l_hubbard**2+im1,l_hubbard**2+im2,1) = occ_off(l_hubbard,l_hubbard**2+im2,l_hubbard**2+im1,1)
                 occ_off(l_hubbard,l_hubbard**2+im1,l_hubbard**2+im2,2) = occ_off(l_hubbard,l_hubbard**2+im2,l_hubbard**2+im1,2)
              end do
           end do
               
           if(i_opt.eq.1.and.iph.eq.1) then
                do im1=1,(2*l_hubbard+1)
                 do im2=1,(2*l_hubbard+1)
                    a_spin1(im1,im2) = real(occ_off(l_hubbard,l_hubbard**2+im1,l_hubbard**2+im2,1))
                    a_spin2(im1,im2) = real(occ_off(l_hubbard,l_hubbard**2+im1,l_hubbard**2+im2,2))
                 enddo
                enddo
                CALL SSYEV(JOBZ, UPLO, N, a_spin1, LDA, evr_spin1, WORK, LWORK,INFO)
                CALL SSYEV(JOBZ, UPLO, N, a_spin2, LDA, evr_spin2, WORK, LWORK,INFO)
                do ip = 0,nphx
                   do il = 0,lx
                     if(il.eq.l_hubbard.and.ip.eq.1) then
                        do im1=1,(2*l_hubbard+1)
                           do im2=1,(2*l_hubbard+1)
                              re1=a_spin1(im2,im1)
                              dim1=veci_spin1(im2,im1)             
                              re2=a_spin2(im2,im1)
                              dim2=veci_spin2(im2,im1)             
                              TFrm(im1,im2,1,il,ip)=cmplx(re1,0.0) 
                              TFrm(im1,im2,2,il,ip)=cmplx(re2,0.0) 
                              TFrmInv(im2,im1,1,il,ip)=cmplx(re1,0.0) 
                              TFrmInv(im2,im1,2,il,ip)=cmplx(re2,0.0)
                           enddo
                        enddo
                     endif
                   enddo
                enddo

                ! Check tranformation matrices.
                do im1 = 1, (2*l_hubbard+1)
                  do im2 = 1, (2*l_hubbard+1)
                     tmp = 0.0
                     do im = 1, (2*l_hubbard+1)
                        tmp = tmp + a_spin1(im, im1)* a_spin1(im,im2)
                     enddo
                     IF(im1.EQ.im2) THEN
                      IF(ABS(tmp-1.d0).GT.1.d-4) PRINT*, REAL(tmp),im1,is, &
                         'Problem with transformation matrix in ff2rho'
                     ELSE
                        IF(ABS(tmp).GT.1.d-4) PRINT*, tmp, im1, im2, is, & 
                         'Problem with transformation matrix in ff2rho'
                     END IF
                  enddo
                enddo

                open(file='transf.dat',unit=63,status='unknown',iostat=ios)
                open(file='Invtransf.dat',unit=64,status='unknown',iostat=ios)
                do is=1,2
                 do ip=0,nphx
                  do il=0,lx
                   do im1=1,(2*l_hubbard+1)
                      write(63,*)(TFrm(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
                      write(64,*)(TFrmInv(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
                   enddo
                  enddo
                 enddo
                enddo
                close(63)
                close(64)
           endif

           ! assigning the diagonalized density matrix,n_mm' to the d-states of metal atoms.
           xsum=0.0d0
           do is=1,2
             do il=0,lx
              do im=(il**2)+1 , (il+1)**2
                 if(iph.eq.1.and.i_opt.eq.1.and.il.eq.l_hubbard) then
                    if (is.eq.1) occup(iph,il,im,is)=evr_spin1(im-il**2)
                    if (is.eq.2) occup(iph,il,im,is)=evr_spin2(im-il**2)
                 else
                    occup(iph,il,im,is)=occ(il,im,is)
                 endif
                 if(iph.eq.0.or.iph.eq.1) then
                    if(il.eq.l_hubbard) xsum=xsum+occup(iph,il,im,is)
                 endif
              enddo
             enddo
           enddo

           U_hubb=U_hubbard/hart
           J_stoner=J_hubbard/hart
           aver_n0=xsum/real(4*l_hubbard+2)
           do is=1,2
             do il=0,lx
              do im=(il**2)+1,(il+1)**2
                 if(iph.eq.0.or.iph.eq.1) then
                    if (il.eq.l_hubbard) then
                       if(is.eq.1) then
                          xsum1=0.0d0
                          do imm=(il**2)+1, (il+1)**2
                             xsum1=xsum1+U_hubb*(occup(iph,il,imm,is+1)-aver_n0)
                          enddo
                       elseif(is.eq.2) then
                          xsum1=0.0d0
                          do imm=(il**2)+1, (il+1)**2
                             xsum1=xsum1+U_hubb*(occup(iph,il,imm,is-1)-aver_n0)
                          enddo
                       endif
                       xsum2=0.0d0
                       do imm=(il**2)+1,(il+1)**2
                          if(imm.ne.im) then
                             xsum2=xsum2+(U_hubb-J_stoner)* (occup(iph,il,imm,is)-aver_n0)
                          endif
                       enddo
                       Vnlm(il,im,is,iph)=xsum1+xsum2
                    else
                       Vnlm(il,im,is,iph)=0.0d0
                    endif
                 else
                    Vnlm(il,im,is,iph)=0.0d0
                 endif
              enddo
             enddo
           enddo

           if(i_opt.eq.1) then
               ! write hubbard_iph.dat
               write(fname,2001) iph
               2001  format('hubbard', i2.2, '.dat')
               open (unit=24, file=fname, status='unknown', iostat=ios)
               call chopen (ios, 'hubbard.dat', 'ff2rho')
               write(24,2825) comment
               2825  format(a,'     l        m      s     n(lms)         V(lms)')
               do is = 1,2
                  do il = 0,lx
                     do im=(il**2)+1, (il+1)**2
                        write(24,*) '  ',il,'  ',im-(il*(il+1)+1),'   ',  &
                         is,'   ',occup(iph,il,im,is),'  ',Vnlm(il,im,is,iph)
                     enddo
                  enddo
               enddo
               close (unit=24)
           endif
      endif  ! if master

      return
      end

