
   subroutine ff2rho_h_step1(xrhoce, xmrhoce, xrhole, &
                        ne, iph, em, qnrm, xnmues, xmu, inclus, Vnlm, TFrm, TFrmInv, &
                         gtr, gtr_m, gtr_off )
    ! Calculates Vnlm, TFrm, and TFrmInv
      use par
      use constants
      use DimsMod, only: lx, nex, nspx=>nspu, nphx=>nphu
      use hubbard_inp, only: l_hubbard, u_hubbard, J_hubbard, fermi_shift

      implicit none
      integer, intent(in) :: iph
      real*8, intent(in) :: xmu
      real*8, intent(in) ::   xrhoce(0:lx,2, nex)
      complex*16, intent(in) ::  xrhole(0:lx, 2, nex)
      complex*16, intent(in) :: em(nex)
      complex, intent(in) :: gtr(0:lx, 2, 0:nphx,nex)
      complex, intent(in) :: gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex)
      complex, intent(in) :: gtr_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,0:nphx,nex)
      real*8, intent(in) :: qnrm (0:nphx), xnmues (0:lx,0:nphx)
      complex, intent(out) ::TFrm(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
      complex, intent(out) ::TFrmInv(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
      real*8, intent(out) :: Vnlm(0:lx,(lx+1)**2,2,0:nphx)
      integer, intent(in) :: ne
      real*8, intent(out) :: xmrhoce(0:lx,nspx*(lx+1)**2,2,nex)

      real*8 xrhoce_sc(0:lx,2, nex)
      real*8 xmrhoce_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,nex)
      real*8 occup(0:nphx,0:lx,(lx+1)**2,2)
      real*8 rho_lm(nex), enm(nex),onlm(600),en(600)
      real*8 occ(0:lx,nspx*(lx+1)**2,2)
      real*8 occ_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2)
      integer inclus(0:nphx)
      character*30 fname
      character*2, parameter :: comment='# '
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
      logical :: verbose = .false.

      !KJ I don't think anything useful happens here for slave nodes
      if (.not.master) return

      N = 2*l_hubbard+1
      LWORK = 3*N
      LDA = 2*l_hubbard+1
      JOBZ = 'V'
      UPLO = 'U'
      xmuO=xmu-fermi_shift/hart
      t1 = 23.0E+00
      t2 = 23.0E+00

      if(iph.eq.0) then
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
600   continue
      endif

      if(verbose) then
!          Open ldos.dat and start header
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
           write(3,1805) comment, dimag(em(1))*hart
           1805    format (a,' Lorentzian broadening with HWHH ', f10.4, ' eV')
           write(3,1815) comment
           1815  format (a, 71('-'))
           write(3,1825) comment
           1825 format(a,'     e        sDOS(up)   pDOS(up)      dDOS(up) ', &
                       '   fDOS(up)   sDOS(down)    pDOS(down)   dDOS9(down)',&
                       '   fDOS(down)    @#')

!          Open lmdos.dat and start header
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
    !KJ       The original instruction writes xnnmues rather than xnmues.  However, this is an uninitialized variable.
    !         I am guessing it's a typo for "xnmues".  It's only a header anyway.
    !         write (15,978) comment, i,im,xnnmues(i,im,iph)
           enddo
           enddo
           write (15,980) comment, inclus(iph)
           980   format (a, ' Number of atoms in cluster: ', i3)
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

!          Open rhoc.dat
           write(fname,985) iph
           985   format('rhoc', i2.2, '.dat')
           open (unit=4, file=fname, status='unknown', iostat=ios)
           call chopen (ios, 'rhoc.dat', 'ff2rho')

!          Open rhocm.dat(output)
           write(fname,986) iph
           986   format('rhocm', i2.2, '.dat')
           open (unit=14, file=fname,status='unknown', iostat=ios)
           call chopen (ios, 'rhocm.dat', 'ff2rho')

      endif  ! if verbose, write files


      cout=0
      do ik = 1, ne
          cout=cout+1
          if (verbose) write(4,1935)  dble(em(ik))*hart, ((xrhoce(il,is,ik), il=0, lx),is=1,2)
          do is = 1,2
          do il = 0,lx
             do im = (il**2)+1, (il+1)**2
                xmrhoce(il,im,is,ik)=xrhoce(il,is,ik)/(2*il+1)+ &
                  dimag(gtr_m(il,im,is,iph,ik)*xrhole(il,is,ik))
             enddo
             if(il.eq.l_hubbard) then
                xmrhoce_off(:,:,:,is,ik)=0.0d0
                do im1=1,(lx+1)**2
                   xmrhoce_off(il,im1,im1,is,ik)=xrhoce(il,is,ik)/(2*il+1)
                   do im2=1,(lx+1)**2
                      xmrhoce_off(il,im1,im2,is,ik)=xmrhoce_off(il,im1,im2,is,ik) &
                           + dimag(gtr_off(il,im1,im2,is,iph,ik)*xrhole(il,is,ik))
                   enddo
                enddo
             endif
             xrhoce_sc(il,is,ik)=xrhoce(il,is,ik)+ dimag(gtr(il,is,iph,ik)*xrhole(il,is,ik))
          enddo
          enddo
          if(verbose) write(3,1935)  dble(em(ik))*hart, ((xrhoce_sc(il,is,ik),il=0, lx),is=1,2)
          1935    format (1x, f10.4, 1x, 9(1pe13.6,1x))
          if(verbose) write(15,1936) dble(em(ik))*hart,(((xmrhoce(il,im,is,ik), &
          im=(il**2)+1,(il+1)**2),il=0,lx),is=1,2)
          1936    format (1x, f10.4, 1x, 33(1pe13.6,1x))
      enddo !ik

      if(verbose) then
           close(4)
           close(15)
           close(14)
           close(3)
      endif
      ! make occ and occ_off
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
                    !IF((im1.GE.10).AND.(im2.GE.10)) write(97,*) enm(ik),rho_lm(ik)
                 enddo
                 en(1)=dble(em(1))
                 do j=1,600
                    en(j)=en(1)+(j-1)*del3
                    x1=en(j)
                    call terp(enm,rho_lm,ne,3,x1,y1)
                    onlm(j)=y1
                    !if((im1.GE.10).AND.(im2.GE.10))   write(98,*) en(j), onlm(j)
                 enddo
                 !if((im1.GE.10).AND.(im2.GE.10))   write(98,*)
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

!      if(iph.eq.1)write(*,*) 'OCC_OFF',occ_off(l_hubbard,10:16,10:16,1)

      !KJ IMPORTANT!  This instruction "if iph=1" hardwires that the Hubbard term must go to atom of potential type 1.
      !This should be replaced by a more general instruction with a variable like "iat_hubbard".
      !For the moment, we cannot add the Hubbard term to iph=0 because then the central atom term is involved, and
      !we haven't derived/implemented the necessary theory.
      !But it is retarded that we now have to have the Hubbard term at iph=1 whereas any iph>0 should be okay, of course.

      ! calculate transformation matrices TFrm, TFrmInv:
      if(iph.eq.1) then
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
            !Write transformation matrices to file:
            if(verbose) then
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
      endif

      ! assigning the diagonalized density matrix,n_mm' to the d/f-states of metal atoms.
      ! calculate occupation numbers in occup
      xsum=0.0d0
      do is=1,2
         do il=0,lx
          do im=(il**2)+1 , (il+1)**2
             if(iph.eq.1.and.il.eq.l_hubbard) then
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

      ! Finally, calculate Hubbard potential Vnlm:
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

      ! write hubbard_iph.dat
      if(verbose) then
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

   return
   end



