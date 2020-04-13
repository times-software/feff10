
  subroutine ff2rho_h_step2( xmrhoce, xmrhole, ne, iph, em, qnrm, xnmues, xmu, inclus, gtr_m )

      use par, only: master
      use constants, only: hart
      use DimsMod, only: lx, nex, nspx=>nspu, nphx=>nphu

      implicit none
!     This routine simply writes the lm-dos to file.  That's all it does.
      integer, intent(in) :: iph
      real*8, intent(in) ::  xmu
      real*8, intent(out) ::   xmrhoce(0:lx,nspx*(lx+1)**2,2,nex)
      complex*16, intent(in) ::  xmrhole(0:lx, nspx*(lx+1)**2, 2, nex)
      complex*16, intent(in) :: em(nex)
      real*8, intent(in) :: qnrm (0:nphx)
      complex, intent(in) :: gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex)
      integer, intent(in) :: inclus(0:nphx)
      real*8, intent(in) :: xnmues (0:lx,0:nphx)
      integer, intent(in) :: ne
      character*30 fname
      character*2, parameter :: comment='# '
      integer i, j, k1, k2
      integer ios, ipp, il, ie, im, ne1, ne3, im1, im2, is, ip, ik

      
      !KJ Only master writes this output
      if (.not.master) return

!     Open lmdos.dat (output) and start header
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

!     Open rhocm.dat(output)
       write(fname,986) iph
       986   format('rhocm', i2.2, '.dat')
       open (unit=14, file=fname,status='unknown', iostat=ios)
       call chopen (ios, 'rhocm.dat', 'ff2rho')

       ! write  l-dos to 'ldosNN.dat'
       do ik = 1, ne
          write(14,1936) dble(em(ik))*hart, (((xmrhoce(il,im,is,ik), im=(il**2)+1,(il+1)**2),il=0,lx),is=1,2)
             do is = 1,2
             do il = 0,lx
                do im = (il**2)+1, (il+1)**2
                   xmrhoce(il,im,is,ik)=xmrhoce(il,im,is,ik)/(2*il+1)+ dimag(gtr_m(il,im,is,iph,ik)*xmrhole(il,im,is,ik))
                enddo
             enddo
             enddo
          write(15,1936) dble(em(ik))*hart,(((xmrhoce(il,im,is,ik), im=(il**2)+1,(il+1)**2),il=0,lx),is=1,2)
          1936    format (1x, f10.4, 1x, 33(1pe13.6,1x))
       enddo !ik

       close (unit=15)
       close (unit=14)


      return
      end

