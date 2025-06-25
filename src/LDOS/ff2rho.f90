!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2rho.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ff2rho (critcw, ne, xrhoce, xrhole, iph, msapp, em, lfms, qnrm, xnmues, xmu, inclus, ldostype)

  use DimsMod, only: nphx=>nphu, nex, lx
  use controls,only : fullpot
  use constants
  use par
  implicit none

!  implicit double precision (a-h, o-z)
  !     the  output is l-dos in xrhoce

  ! Inputs
  integer, intent(in) :: iph,msapp,lfms
  real*8,  intent(in) :: critcw ! not used
  real*8,  intent(in) :: qnrm(0:nphx), xnmues (0:lx,0:nphx),xmu
  integer, intent(in) :: inclus(0:nphx)
  integer, intent(in) :: ldostype
  complex*16, intent(in) :: xrhole(0:lx, nex)
  complex*16, intent(in) :: em(nex)

  ! In/Out
  real*8, intent(inout) :: xrhoce(0:lx,nex)

  ! Outputs
  integer, intent(out) :: ne

  ! Local variables
  integer nph
  integer nlgtr !KJ for m-DOS

  character*30 fname
  character*2,parameter :: coment='# '

  ! Added to satisfy implicit none
  integer :: i,ie,ik,il,im,ipp,ix, ios,iph0,ifms
  integer :: j,ne1,ne3

  ! Allocated variables
  complex,    allocatable :: gtr(:,:,:) !KJ  gtr(0:lx, 0:nphx,nex)
  real*8,     allocatable :: lmdos(:,:)
  complex*16, allocatable :: cchi(:,:) !KJ cchi(0:lx,nex)

  !KJ allocate gtr dynamically ; for mt-potential, sum over m
  ! for full potential, save m-dependent information.

  ! CC, adding ldostype for controling l,m projected DOS
  if(fullpot.or.(ldostype > 0)) then
     nlgtr=(lx+1)**2
     allocate(lmdos(nlgtr,nex))
     lmdos=dble(0)
  else
     nlgtr=lx
  endif

  ! Allocate local variables:
  allocate(gtr(0:nlgtr,0:nphx,nex),cchi(0:nlgtr,nex))

  gtr = cmplx(0,0)
  cchi = dcmplx(0,0)

  if (msapp.eq.0 .or. msapp.eq.2) then
     !        read gtrNN.bin
     iph0 = 0
     if (lfms.le.0) iph0 = iph
     write(fname,705)  iph0
705  format('gtr', i2.2, '.bin')
     open (unit=3, file=fname, status='old',                        &
          &        access='sequential', form='unformatted', iostat=ios)
     read (3) ne, ne1, ne3, nph, ifms
     read (3) (((gtr(il,ipp,ie), il=0,nlgtr), ipp=0,nph), ie=1,ne) !KJ
     !         read (3) (((gtr(il,ipp,ie), il=0,lx), ipp=0,nph), ie=1,ne)
     if (iph0.eq.0) then
        close (unit=3)
     else
        !          clean up disk space, keep gtr00.bin since it may be used
        !          by ff2chi
        !          close (unit=3,status='delete')
        close (unit=3)
     endif
  else
     !        gtr.bin was not written. set gtr to zero
     gtr=cmplx(0,0)
  endif

  !     Open ldos.dat and rhoc.dat (output) and start header
  MASTER_IF: if (master) then
905  format('ldos', i2.2, '.dat')
     write(fname,905) iph
     open (unit=3, file=fname, status='unknown', iostat=ios)
     call chopen (ios, 'ldos.dat', 'ff2rho')

     write (3,915) coment, ' Fermi level (eV): ', xmu*hart
915  format( 2a, f7.3)
     write (3,915) coment, ' Charge transfer : ', qnrm(iph)
     write (3,915)  coment,                                           &
          &           '   Electron counts for each orbital momentum:'
935  format(a,6x,i1,3x,f8.3)
     do i=0,lx
        write (3,935) coment, i, xnmues(i,iph)
     enddo

     write (3,975) coment, inclus(iph)
975  format (a, ' Number of atoms in cluster: ', i3)

     write(fname,985) iph
985  format('rhoc', i2.2, '.dat')
     open (unit=4, file=fname, status='unknown', iostat=ios)
     call chopen (ios, 'rhoc.dat', 'ff2rho')

     !     chi from fms is contained in gtr
     do i = 1, ne
        do j = 0, nlgtr !KJ lx
           cchi(j,i) =  dcmplx(gtr(j,iph,i))
        enddo
     enddo


     !     Write it out
     write(3,1805) coment, dimag(em(1))*hart
1805 format (a,' Lorentzian broadening with HWHH ', f10.4, ' eV')
     write(3,1815) coment
1815 format (a, 71('-'))
     write(3,1825) coment
1825 format(a,'     e        sDOS           pDOS          dDOS    ',   &
          &          '      fDOS    @#')


     if (ldostype > 0) then
  906  format('lmdos', i2.2, '.dat')
       write(fname,906) iph
       open (unit=67, file=fname, status='unknown', iostat=ios)
       call chopen (ios, 'lmdos.dat', 'ff2rho')
       write (67,915) coment, ' Fermi level (eV): ', xmu*hart
       write (67,915) coment, ' Charge transfer : ', qnrm(iph)
       write (67,915)  coment,                                           &
            &           '   Electron counts for each orbital momentum:'
  
       do i=0,lx
          write (67,935) coment, i, xnmues(i,iph)
       enddo
  
       write (67,975) coment, inclus(iph)
  
       ! Write it out
       write(67,1805) coment, dimag(em(1))*hart
       write(67,1815) coment
       if (ldostype == 1) then
         write(67,1826) coment
       elseif (ldostype == 2) then
         write(67,1827) coment
       endif
       1826 format(a, '     e          s            p_y           p_z      ', &
           '     p_x           d_xy          d_yz          d_z2     ', &
           '     d_xz         d_x2y2      f_y(3x2-y2)     f_xyz     ', &
           '    f_yz2          f_z3         f_xz2       f_z(x2-y2)  ', &
           '  f_x(x2-3y2)      @#')

       1827 format(a, '     e        Y_0^0         Y_1^-1        Y_1^0    ', &
           '     Y_1^1         Y_2^-2        Y_2^-1        Y_2^0  ', &
           '       Y_2^1         Y_2^2        Y_3^-3        Y_3^-2', &
           '        Y_3^-1         Y_3^0         Y_3^1         Y_3^2', &
           '         Y_3^3         @#')
     endif


     !      write  l-dos to 'ldosNN.dat'
     do ik = 1, ne
        write(4,1935)  dble(em(ik))*hart, (xrhoce(il,ik),il=0,lx)
        if (msapp.ne.1) then
           if(fullpot.or.(ldostype > 0)) then
              do il = 0,lx
                 do im = 1,2*il+1
                    ix=il**2
                    lmdos(ix+im,ik)= xrhoce(il,ik) / &
                         &  dble(2*il+1)+dimag(cchi(ix+im,ik)*xrhole(il,ik))
                 enddo
              enddo
              xrhoce(:,ik)=dble(0)
              do il=0,lx
                 do im=1,2*il+1
                    ix=il**2
                    xrhoce(il,ik)=xrhoce(il,ik)+lmdos(ix+im,ik)
                 enddo
              enddo
           else
              do il = 0,lx
                 xrhoce(il,ik)=xrhoce(il,ik) + & 
                      &        dimag(cchi(il,ik)*xrhole(il,ik))
              enddo
           endif
        endif
        if(fullpot.or.(ldostype > 0)) then
           write(67,1934) dble(em(ik))*hart, lmdos(:,ik)
1934       format (1x, f10.4, 1x, 50(1pe13.6,1x))
        endif
        write(3,1935)  dble(em(ik))*hart, (xrhoce(il,ik),il=0,lx)
1935    format (1x, f10.4, 1x, 5(1pe13.6,1x))
     enddo
     close (unit=3)
     close (unit=4)
     close (unit=67)
  endif MASTER_IF

  ! Deallocate local variables
  deallocate(gtr,cchi)

  if(fullpot.or.(ldostype > 0)) deallocate(lmdos)
  return
end subroutine ff2rho
