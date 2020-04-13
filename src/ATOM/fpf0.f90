!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fpf0.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fpf0 ( iz, iholep, srho, dr, hx,                       &
     &     dgc0, dpc0, dgc, dpc,                                  &
     &     eatom, xnel, norb, eorb, kappa)
  !      everything is input. output is written in fpf0.dat
  !      to be read by ff2afs.f to get scattering amplitude

  use DimsMod, only: nphx=>nphu
  use constants
  use par
  implicit none

  !     save central atom dirac components, see comments below.

  integer, intent(in) :: iz,iholep, norb
  real*8, intent(in)  :: eatom,hx
  real*8, intent(in), dimension(251) :: dgc0, dpc0, srho,dr
  real*8, intent(in), dimension(30)  :: xnel, eorb !, kappa
  real*8, intent(in), dimension(251, 30, 0:nphx) :: dgc, dpc
  integer,intent(in), dimension(30) :: kappa

  real*8, dimension(251) :: xpc, xqc

  logical open_16

  !     output arrays
  real*8, dimension(13) :: enosc, oscstr !, index
  integer, dimension(13) :: index

  ! Added to satisfy implicit none
  integer :: i,iorb,iq,jkap,kdif,kinit,np,nosc,ios
  real*8  :: xk0,xj0,xirf,twoj,xmult1,xmult2,dq,fpcorr

  if (master) then
     open (unit=16, file='fpf0.dat', status='unknown', iostat=ios)
     fpcorr =  -(iz/82.5)**2.37
     write (16,*)  ' atom Z = ', iz
     write (16,10)  eatom *alphfs**2 *5/3, fpcorr,                   &
          &        ' total energy part of fprime - 5/3*E_tot/mc**2'
10   format (2(1pe19.5), a)
     open_16 = .true.
  else
     open_16 = .false.
  endif

  !     get oscillator strengths
  do i=1,13
     oscstr(i)=0.d0
     enosc(i)=0.d0
  enddo
  enosc(1)= eorb(iholep)
  index(1)= iholep
  kinit = kappa(iholep)
  oscstr(1) = 2*abs(kinit)
  !     always will use first spot to represent initial state
  nosc=1
  np = 251

  do iorb =1, norb
     if (xnel(iorb) .gt.0.d0) then
        !         it is core orbital, check if it satisfies dipole selection
        jkap = kappa(iorb)
        if (jkap+kinit.eq.0 .or. abs(jkap-kinit).eq.1) then
           nosc = nosc+1
           !            calculate reduced dipole matrix element
           kdif= jkap-kinit
           if (abs(kdif).gt.1) kdif=0
           !            xirf = <i |p| f> relativistic version of dipole m.e.
           !            from Grant,Advan.Phys.,v.19,747(1970) eq. 6.30, using
           !            Messiah's "Q.M." appendices to reduce 9j,3j symbols
           !            to simple coefficients xmult1,2. ala 12.12.95
           twoj = 2.0d0*abs(kinit) - 1.0d0
           if (kdif.eq.-1 .and. kinit.gt.0) then
              xmult1 = 0.0d0
              xmult2 = sqrt(2.0d0 * (twoj+1)*(twoj-1)/twoj )
           elseif (kdif.eq.-1 .and. kinit.lt.0) then
              xmult1 = 0.0d0
              xmult2 = - sqrt(2.0d0 * (twoj+1)*(twoj+3)/(twoj+2) )
           elseif (kdif.eq. 0 .and. kinit.gt.0) then
              xmult1 = - sqrt( (twoj+1)*twoj/(twoj+2) )
              xmult2 = - sqrt( (twoj+1)*(twoj+2)/twoj )
           elseif (kdif.eq. 0 .and. kinit.lt.0) then
              xmult1 = sqrt( (twoj+1)*(twoj+2)/twoj )
              xmult2 = sqrt( (twoj+1)*twoj/(twoj+2) )
           elseif (kdif.eq. 1 .and. kinit.gt.0) then
              xmult1 = sqrt(2.0d0 * (twoj+1)*(twoj+3)/(twoj+2) )
              xmult2 = 0.0d0
           elseif (kdif.eq. 1 .and. kinit.lt.0) then
              xmult1 = - sqrt(2.0d0 * (twoj+1)*(twoj-1)/twoj )
              xmult2 = 0.0d0
           endif
           xk0 = abs(eorb(iorb)-eorb(iholep)) * alphfs
           do i = 1, np
              xj0 = sin(xk0*dr(i))/(xk0*dr(i))
              xpc(i) = (xmult1*dgc0(i)*dpc(i,iorb,0)+                 &
                   &            xmult2*dpc0(i)*dgc(i,iorb,0)) * xj0
              xqc(i) = 0.0d0
           enddo
           !            xirf=lfin+linit+2
           xirf=2
           call somm (dr, xpc, xqc, hx, xirf, 0, np)
           oscstr(nosc) = xirf**2/3.0d0 
           enosc(nosc) = eorb(iorb)
           index(nosc) = iorb
        endif
     endif
  enddo

  !     write down information about oscillators
  if(open_16) then
     write(16, *) nosc
     do i=1,nosc
        write(16,220) oscstr(i), enosc(i), index(i)
220     format ( f9.5, f12.3, i4)
     enddo
  endif

  !     calculate and write out f0(Q) on grid delq=0.5 Angstorm**(-1)
  dq=0.5*bohr 
  do iq = 1,81
     xk0 = dq*(iq-1)
     !        srho is 4*pi*density 
     do i = 1, np
        xj0 = 1.d0
        if(iq.gt.1) xj0 = sin(xk0*dr(i))/(xk0*dr(i))
        xpc(i) = srho(i) * (dr(i)**2) *xj0
        xqc(i) = 0.d0
     enddo
     xirf = 2.d0
     call somm (dr, xpc, xqc, hx, xirf, 0, np)
     if (open_16) write (16, 570) 0.5*(iq-1), xirf
570  format ( f5.1, 1x, f9.4)
  enddo

  if (open_16) close(unit=16)

  return
end subroutine fpf0
