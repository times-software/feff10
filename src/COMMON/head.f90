!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: head.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/03/22 18:59:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sthead (ntitle, title, nph, iz, rmt, rnrm,             &
     &                  xion, ihole, ixc,                               &
     &                  vr0, vi0, gamach, xmu, xf, vint, rs,            &
     &                  nohole, lreal,  rgrd)

!     SeT HEAD
!     This routine makes the file header, returned in head array.
!     header lines do not include a leading blank.
!     Last header line is not --------- end-of-header line

!     title lines coming into sthead include carriage control, since
!     they were read from potph.bin


      use dimsmod, only: nphx=>nphu, nheadx
      use constants
      implicit double precision (a-h, o-z)
! Josh Kas - Changed array dimensions from 30 to 40 (and others) for high Z elements
! according to Pavlo Baranov's changes.
      include '../HEADERS/vers.h'

      dimension xion(0:nphx)
      dimension iz(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(nheadx), store
      character*16 s1, s2

      integer istrln
      external istrln
      character*10 shole(0:40)
      character*8  sout(0:7)
      data shole /'no hole',   'K  shell',  'L1 shell',  'L2 shell',    &
     &            'L3 shell',  'M1 shell',  'M2 shell',  'M3 shell',    &
     &            'M4 shell',  'M5 shell',  'N1 shell',  'N2 shell',    &
     &            'N3 shell',  'N4 shell',  'N5 shell',  'N6 shell',    &
     &            'N7 shell',  'O1 shell',  'O2 shell',  'O3 shell',    &
     &            'O4 shell',  'O5 shell',  'O6 shell',  'O7 shell',    &
     &            'O8 shell',  'O9 shell',  'P1 shell',  'P2 shell',    &
     &            'P3 shell',  'P4 shell',  'P5 shell',  'P6 shell',    &
     &            'P7 shell',  'R1 shell',  'R2 shell',  'R3 shell',    &
     &            'R4 shell',  'R5 shell',  'S1 shell',  'S2 shell',    &
     &            'S3 shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 'DH - HL ',        &
     &           'DH + HL ', 'val=s+d ', 'sigmd(r)', 'sigmd=c '/


!     Fills head arrray, n = number of lines used.
!     Does not include line of dashes at the end.

      if (ntitle .ge. 1 ) then
         ii = istrln(title(1)) 
         if (ii.gt.1)  then
            write(store,100)  title(1)(1:), vfeff
         else
            write(store,102)  vfeff
         endif
      else
         write(store,102)   vfeff
      endif
  100 format( a55, t66, a12)
  102 format( t66, a12)
      title(1) = store
      nstor = 1

!     remove empty title lines
      do 120  ititle = 2, ntitle
         ii = istrln ( title (ititle) ) 
         if (ii.le.1)  goto 120
         nstor = nstor+1
         title(nstor) = title (ititle)
  120 continue
      ntitle = nstor

!     add more title lines
      if (xion(0) .ne. 0)  then
         ntitle = ntitle + 1
         write(title(ntitle),130)  iz(0), rmt(0)*bohr,                  &
     &                    rnrm(0)*bohr, xion(0), shole(ihole)
      else
         ntitle = ntitle + 1
         write(title(ntitle),140)  iz(0), rmt(0)*bohr,                  &
     &                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2, ' Rmt=',f6.3, ' Rnm=',f6.3,                 &
     &       ' Ion=',f5.2,  1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)
!     if (nohole.ge.0)  then
!        ntitle = ntitle + 1
!        write(title(ntitle),142)
! 142    format ('Calculations done with no core hole.')
!     endif
      if (lreal.ge.1 .or. (abs(rgrd - 0.05) .gt. 1.0e-5)) then
        ntitle = ntitle + 1
        s1 = ' '
        if (lreal.gt.1)  then
!        write(title(ntitle),144)
! 144    format ('Calculations done using only real phase shifts.')
         s1 = 'RPHASES'
        elseif (lreal.eq.1) then
!        ntitle = ntitle + 1
!        write(title(ntitle),145)
! 145    format ('Calculations done using only real self energy.')
         s1 = 'RSIGMA'
        endif
        s2 = '  '
        if (abs(rgrd - 0.05) .gt. 1.0e-5)  then
         write(s2,146)  rgrd
  146    format ('  RGRID', f7.4)
        endif
        ilen = istrln(s1)
        title(ntitle) = s1(1:ilen) // s2
      endif

      do 150  iph = 1, nph
         if (xion(iph) .ne. 0)  then
            ntitle = ntitle + 1
            write(title(ntitle),160)  iph, iz(iph),  rmt(iph)*bohr,     &
     &           rnrm(iph)*bohr, xion(iph)
         else
            ntitle = ntitle + 1
            write(title(ntitle),170)  iph, iz(iph),  rmt(iph)*bohr,     &
     &           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',f5.2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
  
  !KJ 3-2012 modified so that vr0,vi0 are always written to file.  It just seems that will be clearer for everyone.
  
 !     if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         ntitle = ntitle + 1
         write(title(ntitle),180)  gamach*hart, sout(ixc), vi0*hart,    &
     &                           vr0*hart
  !    else
  !       ntitle = ntitle + 1
  !       write(title(ntitle),190)  gamach*hart, sout(ixc)
  !    endif
      ntitle = ntitle + 1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  !190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, 'eV kf=',1pe9.3, ' Vint=',1pe10.3,            &
     &        'eV Rs_int=',0pf6.3)
      write(title(ntitle),200)  xmu*hart, xf/bohr, vint*hart, rs

      return
      end

      subroutine wthead (io, ntitle, title)
!     Dump title lines to unit io, which must be open. 
      integer io, i, ll
      character*80 title(ntitle)

!     nice for UNIX to use with gnuplot etc.,
      do 310 i = 1, ntitle
         ll = istrln(title(i))
         write(io,300)  title(i)(1:ll)
  300    format ('# ',a)  !KJ added this 7-09 to make xsect.dat gnu-plottable 
  310 continue

      return
      end
