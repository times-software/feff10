!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: axafs.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine axafs(em, emu, xsec,ne1,ik0)
!     extract axafs from xsec
!     written by a.l.ankudinov Dec. 1998

!     the file axafs.dat (format as in xmu.dat) will be written if
!     you use PRINT 0 1 0 0 0 0 (ipr2 > 0), and ran the second module.

!     the code draws a parabola using least mean square method
!     through xsec(i) * ee (i)**xn 
!     the weight for each point i, is defined as (ee(i)-E_F)**mm*
!     (ee(i+1)- ee(i-1)), where the last multiplier is used since the 
!     grid is not regular in energy.
!     E_F - energy that corresponds to Fermi level.

      use constants
      use DimsMod, only: nex

      implicit double precision (a-h, o-z)

      complex*16 em(nex), xsec(nex)
      dimension ee(nex), xmu(nex), wt(nex)
      dimension xx(0:4), yy(0:2), xm(3,3)
      external determ

!     empirically I found that the best curve is drawn if xn=0 and mm=1
!     alex ankudinov, january 1999.
      xn = 0
      mm = 1
      np = ne1 - ik0
      ef = emu

      do 10 ie = 1, np
        ee(ie) = dble(em(ik0+ie)-em(ik0)) +emu
        xmu(ie) = dimag(xsec(ik0+ie)) * ee(ie)**xn
  10  continue
      do 20 ie = 1, np
        if (ie.eq.1) then
          wt(ie) = (ee(ie+1)-ef) * (abs(ee(ie)-ef))**mm
        elseif (ie.eq.np) then
          wt(ie) = (ee(ie)-ee(ie-1)) * (ee(ie)-ef)**mm
        else
          wt(ie) = (ee(ie+1)-ee(ie-1)) * (ee(ie)-ef)**mm
        endif
  20  continue
      do 30 i = 0, 4
  30  xx(i) = 0
      do 40 i = 0, 2
  40  yy(i) = 0

      do 100 ie = 1, np
         do 80 i = 0,4
  80     xx(i) = xx(i) + wt(ie)*ee(ie)**i
         do 90 i = 0,2
  90     yy(i) = yy(i) + wt(ie)*xmu(ie)*ee(ie)**i
 100  continue

      do 105 i=1,3
      do 105 j=1,3
 105  xm(i,j) = xx(i+j-2)
      denom = determ (xm, 3, 3)

      do 110 i=1,3
      do 110 j=1,3
 110  xm(i,j) = xx(i+j-2)
      do 120 i=1,3
 120  xm(i,1) = yy (i-1)
      aa = determ (xm,3,3)
      aa = aa / denom

      do 210 i=1,3
      do 210 j=1,3
 210  xm(i,j) = xx(i+j-2)
      do 220 i=1,3
 220  xm(i,2) = yy (i-1)
      bb = determ (xm,3,3)
      bb = bb / denom

      do 310 i=1,3
      do 310 j=1,3
 310  xm(i,j) = xx(i+j-2)
      do 320 i=1,3
 320  xm(i,3) = yy (i-1)
      cc = determ (xm,3,3)
      cc = cc / denom

!     find normalization at edge+100 eV
      eee = ee(1) + 100/hart
      xnorm = (aa+bb*eee+cc*eee**2) / eee**xn

      open (unit=1,file='axafs.dat', status='unknown')
      write (1,*) '# File contains AXAFS. See manual for details.'
      write (1,*)                                                       &
     & '#--------------------------------------------------------------'
      write(1,*) '#  e, e(wrt edge), k,',                               &
     &           ' mu_at=(1+chi_at)*mu0_at, mu0_at, chi_at @#'
      do 400 ie = 1, np
        xmu(ie) = dimag(xsec(ie+ik0))
        xmu0 = (aa+bb*ee(ie)+cc*ee(ie)**2) / ee(ie)**xn
        chiat = (xmu(ie) - xmu0) / xmu0
        eee = ee(ie) -ef
        if (eee.ge.0.d0) then
           xk = sqrt(2*eee) /bohr
        else
           xk = -sqrt(-2*eee) /bohr
        endif
        write (1, 410) ee(ie)*hart, (ee(ie)-emu)*hart, xk,              &
     &              xmu(ie)/xnorm, xmu0/xnorm, chiat
 410    format (1x, 2f11.3, f8.3, 1p, 3e13.5)
 400  continue
      close (unit=1)

      return
      end
         

