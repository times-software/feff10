!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sumax.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE SUMAX (RN, ANN, AA2, AASUM)
! This is a version of the subroutine sumax found on page 110 of
! Louck's book.  It performs eq 3.22, using simpson's rule and
! taking advantage of the logarithmic grid so that sum f(r)*dr becomes
! sum over f(r)*r*(0.05).  Linear interpolation is used at the end
! caps.  This version does not sum over 14 shells of identical
! atoms, instead it averages the contribution of one or more atoms
! of type 2 at the location of atom 1.  Louck's description (except
! for his integration algorithm) is very clear.
!
! input:  
!         rn        distance from atom 1 to atom 2 in au
!         ann       number of type 2 atoms to add to atom 1, can
!                   be fractional
!         aa2(i)    potential or density at atom 2
! output: aasum(i)  spherically summed contribution added into this
!                   array so that sumax can be called repeatedly
!                   and the overlapped values summed into aasum
!
! Note that this routine requires that all position data be on a
! grid  rr(j) = exp (-8.8d0 + (j-1)*0.05d0), which is the grid
! used by Louck, and also used by ATOM if nuclear options not used.
!
! Coded by Steven Zabinsky, December 1989
! Modified for FEFF cluster code, August 1990, siz
! Bug fixed, May 1991, SIZ
! Another bug fixed, Mar 1992, SIZ
!
! T.L.Louck, "Augmented Plane Wave Method", W.A.Benjamin, Inc., 1967

      subroutine sumax (rn, ann, aa2, aasum)
      implicit double precision (a-h, o-z)
      parameter (nptx=250)
      dimension aa2(nptx), aasum(nptx)
      dimension stor(nptx)
!#mn
       external ii, xx

!     jjchi     index beyond which aa2 is zero
!     jtop      index just below distance to neighbor
!               aasum is calculated only up to index jtop

!     Wigner-Seitz radius is set to 15 in ATOM.
      rws = 15
      jjchi = ii(rws)
      jtop  = ii(rn)

      topx = xx(jjchi)

      do 120  i = 1, jtop
         x = xx(i)
         xint = 0.0
         et = exp(x)
         blx = log(rn-et)
         if (blx .ge. topx)  goto 119
         jbl = 2.0+20.0*(blx+8.8)
         if (jbl .lt. 1)  jbl=1
         if (jbl .ge. 2)  then
!           use linear interp to make end cap near center of neighbor
            xjbl = jbl
            xbl = 0.05 * (xjbl-1.0) - 8.8
            g = xbl-blx
            xint = xint+0.5*g*(aa2(jbl)*(2.0-20.0*g)*exp(2.0*xbl)       &
     &             +20.0*g*aa2(jbl-1)*exp(2.0*(xbl-0.05)))
         endif
         tlx = log(rn+et)
         if (tlx .ge. topx)  then
            jtl = jjchi
            go to 90
         endif
         jtl = 1.0 + 20.0*(tlx+8.8)
         if (jtl .lt. jbl)  then
!           handle peculiar special case at center of atom 1
            fzn = aa2(jtl)*exp(2.0*(xbl-0.05))
            fz3 = aa2(jbl)*exp(2.0*xbl)
            fz2 = fzn+20.0*(fz3-fzn)*(tlx-xbl+0.05)
            fz1 = fzn+20.0*(fz3-fzn)*(blx-xbl+0.05)
            xint = 0.5*(fz1+fz2)*(tlx-blx)
            go to 119
         endif
         xjtl = jtl
         xtl = 0.05*(xjtl-1.0)-8.8
         c = tlx-xtl
         xint = xint+0.5*c*(aa2(jtl)*(2.0-20.0*c)                       &
     &         *exp(2.0*xtl)+aa2(jtl+1)*20.0*c                          &
     &         *exp(2.0*(xtl+0.05)))

   90    if (jtl .gt. jbl)  then
  100       xint = xint+0.5*(aa2(jbl)*exp(2.0*xbl)+aa2(jbl+1)           &
     &             *exp(2.0*(xbl+0.05)))*0.05
            jbl = jbl+1
            if (jbl .lt. jtl) then
               xbl = xbl+0.05
               go to 100
            endif
         endif
  119    stor(i) = 0.5*xint*ann/(rn*et)
  120 continue

      do 190  i = 1, jtop
         aasum(i) = aasum(i) + stor(i)
  190 continue

      return
      end
