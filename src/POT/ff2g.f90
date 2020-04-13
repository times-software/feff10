!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2g.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2012/02/03 07:17:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ff2g (gtr, iph, ie, ilast, xrhoce, xrhole, xrhocp,     &
     &             ee, ep, yrhole, yrhoce,  yrhocp, rhoval, xnmues,     &
     &             xnatph, xntot, iflr, iflrp, fl, fr, iunf)

      use constants
      use DimsMod, only: lx, nphx=>nphu
      implicit double precision (a-h, o-z)
!     the main output is l-dos in xrhoce, and valence density of states at distance r
!     in yrhoce, which at the input are only embedded atom quantities

      ! At input:
      !   gtr    - scattering matrix g (with appropriate trace already taken)
      !   iph    - current potential index
      !   ie     - current energy index (if first point, xrhocp/yrhocp need to be initialized)
      !   ilast  - index of last radial point to include in yrhoXe
      !   yrhoce - embedded atom valence density rho(E,r)
      !   xrhoce - embedded atom valence l-dos rho_L(E) = \int_{norman sphere} rho_L(E,r)
      !   yrhole - scattering contribution to rho(E,r) (except for scattering matrix)
      !   xrhole - scattering contribution to rho_L(E) (except for scattering matrix)
      !
      ! These are combined to form total valence l-dos rho = rhoce + gtr * rhole
      !
      ! At output:
      !   yrhoce - valence density (full, including scattering)
      !   xrhoce - valence l-dos   (full, including scattering)

      ! Additionally, this routine updates the energy-integrated valence l-dos and density
      !
      ! More input:
      !   yrhocp - previous valence density
      !   xrhocp - previous valence l-dos
      !   rhoval - energy-integrated valence density rho(r) (up to previous energy point)
      !   xnmues - energy-integrated number of electrons/atom for each l
      !   xntot  - energy-integrated total number of electrons in cluster
      !   xnatph - number of atoms of current potential type in cluster
      !   iflr   - current "floor" (imaginary part of energy)
      !   iflrp  - previous "floor" (imaginary part of energy)
      !   fl     - \sum_iph 2 rho(E,iph) xnatph at current energy
      !   fr     - \sum_iph 2 rho(E,iph) xnatph at previous energy
      !   iunf   - "unfreeze f" flag specifying whether to include f and higher states in valence
      !
      ! Output:
      !   rhoval - energy-integrated valence density rho(r) (including this point)
      !   xnmues - number of electrons/atom for each l
      !   xntot  - total number of electrons in cluster

      complex*16, intent(in) :: xrhole(0:lx)
      complex*16, intent(in) :: yrhole(251,0:lx)
      complex, intent(in) :: gtr(0:lx)

      complex*16, intent(inout) :: xrhoce(0:lx,0:nphx), xrhocp(0:lx,0:nphx)
      complex*16, intent(inout) :: yrhoce(251), yrhocp(251)
      real*8, intent(inout) :: xnmues(0:lx)

      complex*16 ee, ep, del, der, fl, fr
      dimension rhoval(251)

      do 730 il = 0,lx
        xrhoce(il, iph)=xrhoce(il, iph)+ gtr(il)*xrhole(il)
        if (ie.eq.1) xrhocp(il,iph) = xrhoce(il,iph)
  730 continue

      del = ee-ep
      der = del
!     if iflr=1 add/subtract integral from point to real axis
!     factor 2 below comes from spin degeneracy
      if (iflr.eq.1) der = der - coni * 2 * dimag(ee)
      if (iflrp.eq.1) del = del + coni * 2 * dimag(ep)
      do 750 il = 0, lx
        if (il.le.2 .or. iunf.ne.0) then
         fl = fl + 2 * xrhocp(il,iph) * xnatph
         fr = fr + 2 * xrhoce(il,iph) * xnatph
         xnmues(il) = xnmues(il) + dimag( xrhoce(il,iph) * der + xrhocp(il,iph) * del )
         xntot = xntot + xnmues(il) * xnatph
        endif
  750 continue

!c    calculate r-dependent l-dos for later use
      do 840 il = 0,lx
      do 840 ir = 1,ilast
       if (il.le.2 .or. iunf.ne.0) then
        yrhoce(ir) = yrhoce(ir) + gtr(il)*yrhole(ir,il)
        if (ie.eq.1) yrhocp(ir) = yrhoce(ir)
       endif
  840 continue

      do 850 ir = 1, ilast
         rhoval(ir) = rhoval(ir) + dimag(yrhoce(ir)*der+yrhocp(ir)*del)
  850 continue

      return
      end
