!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xscorratan.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xscorratan(ispec, emxs ,ne1, ne, ik0, xsec, xsnorm, chia, vrcorr, vicorr, cchi)
!     calculate the correction to xsec due to convolution with
!     lorentzian, based on integrtion in complex energy plane
!     Brouder et al PRB ???
!     the output correction is returned via cchi. The rest is input
!      mu(omega) = xsec + xsnorm*chia  + (cchi)

      use constants
	  use dimsmod, only: nex
      implicit double precision (a-h, o-z)

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex) 
      complex*16 xmu(nex), aa, bb, c1, f1, f2, ec
      parameter (eps4 = 1.0d-4)
      complex*16, external :: lorenz
    

      ifp = 0
      ne2 = ne-ne1
      efermi = dble(emxs(ne)) 
      xloss = dimag(emxs(1))
      gamach=xloss/2.0d0
      pii=4.0d0*atan(1.0d0)
      xvert = max(xloss, 0.02/hart)
!
!     Instead of using a fancy contour integral (becuase of the 
!     problems with it ) we will use simple atan
!

      do 10 ie = 1,ne
        cchi(ie) = 0
        write(88,*) ie,dble(emxs(ie)),imag(xsec(ie)),imag(chia(ie))
   10   xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)

      bb = xmu(ik0)
      if (abs(vrcorr).gt.eps4) then
         efermi = efermi - vrcorr
         do 20 ie = 1,ne1
   20    omega(ie) = dble(emxs(ie))
         call terpc(omega, xmu ,ne1, 1, efermi, bb)
         do 30 ie = 1, ne2
   30    emxs(ne1+ie) = emxs(ne1+ie) + vrcorr
      endif
              
      bb = bb/xmu(ik0)
!     rescale values on vertical axis


      do 200 ie = 1, ne1
!        cycle over energy points on horizontal grid
!        make cchi what ever needs to be subtracted
         dele = dble(emxs(ie)) - efermi
         alorstep=-0.50d0+atan2(dele,gamach)/pii
         cchi(ie)=xmu(ie)*alorstep
         if (ispec.eq.2) cchi(ie) = -cchi(ie) - xmu(ie)
         write(87,*) dele,dble(emxs(ie)),imag(xmu(ie)),imag(cchi(ie)),alorstep

  200 continue

!     restore the input energy mesh
      if (vicorr.gt.eps4) then
         do 250 ie = 1, ne1
  250    emxs(ie) = emxs(ie) - coni*vicorr
      endif
      if (abs(vrcorr).gt.eps4) then
         do 260 ie = ne1+1, ne
  260    emxs(ie) = emxs(ie) - vrcorr
      endif

      return
      end


