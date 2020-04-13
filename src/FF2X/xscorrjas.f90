!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xscorrjas.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xscorrjas(ispec, emxs ,ne1, ne, ik0, xsec, xsnorm, chia, vrcorr, vicorr, cchi)
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
      complex*16 lorenz
      external lorenz

      ifp = 0
      ne2 = ne-ne1
      efermi = dble(emxs(ne)) 
      xloss = dimag(emxs(1))
      xvert = max(xloss, 0.02/hart)

      do 10 ie = 1,ne
        cchi(ie) = 0
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
      do 60 ie = ne1+1, ne
   60 xmu(ie) = xmu (ie) * bb 

      if (vicorr.gt.eps4 .and. xloss.eq.xvert) then
         xloss = xloss + vicorr
         xvert = xloss
         do 40 ie=1,ne2
   40    omega(ie) = dimag(emxs(ne1+ie))
         call terpc(omega, xmu(ne1+1) ,ne2, 1, xloss, aa)
         do 50 ie = 1, ne1
            xx = vicorr**2 /(vicorr**2 + (dble(emxs(ie))-efermi)**2)
            xmu(ie) = xmu(ie) + (aa - xmu(ik0)) * xx
            emxs(ie) = emxs(ie) + coni*vicorr
   50    continue
      endif

      do 200 ie = 1, ne1
!        cycle over energy points on horizontal grid

         dele = dble(emxs(ie)) - efermi
         if (abs(dele).lt.eps4) dele = 0.0d0
         w1 = dimag(emxs(ne1+1))
         w2 = dimag(emxs(ne1+2))
         w3 = dimag(emxs(ne1+3))

         if (xloss.lt.xvert) goto 120
!        matsubara pole and  sommerfeld correction
         cchi(ie) = cchi(ie) + lorenz(ifp,xloss,w1,dele)*xmu(ne1+1) *2*coni*w1 &
          + coni * w1**2 / 6 * (lorenz(ifp,xloss,w3,dele)*xmu(ne1+3)- lorenz(ifp,xloss,w2,dele)*xmu(ne1+2)) / (w3-w2) 

  120    continue
         if (dele .le. -eps4) cchi(ie) = cchi(ie) - xmu(ie)
         if (abs(dele) .lt. eps4 .and. xvert.eq.xloss) cchi(ie) = cchi(ie) - xmu(ie)/2

         if (xloss.lt.xvert) goto 130
!        integration over vertical axis to final point
!        use linear interpolation for xmu(ifp=0) or xmu*(i*w-de)(ifp=1)
         do 100 iv = ne1+2,ne-1
           w1 = dimag(emxs(iv))
           w2 = dimag(emxs(iv+1))
           bb = (xmu(iv+1)-xmu(iv))/ (w2-w1)
           aa = xmu(iv) - bb * w1
           c1 = (bb + (aa-coni*dele*bb)/xloss ) / 2
           if (abs(dele).lt.eps4) then
              cchi(ie) = cchi(ie) -  coni*xloss/pi *c1* log( abs((w2-xloss)/(w1-xloss)) )
           else
              cchi(ie) = cchi(ie) -  coni*xloss/pi *c1* log((w2+coni*dele-xloss)/(w1+coni*dele-xloss))
           endif
           c1 = (bb - (aa-coni*dele*bb)/xloss ) / 2
           cchi(ie) = cchi(ie) -  coni*xloss/pi *c1* log((w2+coni*dele+xloss)/(w1+coni*dele+xloss))
  100    continue

!        add the correction from the tail to infinity
!        assume xmu = aa/(w+ec) at high w - like single pole.
         f1 = xmu(ne-1)
!        f1 = xsec(ne-1)
         w1 = dimag(emxs(ne-1))
         f2 = xmu(ne)
!        f2 = xsec(ne)
         w2 = dimag(emxs(ne))
         ec = 0.01*(f2-f1*w1/w2)
         if (abs(ec).gt.abs(f2-f1)) then
!          want be safe if f2=f1
           ec=0
         else
           ec =100*w2*ec/(f1-f2)
         endif
!        do not allow the pole be higher than w2/2 from real axis
         if (dble(-2*ec).gt.w2) ec = -w2/2 + coni * dimag(ec)
         aa = f2*(w2+ec)

!        can obtain analytical results for f' and f".
         ec = ec - coni*dele
         c1 = 1/(w2+coni*dele)
         if (abs(ec/w2).lt.0.1) then
            c1 = c1**2 / 2 - ec * c1**3 / 3
         else
            c1 = c1/ec - log(c1*(w2+ec+coni*dele)) /ec**2
         endif
         cchi(ie) = cchi(ie) - coni*xloss/pi*aa*c1

  130    continue
         if (ispec.eq.2) cchi(ie) = -cchi(ie) - xmu(ie)
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

