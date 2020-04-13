!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fprime.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fprime( ei, emxs ,ne1, ne3, ne, ik0, xsec, xsnorm,chia,&
     &                  vrcorr, vicorr, cchi)
!     calculate f' including solid state and lifetime effects.
!     using algorithm in Ankudinov, Rehr danes paper.
!     the output correction is returned via cchi. The rest is input
!      mu(omega) = xsec + xsnorm*chia  + (cchi)

      use dimsmod, only: nex
	  use constants
      implicit double precision (a-h, o-z)

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex) 
      complex*16 xmu(nex), aa, bb, c1, x1, x2, ec, temp
      complex*16 xmup(nex)
      dimension emp(nex)
      parameter (eps4 = 1.0d-4)
      complex*16 lorenz, funlog, value
      external lorenz, funlog
      dimension dout(7,nex)
      character*72 string
      dimension oscstr(14), enosc(14)
      integer ient
      data ient /0/

!     read data from fpf0.dat
      open (unit=16, file='fpf0.dat', status='old', iostat=ios)
      read  (16,*)  string
      read  (16,*)  eatom
      read  (16,*)  nosc
      do 5 i=1, nosc
        read (16,*) oscstr(i), enosc(i)
   5  continue
!     the rest is f0(Q) and is not currently needed
      close (unit=16)

      ient = ient+1
      ifp = 1
      efermi = dble(emxs(ne1+1)) 
      xloss = dimag(emxs(1))
      ne2 = ne-ne1-ne3
      if (ne2.gt.0) then
!        DANES
         do 10 ie = 1,ne1
   10    xmu(ie) = coni*xsnorm(ie) +  xsnorm(ie)*chia(ie)
         do 11 ie = ne1+1,ne1+ne2
   11    xmu (ie) = xsnorm(ie)*chia(ie)
         do 12 ie = ne-ne3+1, ne
   12    xmu (ie) =  coni*xsnorm(ie)
      else
!        FPRIME
         do 13 ie = 1,ne
   13    xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)
      endif

      if (abs(vrcorr).gt.eps4.and.ne2.gt.0) then
         bb = xmu(ik0)
         efermi = efermi - vrcorr
         do 20 ie = 1,ne1
   20    omega(ie) = dble(emxs(ie))
         call terpc(omega, xmu ,ne1, 1, efermi, bb)
         do 30 ie = 1, ne2
   30    emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
         if (abs(xmu(ik0)).gt. eps4) bb = bb/xmu(ik0)
!        rescale values on vertical axis
         do 60 ie = ne1+1, ne-ne3
   60    xmu(ie) = xmu (ie) * bb 
      endif
              

      if (vicorr.gt.eps4.and.ne2.gt.0) then
         xloss = xloss + vicorr
         do 40 ie=1,ne2
   40    omega(ie) = dimag(emxs(ne1+ie))
         call terpc(omega, xmu(ne1+1) ,ne2, 1, xloss, aa)
         do 50 ie = 1, ne1
            xx = vicorr**2 /(vicorr**2 + (dble(emxs(ie))-efermi)**2)
            xmu(ie) = xmu(ie)*(1.0d0 - xx) + aa * xx
            emxs(ie) = emxs(ie) + coni*vicorr
   50    continue
      endif

      do 200 ie = 1, ne1
!        cycle over energy points on horizontal grid

         dout(1,ie) = dble(emxs(ie)) * hart
         dele = dble(emxs(ie)) - efermi
!        delp correspond to pole with negative frequency
!        see Sakurai for details

         delp = -dele - 2*ei
!        delp = dele
!        dele = delp

         cchi(ie) = 0
         if (ne2.gt.0) then
            if (abs(dele).lt.eps4) dele = 0.0d0
            w1 = dimag(emxs(ne1+1))
            w2 = dimag(emxs(ne1+2))
            w3 = dimag(emxs(ne1+3))

!           matsubara pole
            temp = lorenz(ifp,xloss,w1,dele)*xmu(ne1+1)*2*coni*w1
            temp = temp + lorenz(ifp,xloss,w1,delp)*xmu(ne1+1)*2*coni*w1
            dout(2,ie)=dble(temp)
!           sommerfeld correction
            temp = coni*w1**2/ 6*(lorenz(ifp,xloss,w3,dele)*xmu(ne1+3)- &
     &      lorenz(ifp,xloss,w2,dele)*xmu(ne1+2)) / (w3-w2) 
            dout(3,ie)=dble(temp)

            cchi(ie) = lorenz(ifp,xloss,w1,dele)*xmu(ne1+1) *2*coni*w1  &
     &      + coni * w1**2 / 6 * (lorenz(ifp,xloss,w3,dele)*xmu(ne1+3)- &
     &      lorenz(ifp,xloss,w2,dele)*xmu(ne1+2)) / (w3-w2) 
!           from negative pole has additional minus sign
            cchi(ie) = cchi(ie) +                                       &
     &      lorenz(ifp,xloss,w1,delp)*xmu(ne1+1) *2*coni*w1             &
     &      + coni * w1**2 / 6 * (lorenz(ifp,xloss,w3,delp)*xmu(ne1+3)- &
     &      lorenz(ifp,xloss,w2,delp)*xmu(ne1+2)) / (w3-w2) 

!           theta funcion contribution only for positive pole
            if (dele .lt. eps4)    cchi(ie) = cchi(ie) - xmu(ie)
            if (abs(dele).lt.eps4) cchi(ie) = cchi(ie) + xmu(ie)/2

!           anomalous contribution
            temp = 0
            wp = 2*ei
            if (dele.ge.eps4) temp = xmu(ie)
            if (abs(dele).lt.eps4) temp = xmu(ie)/2
            temp = temp + xmu(ik0)*  funlog(1,xloss,wp,dele)
!               xmu(iko) + xsec(ik0)  if n3 >0
            dout(4,ie)=dble(temp) 

!           integration over vertical axis to final point
            n1 = ne1+2
            n2 = ne-ne3
            call fpint (emxs, xmu, n1, n2, dele, xloss, eps4, efermi,   &
     &                  value)
            cchi(ie) = cchi(ie) + value
!           add contribution from other pole
            call fpint (emxs, xmu, n1, n2, delp, xloss, eps4, efermi,   &
     &                  value)
            cchi(ie) = cchi(ie) + value
         endif 

!        integration over horizontal axis to final point
         temp = 0
         if (ne2.gt.0) then
!           DANES
            n1 = ne1-ik0 + 1
            do 120 i = ik0, ne1
              emp(i-ik0+1) = dble(emxs(i))
              xmup(i-ik0+1) = coni*xsnorm(i)
  120       continue
            do 130 i = 1, ne3
              emp(i+n1) = dble(emxs(i+ne-ne3))
              xmup(i+n1) = xmu(i+ne-ne3)
  130       continue
            n2 = n1 + ne3
         else
!           FPRIME
            n1 = 0
            do 140 i = 1, ne1
              if (n1.eq.0 .and. dble(emxs(i)).gt. dble(emxs(ne1+1)))    &
     &            n1 = i
  140       continue
            do 150 i = 1, ne3
               emp(i) =  dble(emxs(ne1+i))
               xmup(i) =  xmu(ne1+i)
  150       continue
            n2 = ne3
         endif
         call fpintp (emp, xmup , n2, dele, xloss, efermi, value)
         temp  = temp + value
!        add contribution from other pole
         call fpintp (emp, xmup , n2, delp, xloss, efermi, value)
         temp  = temp + value

!         was used before
!c          contribution to fp from poles of the core states
!           temp=0
!           do 110  i=2, nosc
!c             eif = E_f- E_i  in hartrees
!c             eif = enosc(i)-enosc(1) 
!c             deltaf = deltaf - oscstr(i)*2*alpinv**2/eif
!              temp = temp + alpinv**2 * oscstr(i)* (dele -
!    1      enosc(i)+efermi-1)/ ((dele-enosc(i)+efermi-1)**2+xloss**2)
!              temp = temp + alpinv**2 * oscstr(i)* (delp -
!    1      enosc(i)+efermi-1)/ ((delp-enosc(i)+efermi-1)**2+xloss**2)
! 110       continue

         dout(5,ie) = dble(temp)
         cchi(ie) = cchi(ie) + temp

!        total contribution (not normalized)
         temp = xmu(ie) + cchi(ie)
         dout(6,ie) = dble(temp)
!        (integral w2 to wmax) minus (cusp formula)
         dout (7,ie) = dout(6,ie)-dout(4,ie)
  200 continue

!     restore the input energy mesh
      if (vicorr.gt.eps4) then
         do 250 ie = 1, ne1
  250    emxs(ie) = emxs(ie) - coni*vicorr
      endif
      if (abs(vrcorr).gt.eps4) then
         do 260 ie = 1, ne2
  260    emxs(ne1+ie) = emxs(ne1+ie) + vrcorr
      endif

!     if (ient.eq.1) then
      open(unit=3,file='danes.dat', status='unknown', iostat=ios)
      write(3,310) '# E  matsub. sommerf. anomal. tale, total, differ.'
  310 format (a)
      do 300 ie = 1, ne1
         write(3,320) (dout(i,ie), i=1,7)
  320    format ( 7(1x,1pe11.4))
  300 continue
      close(unit=3)
!     endif

      return
      end

      complex*16 function funlog (icase, xloss, w, dele)
!     anomalous fp should have all main features of total fp
!     except smooth difference 
!     analytic expression for anomalous fp (without integral)
!     is obtained by adding and subtracting G(Ef + i*Gamma) / E-w
!     and performing integral for Im axis analytically
!     icase = 1 simplified expression (compared to 2) 
!     icase=2  use real w 
!     icase=3  pure imaginary w (absolute value is input)
      use par
	  use constants,only: pi,coni !KJ 7-09 surprised to see these missing ???
      implicit double precision (a-h, o-z)
      parameter (eps4 = 1.0d-4)

      if (icase.eq.1) then 
         if (abs(dele).ge.eps4) then 
            funlog= coni/2/pi*                                          &
     &      (log((-xloss+coni*dele)/w)+ log((xloss+coni*dele)/w))

         else
            funlog= coni/pi*log(abs(xloss/w))
         endif

      elseif (icase.eq.2) then
        if (abs(dele).ge.eps4) then
          funlog= coni/2/pi* (w+coni*xloss) * (                         &
     &    ( log((-xloss+coni*dele)/w)) / (w+dele+coni*xloss) +          &
     &    ( log(( xloss+coni*dele)/w)) / (w+dele-coni*xloss))
        else
          funlog= coni/pi*(log(abs(xloss/w)))*                          &
     &    (1 + coni*xloss/(w-coni*xloss))
        endif

      elseif (icase.eq.3) then
        if (abs(dele).ge.eps4) then
          funlog= -(w+xloss)/2/pi* (                                    &
     &    log((-xloss+coni*dele)/w) / (dele+coni*(w+xloss)) +           &
     &    log(( xloss+coni*dele)/w) / (dele+coni*(w-xloss)) )
        else
          funlog= coni/pi* log(abs(xloss/w))*                           &
     &    (1 + xloss/(w-xloss))
        endif
      
      endif

      return
      end

      subroutine fpint (emxs, xmu, n1, n2, dele, xloss, eps4, efermi,   &
     &                  value)
!     performs integral for fp calculations between points n1 and n2.

      use dimsmod, only: nex
	  use constants
      implicit double precision (a-h, o-z)

      complex*16 emxs(nex), xmu(nex), value
      complex*16  z1, z2, aa, bb, c1

!     last interval - similar to Matsubara pole ( shift and - sign)
!     notice that this also works for horizontal axis if last value
!     is small
      z1 = emxs(n2)-efermi
      z2 = emxs(n2-1)-efermi
      value =  - coni/pi * (z1-dele) / (xloss**2+(z1-dele)**2)          &
     &          *xmu(n2) * (2 * (z1-z2))
!     all other intervals
      do  300 i = n1, n2-2
         z1 = emxs(i) - efermi
         z2 = emxs(i+1) - efermi
         bb=(xmu(i+1)*(z2-dele) - xmu(i)*(z1-dele)) / xloss / (z2-z1)
         aa = xmu(i)*(z1-dele)/xloss - bb * z1
         c1 = (aa+bb*(dele+coni*xloss )) / 2 /coni
         if (abs(dele-dble(z1)).lt.eps4 .and.                           &
     &       abs(dele-dble(z2)).lt.eps4) then
            value = value  -  coni/pi *c1*                              &
     &      log( abs((z2-dele-coni*xloss)/(z1-dele-coni*xloss)) )
         else
            value    = value   -  coni/pi *c1*                          &
     &      log((z2-dele-coni*xloss)/(z1-dele-coni*xloss))
         endif
         c1 = -(aa+bb*(dele-coni*xloss )) / 2 /coni
         value    = value    -  coni/pi *c1*                            &
     &   log((z2-dele+coni*xloss)/(z1-dele+coni*xloss))
  300  continue

      return
      end

      subroutine fpintp (em, xmu, n2, dele, xloss, efermi, value)
!     performs integral for fp calculations between points 1 and n2.
!     and adds tail to infinity
      use dimsmod, only: nex
	  use constants
      implicit double precision (a-h, o-z)

      dimension em(nex)
      complex*16 xmu(nex), value
      complex*16  z1, z2, aa, bb, cc

      value = 0
!     all intervals 
      do  300 i = 1, n2-1
         x1 = em(i) - efermi
         x2 = em(i+1) - efermi
         de = (x2-x1)/2
         x0 = (em(i) + em(i+1)) / 2
         call terpc(em, xmu, n2, 3, x0, aa)
         bb=(xmu(i+1) - xmu(i)) / (x2-x1)
         cc = (xmu(i+1) - aa - bb * de) / de**2
         z1 =  dele - x0 + efermi - coni*xloss
         z2 =  dele - x0 + efermi + coni*xloss
         value    = value  + 2*de*bb + 2*z1*de*cc +                     &
     &    log((de-z1)/(-de-z1)) * (aa+bb*z1+cc*z1**2)
         value    = value  + 2*de*bb + 2*z2*de*cc +                     &
     &    log((de-z2)/(-de-z2)) * (aa+bb*z2+cc*z2**2)
  300 continue

!     tail of xmu to infinity approximated by aa/(w-bb)**2
      x1 = em(n2-1)
      x2 = em(n2)
      a = sqrt ( dble(xmu(n2-1)/xmu(n2)) )
      b = ( a*x1 - x2) / (a-1)
      if (b.gt. x1) b = 0
      aa = xmu(n2) * (x2-b)**2
      z1 = dele -coni*xloss - b
      z2 = dele +coni*xloss - b
      x0 = x2 - b
      value = value + log( x0/(x0-z1) ) *aa/z1**2 - aa/z1/x0
      value = value + log( x0/(x0-z2) ) *aa/z2**2 - aa/z2/x0

!     multiply by constant factor
      value = - coni /2 /pi *value

      return
      end
