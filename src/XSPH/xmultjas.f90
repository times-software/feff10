!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xmultjas.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xmultjas (k, kp, ls, xm)

      !use constants
      implicit double precision (a-h, o-z)
      complex*16 xm, alslb
!
!     Aleksi removed i**ls and moved to other parts of the code
!
!     xm is a prefactor in
!     see Grant eq. 6.2 calculate the factors 
!     <k|exp(i*q*z)|k'> = (-)**(j-m) * 3j( j l j'; -m 0 m')*R_k,k'
!     R_k,k'(l) = xm * \int dr (P_k*P_k'+ Q_k*Q_k') * j_l(qr)
!     xm = i**ls * (-)**(j+1/2) * (2*l+1) * sqrt((2j+1)(2j'+1)) *
!         3j(j l j'; 1/2  0 -1/2)
!     Additional condition: j+j'+l is even for a.ne.a' and odd for a=a'
      integer a, ap


!     set all angular momenta (multiplied by 2 for call to cwig3j)
      j2 = 2*abs(k) -1
      a = 1
      if (k.gt.0) a=-1
      jp2 = 2*abs(kp) -1
      ap = 1
      if (kp.gt.0) ap=-1
      ls2 = 2*ls

      i4 = j2+jp2+ls2+ a-ap
!     check additional condition
      if (mod(i4,4).eq.0 .or.ls2.lt.abs(j2-jp2) .or.ls2.gt.j2+jp2) then
!       j+j'+ls +(a-a')/2 is even
        xm = 0
      else
!       j+j'+ls +(a-a')/2 is odd
        aa = (j2+1)*(jp2+1)
        alslb =  sqrt(aa) * (ls2+1)
        xm = alslb * cwig3j(j2, ls2, jp2, 1, 0, 2) 
!        if (mod(j2+1,4).ne.0) xm = - xm
      endif

      return
      end
