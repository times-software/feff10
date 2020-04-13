!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yzktec.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine yzktec (f,af,g,ag,dr,ap,h,k,nd,np,idim, dyzk)
! calculation of yk(r)=zk(r)+ r**(k+1) * integral from r to 
!   infinity of  f(u) * u**(-k-1)
! zk(r) = r**(-k) * integral from 0 to r of f(u) * u**k

! at the origin f(r)=sum from i=1 to nd of af(i)*r**(ap+i-1)
! dr tabulation points   h exponential step
! np number of tabulation points for f
! idim dimension of the blocks f,g and dr

! at the origin yk=cte*r**(k+1)-developement limit
! the constant for yk lies in ap
! output functions yk and zk lie in f and g, and their
! development coefficients at the origin in af and ag.

! integration from point to point by a 4 points method.
! integral from r to r+h = h*(-f(r-h)+13*f(r)+13*f(r+h)-f(r+h+h))/24
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)

      complex*16 f,af,g,ag,ap, dyzk
      dimension f(nrptx),af(10),g(nrptx),ag(10),dr(nrptx)
 
!    initialisation and development coefficients of yk
      np= min(np,idim-1)
      f(np+1)=0.0d0
      b = dble(ap)
      ap=0.0d00
      g(1)=0.0d00
      do 15 i=1,nd
         b=b+1.0d00
         ag(i)=af(i)/(b+k)
         if (af(i).ne.0.0d00) then
            c=dr(1)**b
            g(1)=g(1)+ag(i)*c
!         for irregular solution b-k-1 can become zero
            if (abs(b-k-1) .le. 0.00001) then
               af(i) = 0.0
               b = b - 1.0d0
            else
               af(i)=(k+k+1)*ag(i)/(b-k-1)
            endif
            ap=ap+af(i)*c
         endif
 15   continue
      do 21 i=1,np
 21   f(i)=f(i)*dr(i)

!     calcualation of zk
      hk=h*k
      e = exp(-h)
      ehk = e**k 

      if (k.ne.0)then
       b1 = (ehk-1.0d0 +hk) / (hk*k)
      else
       b1=h/2.0
      endif

      b0 = h-(1.0+hk)*b1
      do 51 i=1,np
 51      g(i+1)=g(i)*ehk+b0*f(i)+f(i+1)*b1
 
!     calculation of yk
      f(np+1)=g(np+1) + dyzk
      ehk=ehk*e
      i=k+k+1
      hk=hk+h
      b1 = i*(ehk-1.0d0 +hk) / (hk*(k+1))
      b0 = i*h-(1.0+hk)*b1
      do 75  i=np,1,-1
 75      f(i) = f(i+1)*ehk+b0*g(i+1)+b1*g(i)

      ap=(ap+f(1))/(dr(1)**(k+1))
      return
      end
