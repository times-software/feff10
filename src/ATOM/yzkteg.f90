!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yzkteg.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine yzkteg (f,af,g,ag,dr,ap,h,k,nd,np,idim)
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

      implicit double precision (a-h,o-z)
      dimension f(251),af(10),g(251),ag(10),dr(251)
 
!    initialisation and development coefficients of yk
      np= min(np,idim-2)
      b=ap
      ap=0.0d00
      g(1)=0.0d00
      g(2)=0.0d00
      do 15 i=1,nd
         b=b+1.0d00
         ag(i)=af(i)/(b+k)
         if (af(i).ne.0.0d00) then
            c=dr(1)**b
            g(1)=g(1)+ag(i)*c
            g(2)=g(2)+ag(i)*(dr(2)**b)
            af(i)=(k+k+1)*ag(i)/(b-k-1)
            ap=ap+af(i)*c
         endif
 15   continue
      do 21 i=1,np
 21   f(i)=f(i)*dr(i)
      np1=np+1
      f(np1)=0.0d00
      f(np1+1)=0.0d00

!     calcualation of zk
      eh= exp(h)
      e=eh**(-k)
      b=h/2.4d+01
      c=1.3d+01*b
      ee=e*e*b
      b=b/e
      do 51 i=3,np1
 51   g(i)=g(i-1)*e+(c*(f(i)+f(i-1)*e)-(f(i-2)*ee+f(i+1)*b))
 
!     calcualation of yk
      f(np)=g(np)
      do 61 i=np1,idim
 61   f(i)=f(i-1)*e
      i=k+k+1
      b=i*b*eh
      ee=i*ee/(eh*eh)
      e=e/eh
      c=i*c
      do 71  i=np-1,2,-1
 71   f(i)=f(i+1)*e+(c*(g(i)+g(i+1)*e)-(g(i+2)*ee+g(i-1)*b))
      ee=e*e
      c=8.0d00*c/1.3d+01
      f(1)=f(3)*ee+c*(g(3)*ee+4.0d00*e*g(2)+g(1))
      ap=(ap+f(1))/(dr(1)**(k+1))
      return
      end
