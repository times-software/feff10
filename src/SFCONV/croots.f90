!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: croots.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine croots (a,b,c,d,x1,x2,x3,nrroots)
! Solves the cubic polynomial ax^3+bx^2+cx+d=0.
! Returns the number of real roots nrroots and the three roots
! of the cubic polynomial.
      implicit none
      integer nrroots
      double precision a,b,c,d,p,q,r,ar,br,disc,disc2
      complex*16 AA,BB,yr,yi,y1,y2,y3,x1,x2,x3
      x1=0.d0
      x2=0.d0
      x3=0.d0
      if (a.eq.0.d0) then
        if (b.eq.0.d0) then
          if (c.eq.0.d0) then
            nrroots=0
            return
          endif
          nrroots=1
          x1=-d/c
          return
        endif
        disc=c**2-4.d0*d*b
        if (disc.ge.0.d0) then
          nrroots=2
        else
          nrroots=0
        endif
        x1=(-c+sqrt(dcmplx(disc,0.d0)))/(2.d0*b)
        x2=(-c-sqrt(dcmplx(disc,0.d0)))/(2.d0*b)
        return
      endif
      p=b/a
      q=c/a
      r=d/a
      ar=q-p**2/3.d0
      br=(2.d0*p**3-9.d0*p*q)/27.d0+r
      disc=br**2/4.d0+ar**3/27.d0
      if (disc.gt.0.d0) then
        nrroots=1
        disc2=-br/2.d0+sqrt(disc)
        if (disc2.ge.0.d0) then
          AA=dcmplx(disc2**(1.d0/3.d0),0.d0)
        else
          AA=dcmplx(-((-disc2)**(1.d0/3.d0)),0.d0)
        endif
        disc2=-br/2.d0-sqrt(disc)
        if (disc2.ge.0.d0) then
          BB=dcmplx(disc2**(1.d0/3.d0),0.d0)
        else
          BB=dcmplx(-((-disc2)**(1.d0/3.d0)),0.d0)
        endif
      else
        nrroots=3
        if (br.lt.0.d0) then
          AA=dcmplx(-br/2.d0,sqrt(-disc))**(1.d0/3.d0)
          BB=dcmplx(dble(AA),-dimag(AA))
        else
          AA=-dcmplx(br/2.d0,-sqrt(-disc))**(1.d0/3.d0)
          BB=dcmplx(dble(AA),-dimag(AA))
        endif
      endif
      yr=-(AA+BB)/2.d0
      yi=(AA-BB)*sqrt((-3.d0,0.d0))/2.d0
      y1=yr+yi
      y2=yr-yi
      y3=AA+BB
      x1=y1-p/3.d0
      x2=y2-p/3.d0
      x3=y3-p/3.d0
      return
      end
