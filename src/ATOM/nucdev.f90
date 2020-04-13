!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: nucdev.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine nucdev (av,dr,dv,dz,hx,nuc,np,ndor,dr1)
!        * construction of nuclear potential *
! av coefficients of the development at the origin of nuclear potential
! dr  tabulation points
! dv  nuclear potential 
! dz  nuclear charge 
! hx  exponential step
! nuc index of the nuclear radius
! np  number of tabulation points
! ndor number of the coefficients for development at the origin
! the declared below arguments are saved, dr1 is the first
 
      implicit double precision (a-h,o-z)
      dimension av(10),dr(251),dv(251),at(251)

!    specify atomic mass and thickness of nuclear shell
! a atomic mass (negative or null for the point charge)
! epai parameter of the fermi density distribution
! (negative or null for uniform distribution), which is
!       cte / (1. + exp((r-rn)/epai) )
! with nuclear radius rn= 2.2677e-05 * (a**(1/3))

! calculate radial mesh
      a = 0.0
      epai = 0.0

      if (a.le.1.0d-01) then
         nuc=1
      else
         a=dz*(a**(1./3.))*2.2677d-05
         b=a/ exp(hx*(nuc-1))
         if (b.le.dr1) then
            dr1=b
         else
            b=log(a/dr1)/hx
            nuc=3+2*int(b/2.0)
            if (nuc.ge.np) call par_stop('dr1 too small')
!           index of atomic radius larger than dimension of dr
            dr1=a*exp(-(nuc-1)*hx)
         endif
      endif

      dr(1)=dr1/dz
      do 181 l=2,np
 181  dr(l)=dr(1)* exp(hx*(l-1))

      if (ndor.lt.5) then
!       * it should be at least 5 development coefficients
         call par_stop('stopped in programm nucdev, ndor should be > 4.')
!        stop
      endif
!  calculate nuclear potential on calculated radial mesh
      do 11 i=1,ndor
 11      av(i)=0.0d00
      if (epai.le.0.0) then
         do 15 i=1,np
 15         dv(i)=-dz/dr(i)
         if (nuc.le.1) then
            av(1)=-dz
         else
            av(2)=-3.0d00*dz/(dr(nuc)+dr(nuc))
            av(4)=-av(2)/(3.0d00*dr(nuc)*dr(nuc))
            l=nuc-1
            do 25 i=1,l
 25            dv(i)=av(2)+av(4)*dr(i)*dr(i)
         endif
      else
         b= exp(-dr(nuc)/epai)
         b=1.0d00/(1.0d00+b)
         av(4)=b
         av(5)=epai*b*(b-1.0d00)
         if (ndor.le.5) go to 45
         at(1)=1.0d00
         at(2)=1.0d00
         nf=1
         do 41 i=6,ndor
            n=i-4
            nf=n*nf
            dv(1)=n*at(1)
            n1=n+1
            dv(n1)=1.0d00
            do 35 j=2,n
 35         dv(j)=(n-j+2)*at(j-1)+(n-j+1)*at(j)
            do 37 j=1,n1
               m=n+1-j
               l=1
               if (mod(j,2).eq.0) l=-l
               av(i)=av(i)+l*dv(j)*(b**m)
 37            at(j)=dv(j)
 41         av(i)=b*av(i)*(epai**n)/nf
 45      do 47 i=1,np
            b=1.0d00+ exp((dr(i)-dr(nuc))/epai)
            if ((b*av(4)).gt.1.0d+15) go to 51
            dv(i)=dr(i)*dr(i)*dr(i)/b
 47         l=i
 51      if (l.ge.(np-1)) l=np-2
         k=l+1
         do 55 i=k,np
 55         dv(i)=0.0d00
         at(1)=0.0d00
         at(2)=0.0d00
         k=2
         do 61 i=4,ndor
            k=k+1
            do 58 j=1,2
 58         at(j)=at(j)+av(i)*(dr(j)**k)/k
            av(i)=av(i)/(k*(k-1))
 61         av(2)=av(2)+av(i)*(dr(1)**k)
         a=hx/2.4d+01
         b=a*1.3d+01
         k=l+1
         do 71 i=3,k
 71      at(i)=at(i-1)+b*(dv(i-1)+dv(i))-a*(dv(i-2)+dv(i+1))
         dv(l)=at(l)
         do 75 i=k,np
 75      dv(i)=dv(l)
         e= exp(hx)
         c=1.0d00/(e*e)
         i=l-1
 83      dv(i)=dv(i+1)/e+b*(at(i+1)/e+at(i))-a*(at(i+2)*c+at(i-1)*e)
         i=i-1
         if (i-1) 85,85,83
 85      dv(1)=dv(3)*c+hx*(at(1)+4.0d00*at(2)/e+at(3)*c)/3.0d00
         av(2)=(av(2)+dv(1))/dr(1)
         a=-dz/dv(l)
         do 95 i=4,ndor
 95      av(i)=-a*av(i)
         av(2)=a*av(2)
         do 97 i=1,np
 97      dv(i)=a*dv(i)/dr(i)
      endif

      return
      end
