!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: lipman.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lipman( j1, ph, qh, pir, qir, fxc, jri, imx0,  chik)
!     calculate K*chi_0 matrix, by perorming radial integration
!     chik(r,r") = \int d r' K(r,r') \chi0(r', r") * (hx*r(i))
!     K(r,r') = r< / r>**2  + fxc(r)*\delta(r-r')
!     chi0(r',r") = \phi_j1(r') \phi_j1(r") R_l(r<) H_l(r>)
!     (hh*r(i)) - factor needed on exponential grid integration
!      and is required for matrix inversion in r-space
!     imx0 the last point of tabulation of the wave function
      use dimsmod, only: nrptx
	  use constants 
      implicit double precision (a-h,o-z)
      parameter (npi=6, test=1.0d+5)
      parameter (ccl=2*alpinv, csq=ccl**2 )
      complex*16 en,agi,api
      complex*16 chik (251,251) 
      complex*16 gg,ag,gp,ap,dv,av,eg,ceg,ep,cep, vm(nrptx)
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),  &
     &   av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
      common/dff/cg(nrptx,30),cp(nrptx,30),bg(10,30),bp(10,30),         &
     &             fl(30), fix(30), ibgp

      complex*16 ec,eph,f,g, ph(nrptx), qh(nrptx)
      complex*16 pir(nrptx), qir(nrptx)
      dimension fxc(nrptx)
!
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idm
      complex*16 f1(nrptx), f2(nrptx), f3(nrptx), f4(nrptx)
! hx exponential step
! dr radial mesh
! idm dimension of the block dr

!     initialize
      do 20 i0 = 1,251
      do 20 j0 = 1,251
        chik(i0,j0) = 0
 20   continue

      do 205 i = 1, imx0
        eg(i) = (cg(i,j1) * ph(i) + cp(i,j1) *qh(i)) / dr(i)**2
 205  continue
      f1(1) = 0
      do 210 i = 1, imx0-1
        f1(i+1) = f1(i) + (dr(i+1)-dr(i)) / 2 * (eg(i)+eg(i+1))
 210  continue

      do 215 i = 1, imx0
        eg(i) = (cg(i,j1) * ph(i) + cp(i,j1) *qh(i)) * dr(i)
 215  continue
      f2(1) = 0
      do 220 i = 1, imx0-1
        f2(i+1) = f2(i) + (dr(i+1)-dr(i)) / 2 * (eg(i)+eg(i+1))
 220  continue

      do 225 i = 1, imx0
        eg(i) = (cg(i,j1) * pir(i) + cp(i,j1) *qir(i)) / dr(i)**2
 225  continue
      do 230 i = jri, imx0
 230  f3(i) = 0
      do 235 i =  jri-1, 1, -1
        g = eg(i)*pir(i) + ep(i)*qir(i)
        f3(i) = f3(i+1) + (dr(i+1)-dr(i)) / 2 * (eg(i) + eg(i+1))
 235  continue

      do 245 i = 1, imx0
        eg(i) = (cg(i,j1) * pir(i) + cp(i,j1) *qir(i)) * dr(i)
 245  continue
      do 250 i = jri, imx0
 250  f4(i) = 0
      do 255 i =  jri-1, 1, -1
        f4(i) = f4(i+1) + (dr(i+1)-dr(i)) / 2 * (eg(i) + eg(i+1))
 255  continue

      hh = 0.05d0
      do 265 i = 1, imx0
        gg(i) = cg(i,j1) * ph(i) + cp(i,j1) * qh(i)
        eg(i) = gg(i) * hh * dr(i)
        gp(i) = cg(i,j1) * pir(i) + cp(i,j1) * qir(i)
        ep(i) = gp(i) * hh * dr(i)
 265  continue
      do 280 i0 = 1,251
         i = 1 + 5*(i0-1)
         do 275 j0 = 1,251
           j = 1 + 5*(j0-1)
           if (i.gt.imx0 .or.j.gt.imx0) goto 275

           if (i0.le.j0) then
             chik(i0,j0) = chik(i0,j0) +                                &
     &       (f2(i)/dr(i)**2 + (f1(j)-f1(i))*dr(i)) * ep(j) +           &
     &       f3(j)*dr(i) * eg(j)
!  add xc part
             chik(i0,j0) = chik(i0,j0) + fxc(i)* gg(i) * ep(j)
           else
             chik(i0,j0) = chik(i0,j0) +                                &
     &       f2(j)/dr(i)**2 * ep(j) +                                   &
     &       ((f4(j)-f4(i))/dr(i)**2 + f3(i)*dr(i)) * eg(j)
!  add xc part
             chik(i0,j0) = chik(i0,j0) + fxc(i)* gp(i) * eg(j)
           endif
 275     continue
 280  continue

      return
      end
