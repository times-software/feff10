!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: radjas.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/12/17 01:12:33 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine radjas (ifl, kinit, dgc0, dpc0, dgcx0, dpcx0,  &
          ikap, p, q, pn, qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold, xirf, &
          ljmax,ljneeded,jbess,iorthg,ortcor  )
!     performs radial integration for multipole matrix element
!     or central atom absorption depending on flag 'ifl'.

!     INPUT
!     ifl - number corresponds to the calling order in xsect.f
!     - 1 - calculate matrix element (rkk)
!     - 2 - calculate cross section (xsec)
!     cross term needed for spin-dependent potential only
!     - 3 - cross term (xsec) with irregular part for current kappa
!     - 4 - cross term (xsec) with regular part for current kappa
!     ls - specifies multipole transition 
!     bf - bessel functions for x-ray k-vector for l=0,1,2
!     kinit - initial kappa
!     dgc0,dpc0 - large (small) dirac components for initial orbital 
!     ikap  = final state kappa
!     p,q   Dirac components for regular (R) final state solution
!     pn,qn  Dirac components for irregular(N) final state solution
!     ri,dx - radial grid
!     ilast - last integration point
!     iold  - 0 - do nothing to xrcold, xncold (ic3=0 case)
!     1 - store intermediate results in xrcold, xncold (ic3=1)
!     2 - use intermediate results in xrcold, xncold (ic3=1)
!     
!     OUTPUT
!     xrcold,xncold - coupling to regular (R) and irregular(N) solutions
!     both output and input
!     xirf  - value of the radial integral

      use dimsmod, only :nrptx
	  use global_inp,only: nq
      implicit double precision (a-h, o-z)

      dimension ri(nrptx), dgc0(nrptx), dpc0(nrptx)
      dimension dgcx0(nrptx), dpcx0(nrptx),dgc1(nrptx), dpc1(nrptx)
      dimension bf(0:2, nrptx)
!     dimension bf(0:ltot+4, nrptx)
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx)
      complex*16 p0(nrptx), q0(nrptx),px0(nrptx), qx0(nrptx)
!     storage for calculation of cross term (SPIN 1 only)
!      complex*16 xrcold(nrptx,0:ljmax) , xncold(nrptx,0:ljmax)
      complex*16 xrcold(nrptx,0:1) , xncold(nrptx,0:1)
      
      complex*16  xirf(0:ljmax)
      complex*16  ortcor(0:ljmax)

      integer ljmax,ll
      integer ljneeded(0:ljmax)
      double precision jbess(nrptx,0:ljmax)
      integer iorthg
      
!     local staff
      complex*16  xm(0:ljmax)
      complex*16 xrc(nrptx,0:ljmax), xnc(nrptx,0:ljmax)
      complex*16 ap,aq,apx,aqx
      complex*16 cnormx,cnorm,cS0
      

      linit = kinit
      if (kinit.lt.0) linit = - kinit - 1
      
      lfin = ikap
      if (ikap.lt.0) lfin = - ikap - 1
!     set multipliers  from Grant,Advan.Phys.,v.19,747(1970) eq. 6.2,
!     using Messiah's "Q.M." appendices to calculate 3j symbols

      do i=1,nrptx
         p0(i)=dcmplx(dgc0(i),0.0d0)
         q0(i)=dcmplx(dpc0(i),0.0d0)
      end do

      call getorthg(nrptx, ri, linit,dgc0,dpc0,linit,p0,q0, ap,aq,dx, 0, ilast)
      cnorm=ap+aq
      cnorm=dcmplx(sqrt(dble(cnorm)),0.0d0)
      
      call getorthg(nrptx, ri, linit,dgcx0,dpcx0,linit,p0,q0, apx,aqx,dx, 0, ilast)
      cS0=apx+aqx
      do i=1,nrptx
         dgc1(i)=dgc0(i)
         dpc1(i)=dpc0(i)
      end do

      if (kinit.eq.ikap) then 
!
!     test orthogonality         
!
         if (ifl.eq.1) then
            call getorthg(nrptx, ri, linit,dgc0,dpc0,linit,p,q, ap,aq,dx, 0, ilast)
            call getorthg(nrptx, ri, linit,dgcx0,dpcx0,linit,p,q, apx,aqx,dx, 0, ilast)
         else

            call getorthg(nrptx, ri, linit,dgc0,dpc0,-(linit+1),pn,qn, ap,aq,dx, 0, ilast)
            call getorthg(nrptx, ri, linit,dgcx0,dpcx0,-(linit+1),pn,qn, apx,aqx,dx, 0, ilast)

         end if

      end if

      do ll=0,ljmax

         if (ljneeded(ll).gt.0) then 
            call xmultjas( kinit,ikap, ll, xm(ll))
         end if
      end do
!     radial integrals depending on case
      if (ifl.eq.1) then
!     single radial integral for rkk - reduced matrix elements
!     xirf = <f |p| i> relativistic version of dipole m.e.
         do ll=0,ljmax
            if (ljneeded(ll).gt.0) then 
               if (kinit.eq.ikap) then 
                  do  i = 1, ilast
                     xrc(i,ll) = (jbess(i,ll)-ortcor(ll)) * (dgc1(i)*p(i)+dpc1(i)*q(i))
                     xnc(i,ll) = 0.0d0
                  end do
               else
                  do  i = 1, ilast
                     xrc(i,ll) = (jbess(i,ll)) * (dgc1(i)*p(i) + dpc1(i)*q(i))
                     xnc(i,ll) = 0.0d0
                  end do
               end if
!     store xrc if needed
                  if (iold.eq.1) then 
                     stop 'should not be here'
                  end if
               xirf(ll) = lfin + linit + ll +1
               call csommjas (ri, xrc(1,ll), xnc(1,ll), dx, xirf(ll), 0, ilast)
            end if
         end do
      else
!     need to perform double radial integral in all cases below
         if (ifl.eq.2) then
!     combine regular(kdif) and irregular(kdif) solution into
!     the central atom absorption coefficient xsec (mu = dimag(xsec))
!     thus for real energy dimag(xsec)=xsnorm
            do ll=0,ljmax
               if (ljneeded(ll).gt.0) then
                  if (ikap.eq.kinit) then
                     do i = 1, ilast
                        xnc(i,ll) = (jbess(i,ll)-ortcor(ll)) * (dgc1(i)*pn(i) + dpc1(i)*qn(i))
                     end do
                  else
                     do i = 1, ilast
                        xnc(i,ll) = (jbess(i,ll)) * (dgc1(i)*pn(i) + dpc1(i)*qn(i))
                     end do
                  end if
!     store irregular contribution for later use
                  if (iold.eq.1) then 
                     stop 'should not be here'
                  end if
               end if
            end do
         elseif (ifl.eq.3 .and. iold.eq.2) then
            stop 'should not be here'
         elseif(ifl.eq.4 .and. iold.eq.2) then
            stop 'should not be here'
         endif
         
!     same staff for all double integrals
         if ((iold.eq.0.and.ifl.eq.2) .or. (ifl.gt.2.and.iold.eq.2)) then
!     do radial integration for r'>r first
!     power of xrc near zero
            do ll=0,ljmax
               if (ljneeded(ll).gt.0) then  
                  
                  lpwr = lfin + linit +2
!     factor 2 since integral(r<r')=integral(r>r')
                  xirf(ll) = 2 * xrc(1,ll) * ri(1) /(lpwr+1)
                  xnc(1,ll) = xnc(1,ll) * xirf(ll)
                  do 70 i = 2, ilast
                     xirf(ll) = xirf(ll) + (xrc(i-1,ll)+xrc(i,ll)) * (ri(i)-ri(i-1))
                     xnc(i,ll) = xnc(i,ll) * xirf(ll)
 70               continue
                  do 80 i = 1,ilast
 80                  xrc(i,ll) = 0
                     xirf(ll) = lpwr+1+linit+1-lfin
!     ready for second integral over r from 0 to \infty
                     call csommjas (ri, xrc(1,ll), xnc(1,ll), dx, xirf(ll), 0, ilast)
               end if
            end do
         end if
      end if
      do ll=0,ljmax   
         if (ifl.eq.1) then 
            if (ljneeded(ll).gt.0) then  
               xirf(ll) = xirf(ll) * xm(ll)
            end if
         else
            if (ljneeded(ll).gt.0) then  
               xirf(ll) = xirf(ll) * xm(ll)*dconjg(xm(ll))
            end if
         end if
      end do
      return
      end


      subroutine getorthg(nrptx, ri, linit,dgc0,dpc0,lfin,p,q, ap,aq,dx, m, np)
      implicit none
      integer nrptx
      double precision ri(nrptx)
      integer linit,lfin,np
      double precision dgc0(nrptx),dpc0(nrptx)
      complex*16 p(nrptx),q(nrptx)
      complex*16 ap,aq
      double precision dx
      integer m,mm
      complex*16 da,dc
      double precision d1,dpas,db,dd
      integer k,j,l

      mm=m+1
      dpas=dx

      
      d1=mm+lfin+linit+1.0d0
      
      ap=dcmplx(0.0d0,0.0d0)
      aq=dcmplx(0.0d0,0.0d0)
      
      k=np

      
      do while (k.gt.0)
         if (k.eq.np .or. k.lt.5) then
            ap = ap + 14.d0*dgc0(k)*p(k)*ri(k)**mm
            aq = aq + 14.d0*dpc0(k)*q(k)*ri(k)**mm
         else
            ap = ap + 28.d0*dgc0(k)*p(k)*ri(k)**mm
            aq = aq + 28.d0*dpc0(k)*q(k)*ri(k)**mm
         end if
         k = k - 4
      end do
      k=k+4

      j=np-1
      do while (j.gt.k)
         ap = ap + 64.d0*dgc0(j)*p(j)*ri(j)**mm
         aq = aq + 64.d0*dpc0(j)*q(j)*ri(j)**mm
         j = j - 2
      end do
!

      l=np-2
      do while (l.gt.k)
         ap = ap + 24.d0*dgc0(l)*p(l)*ri(l)**mm
         aq = aq + 24.d0*dpc0(l)*q(l)*ri(l)**mm
         l = l - 4
      end do

      ap = ap * dx / 45.d0
      aq = aq * dx / 45.d0
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=ri(1)*(ri(2)**m)/db
      dd=(ri(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      ap=ap+dd*dgc0(1)*p(1)-db*dgc0(2)*p(2)
      aq=aq+dd*dpc0(1)*q(1)-db*dpc0(2)*q(2)
      return
      end

      subroutine getorthg2(nrptx, ri, linit,dgc0,dpc0,lfin,p,q, ap,aq,dx, m, np)
      implicit none
      integer nrptx
      double precision ri(nrptx)
      integer linit,lfin,np
      double precision dgc0(nrptx),dpc0(nrptx)
      complex*16 p(nrptx),q(nrptx)
      complex*16 ap,aq
      double precision dx
      integer m,mm
      complex*16 da,dc
      double precision d1,dpas,db,dd
      integer k,j,l




      
      ap=dcmplx(0.0d0,0.0d0)
      aq=dcmplx(0.0d0,0.0d0)
      
      do j=2,np
         ap=ap+0.50d0*(ri(j)-ri(j-1)) *(dgc0(j)*p(j)+dgc0(j-1)*p(j-1))
         aq=aq+0.50d0*(ri(j)-ri(j-1)) *(dpc0(j)*q(j)+dpc0(j-1)*q(j-1))
      end do


      return
      end


      subroutine getcorrection(nrptx,ri,dx,jinit,linit,dpc0,dgc0,ljmax,qjbess,ilast, ortcor)
      use global_inp,only: nq
      implicit none 
      integer nrptx
      double precision ri(nrptx)
      double precision dx
      integer jinit,linit
      double precision dpc0(nrptx),dgc0(nrptx)
      integer ljmax
      double precision qjbess(nrptx,0:ljmax,nq)
      integer ilast
      complex*16 ortcor(0:ljmax,nq)
      integer ir,ll,mm,iq
      complex*16 xrc(nrptx), xnc(nrptx)
      complex*16 norml
      
      do ir=1,nrptx
         xrc(ir)=dcmplx(0.0d0,0.0d0)
         xnc(ir)=dcmplx(0.0d0,0.0d0)
      end do
      do ir=1,ilast
         xrc(ir)=dcmplx(dpc0(ir)*dpc0(ir),0.0d0)
         xrc(ir)=xrc(ir)+dcmplx(dgc0(ir)*dgc0(ir),0.0d0)
         xnc(ir)=dcmplx(0.0d0,0.0d0)
      end do
      norml=2*linit+1
      call csommjas(ri, xrc, xnc, dx, norml, 0, ilast)
	  do iq=1,nq
      do ll=0,ljmax
         if (ll.le.jinit) then 
            do ir=1,ilast
               xrc(ir)=dcmplx(dpc0(ir)*qjbess(ir,ll,iq)*dpc0(ir),0.0d0)
               xrc(ir)=xrc(ir)+dcmplx(dgc0(ir)*qjbess(ir,ll,iq)*dgc0(ir),0.0d0)
               xnc(ir)=dcmplx(0.0d0,0.0d0)
            end do
            ortcor(ll,iq)=2*linit+ll+1
            call csommjas(ri, xrc, xnc, dx, ortcor(ll,iq), 0, ilast)
         else
            ortcor(ll,iq)=dcmplx(0.0d0,0.0d0)
         end if
         ortcor(ll,iq)=ortcor(ll,iq)/norml
      end do
      enddo

      return 
      end
