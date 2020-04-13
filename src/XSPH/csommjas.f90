!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: csommjas.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine csommjas (dr,dp,dq,dpas,da,m,np)
! Modified to use complex p and q.  SIZ 4/91
! integration by the method of simpson of (dp+dq)*dr**m from 
! 0 to r=dr(np)
! dpas=exponential step;
! for r in the neighborhood of zero (dp+dq)=cte*r**da
!
! Aleksi used a bintegrate routine from Eric Shirley Et. al.
! as a starting point and rewrote this 
!
!
! **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),dq(*),da,dc
      integer k,j,l
      mm=m+1
      d1=da+mm

      da=dcmplx(0.d0,0.0d0)
!
      k=np
      do while (k.gt.0)
         if (k.eq.np .or. k.lt.5) then
            da = da + 14.d0*dp(k)*dr(k)**mm
            da = da + 14.d0*dq(k)*dr(k)**mm
         else
            da = da + 28.d0*dp(k)*dr(k)**mm
            da = da + 28.d0*dq(k)*dr(k)**mm
         end if
         k = k - 4
      end do
      k = k + 4
!
      j=np-1
      do while (j.gt.k)
         da = da + 64.d0*dp(j)*dr(j)**mm
         da = da + 64.d0*dq(j)*dr(j)**mm
         j = j - 2
      end do
!
      l=np-2
      do while (l.gt.k)
         da = da + 24.d0*dp(l)*dr(l)**mm
         da = da + 24.d0*dq(l)*dr(l)**mm
         l = l - 4
      end do
!
      da = da * dpas / 45.d0


!
!     initial point correction from the initial csomm routine
!e
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
