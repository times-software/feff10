!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: somm2.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine somm2 (dr,dp,dpas,da,rnrm,m,np)
! Modified to use complex p and q.  SIZ 4/91
! Modified to use double simpson integration ALA 3/97
! integration by the method of simpson of dp*dr from 
! 0 to r=rnrm  with proper end corrections
! dpas=exponential step;
! for r in the neighborhood of zero dp=cte*r**da
! **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      dimension  dp(*)

      mm = m + 1
      d1=dble(da)+mm
      da=0.0
      db=0.0
!      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)**mm*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)**mm*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)**mm*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)**mm*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)**mm*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)**mm*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)**mm*(-a2+a3)
         else
!           like trapesoidal rule
            dc=dp(i)*dr(i)**mm
         endif
         da=da+dc
   70 continue
      da=dpas*da

!     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
