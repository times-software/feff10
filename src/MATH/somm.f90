!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: somm.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine somm (dr,dp,dq,dpas,da,m,np)
!
! integration by the method of simpson of (dp+dq)*dr**m from
! 0 to r=dr(np)
! dpas=exponential step;
! for r in the neighborhood of zero (dp+dq)=cte*r**da
! **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(np), dp(np), dq(np)
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      if (dc) 20,40,30
   20 db=db+dc
      go to 40
   30 da=da+dc
   40 dc=dq(i)*dl
      if (dc) 50,70,60
   50 db=db+dc
      go to 70
   60 da=da+dc
   70 continue
      da = dpas * (da + db) / 3.0
      dc=exp(dpas)-1.0
      db=d1*(d1+1.0)*dc*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.0+1.0/(dc*(d1+1.0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
