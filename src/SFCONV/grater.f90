!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: grater.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************
!   This is Steve White's rewrite of Mike Teter's integration routine.
!   Modified by J. Rehr for complex integration.
!   The following is a listing of the arguments in the initial function
!   statement:
!      fn    -- routine requires external function statement in MAIN
!      xmin  -- lower limit
!      xmax  -- upper limit
!      abr   -- absolute tolerable error
!      rlr   -- relative tolerable error
!      nsing -- number of singularities or regions requiring
!                   special attention
!      xsing -- array of locations of singularities or endpoints
!                   of special regions
!      error -- output for routine error messages
!      numcal-- the number of times fn was called
!      maxns -- the maximum number of regions being considered simultaneously
!       function grater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
!       fn declared double precision
!       double precision function grater(fn,xmin,xmax,abr,rlr,
!       fn declared complex*16
!      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1

       double precision function grater(fn,xmin,xmax,abr,rlr,           &
     & nsing,xsing,error,numcal,maxns)

       implicit double precision (a-h,o-z)
       parameter (mx=1500)
       dimension xleft(mx),fval(3,mx),dx(3),wt(3)
       dimension wt9(9), xsing(20)
       double precision,external :: fn
        logical atsing
        save dx,wt,wt9
        data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
        data wt/0.277777777777777778  ,0.4444444444444444444  ,         &
     &                               0.2777777777777777778  /
        data wt9/0.0616938806304841571  ,0.108384229110206161  ,        &
     &           0.0398463603260281088  ,0.175209035316976464  ,        &
     &           0.229732989232610220  ,0.175209035316976464  ,         &
     &           0.0398463603260281088  ,0.108384229110206161  ,        &
     &           0.0616938806304841571  /
! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater) and the region
! is removed from consideration (and a new region is the right-most).
        nstack=nsing+1
        maxns = nstack
        error=0.
        grater=0.
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
        xleft(1)=xmin
        xleft(nsing+2)=xmax
        if(nsing.gt.0) then
          do 9 j=1,nsing
9           xleft(j+1)=xsing(j)
        endif
! For each region, calculate the function and store at three selected points.
        do 1 jj=1,nstack
          del=xleft(jj+1)-xleft(jj)
!         print*, 'fn call j= ,'
          do 1 j=1,3
!         print*, 'fn call in grater j= ',j
1           fval(j,jj)=fn(xleft(jj)+del*dx(j))
!         print*, 'output of fn call, fval(j,jj)',fval(j,jj)
        numcal = nstack * 3
6       continue
          if(nstack+3.ge.mx) then
            write(*,*) 'TOO MANY REGIONS'
            stop 0006
          endif
! Divide the rightmost region into three subregions.
          del=xleft(nstack+1)-xleft(nstack)
          xleft(nstack+3)=xleft(nstack+1)
          xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
          xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
          fval(2,nstack+2)=fval(3,nstack)
          fval(2,nstack+1)=fval(2,nstack)
          fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
          icount=0
          value=0.
          valu=0.
          do 3 j=nstack,nstack+2
            del1=xleft(j+1)-xleft(j)
!         print*, 'fn call 2'
            fval(1,j)=fn(xleft(j)+dx(1)*del1)
            fval(3,j)=fn(xleft(j)+dx(3)*del1)
!         print*, 'fn call 2'
            numcal = numcal + 2
            do 5 k=1,3
              icount=icount+1
              value=value+wt9(icount)*fval(k,j)*del
5             valu=valu+fval(k,j)*wt(k)*del1
3         continue
          dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
          frac = del / (xmax - xmin)
          atsing = .false.
          if(frac .le. 1.0e-8) atsing = .true.
          if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.          &
     &       (atsing .and.                                              &
     &     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
! The following commented out line is Teeter's old error criterion.
!          if(dif.le.abr.or.dif.le.rlr*abs(value))then
            grater=grater+value
            error=error+abs(dif)
            nstack=nstack-1
! If no more regions, we are done.
            if(nstack.le.0) return
          else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
            nstack=nstack+2
            maxns = max(maxns,nstack)
          endif
        go to 6
        end
