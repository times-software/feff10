!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: grids.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  grids ( ecv, xmu, negx, neg, emg , step, nflrx)
!     makes a grid in complex e-plane for scmt calculation
!     add complications for complex cases later. ala
!     emg is comlex energy in hartrees
      use constants
      implicit double precision (a-h, o-z)

      complex*16 emg(negx), eim, eimmin
      dimension step(nflrx)
!     the choice of e_cv should be automated later
!     all l-dos should be zero at ecv
!     fix it by hand if needed below
!     for some complicated materials may need multiple e_cv
!     it may also depend on core-valence separation



!     eimmin = the lowest im energy to search for fermi level
!     may simulate Fermi distr for occ numbers, thus may want
!     to lower eimmin for low temperatures.
      eimmin = coni*0.05/hart
      neg1 = (nflrx+1)/2
      neg3 = nflrx - 1
      neg2mx = negx-neg1-neg3
!     never do calculations on real axis.
      eim = eimmin*neg1**2
      eim = eimmin 
      de = dimag(eim)/4

      do 10 i =1, neg1
!        step linearly increases as one get farther from real axis
         eim = eimmin *i**2
         emg(i) = ecv +eim
  10  continue
      step(nflrx) = dimag(eim)/4

!     set energy step for integration eim above real axis
      de = dimag(emg(neg1))/4
      neg2= nint((xmu-ecv)/de)
      if (neg2.gt.neg2mx) neg2=neg2mx
      if (neg2.lt.neg1) neg2 = neg1
      de = (xmu-ecv) / neg2
      do 20 i = neg1+1,neg1+neg2
  20  emg(i) = emg(i-1) + de

      neg = neg1 + neg2 + neg3
      do 30 i =1, neg3
!        step linearly increases as one get farther from real axis
         eim = eimmin *(i+1)**2 /4.d0
         if (i.le.nflrx) step(i) = dimag(eim)/4
         emg(neg-i+1) = xmu + eim
  30  continue

	  
      return
      end
