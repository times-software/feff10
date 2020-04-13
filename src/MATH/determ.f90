!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: determ.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function determ(array,nord,nrows)
!
!     calculate determinate of a square matrix
!        (from bevington "data reduction and error analysis
!         for the physical sciences" pg 294)
!     array: matrix to be analyzed
!     nord: order of matrix
!     nrows:  first dimension of matrix in calling routine
!
      double precision array(nrows,nrows)
      determ = 1.
      do 150 k=1,nord
!
!
        if (array(k,k).ne.0) go to 130
        do 100 j=k,nord
          if (array(k,j).ne.0) go to 110
  100   continue
        determ = 0.
        go to 160
!
  110   do 120 i=k,nord
          saved = array(i,j)
          array(i,j) = array(i,k)
  120   array(i,k) = saved
        determ = -determ
!
  130   determ = determ*array(k,k)
        if (k.ge.nord) go to 150
        k1 = k+1
        do 140 i=k1,nord
          do 140 j=k1,nord
  140   array(i,j) = array(i,j)-array(i,k)*array(k,j)/array(k,k)
  150 continue
  160 return
! end double precision function determ
      end
