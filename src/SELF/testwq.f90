!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: testwq.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Test
  INTEGER i1
  DOUBLE PRECISION q, wp, wq, Gamq
  EXTERNAL wq, Gamq

  wp = 1.d0
  DO i1 = 1, 3000
     q = DBLE(i1)/DBLE(100)
     PRINT*, Omegaq(wp,q), Gamq(wp,q)
  END DO

END PROGRAM Test
  
