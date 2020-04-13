!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: TestPadi.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM TestPADi
  USE PADIO
  INTEGER i, j
  CHARACTER(1000) str
  j = 0
  i = 2345
  print*, i
  CALL wrpadisc(str,i)
  PRINT*, TRIM(str)
  CALL rdpadisc(str,j)
  print*, j
END PROGRAM TestPADi
