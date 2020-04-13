!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: logi.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION Logi(a1, isgn)
      COMPLEX*16 Logi
      COMPLEX*16 a1
      DOUBLE PRECISION pi
      INTEGER isgn
      PARAMETER(pi = 3.1415926535897932384626433d0)
      
      Logi = 0.d0
      IF(DBLE(a1).lt.0.d0) THEN
         Logi = pi*DBLE(isgn)
      END IF
      Logi = (0.d0,1.d0)*Logi
      
      RETURN
      END
