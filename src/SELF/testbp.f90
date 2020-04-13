!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: testbp.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM testbp

      COMPLEX*16 q, cpar(2), bpr1, r1, bpr2, r2, bpr3, r3
      DOUBLE PRECISION dppar(4), x
      INTEGER ie
      EXTERNAL bpr1, bpr2, bpr3

      PRINT*, "# Enter wp, gamma:"
      READ*, dppar(1), dppar(2)
      PRINT*, 'Enter q: '
      READ*, x
      q = x
      dppar(1) = dppar(1)
      dppar(2) = dppar(2)
      dppar(3) = 1.d0 
      OPEN(unit=31,file='bpr1.dat',status='replace')
      OPEN(unit=32,file='bpr2.dat',status='replace')
      OPEN(unit=33,file='bpr3.dat',status='replace')
      DO ie = 1, 4000
         q = ie*0.01d0
!        Set k
         cpar(1) = 2.d0
!        Set E
         cpar(2) = cpar(1)**2/2.d0
         r1 = bpr1(q,dppar,cpar)
         r2 = bpr2(q,dppar,cpar)
         r3 = bpr3(q,dppar,cpar)
         WRITE(31,*) DBLE(q), DBLE(r1), DIMAG(r1)
         WRITE(32,*) DBLE(q), DBLE(r2), DIMAG(r2)
         WRITE(33,*) DBLE(q), DBLE(r3), DIMAG(r3)
      END DO

      END
