!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: testrs.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM testrs
      DOUBLE PRECISION dppar(4)
      COMPLEX*16 q, cpar(2), bpr1, bpr2, bpr3
      INTEGER iE, iq
      EXTERNAL bpr1, bpr2, bpr3
      
      PRINT*, "Enter w0: "
      READ*, dppar(1)
      PRINT*, "Enter gam0: "
      READ*, dppar(2)
      dppar(3) = 0.d0
      dppar(4) = 0.d0

      DO iE = 1, 100
         cpar(2) = DBLE(iE)
         cpar(1) = SQRT(2.d0*cpar(2))
         DO iq = 1, 100
            q = DBLE(iq)
            WRITE(33,'(I4,I4,4F30.10)') iE, iq,                         &
     &           DBLE(bpr1(q,dppar,cpar)), DBLE(bpr2(q,dppar,cpar)),    &
     &           DBLE(bpr3(q,dppar,cpar))
         END DO
         WRITE(33,*) '# New Block'
         WRITE(33,*) ' '
      END DO

      END 
