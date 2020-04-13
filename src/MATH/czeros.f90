!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: czeros.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CQdrtc(Coef,Sol,NSol)
!     Combutes the zeros of a quadratic polynomial
!cccccccccccccccccccccccccccccccccccccccccccccccc            
!     Input
!     Coef - array of coefficients
      COMPLEX*16 Coef(3)
!cccccccccccccccccccccccccccccccccccccccccccccccc
!     Output
!     Sol  - Array of solutions
!     NSol - # of solutions (only one if Coef(1) = 0 etc.)
!     NSol = -1 means a and b are zero
      COMPLEX*16 Sol(2)
      INTEGER NSol
!cccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables
      COMPLEX*16 q, Sqrt
      DOUBLE PRECISION Sgn

      IF(Coef(1).eq.0.d0) THEN
         IF(Coef(2).eq.0.d0) THEN
            NSol = -1
            RETURN
         ELSE
            NSol = 1
            Sol(1) = -Coef(3)/Coef(2)
         END IF
      ELSE
         NSol = 2
         Root = Sqrt(Coef(2)**2-4.d0*Coef(1)*Coef(3))
         Sgn  = SIGN(DBLE(CONJG(Coef(2))*Root),1.d0)
         q    = -0.5d0*(Coef(2) + Sgn*Root)
         
         Sol(1) = q/Coef(1)
         Sol(2) = Coef(3)/q
      END IF

      RETURN
      END


      SUBROUTINE CCubic(Coef,Sol,NSol)
!     Combutes the zeros of a cubic polynomial
!cccccccccccccccccccccccccccccccccccccccccccccccc            
!     Input
!     Coef - array of coefficients
      COMPLEX*16 Coef(4)
!cccccccccccccccccccccccccccccccccccccccccccccccc
!     Output
!     Sol  - Array of solutions
!     NSol - # of solutions (only one if Coef(1) = 0 etc.)
!     NSol = -1 means a, b, and c are zero
      COMPLEX*16 Sol(4)
      INTEGER NSol
!cccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables
      COMPLEX*16 P1, P2, Q, R, Coef2(3), a, b, c
      DOUBLE PRECISION Sgn, Theta
!     PARAMETERS
      COMPLEX*16 I
      PARAMETER(I = (0.d0, 1.d0))
      DOUBLE PRECISION Pi
      PARAMETER(Pi = 3.141592653589793238462643d0)

      IF(Coef(1).eq.0.d0) THEN
         Coef2(1) = Coef(2)
         Coef2(2) = Coef(3)
         Coef2(3) = Coef(4)         
         CALL CQdrtc(Coef2,Sol,NSol)
      ELSE
         a = Coef(2)/Coef(1)
         b = Coef(3)/Coef(1)
         c = Coef(4)/Coef(1)
         NSol = 3
         Q = (a**2 - 3.d0*b)/9.d0
         R = (2.d0*a**3 - 9.d0*a*b + 27.d0*c)/54.d0

         IF(((DIMAG(Q).eq.0.d0).and.(DIMAG(R).eq.0.d0)).and.            &
     &        (DIMAG(R**2).lt.DIMAG(Q**3))) THEN
            Theta = ACOS (DBLE(R/SQRT(Q**3)))
            Sol(1) = -2*SQRT(Q)*Cos(Theta/3.d0) - a/3.d0
            Sol(2) = -2*SQRT(Q)*Cos((Theta+2.d0*Pi)/3.d0) - a/3.d0
            Sol(3) = -2*SQRT(Q)*Cos((Theta-2.d0*Pi)/3.d0) - a/3.d0
         ELSE
            Sgn = SIGN(1.d0, DBLE(CONJG(R)*SQRT(R**2-Q**3)))
            P1 = -(R + Sgn*SQRT(R**2-Q**3))**(1.d0/3.d0)
            IF(P1.eq.0.d0) THEN
               P2 = 0.d0
            ELSE
               P2 = Q/P1
            END IF
            Sol(1) = (P1 + P2) - a/3.d0
            Sol(2) = -0.5d0*(P1 + P2) - a/3.d0 +                        &
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
            Sol(3) = -0.5d0*(P1 + P2) - a/3.d0 -                        &
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
         END IF
      END IF

      RETURN
      END
      
