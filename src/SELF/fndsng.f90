!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fndsng.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
!     Josh Kas
!     This subroutine finds the singularities in the integrands of eq. 13
!     in
!     Single-particle Spectrum of the Degenerate Electron Gas
!     II. Numerical Results for Electrons Coupled to Plasmons
!     Phys. kondens. Materie, Bd. 6 (1967)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     In practice this routine only solves for the singularities of one
!     of the three integrands, then checks to see that the singularity
!     is within the limits of integration, and throws out singularities
!     that are not.
!     In units of the Fermi energy the equations to solve are:
!     1) +/- k*q**3 + 2*(3*k**2 - E - 2/3)*q**2 +/- 4*k*(k**2 - E)*q +
!                     [(k**2 - E)**2 - Wp**2] = 0 
!     2) q**4 + 4/3*q**2 + Wp**2 - (1 - E)**2 = 0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     Input:
!     Limit1 - Lower limit of integration
!     Limit2 - Upper limit of integration
!     CPar   - Array of complex parameters passed to function
!              CPar(1) = ck/kFermi
!              CPar(2) = Energy/EFermi + i*Gamma/EFermi
!     DPPar  - Array of double precision parameters passed to function
!              DPPar(1) = Wp/EFermi
!              DPPar(2) = Gamma/EFermi
!              DPPar(3) = Energy/EFermi
!              DPPar(4) = xeg (gap energy)
!     iFcn   - Integer denoting which function is the integrand
!              iFcn = 1: solve eqs 1 and 2 for q
!              iFcn = 2: solve eq 1 for q      
      COMPLEX*16 Limit1, Limit2, CPar(10)
      DOUBLE PRECISION DPPar(10)
      INTEGER iFcn
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output:
!     XSing  - Array of singularities
!     NSing  - Number of singularities
      DOUBLE PRECISION XSing(20)
      INTEGER NSing
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables:
!     Coef   - Coefficients of q**n of eq. to solve.
!              eq = Coef(1)*q**n + Coef(2)*q**(n-1)...
!     Sol(4) - Array of solutions to the equation.
!     XSing2 - Temp XSing      
!     Test   - Used to test solution of equation.
!     Zero   - Tolerance for testing solution to eqs.
!     NSol   - Number of solutions to eq.
!     Order  - Used to order singulaties from smallest to largest
      COMPLEX*16 Coef(4), Sol(4)
      DOUBLE PRECISION Test, Zero, XSing2(4)
      INTEGER NSol, Order(100)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Loop variables
      INTEGER i1, i2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Initialization
      NSing = 0
      Zero=1.d-4

!     Solve eq 1 for q with + sign
      Coef(1) = 4.d0*CPar(1)
      Coef(2) = 2.d0*(3.d0*CPar(1)**2 - DPPar(3) - 2.d0/3.d0)
      Coef(3) = 4.d0*CPar(1)*(CPar(1)**2 - DPPar(3))
      Coef(4) = (CPar(1)**2 - DPPar(3))**2 - DPPar(1)**2
      
      CALL CCubic(Coef, Sol, NSol)
      
!     Test solutions. If Sol is a solution and it is real
!     and it lies between Limit1 and Limit2, add it to list
!     of singularities.         
      DO i1 = 1, NSol
         Test = ABS((CPar(1)+Sol(i1))**2 - DPPar(3) +                   &
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN               
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.                     &
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.                   &
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN               
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO
      
!     Now solve eq. 1 for q with - sign
      Coef(1) = -Coef(1)
      Coef(3) = -Coef(3)
      
      CALL CCubic(Coef, Sol, NSol)
      
!     Test solutions as before.
      DO i1 = 1, NSol
         Test = ABS((CPar(1)-Sol(i1))**2 - DPPar(3) -                   &
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.                     &
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.                   &
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO

!     If iFcn = 1 (Solving for singularities of r1(q))
      IF(iFcn.eq.1) THEN         
!        Solve eq. 2 for q
         Coef(1) = 1.d0
         Coef(2) = 4.d0/3.d0
         Coef(3) = DPPar(1)**2

         CALL CQdrtc(Coef,Sol,NSol)
         DO i1 = 1, NSol
            XSing2(2*i1-1) =  DBLE(SQRT(Sol(i1)))
            XSing2(2*i1)   = -DBLE(SQRT(Sol(i1)))
         END DO

!        Test Solutions
         DO i1 = 1, 2*NSol
            IF((DBLE(XSing2(i1)).ge.DBLE(Limit1)).and.                  &
     &           (DBLE(XSing2(i1)).le.DBLE(Limit2)).and.                &
     &              (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = XSing2(i1)
            END IF
         END DO
      END IF
      
!     Sort singularities
      CALL QSORTI(Order,NSing,Xsing)
      DO i1 = 1, NSing
         XSing2(i1) = XSing(i1)
      END DO
      DO i1 = 1, NSing
         XSing(i1) = XSing2(Order(i1))
      END DO

      RETURN
      END
