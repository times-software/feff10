!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ylm.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !ROUTINE: ylm
! !INTERFACE:
      subroutine ylm(v,lmax,y)
! !USES:
   !!!   use constants,only : pi
! !INPUT/OUTPUT PARAMETERS:
!   lmax   : spherical harmonics are calculated for l = 0 to lmax
!   v      : vector, argument of the spherical harmonics (we calculate Ylm(v/norm(v))
!   y      : array containing Ylm(v) for several l,m

! !DESCRIPTION:
!   1.  PURPOSE
!        The spherical harmonics (Condon and Shortley convention)
!          Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!        for vector V (given in Cartesian coordinates)
!        are calculated. In the Condon Shortley convention the
!        spherical harmonics are defined as
!        $$ Y(l,m) = (-1)^m \sqrt{\frac{1}{\pi}} P_{lm}(\cos{\theta}) \rm
!        e^{\rm i m \phi} $$
!                        
!        where  $P_{lm}(\cos{\theta})$ is the normalized Associated Legendre
!                  
!        function. Thus,
!                                             
!                     $$  Y(l,-m) = (-1)^m Y^*(l,m) $$
!           
!
!   2.  USAGE
!        DOUBLE PRECISION V(3), Y(5*5)
!        V(1) = ...
!        V(2) = ...
!        V(3) = ...
!        CALL YLM(V,4,Y)
!
!       ARGUMENT-DESCRIPTION
!          V      - DOUBLE PRECISION vector, dimension 3        (input)
!                   Must be given in Cartesian coordinates.
!                   Conversion of V to polar coordinates gives the
!                   angles Theta and Phi necessary for the calculation
!                   of the spherical harmonics.
!          LMAX   - INTEGER value                               (input)
!                   upper bound of L for which spherical harmonics
!                   will be calculated
!                   constraint:
!                      LMAX >= 0
!          Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                   contains the calculated spherical harmonics
!                   Y(1)                   for L .EQ. 0 (M = 0)
!                   Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                   ...
!                   Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                          for L .EQ. LMAX
!                                              (M = -L,...,L)
!                   constraint:
!                      Dimension of Y .GE. (LMAX+1)**2 (not checked)
!!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
!   3.     METHOD
!           The basic algorithm used to calculate the spherical
!           harmonics for vector V is as follows:
!
!           Y(0,0)
!           Y(1,0)
!           Y(1,1)
!           Y(1,-1) = -Y(1,1)
!           DO L = 2, LMAX
!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!              DO M = L-2, 0, -1
!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!                 Y(L,-M)= (-1)**M*Y(L,M)
!              ENDDO
!           ENDDO
!
!           In the following the necessary recursion formulas and
!           starting values are given:
!
!        Start:
!%                        +------+
!%                        |   1     
!%           Y(0,0) =  -+ | -----  
!%                       \| 4(Pi)  
!%
!%                                   +------+
!%                                   |   3     
!%           Y(1,0) =  cos(Theta) -+ | -----  
!%                                  \| 4(Pi)  
!%
!%                                     +------+
!%                                     |   3    i(Phi)
!%           Y(1,1) =  - sin(Theta) -+ | ----- e
!%                                    \| 8(Pi)  
!%
!%        Formula 1:
!%
!%           Y(l,l) =
!%                           +--------+
!%                           | (2l+1)   i(Phi)
!%            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!%                          \|   2l  
!%
!%        Formula 2:
!%                                  +---------------+  
!%                                  |  (2l-1)(2l+1)   
!%           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!%                                 \|   (l-m)(l+m)       
!%
!%                                    +--------------------+  
!%                                    |(l-1+m)(l-1-m)(2l+1)
!%                              -  -+ |-------------------- Y(l-2,m)
!%                                   \|  (2l-3)(l-m)(l+m)                 
!%
!%        Formula 3: (not used in the algorithm because of the division
!%                    by sin(Theta) which may be zero)
!%
!%                                    +--------------+  
!%                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!%           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!%                      sin(Theta)   \| (l+m+1)(l-m)       
!%
!%                                    +--------------+  
!%                                    |(l-m-1)(l+m+2)  -2i(Phi)
!%                              -  -+ |-------------- e        Y(l,m+2)
!%                                   \| (l-m)(l+m+1)                         
!%                                  
!%
! !REVISION HISTORY:
!   26. April 1994                                   Version 1.2
!   Taken 8 1 98 from SRC_lapw2 to SRC_telnes
!   Updated November 2004 (Kevin Jorissen)
!   cosmetics March 2005 (Kevin Jorissen)
!EOP

      implicit none
!
!   In/Output :

      integer,intent(in) ::     LMAX
      real*8,intent(in)  ::     V(3)
      complex*16,intent(out) :: Y(*)
!   Local variables :
      real*8,parameter :: pi = 3.1415926535897932384626433

      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
!
!        Y(0,0)
!
      YLLR = dble(1)/dsqrt(dble(4)*PI)
      YLLI = dble(0)
      Y(1) = DCMPLX(YLLR,YLLI)
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0) GOTO 999
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. dble(0)) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = SQRT(A*A+B*B)
         COSPH = A/AB
         SINPH = B/AB
      ELSE
         COSPH = dble(1)
         SINPH = dble(0)
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. dble(0)) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = SQRT(AB + C*C)
         COSTH = C/ABC
         SINTH = SQRT(AB)/ABC
      ELSE
         COSTH = dble(1)
         SINTH = dble(0)
      ENDIF
!
!        Y(1,0)
!
      Y(3) = DCMPLX(dsqrt(dble(3))*YLLR*COSTH,dble(0))
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -SQRT(dble(1.5))*YLLR*SINTH
      Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
      Y(2) = -DCONJG(Y(4))
!
      DO L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = DBLE(Y(INDEX-1))
         YL1L1I = DIMAG(Y(INDEX-1))
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
         YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
         Y(INDEX2) = DCMPLX(YLLR,YLLI)
         Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = dSQRT(DBLE(2*L+1))*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
         Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX2 - 4*L + 2
         I2L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = dble(1)/dSQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*dSQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
            YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
            Y(INDEX2) = DCMPLX(YLMR,YLMI)
            Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
!
            MSIGN  = -MSIGN
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I4L2   = I4L2   - 1
            I2L    = I2L    - 1
         ENDDO
      ENDDO
!
  999 RETURN
!
!        End of 'YLM'
!
      END
