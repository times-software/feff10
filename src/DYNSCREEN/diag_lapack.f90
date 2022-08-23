!> \brief \b DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
!>          Input to test for NaN.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END FUNCTION DISNAN
!> \brief \b DLABAD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLABAD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabad.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabad.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabad.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLABAD( SMALL, LARGE )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   LARGE, SMALL
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLABAD takes as input the values computed by DLAMCH for underflow and
!> overflow, and returns the square root of each of these values if the
!> log of LARGE is sufficiently large.  This subroutine is intended to
!> identify machines with a large exponent range, such as the Crays, and
!> redefine the underflow and overflow limits to be the square roots of
!> the values computed by DLAMCH.  This subroutine is needed because
!> DLAMCH does not compensate for poor arithmetic in the upper half of
!> the exponent range, as is found on a Cray.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SMALL
!> \verbatim
!>          SMALL is DOUBLE PRECISION
!>          On entry, the underflow threshold as computed by DLAMCH.
!>          On exit, if LOG10(LARGE) is sufficiently large, the square
!>          root of SMALL, otherwise unchanged.
!> \endverbatim
!>
!> \param[in,out] LARGE
!> \verbatim
!>          LARGE is DOUBLE PRECISION
!>          On entry, the overflow threshold as computed by DLAMCH.
!>          On exit, if LOG10(LARGE) is sufficiently large, the square
!>          root of LARGE, otherwise unchanged.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of DLABAD
!
      END SUBROUTINE DLABAD
!> \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, D, P, Q
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLADIV performs complex division in  real arithmetic
!>
!>                       a + i*b
!>            p + i*q = ---------
!>                       c + i*d
!>
!> The algorithm is due to Michael Baudin and Robert L. Smith
!> and can be found in the paper
!> "A Robust Complex Division in Scilab"
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION
!>          The scalars p and q in the above expression.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   BS
      PARAMETER          ( BS = 2.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLADIV1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0D0

      OV = DLAMCH( 'Overflow threshold' )
      UN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Epsilon' )
      BE = BS / (EPS*EPS)

      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL DLADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL DLADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
!
      RETURN
!
!     End of DLADIV
!
      END SUBROUTINE DLADIV

!> \ingroup doubleOTHERauxiliary


      SUBROUTINE DLADIV1( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   R, T
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLADIV2
      EXTERNAL           DLADIV2
!     ..
!     .. Executable Statements ..
!
      R = D / C
      T = ONE / (C + D * R)
      P = DLADIV2(A, B, C, D, R, T)
      A = -A
      Q = DLADIV2(B, A, C, D, R, T)
!
      RETURN
!
!     End of DLADIV1
!
      END SUBROUTINE DLADIV1

!> \ingroup doubleOTHERauxiliary

      DOUBLE PRECISION FUNCTION DLADIV2( A, B, C, D, R, T )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, R, T
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   BR
!     ..
!     .. Executable Statements ..
!
      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            DLADIV2 = (A + BR) * T
         ELSE
            DLADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         DLADIV2 = (A + D * (B / C)) * T
      END IF
!
      RETURN
!
!     End of DLADIV2
!
      END FUNCTION DLADIV2
!> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in DISNAN.
!>
!> DLAISNAN checks for NaNs by comparing its two arguments for
!> inequality.  NaN is the only floating-point value where NaN != NaN
!> returns .TRUE.  To check for NaNs, pass the same variable as both
!> arguments.
!>
!> A compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> Interprocedural or whole-program optimization may delete this
!> test.  The ISNAN functions will be replaced by the correct
!> Fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN1
!> \verbatim
!>          DIN1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DIN2
!> \verbatim
!>          DIN2 is DOUBLE PRECISION
!>          Two numbers to compare for inequality.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!     ..
!
!  =====================================================================
!
!  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END FUNCTION DLAISNAN
!> \brief \b DLAPY3 returns sqrt(x2+y2+z2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAPY3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow and unnecessary underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION
!>          X, Y and Z specify the values x, y and z.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS, HUGEVAL
!     ..
!     .. External Subroutines ..
      DOUBLE PRECISION   DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      HUGEVAL = DLAMCH( 'Overflow' )
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
!     W can be zero for max(0,nan,0)
!     adding all three entries together will make sure
!     NaN will not disappear.
         DLAPY3 =  XABS + YABS + ZABS
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ &
      &            ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY3
!
      END FUNCTION DLAPY3
      
!> \brief \b DZASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
!>    returns a double precision result.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION STEMP
      INTEGER I,NINCX
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
!     ..
      DZASUM = 0.0d0
      STEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DO I = 1,N
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      END IF
      DZASUM = STEMP
      RETURN
!
!     End of DZASUM
!
      END FUNCTION DZASUM
!> \brief \b ILAZLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAZLC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAZLC( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAZLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILAZLC( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16       ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
!     ..
!     .. Local Scalars ..
      INTEGER I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILAZLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILAZLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILAZLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END FUNCTION ILAZLC
!> \brief \b ILAZLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAZLR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAZLR( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAZLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILAZLR( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16       ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILAZLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILAZLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILAZLR = MAX( ILAZLR, I )
         END DO
      END IF
      RETURN
      END FUNCTION ILAZLR
!> \brief \b ZAXPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ZA
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*),ZY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZAXPY constant times a vector plus a vector.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] ZA
!> \verbatim
!>          ZA is COMPLEX*16
!>           On entry, ZA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!>
!> \param[in,out] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of ZY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
!     ..
      IF (N.LE.0) RETURN
      IF (DCABS1(ZA).EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            ZY(I) = ZY(I) + ZA*ZX(I)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZY(IY) = ZY(IY) + ZA*ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
!
      RETURN
!
!     End of ZAXPY
!
      END SUBROUTINE ZAXPY
!> \brief \b ZDOTC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*),ZY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDOTC forms the dot product of two complex vectors
!>      ZDOTC = X^H * Y
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!>
!> \param[in] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of ZY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG
!     ..
      ZTEMP = (0.0d0,0.0d0)
      ZDOTC = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTC = ZTEMP
      RETURN
!
!     End of ZDOTC
!
      END FUNCTION ZDOTC
!> \brief \b ZDOTU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*),ZY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDOTU forms the dot product of two complex vectors
!>      ZDOTU = X^T * Y
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!>
!> \param[in] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of ZY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
!     ..
      ZTEMP = (0.0d0,0.0d0)
      ZDOTU = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(I)*ZY(I)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTU = ZTEMP
      RETURN
!
!     End of ZDOTU
!
      END FUNCTION ZDOTU
!> \brief \b ZDSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZDSCAL scales a vector by a constant.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 ZX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DO I = 1,N
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      END IF
      RETURN
!
!     End of ZDSCAL
!
      END SUBROUTINE ZDSCAL
!> \brief \b ZGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebak.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   SCALE( * )
!       COMPLEX*16         V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEBAK forms the right or left eigenvectors of a complex general
!> matrix by backward transformation on the computed eigenvectors of the
!> balanced matrix output by ZGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N': do nothing, return immediately;
!>          = 'P': do backward transformation for permutation only;
!>          = 'S': do backward transformation for scaling only;
!>          = 'B': do backward transformations for both permutation and
!>                 scaling.
!>          JOB must be the same as the argument JOB supplied to ZGEBAL.
!> \endverbatim
!>
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  V contains right eigenvectors;
!>          = 'L':  V contains left eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrix V.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          The integers ILO and IHI determined by ZGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by ZGEBAL.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix V.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by ZHSEIN or ZTREVC.
!>          On exit, V is overwritten by the transformed eigenvectors.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
      &                   INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
      &    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEBAK', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
      &   RETURN
      IF( M.EQ.0 ) &
      &   RETURN
      IF( LSAME( JOB, 'N' ) ) &
      &   RETURN
!
      IF( ILO.EQ.IHI ) &
      &   GO TO 30
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL ZDSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL ZDSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         END IF
!
      END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
         IF( RIGHTV ) THEN
            DO 40 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
      &            GO TO 40
               IF( I.LT.ILO ) &
      &            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
      &            GO TO 40
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
      &            GO TO 50
               IF( I.LT.ILO ) &
      &            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
      &            GO TO 50
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of ZGEBAK
!
      END SUBROUTINE ZGEBAK
!> \brief \b ZGEBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebal.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   SCALE( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEBAL balances a general complex matrix A.  This involves, first,
!> permuting A by a similarity transformation to isolate eigenvalues
!> in the first 1 to ILO-1 and last IHI+1 to N elements on the
!> diagonal; and second, applying a diagonal similarity transformation
!> to rows and columns ILO to IHI to make the rows and columns as
!> close in norm as possible.  Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrix, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!>                  for i = 1,...,N;
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to INTEGER such that on exit
!>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied to
!>          A.  If P(j) is the index of the row and column interchanged
!>          with row and column j and D(j) is the scaling factor
!>          applied to row and column j, then
!>          SCALE(j) = P(j)    for j = 1,...,ILO-1
!>                   = D(j)    for j = ILO,...,IHI
!>                   = P(j)    for j = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The permutations consist of row and column interchanges which put
!>  the matrix in the form
!>
!>             ( T1   X   Y  )
!>     P A P = (  0   B   Z  )
!>             (  0   0   T2 )
!>
!>  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!>  along the diagonal.  The column indices ILO and IHI mark the starting
!>  and ending columns of the submatrix B. Balancing consists of applying
!>  a diagonal similarity transformation inv(D) * B * D to make the
!>  1-norms of each row of B and its corresponding column nearly equal.
!>  The output matrix is
!>
!>     ( T1     X*D          Y    )
!>     (  0  inv(D)*B*D  inv(D)*Z ).
!>     (  0      0           T2   )
!>
!>  Information about the permutations P and the diagonal matrix D is
!>  returned in the vector SCALE.
!>
!>  This subroutine is based on the EISPACK routine CBAL.
!>
!>  Modified by Tzu-Yi Chen, Computer Science Division, University of
!>    California at Berkeley, USA
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 2.0D+0 )
      DOUBLE PRECISION   FACTOR
      PARAMETER          ( FACTOR = 0.95D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
      &                   SFMIN2
!     ..
!     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZNRM2
      EXTERNAL           DISNAN, LSAME, IZAMAX, DLAMCH, DZNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
      &    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEBAL', -INFO )
         RETURN
      END IF
!
      K = 1
      L = N
!
      IF( N.EQ.0 ) &
      &   GO TO 210
!
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
!
      IF( LSAME( JOB, 'S' ) ) &
      &   GO TO 120
!
!     Permutation to isolate eigenvalues if possible
!
      GO TO 50
!
!     Row and column exchange.
!
   20 CONTINUE
      SCALE( M ) = J
      IF( J.EQ.M ) &
      &   GO TO 30
!
      CALL ZSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL ZSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
!
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
!
!     Search for rows isolating an eigenvalue and push them down.
!
   40 CONTINUE
      IF( L.EQ.1 ) &
      &   GO TO 210
      L = L - 1
!
   50 CONTINUE
      DO 70 J = L, 1, -1
!
         DO 60 I = 1, L
            IF( I.EQ.J ) &
      &         GO TO 60
            IF( DBLE( A( J, I ) ).NE.ZERO .OR. DIMAG( A( J, I ) ).NE. &
      &          ZERO )GO TO 70
   60    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
   70 CONTINUE
!
      GO TO 90
!
!     Search for columns isolating an eigenvalue and push them left.
!
   80 CONTINUE
      K = K + 1
!
   90 CONTINUE
      DO 110 J = K, L
!
         DO 100 I = K, L
            IF( I.EQ.J ) &
      &         GO TO 100
            IF( DBLE( A( I, J ) ).NE.ZERO .OR. DIMAG( A( I, J ) ).NE. &
      &          ZERO )GO TO 110
  100    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  110 CONTINUE
!
  120 CONTINUE
      DO 130 I = K, L
         SCALE( I ) = ONE
  130 CONTINUE
!
      IF( LSAME( JOB, 'P' ) ) &
      &   GO TO 210
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
!
      DO 200 I = K, L
!
         C = DZNRM2( L-K+1, A( K, I ), 1 )
         R = DZNRM2( L-K+1, A( I, K ), LDA )
         ICA = IZAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IZAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!        Guard against zero C or R due to underflow.
!
         IF( C.EQ.ZERO .OR. R.EQ.ZERO ) &
      &      GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR. &
      &       MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
            IF( DISNAN( C+F+CA+R+G+RA ) ) THEN
!
!           Exit if NaN to avoid infinite loop
!
            INFO = -3
            CALL XERBLA( 'ZGEBAL', -INFO )
            RETURN
         END IF
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 160
!
  170    CONTINUE
         G = C / SCLFAC
  180    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR. &
      &       MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 180
!
!        Now balance.
!
  190    CONTINUE
         IF( ( C+R ).GE.FACTOR*S ) &
      &      GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 ) &
      &         GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F ) &
      &         GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL ZDSCAL( N-K+1, G, A( I, K ), LDA )
         CALL ZDSCAL( L, F, A( 1, I ), 1 )
!
  200 CONTINUE
!
      IF( NOCONV ) &
      &   GO TO 140
!
  210 CONTINUE
      ILO = K
      IHI = L
!
      RETURN
!
!     End of ZGEBAL
!
      END SUBROUTINE ZGEBAL
!> \brief <b> ZGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
!                         WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
!> eigenvalues and, optionally, the left and/or right eigenvectors.
!>
!> The right eigenvector v(j) of A satisfies
!>                  A * v(j) = lambda(j) * v(j)
!> where lambda(j) is its eigenvalue.
!> The left eigenvector u(j) of A satisfies
!>               u(j)**H * A = lambda(j) * u(j)**H
!> where u(j)**H denotes the conjugate transpose of u(j).
!>
!> The computed eigenvectors are normalized to have Euclidean norm
!> equal to 1 and largest component real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N': left eigenvectors of A are not computed;
!>          = 'V': left eigenvectors of are computed.
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N': right eigenvectors of A are not computed;
!>          = 'V': right eigenvectors of A are computed.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A has been overwritten.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          W contains the computed eigenvalues.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is COMPLEX*16 array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!>          after another in the columns of VL, in the same order
!>          as their eigenvalues.
!>          If JOBVL = 'N', VL is not referenced.
!>          u(j) = VL(:,j), the j-th column of VL.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1; if
!>          JOBVL = 'V', LDVL >= N.
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is COMPLEX*16 array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!>          after another in the columns of VR, in the same order
!>          as their eigenvalues.
!>          If JOBVR = 'N', VR is not referenced.
!>          v(j) = VR(:,j), the j-th column of VR.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1; if
!>          JOBVR = 'V', LDVR >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,2*N).
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the QR algorithm failed to compute all the
!>                eigenvalues, and no eigenvectors have been computed;
!>                elements i+1:N of W contain eigenvalues which have
!>                converged.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!
!  @precisions fortran z -> c
!
!> \ingroup complex16GEeigen
!
!  =====================================================================
      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
      &                  WORK, LWORK, RWORK, INFO )
      implicit none
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
      &                   W( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU, &
      &                   IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX*16         TMP
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, ZGEHRD, &
      &                   ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC3, ZUNGHR
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DZNRM2, ZLANGE
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DZNRM2, ZLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, CONJG, AIMAG, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       CWorkspace refers to complex workspace, and RWorkspace to real
!       workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by ZHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N
            IF( WANTVL ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', &
      &                       ' ', N, 1, N, -1 ) )
               CALL ZTREVC3( 'L', 'B', SELECT, N, A, LDA, &
      &                       VL, LDVL, VR, LDVR, &
      &                       N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VL, LDVL, &
      &                      WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', &
      &                       ' ', N, 1, N, -1 ) )
               CALL ZTREVC3( 'R', 'B', SELECT, N, A, LDA, &
      &                       VL, LDVL, VR, LDVR, &
      &                       N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VR, LDVR, &
      &                      WORK, -1, INFO )
            ELSE
               CALL ZHSEQR( 'E', 'N', N, 1, N, A, LDA, W, VR, LDVR, &
      &                      WORK, -1, INFO )
            END IF
            HSWORK = INT( WORK(1) )
            MAXWRK = MAX( MAXWRK, HSWORK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
!
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
      &   RETURN
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) &
      &   CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Balance the matrix
!     (CWorkspace: none)
!     (RWorkspace: need N)
!
      IBAL = 1
      CALL ZGEBAL( 'B', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
      ITAU = 1
      IWRK = ITAU + N
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
      &             LWORK-IWRK+1, IERR )
!
      IF( WANTVL ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL ZLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate unitary matrix in VL
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL ZUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
      &                LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, &
      &                WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         IF( WANTVR ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL ZLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
!
      ELSE IF( WANTVR ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL ZLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate unitary matrix in VR
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL ZUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
      &                LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, &
      &                WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, &
      &                WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
!
!     If INFO .NE. 0 from ZHSEQR, then quit
!
      IF( INFO.NE.0 ) &
      &   GO TO 50
!
      IF( WANTVL .OR. WANTVR ) THEN
!
!        Compute left and/or right eigenvectors
!        (CWorkspace: need 2*N, prefer N + 2*N*NB)
!        (RWorkspace: need 2*N)
!
         IRWORK = IBAL + N
         CALL ZTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
      &                 N, NOUT, WORK( IWRK ), LWORK-IWRK+1, &
      &                 RWORK( IRWORK ), N, IERR )
      END IF
!
      IF( WANTVL ) THEN
!
!        Undo balancing of left eigenvectors
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
         CALL ZGEBAK( 'B', 'L', N, ILO, IHI, RWORK( IBAL ), N, VL, LDVL, &
      &                IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            SCL = ONE / DZNRM2( N, VL( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VL( 1, I ), 1 )
            DO 10 K = 1, N
               RWORK( IRWORK+K-1 ) = DBLE( VL( K, I ) )**2 + &
      &                               AIMAG( VL( K, I ) )**2
   10       CONTINUE
            K = IDAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VL( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL ZSCAL( N, TMP, VL( 1, I ), 1 )
            VL( K, I ) = DCMPLX( DBLE( VL( K, I ) ), ZERO )
   20    CONTINUE
      END IF
!
      IF( WANTVR ) THEN
!
!        Undo balancing of right eigenvectors
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
         CALL ZGEBAK( 'B', 'R', N, ILO, IHI, RWORK( IBAL ), N, VR, LDVR, &
      &                IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            SCL = ONE / DZNRM2( N, VR( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VR( 1, I ), 1 )
            DO 30 K = 1, N
               RWORK( IRWORK+K-1 ) = DBLE( VR( K, I ) )**2 + &
      &                               AIMAG( VR( K, I ) )**2
   30       CONTINUE
            K = IDAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VR( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL ZSCAL( N, TMP, VR( 1, I ), 1 )
            VR( K, I ) = DCMPLX( DBLE( VR( K, I ) ), ZERO )
   40    CONTINUE
      END IF
!
!     Undo scaling if necessary
!
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), &
      &                MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
         END IF
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of ZGEEV
!
      END SUBROUTINE ZGEEV
!> \brief \b ZGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEHD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H
!> by a unitary similarity transformation:  Q**H * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to ZGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= max(1,N).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the n by n general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the unitary matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZLARFG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEHD2', -INFO )
         RETURN
      END IF
!
      DO 10 I = ILO, IHI - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         ALPHA = A( I+1, I )
         CALL ZLARFG( IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
         A( I+1, I ) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL ZLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), &
      &               A( 1, I+1 ), LDA, WORK )
!
!        Apply H(i)**H to A(i+1:ihi,i+1:n) from the left
!
         CALL ZLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, &
      &               DCONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK )
!
         A( I+1, I ) = ALPHA
   10 CONTINUE
!
      RETURN
!
!     End of ZGEHD2
!
      END SUBROUTINE ZGEHD2
!> \brief \b ZGEHRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEHRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16        A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEHRD reduces a complex general matrix A to upper Hessenberg form H by
!> an unitary similarity transformation:  Q**H * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to ZGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the N-by-N general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the unitary matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!>          zero.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,N).
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This file is a slight modification of LAPACK-3.0's ZGEHRD
!>  subroutine incorporating improvements proposed by Quintana-Orti and
!>  Van de Geijn (2006). (See ZLAHR2.)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
      &                     TSIZE = LDT*NBMAX )
      COMPLEX*16        ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
      &                     ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB, &
      &                   NBMIN, NH, NX
      COMPLEX*16        EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEHD2, ZGEMM, ZLAHR2, ZLARFB, ZTRMM, &
      &                   XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         NB = MIN( NBMAX, ILAENV( 1, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
         LWKOPT = N*NB + TSIZE
         WORK( 1 ) = LWKOPT
      ENDIF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
!
!     Quick return if possible
!
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size
!
      NB = MIN( NBMAX, ILAENV( 1, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code)
!
         NX = MAX( NB, ILAENV( 3, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
!
!           Determine if workspace is large enough for blocked code
!
            IF( LWORK.LT.N*NB+TSIZE ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code
!
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEHRD', ' ', N, ILO, IHI, &
      &                 -1 ) )
               IF( LWORK.GE.(N*NBMIN + TSIZE) ) THEN
                  NB = (LWORK-TSIZE) / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
!
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
!
!        Use unblocked code below
!
         I = ILO
!
      ELSE
!
!        Use blocked code
!
         IWT = 1 + N*NB
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V**H
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL ZLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), &
      &                   WORK( IWT ), LDT, WORK, LDWORK )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
!           to 1
!
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
      &                  IHI, IHI-I-IB+1, &
      &                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, &
      &                  A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
!
!           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
!           right
!
            CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
      &                  'Unit', I, IB-1, &
      &                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL ZAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1, &
      &                     A( 1, I+J+1 ), 1 )
   30       CONTINUE
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL ZLARFB( 'Left', 'Conjugate transpose', 'Forward', &
      &                   'Columnwise', &
      &                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, &
      &                   WORK( IWT ), LDT, A( I+1, I+IB ), LDA, &
      &                   WORK, LDWORK )
   40    CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL ZGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of ZGEHRD
!
      END SUBROUTINE ZGEHRD
!> \brief \b ZGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
!>
!>    y := alpha*A**H*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
      &    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
      &    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
      NOCONJ = LSAME(TRANS,'T')
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
!
!        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
!
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of ZGEMV
!
      END SUBROUTINE ZGEMV
!> \brief \b ZGERC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGERC  performs the rank 1 operation
!>
!>    A := alpha*x*y**H + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERC ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*DCONJG(Y(JY))
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*DCONJG(Y(JY))
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
!
      RETURN
!
!     End of ZGERC
!
      END SUBROUTINE ZGERC
!> \brief \b ZHSEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHSEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhseqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhseqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhseqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!       CHARACTER          COMPZ, JOB
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZHSEQR computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**H, where T is an upper triangular matrix (the
!>    Schur form), and Z is the unitary matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input unitary
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>           = 'E':  compute eigenvalues only;
!>           = 'S':  compute eigenvalues and the Schur form T.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>           = 'N':  no Schur vectors are computed;
!>           = 'I':  Z is initialized to the unit matrix and the matrix Z
!>                   of Schur vectors of H is returned;
!>           = 'V':  Z must contain an unitary matrix Q on entry, and
!>                   the product Q*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>           set by a previous call to ZGEBAL, and then passed to ZGEHRD
!>           when the matrix output by ZGEBAL is reduced to Hessenberg
!>           form. Otherwise ILO and IHI should be set to 1 and N
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and JOB = 'S', H contains the upper
!>           triangular matrix T from the Schur decomposition (the
!>           Schur form). If INFO = 0 and JOB = 'E', the contents of
!>           H are unspecified on exit.  (The output value of H when
!>           INFO > 0 is given under the description of INFO below.)
!>
!>           Unlike earlier versions of ZHSEQR, this subroutine may
!>           explicitly H(i,j) = 0 for i > j and j = 1, 2, ... ILO-1
!>           or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>           The computed eigenvalues. If JOB = 'S', the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>           If COMPZ = 'N', Z is not referenced.
!>           If COMPZ = 'I', on entry Z need not be set and on exit,
!>           if INFO = 0, Z contains the unitary matrix Z of the Schur
!>           vectors of H.  If COMPZ = 'V', on entry Z must contain an
!>           N-by-N matrix Q, which is assumed to be equal to the unit
!>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
!>           if INFO = 0, Z contains Q*Z.
!>           Normally Q is the unitary matrix generated by ZUNGHR
!>           after the call to ZGEHRD which formed the Hessenberg matrix
!>           H. (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if COMPZ = 'I' or
!>           COMPZ = 'V', then LDZ >= MAX(1,N).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!>           On exit, if INFO = 0, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient and delivers very good and sometimes
!>           optimal performance.  However, LWORK as large as 11*N
!>           may be required for optimal performance.  A workspace
!>           query is recommended to determine the optimal workspace
!>           size.
!>
!>           If LWORK = -1, then ZHSEQR does a workspace query.
!>           In this case, ZHSEQR checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             < 0:  if INFO = -i, the i-th argument had an illegal
!>                    value
!>             > 0:  if INFO = i, ZHSEQR failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of W
!>                contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and JOB = 'E', then on exit, the
!>                remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and JOB   = 'S', then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and COMPZ = 'V', then on exit
!>
!>                  (final value of Z)  =  (initial value of Z)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'I', then on exit
!>                      (final value of Z)  = U
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'N', then Z is not
!>                accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>             Default values supplied by
!>             ILAENV(ISPEC,'ZHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
!>             It is suggested that these defaults be adjusted in order
!>             to attain best performance in each particular
!>             computational environment.
!>
!>            ISPEC=12: The ZLAHQR vs ZLAQR0 crossover point.
!>                      Default: 75. (Must be at least 11.)
!>
!>            ISPEC=13: Recommended deflation window size.
!>                      This depends on ILO, IHI and NS.  NS is the
!>                      number of simultaneous shifts returned
!>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
!>                      The default for (IHI-ILO+1) <= 500 is NS.
!>                      The default for (IHI-ILO+1) >  500 is 3*NS/2.
!>
!>            ISPEC=14: Nibble crossover point. (See IPARMQ for
!>                      details.)  Default: 14% of deflation window
!>                      size.
!>
!>            ISPEC=15: Number of simultaneous shifts in a multishift
!>                      QR iteration.
!>
!>                      If IHI-ILO+1 is ...
!>
!>                      greater than      ...but less    ... the
!>                      or equal to ...      than        default is
!>
!>                           1               30          NS =   2(+)
!>                          30               60          NS =   4(+)
!>                          60              150          NS =  10(+)
!>                         150              590          NS =  **
!>                         590             3000          NS =  64
!>                        3000             6000          NS = 128
!>                        6000             infinity      NS = 256
!>
!>                  (+)  By default some or all matrices of this order
!>                       are passed to the implicit double shift routine
!>                       ZLAHQR and this parameter is ignored.  See
!>                       ISPEC=12 above and comments in IPARMQ for
!>                       details.
!>
!>                 (**)  The asterisks (**) indicate an ad-hoc
!>                       function of N increasing from 10 to 64.
!>
!>            ISPEC=16: Select structured matrix multiply.
!>                      If the number of simultaneous shifts (specified
!>                      by ISPEC=15) is less than 14, then the default
!>                      for ISPEC=16 is 0.  Otherwise the default for
!>                      ISPEC=16 is 2.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!  =====================================================================
      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, &
      &                   WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    ZLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 15 )
!
!     ==== NL allocates some local workspace to help small matrices
!     .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
!     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
!     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
!     .    allows up to six simultaneous shifts and a 16-by-16
!     .    deflation window.  ====
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0d0 )
!     ..
!     .. Local Arrays ..
      COMPLEX*16         HL( NL, NL ), WORKL( NL )
!     ..
!     .. Local Scalars ..
      INTEGER            KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ILAENV, LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZLACPY, ZLAHQR, ZLAQR0, ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     ==== Decode and check the input parameters. ====
!
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DCMPLX( DBLE( MAX( 1, N ) ), RZERO )
      LQUERY = LWORK.EQ.-1
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.NE.0 ) THEN
!
!        ==== Quick return in case of invalid argument. ====
!
         CALL XERBLA( 'ZHSEQR', -INFO )
         RETURN
!
      ELSE IF( N.EQ.0 ) THEN
!
!        ==== Quick return in case N = 0; nothing to do. ====
!
         RETURN
!
      ELSE IF( LQUERY ) THEN
!
!        ==== Quick return in case of a workspace query ====
!
         CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, &
      &                LDZ, WORK, LWORK, INFO )
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
         WORK( 1 ) = DCMPLX( MAX( DBLE( WORK( 1 ) ), DBLE( MAX( 1, &
      &               N ) ) ), RZERO )
         RETURN
!
      ELSE
!
!        ==== copy eigenvalues isolated by ZGEBAL ====
!
         IF( ILO.GT.1 ) &
      &      CALL ZCOPY( ILO-1, H, LDH+1, W, 1 )
         IF( IHI.LT.N ) &
      &      CALL ZCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )
!
!        ==== Initialize Z, if requested ====
!
         IF( INITZ ) &
      &      CALL ZLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
!
!        ==== Quick return if possible ====
!
         IF( ILO.EQ.IHI ) THEN
            W( ILO ) = H( ILO, ILO )
            RETURN
         END IF
!
!        ==== ZLAHQR/ZLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'ZHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, &
      &          ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====
!
         IF( N.GT.NMIN ) THEN
            CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, &
      &                   Z, LDZ, WORK, LWORK, INFO )
         ELSE
!
!           ==== Small matrix ====
!
            CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, &
      &                   Z, LDZ, INFO )
!
            IF( INFO.GT.0 ) THEN
!
!              ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
!              .    when ZLAHQR fails. ====
!
               KBOT = INFO
!
               IF( N.GE.NL ) THEN
!
!                 ==== Larger matrices have enough subdiagonal scratch
!                 .    space to call ZLAQR0 directly. ====
!
                  CALL ZLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, &
      &                         ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
!
               ELSE
!
!                 ==== Tiny matrices don't have enough subdiagonal
!                 .    scratch space to benefit from ZLAQR0.  Hence,
!                 .    tiny matrices must be copied into a larger
!                 .    array before calling ZLAQR0. ====
!
                  CALL ZLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL ZLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), &
      &                         NL )
                  CALL ZLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, &
      &                         ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                  IF( WANTT .OR. INFO.NE.0 ) &
      &               CALL ZLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
!
!        ==== Clear out the trash, if necessary. ====
!
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) &
      &      CALL ZLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
!
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
!
         WORK( 1 ) = DCMPLX( MAX( DBLE( MAX( 1, N ) ), &
      &               DBLE( WORK( 1 ) ) ), RZERO )
      END IF
!
!     ==== End of ZHSEQR ====
!
      END SUBROUTINE ZHSEQR
!> \brief \b ZLACGV conjugates a complex vector.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLACGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacgv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacgv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacgv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACGV( N, X, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACGV conjugates a complex vector of length N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the vector X.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
!>                         (1+(N-1)*abs(INCX))
!>          On entry, the vector of length N to be conjugated.
!>          On exit, X is overwritten with conjg(X).
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive elements of X.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACGV( N, X, INCX )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IOFF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 ) &
      &      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZLACGV
!
      END SUBROUTINE ZLACGV
!> \brief \b ZLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacpy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
!>          is accessed; if UPLO = 'L', only the lower trapezium is
!>          accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
!
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLACPY
!
      END SUBROUTINE ZLACPY
!> \brief \b ZLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       COMPLEX*16     FUNCTION ZLADIV( X, Y )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!> will not overflow on an intermediary step unless the results
!> overflows.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX*16
!>          The complex scalars X and Y.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      COMPLEX*16     FUNCTION ZLADIV( X, Y )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16         X, Y
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
!     ..
!     .. Executable Statements ..
!
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR, &
      &             ZI )
      ZLADIV = DCMPLX( ZR, ZI )
!
      RETURN
!
!     End of ZLADIV
!
      END FUNCTION ZLADIV
!> \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAHQR is an auxiliary routine called by CHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by CHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper triangular in rows and
!>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
!>          ZLAHQR works primarily with the Hessenberg submatrix in rows
!>          and columns ILO to IHI, but applies transformations to all of
!>          H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., then H
!>          is upper triangular in rows and columns ILO:IHI.  If INFO
!>          is zero and if WANTT is .FALSE., then the contents of H
!>          are unspecified on exit.  The output state of H in case
!>          INF is positive is below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          The computed eigenvalues ILO to IHI are stored in the
!>          corresponding elements of W. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by CHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:   successful exit
!>           > 0:   if INFO = i, ZLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of W contain
!>                  those eigenvalues which have been successfully
!>                  computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix
!>                  rows and columns ILO through INFO of the final,
!>                  output value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of ZLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
      &                   IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
!     ..
!
!  =========================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE, HALF
      PARAMETER          ( RZERO = 0.0d0, RONE = 1.0d0, HALF = 0.5d0 )
      DOUBLE PRECISION   DAT1
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U, &
      &                   V2, X, Y
      DOUBLE PRECISION   AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX, &
      &                   SAFMIN, SMLNUM, SX, T2, TST, ULP
      INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M, &
      &                   NH, NZ, KDEFL
!     ..
!     .. Local Arrays ..
      COMPLEX*16         V( 2 )
!     ..
!     .. External Functions ..
      COMPLEX*16         ZLADIV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           ZLADIV, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZLARFG, ZSCAL
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
      &   RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
!
!     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 ) &
      &   H( IHI, IHI-2 ) = ZERO
!     ==== ensure that subdiagonal entries are real ====
      IF( WANTT ) THEN
         JLO = 1
         JHI = N
      ELSE
         JLO = ILO
         JHI = IHI
      END IF
      DO 20 I = ILO + 1, IHI
         IF( DIMAG( H( I, I-1 ) ).NE.RZERO ) THEN
!           ==== The following redundant normalization
!           .    avoids problems with both gradual and
!           .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
            SC = DCONJG( SC ) / ABS( SC )
            H( I, I-1 ) = ABS( H( I, I-1 ) )
            CALL ZSCAL( JHI-I+1, SC, H( I, I ), LDH )
            CALL ZSCAL( MIN( JHI, I+1 )-JLO+1, DCONJG( SC ), &
      &                  H( JLO, I ), 1 )
            IF( WANTZ ) &
      &         CALL ZSCAL( IHIZ-ILOZ+1, DCONJG( SC ), Z( ILOZ, I ), 1 )
         END IF
   20 CONTINUE
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITMAX is the total number of QR iterations allowed.
!
      ITMAX = 30 * MAX( 10, NH )
!
!     KDEFL counts the number of iterations since a deflation
!
      KDEFL = 0
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   30 CONTINUE
      IF( I.LT.ILO ) &
      &   GO TO 150
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      L = ILO
      DO 130 ITS = 0, ITMAX
!
!        Look for a single small subdiagonal element.
!
         DO 40 K = I, L + 1, -1
            IF( CABS1( H( K, K-1 ) ).LE.SMLNUM ) &
      &         GO TO 50
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO ) &
      &            TST = TST + ABS( DBLE( H( K-1, K-2 ) ) )
               IF( K+1.LE.IHI ) &
      &            TST = TST + ABS( DBLE( H( K+1, K ) ) )
            END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some examples.  ====
            IF( ABS( DBLE( H( K, K-1 ) ) ).LE.ULP*TST ) THEN
               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               AA = MAX( CABS1( H( K, K ) ), &
      &              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( CABS1( H( K, K ) ), &
      &              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM, &
      &             ULP*( BB*( AA / S ) ) ) )GO TO 50
            END IF
   40    CONTINUE
   50    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 has split off.
!
         IF( L.GE.I ) &
      &      GO TO 140
         KDEFL = KDEFL + 1
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
!
!           Exceptional shift.
!
            S = DAT1*ABS( DBLE( H( I, I-1 ) ) )
            T = S + H( I, I )
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
!
!           Exceptional shift.
!
            S = DAT1*ABS( DBLE( H( L+1, L ) ) )
            T = S + H( L, L )
         ELSE
!
!           Wilkinson's shift.
!
            T = H( I, I )
            U = SQRT( H( I-1, I ) )*SQRT( H( I, I-1 ) )
            S = CABS1( U )
            IF( S.NE.RZERO ) THEN
               X = HALF*( H( I-1, I-1 )-T )
               SX = CABS1( X )
               S = MAX( S, CABS1( X ) )
               Y = S*SQRT( ( X / S )**2+( U / S )**2 )
               IF( SX.GT.RZERO ) THEN
                  IF( DBLE( X / SX )*DBLE( Y )+DIMAG( X / SX )* &
      &                DIMAG( Y ).LT.RZERO )Y = -Y
               END IF
               T = T - U*ZLADIV( U, ( X+Y ) )
            END IF
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 60 M = I - 1, L + 1, -1
!
!           Determine the effect of starting the single-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = DBLE( H( M+1, M ) )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = DBLE( H( M, M-1 ) )
            IF( ABS( H10 )*ABS( H21 ).LE.ULP* &
      &          ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) ) &
      &          GO TO 70
   60    CONTINUE
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = DBLE( H( L+1, L ) )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
   70    CONTINUE
!
!        Single-shift QR step
!
         DO 120 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix.
!
!           V(2) is always real before the call to ZLARFG, and hence
!           after the call T2 ( = T1*V(2) ) is also real.
!
            IF( K.GT.M ) &
      &         CALL ZCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = DBLE( T1*V2 )
!
!           Apply G from the left to transform the rows of the matrix
!           in columns K to I2.
!
            DO 80 J = K, I2
               SUM = DCONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
   80       CONTINUE
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+2,I).
!
            DO 90 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
   90       CONTINUE
!
            IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
               DO 100 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
  100          CONTINUE
            END IF
!
            IF( K.EQ.M .AND. M.GT.L ) THEN
!
!              If the QR step was started at row M > L because two
!              consecutive small subdiagonals were found, then extra
!              scaling must be performed to ensure that H(M,M-1) remains
!              real.
!
               TEMP = ONE - T1
               TEMP = TEMP / ABS( TEMP )
               H( M+1, M ) = H( M+1, M )*DCONJG( TEMP )
               IF( M+2.LE.I ) &
      &            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 110 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J ) &
      &                  CALL ZSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL ZSCAL( J-I1, DCONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL ZSCAL( NZ, DCONJG( TEMP ), Z( ILOZ, J ), &
      &                              1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
!
!        Ensure that H(I,I-1) is real.
!
         TEMP = H( I, I-1 )
         IF( DIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I ) &
      &         CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH )
            CALL ZSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL ZSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
            END IF
         END IF
!
  130 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  140 CONTINUE
!
!     H(I,I-1) is negligible: one eigenvalue has converged.
!
      W( I ) = H( I, I )
!     reset deflation counter
      KDEFL = 0
!
!     return to start of the main loop with new value of I.
!
      I = L - 1
      GO TO 30
!
  150 CONTINUE
      RETURN
!
!     End of ZLAHQR
!
      END SUBROUTINE ZLAHQR
!> \brief \b ZLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAHR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LDT, LDY, N, NB
!       ..
!       .. Array Arguments ..
!       COMPLEX*16        A( LDA, * ), T( LDT, NB ), TAU( NB ),
!      $                   Y( LDY, NB )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)
!> matrix A so that elements below the k-th subdiagonal are zero. The
!> reduction is performed by an unitary similarity transformation
!> Q**H * A * Q. The routine returns the matrices V and T which determine
!> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
!>
!> This is an auxiliary routine called by ZGEHRD.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The offset for the reduction. Elements below the k-th
!>          subdiagonal in the first NB columns are reduced to zero.
!>          K < N.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N-K+1)
!>          On entry, the n-by-(n-k+1) general matrix A.
!>          On exit, the elements on and above the k-th subdiagonal in
!>          the first NB columns are overwritten with the corresponding
!>          elements of the reduced matrix; the elements below the k-th
!>          subdiagonal, with the array TAU, represent the matrix Q as a
!>          product of elementary reflectors. The other columns of A are
!>          unchanged. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (NB)
!>          The scalar factors of the elementary reflectors. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,NB)
!>          The upper triangular matrix T.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension (LDY,NB)
!>          The n-by-nb matrix Y.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y. LDY >= N.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of nb elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!>  A(i+k+1:n,i), and tau in TAU(i).
!>
!>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!>  V which is needed, with T and Y, to apply the transformation to the
!>  unreduced part of the matrix, using an update of the form:
!>  A := (I - V*T*V**H) * (A - Y*V**H).
!>
!>  The contents of A on exit are illustrated by the following example
!>  with n = 7, k = 3 and nb = 2:
!>
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( h   h   a   a   a )
!>     ( v1  h   a   a   a )
!>     ( v1  v2  a   a   a )
!>     ( v1  v2  a   a   a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This subroutine is a slight modification of LAPACK-3.0's ZLAHRD
!>  incorporating improvements proposed by Quintana-Orti and Van de
!>  Gejin. Note that the entries of A(1:K,2:NB) differ from those
!>  returned by the original LAPACK-3.0's ZLAHRD routine. (This
!>  subroutine is not backward compatible with LAPACK-3.0's ZLAHRD.)
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
!>  performance of reduction to Hessenberg form," ACM Transactions on
!>  Mathematical Software, 32(2):180-194, June 2006.
!>
!  =====================================================================
      SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), T( LDT, NB ), TAU( NB ), &
      &                   Y( LDY, NB )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
      &                     ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16        EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZGEMM, ZGEMV, ZLACPY, &
      &                   ZLARFG, ZSCAL, ZTRMM, ZTRMV, ZLACGV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
      &   RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(K+1:N,I)
!
!           Update I-th column of A - Y * V**H
!
            CALL ZLACGV( I-1, A( K+I-1, 1 ), LDA )
            CALL ZGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY, &
      &                  A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
            CALL ZLACGV( I-1, A( K+I-1, 1 ), LDA )
!
!           Apply I - V * T**H * V**H to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1**H * b1
!
            CALL ZCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL ZTRMV( 'Lower', 'Conjugate transpose', 'UNIT', &
      &                  I-1, A( K+1, 1 ), &
      &                  LDA, T( 1, NB ), 1 )
!
!           w := w + V2**H * b2
!
            CALL ZGEMV( 'Conjugate transpose', N-K-I+1, I-1, &
      &                  ONE, A( K+I, 1 ), &
      &                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
!
!           w := T**H * w
!
            CALL ZTRMV( 'Upper', 'Conjugate transpose', 'NON-UNIT', &
      &                  I-1, T, LDT, &
      &                  T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL ZGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE, &
      &                  A( K+I, 1 ), &
      &                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL ZTRMV( 'Lower', 'NO TRANSPOSE', &
      &                  'UNIT', I-1, &
      &                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL ZAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(I) to annihilate
!        A(K+I+1:N,I)
!
         CALL ZLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
      &                TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
!
!        Compute  Y(K+1:N,I)
!
         CALL ZGEMV( 'NO TRANSPOSE', N-K, N-K-I+1, &
      &               ONE, A( K+1, I+1 ), &
      &               LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
         CALL ZGEMV( 'Conjugate transpose', N-K-I+1, I-1, &
      &               ONE, A( K+I, 1 ), LDA, &
      &               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL ZGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, &
      &               Y( K+1, 1 ), LDY, &
      &               T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
         CALL ZSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
!
!        Compute T(1:I,I)
!
         CALL ZSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL ZTRMV( 'Upper', 'No Transpose', 'NON-UNIT', &
      &               I-1, T, LDT, &
      &               T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
!     Compute Y(1:K,1:NB)
!
      CALL ZLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
      CALL ZTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE', &
      &            'UNIT', K, NB, &
      &            ONE, A( K+1, 1 ), LDA, Y, LDY )
      IF( N.GT.K+NB ) &
      &   CALL ZGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, &
      &               NB, N-K-NB, ONE, &
      &               A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y, &
      &               LDY )
      CALL ZTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE', &
      &            'NON-UNIT', K, NB, &
      &            ONE, T, LDT, Y, LDY )
!
      RETURN
!
!     End of ZLAHR2
!
      END SUBROUTINE ZLAHR2
!> \brief \b ZLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlange.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   WORK( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex matrix A.
!> \endverbatim
!>
!> \return ZLANGE
!> \verbatim
!>
!>    ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in ZLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          ZLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          ZLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!>          referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16GEauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               TEMP = ABS( A( I, J ) )
               IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            IF( VALUE.LT.SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL ZLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANGE = VALUE
      RETURN
!
!     End of ZLANGE
!
      END FUNCTION ZLANGE
!> \brief \b ZLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAQR0 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**H, where T is an upper triangular matrix (the
!>    Schur form), and Z is the unitary matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input unitary
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to ZGEBAL, and then passed to ZGEHRD when the
!>           matrix output by ZGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H
!>           contains the upper triangular matrix T from the Schur
!>           decomposition (the Schur form). If INFO = 0 and WANT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO > 0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i > j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
!>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>           Specify the rows of Z to which transformations must be
!>           applied if WANTZ is .TRUE..
!>           1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then ZLAQR0 does a workspace query.
!>           In this case, ZLAQR0 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             > 0:  if INFO = i, ZLAQR0 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!>
!  =====================================================================
      SUBROUTINE ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
      &                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    ZLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 15 )
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by varying the size of the
!     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      PARAMETER          ( KEXNW = 5 )
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    ====
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 6 )
!
!     ==== The constant WILK1 is used to form the exceptional
!     .    shifts. ====
      DOUBLE PRECISION   WILK1
      PARAMETER          ( WILK1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
      &                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, &
      &                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, &
      &                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLACPY, ZLAHQR, ZLAQR3, ZLAQR4, ZLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, INT, MAX, MIN, MOD, &
      &                   SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
      INFO = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
!
      IF( N.LE.NTINY ) THEN
!
!        ==== Tiny matrices must use ZLAHQR. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 ) &
      &      CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
      &                   IHIZ, Z, LDZ, INFO )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         INFO = 0
!
!        ==== Set up job flags for ILAENV. ====
!
         IF( WANTT ) THEN
            JBCMPZ( 1: 1 ) = 'S'
         ELSE
            JBCMPZ( 1: 1 ) = 'E'
         END IF
         IF( WANTZ ) THEN
            JBCMPZ( 2: 2 ) = 'V'
         ELSE
            JBCMPZ( 2: 2 ) = 'N'
         END IF
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 15, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         NWR = ILAENV( 13, 'ZLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'ZLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to ZLAQR3 ====
!
         CALL ZLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, &
      &                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, &
      &                LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR3) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DCMPLX( LWKOPT, 0 )
            RETURN
         END IF
!
!        ==== ZLAHQR/ZLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'ZLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'ZLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'ZLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         NDFL = 1
!
!        ==== ITMAX = iteration limit ====
!
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
!
!        ==== Last row and column in the active block ====
!
         KBOT = IHI
!
!        ==== Main Loop ====
!
         DO 70 IT = 1, ITMAX
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( KBOT.LT.ILO ) &
      &         GO TO 80
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) &
      &            GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
!
!           ==== Select deflation window size:
!           .    Typical Case:
!           .      If possible and advisable, nibble the entire
!           .      active block.  If not, use size MIN(NWR,NWMAX)
!           .      or MIN(NWR+1,NWMAX) depending upon which has
!           .      the smaller corresponding subdiagonal entry
!           .      (a heuristic).
!           .
!           .    Exceptional Case:
!           .      If there have been no deflations in KEXNW or
!           .      more iterations, then vary the deflation window
!           .      size.   At first, because, larger windows are,
!           .      in general, more powerful than smaller ones,
!           .      rapidly increase the window to the maximum possible.
!           .      Then, gradually reduce the window size. ====
!
            NH = KBOT - KTOP + 1
            NWUPBD = MIN( NH, NWMAX )
            IF( NDFL.LT.KEXNW ) THEN
               NW = MIN( NWUPBD, NWR )
            ELSE
               NW = MIN( NWUPBD, 2*NW )
            END IF
            IF( NW.LT.NWMAX ) THEN
               IF( NW.GE.NH-1 ) THEN
                  NW = NH
               ELSE
                  KWTOP = KBOT - NW + 1
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT. &
      &                CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 ) &
      &            NDEC = 0
               NW = NW - NDEC
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
!
!           ==== Aggressive early deflation ====
!
            CALL ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
      &                   IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO, &
      &                   H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK, &
      &                   LWORK )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            KBOT = KBOT - LD
!
!           ==== KS points to the shifts. ====
!
            KS = KBOT - LS + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- &
      &          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if ZLAQR3
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    ZLAQR3 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, KS + 1, -2
                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
                     W( I-1 ) = W( I )
   30             CONTINUE
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use ZLAQR4 or
!                 .    ZLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL ZLACPY( 'A', NS, NS, H( KS, KS ), LDH, &
      &                            H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL ZLAQR4( .false., .false., NS, 1, NS, &
      &                               H( KT, 1 ), LDH, W( KS ), 1, 1, &
      &                               ZDUM, 1, WORK, LWORK, INF )
                     ELSE
                        CALL ZLAHQR( .false., .false., NS, 1, NS, &
      &                               H( KT, 1 ), LDH, W( KS ), 1, 1, &
      &                               ZDUM, 1, INF )
                     END IF
                     KS = KS + INF
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  Scale to avoid
!                    .    overflows, underflows and subnormals.
!                    .    (The scale factor S can not be zero,
!                    .    because H(KBOT,KBOT-1) is nonzero.) ====
!
                     IF( KS.GE.KBOT ) THEN
                        S = CABS1( H( KBOT-1, KBOT-1 ) ) + &
      &                      CABS1( H( KBOT, KBOT-1 ) ) + &
      &                      CABS1( H( KBOT-1, KBOT ) ) + &
      &                      CABS1( H( KBOT, KBOT ) )
                        AA = H( KBOT-1, KBOT-1 ) / S
                        CC = H( KBOT, KBOT-1 ) / S
                        BB = H( KBOT-1, KBOT ) / S
                        DD = H( KBOT, KBOT ) / S
                        TR2 = ( AA+DD ) / TWO
                        DET = ( AA-TR2 )*( DD-TR2 ) - BB*CC
                        RTDISC = SQRT( -DET )
                        W( KBOT-1 ) = ( TR2+RTDISC )*S
                        W( KBOT ) = ( TR2-RTDISC )*S
!
                        KS = KBOT - 1
                     END IF
                  END IF
!
                  IF( KBOT-KS+1.GT.NS ) THEN
!
!                    ==== Sort the shifts (Helps a little) ====
!
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) &
      &                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) ) &
      &                          THEN
                              SORTED = .false.
                              SWAP = W( I )
                              W( I ) = W( I+1 )
                              W( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
!
!              ==== If there are only two shifts, then use
!              .    only one.  ====
!
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT. &
      &                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                     W( KBOT-1 ) = W( KBOT )
                  ELSE
                     W( KBOT ) = W( KBOT-1 )
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnitude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               KDU = 2*NS
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, &
      &                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK, &
      &                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH, &
      &                      NHO, H( KU, KWH ), LDH )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
!
!           ==== End of main loop ====
   70    CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         INFO = KBOT
   80    CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      WORK( 1 ) = DCMPLX( LWKOPT, 0 )
!
!     ==== End of ZLAQR0 ====
!
      END SUBROUTINE ZLAQR0
!> \brief \b ZLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         S1, S2
!       INTEGER            LDH, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a
!>      scalar multiple of the first column of the product
!>
!>      (*)  K = (H - s1*I)*(H - s2*I)
!>
!>      scaling to avoid overflows and most underflows.
!>
!>      This is useful for starting double implicit shift bulges
!>      in the QR algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>              Order of the matrix H. N must be either 2 or 3.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>              The 2-by-2 or 3-by-3 matrix H in (*).
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>              The leading dimension of H as declared in
!>              the calling procedure.  LDH >= N
!> \endverbatim
!>
!> \param[in] S1
!> \verbatim
!>          S1 is COMPLEX*16
!> \endverbatim
!>
!> \param[in] S2
!> \verbatim
!>          S2 is COMPLEX*16
!>
!>          S1 and S2 are the shifts defining K in (*) above.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (N)
!>              A scalar multiple of the first column of the
!>              matrix K in (*).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16         S1, S2
      INTEGER            LDH, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), V( * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         CDUM, H21S, H31S
      DOUBLE PRECISION   S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.NE.2 .AND. N.NE.3 ) THEN
         RETURN
      END IF
!
      IF( N.EQ.2 ) THEN
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) )
         IF( S.EQ.RZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-S1 )* &
      &               ( ( H( 1, 1 )-S2 ) / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 )
         END IF
      ELSE
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) ) + &
      &       CABS1( H( 3, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-S1 )*( ( H( 1, 1 )-S2 ) / S ) + &
      &               H( 1, 2 )*H21S + H( 1, 3 )*H31S
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 ) + H( 2, 3 )*H31S
            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-S1-S2 ) + H21S*H( 3, 2 )
         END IF
      END IF
      END SUBROUTINE ZLAQR1
!> \brief \b ZLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
!                          NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAQR2 is identical to ZLAQR3 except that it avoids
!>    recursion by calling ZLAHQR instead of ZLAQR4.
!>
!>    Aggressive early deflation:
!>
!>    ZLAQR2 accepts as input an upper Hessenberg matrix
!>    H and performs an unitary similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an unitary similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the unitary matrix Z is updated so
!>          so that the unitary Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the unitary matrix Z.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!>          KBOT and KTOP together determine an isolated block
!>          along the diagonal of the Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>          It is assumed without a check that either
!>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!>          determine an isolated block along the diagonal of the
!>          Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          Deflation window size.  1 <= NW <= (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by a unitary
!>          similarity transformation, perturbed, and the returned
!>          to Hessenberg form that (it is to be hoped) has some
!>          zero subdiagonal entries.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N <= LDH
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the unitary
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 <= LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SH
!> \verbatim
!>          SH is COMPLEX*16 array, dimension (KBOT)
!>          On output, approximate eigenvalues that may
!>          be used for shifts are stored in SH(KBOT-ND-NS+1)
!>          through SR(KBOT-ND).  Converged eigenvalues are
!>          stored in SH(KBOT-ND+1) through SH(KBOT).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>          The number of columns of T.  NH >= NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW <= LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>          The number of rows of work array WV available for
!>          workspace.  NV >= NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX*16 array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; ZLAQR2
!>          only estimates the optimal workspace size for the given
!>          values of N, NW, KTOP and KBOT.  The estimate is returned
!>          in WORK(1).  No error message related to LWORK is issued
!>          by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
      &                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, &
      &                   NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
      &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), &
      &                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, RONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN, &
      &                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZGEHRD, ZGEMM, ZLACPY, ZLAHQR, &
      &                   ZLARF, ZLARFG, ZLASET, ZTREXC, ZUNMHR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
!
!        ==== Workspace query call to ZGEHRD ====
!
         CALL ZGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to ZUNMHR ====
!
         CALL ZUNMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, &
      &                WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
!
!        ==== Optimal workspace ====
!
         LWKOPT = JW + MAX( LWK1, LWK2 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DCMPLX( LWKOPT, 0 )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT ) &
      &   RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) &
      &   RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Setup deflation window ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
!
      IF( KBOT.EQ.KWTOP ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         SH( KWTOP ) = H( KWTOP, KWTOP )
         NS = 1
         ND = 0
         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP, &
      &       KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) &
      &         H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         WORK( 1 ) = ONE
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL ZLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL ZCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL ZLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      CALL ZLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1, &
      &             JW, V, LDV, INFQR )
!
!     ==== Deflation detection loop ====
!
      NS = JW
      ILST = INFQR + 1
      DO 10 KNT = INFQR + 1, JW
!
!        ==== Small spike tip deflation test ====
!
         FOO = CABS1( T( NS, NS ) )
         IF( FOO.EQ.RZERO ) &
      &      FOO = CABS1( S )
         IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) &
      &        THEN
!
!           ==== One more converged eigenvalue ====
!
            NS = NS - 1
         ELSE
!
!           ==== One undeflatable eigenvalue.  Move it up out of the
!           .    way.   (ZTREXC can not fail in this case.) ====
!
            IFST = NS
            CALL ZTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
            ILST = ILST + 1
         END IF
   10 CONTINUE
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 ) &
      &   S = ZERO
!
      IF( NS.LT.JW ) THEN
!
!        ==== sorting the diagonal of T improves accuracy for
!        .    graded matrices.  ====
!
         DO 30 I = INFQR + 1, NS
            IFST = I
            DO 20 J = I + 1, NS
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) ) &
      &            IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST ) &
      &         CALL ZTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
   30    CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      DO 40 I = INFQR + 1, JW
         SH( KWTOP+I-1 ) = T( I, I )
   40 CONTINUE
!
!
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL ZCOPY( NS, V, LDV, WORK, 1 )
            DO 50 I = 1, NS
               WORK( I ) = DCONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            CALL ZLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL ZLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL ZLARF( 'L', NS, JW, WORK, 1, DCONJG( TAU ), T, LDT, &
      &                  WORK( JW+1 ) )
            CALL ZLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, &
      &                  WORK( JW+1 ) )
            CALL ZLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, &
      &                  WORK( JW+1 ) )
!
            CALL ZGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
      &                   LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 ) &
      &      H( KWTOP, KWTOP-1 ) = S*DCONJG( V( 1, 1 ) )
         CALL ZLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL ZCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), &
      &               LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO ) &
      &      CALL ZUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, &
      &                   WORK( JW+1 ), LWORK-JW, INFO )
!
!        ==== Update vertical slab in H ====
!
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 60 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL ZGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), &
      &                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL ZLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   60    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL ZGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, &
      &                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL ZLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), &
      &                      LDH )
   70       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL ZGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), &
      &                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL ZLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), &
      &                      LDZ )
   80       CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      ND = JW - NS
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      NS = NS - INFQR
!
!      ==== Return optimal workspace. ====
!
      WORK( 1 ) = DCMPLX( LWKOPT, 0 )
!
!     ==== End of ZLAQR2 ====
!
      END SUBROUTINE ZLAQR2
!> \brief \b ZLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
!                          NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Aggressive early deflation:
!>
!>    ZLAQR3 accepts as input an upper Hessenberg matrix
!>    H and performs an unitary similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an unitary similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the unitary matrix Z is updated so
!>          so that the unitary Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the unitary matrix Z.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!>          KBOT and KTOP together determine an isolated block
!>          along the diagonal of the Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>          It is assumed without a check that either
!>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!>          determine an isolated block along the diagonal of the
!>          Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          Deflation window size.  1 <= NW <= (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by a unitary
!>          similarity transformation, perturbed, and the returned
!>          to Hessenberg form that (it is to be hoped) has some
!>          zero subdiagonal entries.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N <= LDH
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the unitary
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 <= LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SH
!> \verbatim
!>          SH is COMPLEX*16 array, dimension (KBOT)
!>          On output, approximate eigenvalues that may
!>          be used for shifts are stored in SH(KBOT-ND-NS+1)
!>          through SR(KBOT-ND).  Converged eigenvalues are
!>          stored in SH(KBOT-ND+1) through SH(KBOT).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>          The number of columns of T.  NH >= NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW <= LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>          The number of rows of work array WV available for
!>          workspace.  NV >= NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX*16 array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; ZLAQR3
!>          only estimates the optimal workspace size for the given
!>          values of N, NW, KTOP and KBOT.  The estimate is returned
!>          in WORK(1).  No error message related to LWORK is issued
!>          by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
      &                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, &
      &                   NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
      &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), &
      &                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, RONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN, &
      &                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3, &
      &                   LWKOPT, NMIN
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      EXTERNAL           DLAMCH, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZGEHRD, ZGEMM, ZLACPY, ZLAHQR, &
      &                   ZLAQR4, ZLARF, ZLARFG, ZLASET, ZTREXC, ZUNMHR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
!
!        ==== Workspace query call to ZGEHRD ====
!
         CALL ZGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to ZUNMHR ====
!
         CALL ZUNMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, &
      &                WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to ZLAQR4 ====
!
         CALL ZLAQR4( .true., .true., JW, 1, JW, T, LDT, SH, 1, JW, V, &
      &                LDV, WORK, -1, INFQR )
         LWK3 = INT( WORK( 1 ) )
!
!        ==== Optimal workspace ====
!
         LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DCMPLX( LWKOPT, 0 )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT ) &
      &   RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) &
      &   RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Setup deflation window ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
!
      IF( KBOT.EQ.KWTOP ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         SH( KWTOP ) = H( KWTOP, KWTOP )
         NS = 1
         ND = 0
         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP, &
      &       KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) &
      &         H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         WORK( 1 ) = ONE
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL ZLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL ZCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL ZLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      NMIN = ILAENV( 12, 'ZLAQR3', 'SV', JW, 1, JW, LWORK )
      IF( JW.GT.NMIN ) THEN
         CALL ZLAQR4( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1, &
      &                JW, V, LDV, WORK, LWORK, INFQR )
      ELSE
         CALL ZLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1, &
      &                JW, V, LDV, INFQR )
      END IF
!
!     ==== Deflation detection loop ====
!
      NS = JW
      ILST = INFQR + 1
      DO 10 KNT = INFQR + 1, JW
!
!        ==== Small spike tip deflation test ====
!
         FOO = CABS1( T( NS, NS ) )
         IF( FOO.EQ.RZERO ) &
      &      FOO = CABS1( S )
         IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) &
      &        THEN
!
!           ==== One more converged eigenvalue ====
!
            NS = NS - 1
         ELSE
!
!           ==== One undeflatable eigenvalue.  Move it up out of the
!           .    way.   (ZTREXC can not fail in this case.) ====
!
            IFST = NS
            CALL ZTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
            ILST = ILST + 1
         END IF
   10 CONTINUE
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 ) &
      &   S = ZERO
!
      IF( NS.LT.JW ) THEN
!
!        ==== sorting the diagonal of T improves accuracy for
!        .    graded matrices.  ====
!
         DO 30 I = INFQR + 1, NS
            IFST = I
            DO 20 J = I + 1, NS
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) ) &
      &            IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST ) &
      &         CALL ZTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
   30    CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      DO 40 I = INFQR + 1, JW
         SH( KWTOP+I-1 ) = T( I, I )
   40 CONTINUE
!
!
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL ZCOPY( NS, V, LDV, WORK, 1 )
            DO 50 I = 1, NS
               WORK( I ) = DCONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            CALL ZLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL ZLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL ZLARF( 'L', NS, JW, WORK, 1, DCONJG( TAU ), T, LDT, &
      &                  WORK( JW+1 ) )
            CALL ZLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, &
      &                  WORK( JW+1 ) )
            CALL ZLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, &
      &                  WORK( JW+1 ) )
!
            CALL ZGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
      &                   LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 ) &
      &      H( KWTOP, KWTOP-1 ) = S*DCONJG( V( 1, 1 ) )
         CALL ZLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL ZCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), &
      &               LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO ) &
      &      CALL ZUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, &
      &                   WORK( JW+1 ), LWORK-JW, INFO )
!
!        ==== Update vertical slab in H ====
!
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 60 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL ZGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), &
      &                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL ZLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   60    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL ZGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, &
      &                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL ZLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), &
      &                      LDH )
   70       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL ZGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), &
      &                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL ZLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), &
      &                      LDZ )
   80       CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      ND = JW - NS
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      NS = NS - INFQR
!
!      ==== Return optimal workspace. ====
!
      WORK( 1 ) = DCMPLX( LWKOPT, 0 )
!
!     ==== End of ZLAQR3 ====
!
      END SUBROUTINE ZLAQR3
!> \brief \b ZLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAQR4 implements one level of recursion for ZLAQR0.
!>    It is a complete implementation of the small bulge multi-shift
!>    QR algorithm.  It may be called by ZLAQR0 and, for large enough
!>    deflation window size, it may be called by ZLAQR3.  This
!>    subroutine is identical to ZLAQR0 except that it calls ZLAQR2
!>    instead of ZLAQR3.
!>
!>    ZLAQR4 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**H, where T is an upper triangular matrix (the
!>    Schur form), and Z is the unitary matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input unitary
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to ZGEBAL, and then passed to ZGEHRD when the
!>           matrix output by ZGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H
!>           contains the upper triangular matrix T from the Schur
!>           decomposition (the Schur form). If INFO = 0 and WANT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO > 0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i > j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
!>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>           Specify the rows of Z to which transformations must be
!>           applied if WANTZ is .TRUE..
!>           1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then ZLAQR4 does a workspace query.
!>           In this case, ZLAQR4 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             =  0:  successful exit
!>             > 0:  if INFO = i, ZLAQR4 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!>
!  =====================================================================
      SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
      &                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    ZLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 15 )
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by varying the size of the
!     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      PARAMETER          ( KEXNW = 5 )
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    ====
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 6 )
!
!     ==== The constant WILK1 is used to form the exceptional
!     .    shifts. ====
      DOUBLE PRECISION   WILK1
      PARAMETER          ( WILK1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
      &                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, &
      &                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, &
      &                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLACPY, ZLAHQR, ZLAQR2, ZLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, INT, MAX, MIN, MOD, &
      &                   SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
      INFO = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
!
      IF( N.LE.NTINY ) THEN
!
!        ==== Tiny matrices must use ZLAHQR. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 ) &
      &      CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
      &                   IHIZ, Z, LDZ, INFO )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         INFO = 0
!
!        ==== Set up job flags for ILAENV. ====
!
         IF( WANTT ) THEN
            JBCMPZ( 1: 1 ) = 'S'
         ELSE
            JBCMPZ( 1: 1 ) = 'E'
         END IF
         IF( WANTZ ) THEN
            JBCMPZ( 2: 2 ) = 'V'
         ELSE
            JBCMPZ( 2: 2 ) = 'N'
         END IF
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 15, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         NWR = ILAENV( 13, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to ZLAQR2 ====
!
         CALL ZLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, &
      &                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, &
      &                LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR2) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DCMPLX( LWKOPT, 0 )
            RETURN
         END IF
!
!        ==== ZLAHQR/ZLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         NDFL = 1
!
!        ==== ITMAX = iteration limit ====
!
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
!
!        ==== Last row and column in the active block ====
!
         KBOT = IHI
!
!        ==== Main Loop ====
!
         DO 70 IT = 1, ITMAX
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( KBOT.LT.ILO ) &
      &         GO TO 80
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) &
      &            GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
!
!           ==== Select deflation window size:
!           .    Typical Case:
!           .      If possible and advisable, nibble the entire
!           .      active block.  If not, use size MIN(NWR,NWMAX)
!           .      or MIN(NWR+1,NWMAX) depending upon which has
!           .      the smaller corresponding subdiagonal entry
!           .      (a heuristic).
!           .
!           .    Exceptional Case:
!           .      If there have been no deflations in KEXNW or
!           .      more iterations, then vary the deflation window
!           .      size.   At first, because, larger windows are,
!           .      in general, more powerful than smaller ones,
!           .      rapidly increase the window to the maximum possible.
!           .      Then, gradually reduce the window size. ====
!
            NH = KBOT - KTOP + 1
            NWUPBD = MIN( NH, NWMAX )
            IF( NDFL.LT.KEXNW ) THEN
               NW = MIN( NWUPBD, NWR )
            ELSE
               NW = MIN( NWUPBD, 2*NW )
            END IF
            IF( NW.LT.NWMAX ) THEN
               IF( NW.GE.NH-1 ) THEN
                  NW = NH
               ELSE
                  KWTOP = KBOT - NW + 1
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT. &
      &                CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 ) &
      &            NDEC = 0
               NW = NW - NDEC
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
!
!           ==== Aggressive early deflation ====
!
            CALL ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
      &                   IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO, &
      &                   H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK, &
      &                   LWORK )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            KBOT = KBOT - LD
!
!           ==== KS points to the shifts. ====
!
            KS = KBOT - LS + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- &
      &          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if ZLAQR2
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    ZLAQR2 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, KS + 1, -2
                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
                     W( I-1 ) = W( I )
   30             CONTINUE
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use ZLAHQR
!                 .    on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL ZLACPY( 'A', NS, NS, H( KS, KS ), LDH, &
      &                            H( KT, 1 ), LDH )
                     CALL ZLAHQR( .false., .false., NS, 1, NS, &
      &                            H( KT, 1 ), LDH, W( KS ), 1, 1, ZDUM, &
      &                            1, INF )
                     KS = KS + INF
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  Scale to avoid
!                    .    overflows, underflows and subnormals.
!                    .    (The scale factor S can not be zero,
!                    .    because H(KBOT,KBOT-1) is nonzero.) ====
!
                     IF( KS.GE.KBOT ) THEN
                        S = CABS1( H( KBOT-1, KBOT-1 ) ) + &
      &                      CABS1( H( KBOT, KBOT-1 ) ) + &
      &                      CABS1( H( KBOT-1, KBOT ) ) + &
      &                      CABS1( H( KBOT, KBOT ) )
                        AA = H( KBOT-1, KBOT-1 ) / S
                        CC = H( KBOT, KBOT-1 ) / S
                        BB = H( KBOT-1, KBOT ) / S
                        DD = H( KBOT, KBOT ) / S
                        TR2 = ( AA+DD ) / TWO
                        DET = ( AA-TR2 )*( DD-TR2 ) - BB*CC
                        RTDISC = SQRT( -DET )
                        W( KBOT-1 ) = ( TR2+RTDISC )*S
                        W( KBOT ) = ( TR2-RTDISC )*S
!
                        KS = KBOT - 1
                     END IF
                  END IF
!
                  IF( KBOT-KS+1.GT.NS ) THEN
!
!                    ==== Sort the shifts (Helps a little) ====
!
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) &
      &                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) ) &
      &                          THEN
                              SORTED = .false.
                              SWAP = W( I )
                              W( I ) = W( I+1 )
                              W( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
!
!              ==== If there are only two shifts, then use
!              .    only one.  ====
!
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT. &
      &                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                     W( KBOT-1 ) = W( KBOT )
                  ELSE
                     W( KBOT ) = W( KBOT-1 )
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnitude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               KDU = 2*NS
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, &
      &                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK, &
      &                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH, &
      &                      NHO, H( KU, KWH ), LDH )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
!
!           ==== End of main loop ====
   70    CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         INFO = KBOT
   80    CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      WORK( 1 ) = DCMPLX( LWKOPT, 0 )
!
!     ==== End of ZLAQR4 ====
!
      END SUBROUTINE ZLAQR4
!> \brief \b ZLAQR5 performs a single small-bulge multi-shift QR sweep.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
!                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
!                          WV, LDWV, NH, WH, LDWH )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
!      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
!      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAQR5, called by ZLAQR0, performs a
!>    single small-bulge multi-shift QR sweep.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>             WANTT = .true. if the triangular Schur factor
!>             is being computed.  WANTT is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>             WANTZ = .true. if the unitary Schur factor is being
!>             computed.  WANTZ is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] KACC22
!> \verbatim
!>          KACC22 is INTEGER with value 0, 1, or 2.
!>             Specifies the computation mode of far-from-diagonal
!>             orthogonal updates.
!>        = 0: ZLAQR5 does not accumulate reflections and does not
!>             use matrix-matrix multiply to update far-from-diagonal
!>             matrix entries.
!>        = 1: ZLAQR5 accumulates reflections and uses matrix-matrix
!>             multiply to update the far-from-diagonal matrix entries.
!>        = 2: Same as KACC22 = 1. This option used to enable exploiting
!>             the 2-by-2 structure during matrix multiplications, but
!>             this is no longer supported.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>             N is the order of the Hessenberg matrix H upon which this
!>             subroutine operates.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>             These are the first and last rows and columns of an
!>             isolated diagonal block upon which the QR sweep is to be
!>             applied. It is assumed without a check that
!>                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
!>             and
!>                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
!> \endverbatim
!>
!> \param[in] NSHFTS
!> \verbatim
!>          NSHFTS is INTEGER
!>             NSHFTS gives the number of simultaneous shifts.  NSHFTS
!>             must be positive and even.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is COMPLEX*16 array, dimension (NSHFTS)
!>             S contains the shifts of origin that define the multi-
!>             shift QR sweep.  On output S may be reordered.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>             On input H contains a Hessenberg matrix.  On output a
!>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
!>             to the isolated diagonal block in rows and columns KTOP
!>             through KBOT.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>             LDH is the leading dimension of H just as declared in the
!>             calling procedure.  LDH >= MAX(1,N).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>             Specify the rows of Z to which transformations must be
!>             applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,IHIZ)
!>             If WANTZ = .TRUE., then the QR Sweep unitary
!>             similarity transformation is accumulated into
!>             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>             If WANTZ = .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>             LDA is the leading dimension of Z just as declared in
!>             the calling procedure. LDZ >= N.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,NSHFTS/2)
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>             LDV is the leading dimension of V as declared in the
!>             calling procedure.  LDV >= 3.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>             LDU is the leading dimension of U just as declared in the
!>             in the calling subroutine.  LDU >= 2*NSHFTS.
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>             NV is the number of rows in WV agailable for workspace.
!>             NV >= 1.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX*16 array, dimension (LDWV,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>             LDWV is the leading dimension of WV as declared in the
!>             in the calling subroutine.  LDWV >= NV.
!> \endverbatim
!
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>             NH is the number of columns in array WH available for
!>             workspace. NH >= 1.
!> \endverbatim
!>
!> \param[out] WH
!> \verbatim
!>          WH is COMPLEX*16 array, dimension (LDWH,NH)
!> \endverbatim
!>
!> \param[in] LDWH
!> \verbatim
!>          LDWH is INTEGER
!>             Leading dimension of WH just as declared in the
!>             calling procedure.  LDWH >= 2*NSHFTS.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang
!>
!>       Thijs Steel, Department of Computer science,
!>       KU Leuven, Belgium
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang, Optimally packed
!>       chains of bulges in multishift QR algorithms.
!>       ACM Trans. Math. Softw. 40, 2, Article 12 (February 2014).
!>
!  =====================================================================
      SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, &
      &                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, &
      &                   WV, LDWV, NH, WH, LDWH )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
      &                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), &
      &                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ), &
      &                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0d0, RONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      COMPLEX*16         ALPHA, BETA, CDUM, REFSUM
      DOUBLE PRECISION   H11, H12, H21, H22, SAFMAX, SAFMIN, SCL, &
      &                   SMLNUM, TST1, TST2, ULP
      INTEGER            I2, I4, INCOL, J, JBOT, JCOL, JLEN, &
      &                   JROW, JTOP, K, K1, KDU, KMS, KRCOL, &
      &                   M, M22, MBOT, MTOP, NBMPS, NDCOL, &
      &                   NS, NU
      LOGICAL            ACCUM, BMP22
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD
!     ..
!     .. Local Arrays ..
      COMPLEX*16         VT( 3 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZGEMM, ZLACPY, ZLAQR1, ZLARFG, ZLASET, &
      &                   ZTRMM
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     ==== If there are no shifts, then there is nothing to do. ====
!
      IF( NSHFTS.LT.2 ) &
      &   RETURN
!
!     ==== If the active block is empty or 1-by-1, then there
!     .    is nothing to do. ====
!
      IF( KTOP.GE.KBOT ) &
      &   RETURN
!
!     ==== NSHFTS is supposed to be even, but if it is odd,
!     .    then simply reduce it by one.  ====
!
      NS = NSHFTS - MOD( NSHFTS, 2 )
!
!     ==== Machine constants for deflation ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Use accumulated reflections to update far-from-diagonal
!     .    entries ? ====
!
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
!
!     ==== clear trash ====
!
      IF( KTOP+2.LE.KBOT ) &
      &   H( KTOP+2, KTOP ) = ZERO
!
!     ==== NBMPS = number of 2-shift bulges in the chain ====
!
      NBMPS = NS / 2
!
!     ==== KDU = width of slab ====
!
      KDU = 4*NBMPS
!
!     ==== Create and chase chains of NBMPS bulges ====
!
      DO 180 INCOL = KTOP - 2*NBMPS + 1, KBOT - 2, 2*NBMPS
!
!        JTOP = Index from which updates from the right start.
!
         IF( ACCUM ) THEN
            JTOP = MAX( KTOP, INCOL )
         ELSE IF( WANTT ) THEN
            JTOP = 1
         ELSE
            JTOP = KTOP
         END IF
!
         NDCOL = INCOL + KDU
         IF( ACCUM ) &
      &      CALL ZLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
!
!        ==== Near-the-diagonal bulge chase.  The following loop
!        .    performs the near-the-diagonal part of a small bulge
!        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
!        .    chunk extends from column INCOL to column NDCOL
!        .    (including both column INCOL and column NDCOL). The
!        .    following loop chases a 2*NBMPS+1 column long chain of
!        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
!        .    may be less than KTOP and and NDCOL may be greater than
!        .    KBOT indicating phantom columns from which to chase
!        .    bulges before they are actually introduced or to which
!        .    to chase bulges beyond column KBOT.)  ====
!
         DO 145 KRCOL = INCOL, MIN( INCOL+2*NBMPS-1, KBOT-2 )
!
!           ==== Bulges number MTOP to MBOT are active double implicit
!           .    shift bulges.  There may or may not also be small
!           .    2-by-2 bulge, if there is room.  The inactive bulges
!           .    (if any) must wait until the active bulges have moved
!           .    down the diagonal to make room.  The phantom matrix
!           .    paradigm described above helps keep track.  ====
!
            MTOP = MAX( 1, ( KTOP-KRCOL ) / 2+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 2 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+2*( M22-1 ) ).EQ. &
      &              ( KBOT-2 )
!
!           ==== Generate reflections to chase the chain right
!           .    one column.  (The minimum value of K is KTOP-1.) ====
!
            IF ( BMP22 ) THEN
!
!              ==== Special case: 2-by-2 reflection at bottom treated
!              .    separately ====
!
               K = KRCOL + 2*( M22-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL ZLAQR1( 2, H( K+1, K+1 ), LDH, S( 2*M22-1 ), &
      &                         S( 2*M22 ), V( 1, M22 ) )
                  BETA = V( 1, M22 )
                  CALL ZLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  CALL ZLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               END IF

!
!              ==== Perform update from right within 
!              .    computational window. ====
!
               DO 30 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )* &
      &                     H( J, K+2 ) )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
                  H( J, K+2 ) = H( J, K+2 ) - &
      &                          REFSUM*DCONJG( V( 2, M22 ) )
   30          CONTINUE
!
!              ==== Perform update from left within 
!              .    computational window. ====
!
               IF( ACCUM ) THEN
                  JBOT = MIN( NDCOL, KBOT )
               ELSE IF( WANTT ) THEN
                  JBOT = N
               ELSE
                  JBOT = KBOT
               END IF
               DO 40 J = K+1, JBOT
                  REFSUM = DCONJG( V( 1, M22 ) )* &
      &                     ( H( K+1, J )+DCONJG( V( 2, M22 ) )* &
      &                     H( K+2, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
   40          CONTINUE
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF( K.GE.KTOP ) THEN
                  IF( H( K+1, K ).NE.ZERO ) THEN
                     TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
                     IF( TST1.EQ.RZERO ) THEN
                        IF( K.GE.KTOP+1 ) &
      &                     TST1 = TST1 + CABS1( H( K, K-1 ) )
                        IF( K.GE.KTOP+2 ) &
      &                     TST1 = TST1 + CABS1( H( K, K-2 ) )
                        IF( K.GE.KTOP+3 ) &
      &                     TST1 = TST1 + CABS1( H( K, K-3 ) )
                        IF( K.LE.KBOT-2 ) &
      &                     TST1 = TST1 + CABS1( H( K+2, K+1 ) )
                        IF( K.LE.KBOT-3 ) &
      &                     TST1 = TST1 + CABS1( H( K+3, K+1 ) )
                        IF( K.LE.KBOT-4 ) &
      &                     TST1 = TST1 + CABS1( H( K+4, K+1 ) )
                     END IF
                     IF( CABS1( H( K+1, K ) ) &
      &                   .LE.MAX( SMLNUM, ULP*TST1 ) ) THEN
                        H12 = MAX( CABS1( H( K+1, K ) ), &
      &                     CABS1( H( K, K+1 ) ) )
                        H21 = MIN( CABS1( H( K+1, K ) ), &
      &                     CABS1( H( K, K+1 ) ) )
                        H11 = MAX( CABS1( H( K+1, K+1 ) ), &
      &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                        H22 = MIN( CABS1( H( K+1, K+1 ) ), &
      &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                        SCL = H11 + H12
                        TST2 = H22*( H11 / SCL )
!
                        IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE. &
      &                      MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                     END IF
                  END IF
               END IF
!
!              ==== Accumulate orthogonal transformations. ====
!
               IF( ACCUM ) THEN
                  KMS = K - INCOL
                  DO 50 J = MAX( 1, KTOP-INCOL ), KDU
                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+ &
      &                        V( 2, M22 )*U( J, KMS+2 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - &
      &                               REFSUM*DCONJG( V( 2, M22 ) )
  50                 CONTINUE
               ELSE IF( WANTZ ) THEN
                  DO 60 J = ILOZ, IHIZ
                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )* &
      &                        Z( J, K+2 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - &
      &                             REFSUM*DCONJG( V( 2, M22 ) )
  60              CONTINUE
               END IF
            END IF
!
!           ==== Normal case: Chain of 3-by-3 reflections ====
!
            DO 80 M = MBOT, MTOP, -1
               K = KRCOL + 2*( M-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL ZLAQR1( 3, H( KTOP, KTOP ), LDH, S( 2*M-1 ), &
      &                         S( 2*M ), V( 1, M ) )
                  ALPHA = V( 1, M )
                  CALL ZLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
               ELSE
!
!                 ==== Perform delayed transformation of row below
!                 .    Mth bulge. Exploit fact that first two elements
!                 .    of row are actually zero. ====
!
                  REFSUM = V( 1, M )*V( 3, M )*H( K+3, K+2 )
                  H( K+3, K   ) = -REFSUM
                  H( K+3, K+1 ) = -REFSUM*DCONJG( V( 2, M ) )
                  H( K+3, K+2 ) = H( K+3, K+2 ) - &
      &                            REFSUM*DCONJG( V( 3, M ) )
!
!                 ==== Calculate reflection to move
!                 .    Mth bulge one step. ====
!
                  BETA      = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  CALL ZLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
!
!                 ==== A Bulge may collapse because of vigilant
!                 .    deflation or destructive underflow.  In the
!                 .    underflow case, try the two-small-subdiagonals
!                 .    trick to try to reinflate the bulge.  ====
!
                  IF( H( K+3, K ).NE.ZERO .OR. H( K+3, K+1 ).NE. &
      &                ZERO .OR. H( K+3, K+2 ).EQ.ZERO ) THEN
!
!                    ==== Typical case: not collapsed (yet). ====
!
                     H( K+1, K ) = BETA
                     H( K+2, K ) = ZERO
                     H( K+3, K ) = ZERO
                  ELSE
!
!                    ==== Atypical case: collapsed.  Attempt to
!                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
!                    .    If the fill resulting from the new
!                    .    reflector is too large, then abandon it.
!                    .    Otherwise, use the new one. ====
!
                     CALL ZLAQR1( 3, H( K+1, K+1 ), LDH, S( 2*M-1 ), &
      &                            S( 2*M ), VT )
                     ALPHA = VT( 1 )
                     CALL ZLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
                     REFSUM = DCONJG( VT( 1 ) )* &
      &                        ( H( K+1, K )+DCONJG( VT( 2 ) )* &
      &                        H( K+2, K ) )
!
                     IF( CABS1( H( K+2, K )-REFSUM*VT( 2 ) )+ &
      &                   CABS1( REFSUM*VT( 3 ) ).GT.ULP* &
      &                   ( CABS1( H( K, K ) )+CABS1( H( K+1, &
      &                   K+1 ) )+CABS1( H( K+2, K+2 ) ) ) ) THEN
!
!                       ==== Starting a new bulge here would
!                       .    create non-negligible fill.  Use
!                       .    the old one with trepidation. ====
!
                        H( K+1, K ) = BETA
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                     ELSE
!
!                       ==== Starting a new bulge here would
!                       .    create only negligible fill.
!                       .    Replace the old reflector with
!                       .    the new one. ====
!
                        H( K+1, K ) = H( K+1, K ) - REFSUM
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                        V( 1, M ) = VT( 1 )
                        V( 2, M ) = VT( 2 )
                        V( 3, M ) = VT( 3 )
                     END IF
                  END IF
               END IF
!
!              ====  Apply reflection from the right and
!              .     the first column of update from the left.
!              .     These updates are required for the vigilant
!              .     deflation check. We still delay most of the
!              .     updates from the left for efficiency. ====
!
               DO 70 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )* &
      &                     H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
                  H( J, K+2 ) = H( J, K+2 ) - &
      &                          REFSUM*DCONJG( V( 2, M ) )
                  H( J, K+3 ) = H( J, K+3 ) - &
      &                          REFSUM*DCONJG( V( 3, M ) )
   70          CONTINUE
!
!              ==== Perform update from left for subsequent
!              .    column. ====
!
               REFSUM =  DCONJG( V( 1, M ) )*( H( K+1, K+1 ) &
      &                  +DCONJG( V( 2, M ) )*H( K+2, K+1 ) &
      &                  +DCONJG( V( 3, M ) )*H( K+3, K+1 ) )
               H( K+1, K+1 ) = H( K+1, K+1 ) - REFSUM
               H( K+2, K+1 ) = H( K+2, K+1 ) - REFSUM*V( 2, M )
               H( K+3, K+1 ) = H( K+3, K+1 ) - REFSUM*V( 3, M )
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF( K.LT.KTOP) &
      &              CYCLE
               IF( H( K+1, K ).NE.ZERO ) THEN
                  TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
                  IF( TST1.EQ.RZERO ) THEN
                     IF( K.GE.KTOP+1 ) &
      &                  TST1 = TST1 + CABS1( H( K, K-1 ) )
                     IF( K.GE.KTOP+2 ) &
      &                  TST1 = TST1 + CABS1( H( K, K-2 ) )
                     IF( K.GE.KTOP+3 ) &
      &                  TST1 = TST1 + CABS1( H( K, K-3 ) )
                     IF( K.LE.KBOT-2 ) &
      &                  TST1 = TST1 + CABS1( H( K+2, K+1 ) )
                     IF( K.LE.KBOT-3 ) &
      &                  TST1 = TST1 + CABS1( H( K+3, K+1 ) )
                     IF( K.LE.KBOT-4 ) &
      &                  TST1 = TST1 + CABS1( H( K+4, K+1 ) )
                  END IF
                  IF( CABS1( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) ) &
      &                 THEN
                     H12 = MAX( CABS1( H( K+1, K ) ), &
      &                     CABS1( H( K, K+1 ) ) )
                     H21 = MIN( CABS1( H( K+1, K ) ), &
      &                     CABS1( H( K, K+1 ) ) )
                     H11 = MAX( CABS1( H( K+1, K+1 ) ), &
      &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                     H22 = MIN( CABS1( H( K+1, K+1 ) ), &
      &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )
!
                     IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE. &
      &                   MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  END IF
               END IF
   80       CONTINUE
!
!           ==== Multiply H by reflections from the left ====
!
            IF( ACCUM ) THEN
               JBOT = MIN( NDCOL, KBOT )
            ELSE IF( WANTT ) THEN
               JBOT = N
            ELSE
               JBOT = KBOT
            END IF
!
            DO 100 M = MBOT, MTOP, -1
               K = KRCOL + 2*( M-1 )
               DO 90 J = MAX( KTOP, KRCOL + 2*M ), JBOT
                  REFSUM = DCONJG( V( 1, M ) )* &
      &                     ( H( K+1, J )+DCONJG( V( 2, M ) )* &
      &                     H( K+2, J )+DCONJG( V( 3, M ) )*H( K+3, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
   90          CONTINUE
  100       CONTINUE
!
!           ==== Accumulate orthogonal transformations. ====
!
            IF( ACCUM ) THEN
!
!              ==== Accumulate U. (If needed, update Z later
!              .    with an efficient matrix-matrix
!              .    multiply.) ====
!
               DO 120 M = MBOT, MTOP, -1
                  K = KRCOL + 2*( M-1 )
                  KMS = K - INCOL
                  I2 = MAX( 1, KTOP-INCOL )
                  I2 = MAX( I2, KMS-(KRCOL-INCOL)+1 )
                  I4 = MIN( KDU, KRCOL + 2*( MBOT-1 ) - INCOL + 5 )
                  DO 110 J = I2, I4
                     REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )* &
      &                        U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - &
      &                               REFSUM*DCONJG( V( 2, M ) )
                     U( J, KMS+3 ) = U( J, KMS+3 ) - &
      &                               REFSUM*DCONJG( V( 3, M ) )
  110             CONTINUE
  120          CONTINUE
            ELSE IF( WANTZ ) THEN
!
!              ==== U is not accumulated, so update Z
!              .    now by multiplying by reflections
!              .    from the right. ====
!
               DO 140 M = MBOT, MTOP, -1
                  K = KRCOL + 2*( M-1 )
                  DO 130 J = ILOZ, IHIZ
                     REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )* &
      &                        Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - &
      &                             REFSUM*DCONJG( V( 2, M ) )
                     Z( J, K+3 ) = Z( J, K+3 ) - &
      &                             REFSUM*DCONJG( V( 3, M ) )
  130             CONTINUE
  140          CONTINUE
            END IF
!
!           ==== End of near-the-diagonal bulge chase. ====
!
  145    CONTINUE
!
!        ==== Use U (if accumulated) to update far-from-diagonal
!        .    entries in H.  If required, use U to update Z as
!        .    well. ====
!
         IF( ACCUM ) THEN
            IF( WANTT ) THEN
               JTOP = 1
               JBOT = N
            ELSE
               JTOP = KTOP
               JBOT = KBOT
            END IF
            K1 = MAX( 1, KTOP-INCOL )
            NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
!
!           ==== Horizontal Multiply ====
!
            DO 150 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
               JLEN = MIN( NH, JBOT-JCOL+1 )
               CALL ZGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ), &
      &                     LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH, &
      &                     LDWH )
               CALL ZLACPY( 'ALL', NU, JLEN, WH, LDWH, &
      &                      H( INCOL+K1, JCOL ), LDH )
  150       CONTINUE
!
!           ==== Vertical multiply ====
!
            DO 160 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
               JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
               CALL ZGEMM( 'N', 'N', JLEN, NU, NU, ONE, &
      &                     H( JROW, INCOL+K1 ), LDH, U( K1, K1 ), &
      &                     LDU, ZERO, WV, LDWV )
               CALL ZLACPY( 'ALL', JLEN, NU, WV, LDWV, &
      &                      H( JROW, INCOL+K1 ), LDH )
  160       CONTINUE
!
!           ==== Z multiply (also vertical) ====
!
            IF( WANTZ ) THEN
               DO 170 JROW = ILOZ, IHIZ, NV
                  JLEN = MIN( NV, IHIZ-JROW+1 )
                  CALL ZGEMM( 'N', 'N', JLEN, NU, NU, ONE, &
      &                        Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ), &
      &                        LDU, ZERO, WV, LDWV )
                  CALL ZLACPY( 'ALL', JLEN, NU, WV, LDWV, &
      &                         Z( JROW, INCOL+K1 ), LDZ )
  170          CONTINUE
            END IF
         END IF
  180 CONTINUE
!
!     ==== End of ZLAQR5 ====
!
      END SUBROUTINE ZLAQR5
!> \brief \b ZLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFB applies a complex block reflector H or its transpose H**H to a
!> complex M-by-N matrix C, from either the left or the right.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**H from the Left
!>          = 'R': apply H or H**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'C': apply H**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The triangular K-by-K matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
      &                   T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), &
      &                   WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZLACGV, ZTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
      &   RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C1**H
!
               DO 10 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
      &                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2**H * V2
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, &
      &                        K, M-K, ONE, C( K+1, 1 ), LDC, &
      &                        V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
      &                        M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK, &
      &                        LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
      &                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
      &                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K, &
      &                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
      &                        ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
      &                        N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ), &
      &                        LDV, ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
      &                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C2**H
!
               DO 70 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
      &                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**H * V1
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, &
      &                        K, M-K, ONE, C, LDC, V, LDV, ONE, WORK, &
      &                        LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
      &                        M-K, N, K, -ONE, V, LDV, WORK, LDWORK, &
      &                        ONE, C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
      &                     'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK, &
      &                     LDWORK )
!
!              C2 := C2 - W**H
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - &
      &                               DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
      &                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K, &
      &                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
      &                        N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE, &
      &                        C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
      &                     'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, &
      &                     LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C1**H
!
               DO 130 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
      &                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2**H * V2**H
!
                  CALL ZGEMM( 'Conjugate transpose', &
      &                        'Conjugate transpose', N, K, M-K, ONE, &
      &                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
      &                        WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2**H * W**H
!
                  CALL ZGEMM( 'Conjugate transpose', &
      &                        'Conjugate transpose', M-K, N, K, -ONE, &
      &                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
      &                        C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
      &                     K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', &
      &                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
      &                        K, N-K, ONE, C( 1, K+1 ), LDC, &
      &                        V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K, &
      &                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
      &                        C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
      &                     K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C2**H
!
               DO 190 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
      &                     'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK, &
      &                     LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**H * V1**H
!
                  CALL ZGEMM( 'Conjugate transpose', &
      &                        'Conjugate transpose', N, K, M-K, ONE, C, &
      &                        LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1**H * W**H
!
                  CALL ZGEMM( 'Conjugate transpose', &
      &                        'Conjugate transpose', M-K, N, K, -ONE, V, &
      &                        LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
      &                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W**H
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - &
      &                               DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', &
      &                     'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK, &
      &                     LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
      &                        K, N-K, ONE, C, LDC, V, LDV, ONE, WORK, &
      &                        LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
      &                     ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K, &
      &                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
      &                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of ZLARFB
!
      END SUBROUTINE ZLARFB
!> \brief \b ZLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       COMPLEX*16         TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARF applies a complex elementary reflector H to a complex M-by-N
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!>
!> To apply H**H, supply conjg(tau) instead
!> tau.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
      &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
!     ..
!     .. Executable Statements ..
!
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
!
!        Form  H * C
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
!
            CALL ZGEMV( 'Conjugate transpose', LASTV, LASTC, ONE, &
      &           C, LDC, V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
!
            CALL ZGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            CALL ZGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
      &           V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
!
            CALL ZGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of ZLARF
!
      END SUBROUTINE ZLARF
!> \brief \b ZLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       COMPLEX*16         ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFG generates a complex elementary reflector H of order n, such
!> that
!>
!>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
!>              (   x   )   (   0  )
!>
!> where alpha and beta are scalars, with beta real, and x is an
!> (n-1)-element complex vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
!>                     ( v )
!>
!> where tau is a complex scalar and v is a complex (n-1)-element
!> vector. Note that H is not hermitian.
!>
!> If the elements of x are all zero and alpha is real, then tau = 0
!> and H is taken to be the unit matrix.
!>
!> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
!
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
!
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
      &         GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
         CALL ZSCAL( N-1, ALPHA, X, INCX )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      RETURN
!
!     End of ZLARFG
!
      END SUBROUTINE ZLARFG
!> \brief \b ZLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarft.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarft.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarft.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFT forms the triangular factor T of a complex block reflector H
!> of order n, which is defined as a product of k elementary reflectors.
!>
!> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!>
!> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!>
!> If STOREV = 'C', the vector which defines the elementary reflector
!> H(i) is stored in the i-th column of the array V, and
!>
!>    H  =  I - V * T * V**H
!>
!> If STOREV = 'R', the vector which defines the elementary reflector
!> H(i) is stored in the i-th row of the array V, and
!>
!>    H  =  I - V**H * T * V
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies the order in which the elementary reflectors are
!>          multiplied to form the block reflector:
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Specifies how the vectors which define the elementary
!>          reflectors are stored (see also Further Details):
!>          = 'C': columnwise
!>          = 'R': rowwise
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the block reflector H. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the triangular factor T (= the number of
!>          elementary reflectors). K >= 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                               (LDV,K) if STOREV = 'C'
!>                               (LDV,N) if STOREV = 'R'
!>          The matrix V. See further details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The k by k triangular factor T of the block reflector.
!>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!>          lower triangular. The rest of the array is not used.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
      &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZTRMV, ZGEMM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
      &   RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( PREVLASTV, I )
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
!
!              general case
!
               IF( LSAME( STOREV, 'C' ) ) THEN
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * CONJG( V( I , J ) )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
!
                  CALL ZGEMV( 'Conjugate transpose', J-I, I-1, &
      &                        -TAU( I ), V( I+1, 1 ), LDV, &
      &                        V( I+1, I ), 1, ONE, T( 1, I ), 1 )
               ELSE
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
!
                  CALL ZGEMM( 'N', 'C', I-1, 1, J-I, -TAU( I ), &
      &                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, &
      &                        ONE, T( 1, I ), LDT )
               END IF
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
      &                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
             END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * CONJG( V( N-K+I , J ) )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
!
                     CALL ZGEMV( 'Conjugate transpose', N-K+I-J, K-I, &
      &                           -TAU( I ), V( J, I+1 ), LDV, V( J, I ), &
      &                           1, ONE, T( I+1, I ), 1 )
                  ELSE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
!
                     CALL ZGEMM( 'N', 'C', K-I, 1, N-K+I-J, -TAU( I ), &
      &                           V( I+1, J ), LDV, V( I, J ), LDV, &
      &                           ONE, T( I+1, I ), LDT )
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
      &                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
         END DO
      END IF
      RETURN
!
!     End of ZLARFT
!
      END SUBROUTINE ZLARFT
!> \brief \b ZLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlascl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASCL multiplies the M by N complex matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See ZGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
      &         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
      &            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
      &             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
      &            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
      &            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) &
      &   RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
      &   GO TO 10
!
      RETURN
!
!     End of ZLASCL
!
      END SUBROUTINE ZLASCL
!> \brief \b ZLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASET + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaset.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaset.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaset.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       COMPLEX*16         ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASET initializes a 2-D array A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set. The lower triangle
!>                      is unchanged.
!>          = 'L':      Lower triangular part is set. The upper triangle
!>                      is unchanged.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          On entry, M specifies the number of rows of A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          On entry, N specifies the number of columns of A.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>          All the offdiagonal array elements are set to ALPHA.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!>          All the diagonal array elements are set to BETA.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
!>                   A(i,i) = BETA , 1 <= i <= min(m,n)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the diagonal to BETA and the strictly upper triangular
!        part of the array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the diagonal to BETA and the strictly lower triangular
!        part of the array to ALPHA.
!
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
!
      ELSE
!
!        Set the array to BETA on the diagonal and ALPHA on the
!        offdiagonal.
!
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASET
!
      END SUBROUTINE ZLASET
!> \brief \b ZLATRS solves a triangular system of equations with the scale factor set to prevent overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLATRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
!                          CNORM, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORMIN, TRANS, UPLO
!       INTEGER            INFO, LDA, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   CNORM( * )
!       COMPLEX*16         A( LDA, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATRS solves one of the triangular systems
!>
!>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
!>
!> with scaling to prevent overflow.  Here A is an upper or lower
!> triangular matrix, A**T denotes the transpose of A, A**H denotes the
!> conjugate transpose of A, x and b are n-element vectors, and s is a
!> scaling factor, usually less than or equal to 1, chosen so that the
!> components of x will be less than the overflow threshold.  If the
!> unscaled problem will not cause overflow, the Level 2 BLAS routine
!> ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
!> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  Solve A * x = s*b     (No transpose)
!>          = 'T':  Solve A**T * x = s*b  (Transpose)
!>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] NORMIN
!> \verbatim
!>          NORMIN is CHARACTER*1
!>          Specifies whether CNORM has been set or not.
!>          = 'Y':  CNORM contains the column norms on entry
!>          = 'N':  CNORM is not set on entry.  On exit, the norms will
!>                  be computed and stored in CNORM.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max (1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>          On entry, the right hand side b of the triangular system.
!>          On exit, X is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scaling factor s for the triangular system
!>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
!>          If SCALE = 0, the matrix A is singular or badly scaled, and
!>          the vector x is an exact or approximate solution to A*x = 0.
!> \endverbatim
!>
!> \param[in,out] CNORM
!> \verbatim
!>          CNORM is DOUBLE PRECISION array, dimension (N)
!>
!>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!>          contains the norm of the off-diagonal part of the j-th column
!>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!>          must be greater than or equal to the 1-norm.
!>
!>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!>          returns the 1-norm of the offdiagonal part of the j-th column
!>          of A.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -k, the k-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  A rough bound on x is computed; if that is less than overflow, ZTRSV
!>  is called, otherwise, specific code is used which checks for possible
!>  overflow or divide-by-zero at every operation.
!>
!>  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!>  if A is lower triangular is
!>
!>       x[1:n] := b[1:n]
!>       for j = 1, ..., n
!>            x(j) := x(j) / A(j,j)
!>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!>       end
!>
!>  Define bounds on the components of x after j iterations of the loop:
!>     M(j) = bound on x[1:j]
!>     G(j) = bound on x[j+1:n]
!>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!>
!>  Then for iteration j+1 we have
!>     M(j+1) <= G(j) / | A(j+1,j+1) |
!>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!>
!>  where CNORM(j+1) is greater than or equal to the infinity-norm of
!>  column j+1 of A, not counting the diagonal.  Hence
!>
!>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!>                  1<=i<=j
!>  and
!>
!>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!>                                   1<=i< j
!>
!>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTRSV if the
!>  reciprocal of the largest M(j), j=1,..,n, is larger than
!>  max(underflow, 1/overflow).
!>
!>  The bound on x(j) is also used to determine when a step in the
!>  columnwise method can be performed without fear of overflow.  If
!>  the computed bound is greater than a large constant, x is scaled to
!>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!>
!>  Similarly, a row-wise scheme is used to solve A**T *x = b  or
!>  A**H *x = b.  The basic algorithm for A upper triangular is
!>
!>       for j = 1, ..., n
!>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!>       end
!>
!>  We simultaneously compute two bounds
!>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!>       M(j) = bound on x(i), 1<=i<=j
!>
!>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!>  Then the bound on x(j) is
!>
!>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!>
!>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!>                      1<=i<=j
!>
!>  and we can safely call ZTRSV if 1/M(n) and 1/G(n) are both greater
!>  than max(underflow, 1/overflow).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, &
      &                   CNORM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         A( LDA, * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, &
      &                   TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, &
      &                   XBND, XJ, XMAX
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           LSAME, IDAMAX, IZAMAX, DLAMCH, DZASUM, ZDOTC, &
      &                   ZDOTU, ZLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, XERBLA, ZAXPY, ZDSCAL, ZTRSV, DLABAD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1, CABS2
!     ..
!     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( DBLE( ZDUM ) / 2.D0 ) + &
      &                ABS( DIMAG( ZDUM ) / 2.D0 )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
      &         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
      &         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLATRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
      &   RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
!
      IF( LSAME( NORMIN, 'N' ) ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         IF( UPPER ) THEN
!
!           A is upper triangular.
!
            DO 10 J = 1, N
               CNORM( J ) = DZASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
!
!           A is lower triangular.
!
            DO 20 J = 1, N - 1
               CNORM( J ) = DZASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM/2.
!
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM*HALF ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = HALF / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine ZTRSV can be used.
!
      XMAX = ZERO
      DO 30 J = 1, N
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
   30 CONTINUE
      XBND = XMAX
!
      IF( NOTRAN ) THEN
!
!        Compute the growth in A * x = b.
!
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 60
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 40 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
      &            GO TO 60
!
               TJJS = A( J, J )
               TJJ = CABS1( TJJS )
!
               IF( TJJ.GE.SMLNUM ) THEN
!
!                 M(j) = G(j-1) / abs(A(j,j))
!
                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  XBND = ZERO
               END IF
!
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  GROW = ZERO
               END IF
   40       CONTINUE
            GROW = XBND
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
      &            GO TO 60
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
!
      ELSE
!
!        Compute the growth in A**T * x = b  or  A**H * x = b.
!
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 90
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 70 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
      &            GO TO 90
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
!
               TJJS = A( J, J )
               TJJ = CABS1( TJJS )
!
               IF( TJJ.GE.SMLNUM ) THEN
!
!                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
                  IF( XJ.GT.TJJ ) &
      &               XBND = XBND*( TJJ / XJ )
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  XBND = ZERO
               END IF
   70       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 80 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
      &            GO TO 90
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
!
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL ZTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         IF( XMAX.GT.BIGNUM*HALF ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            SCALE = ( BIGNUM*HALF ) / XMAX
            CALL ZDSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         ELSE
            XMAX = XMAX*TWO
         END IF
!
         IF( NOTRAN ) THEN
!
!           Solve A * x = b
!
            DO 120 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) &
      &               GO TO 110
               END IF
               TJJ = CABS1( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by 1/b(j).
!
                        REC = ONE / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     END IF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO 100 I = 1, N
                     X( I ) = ZERO
  100             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  110          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     REC = REC*HALF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL ZDSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
!
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
!
!                    Compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
!
                     CALL ZAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X, &
      &                           1 )
                     I = IZAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
!
!                    Compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
!
                     CALL ZAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1, &
      &                           X( J+1 ), 1 )
                     I = J + IZAMAX( N-J, X( J+1 ), 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
               END IF
  120       CONTINUE
!
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!           Solve A**T * x = b
!
            DO 170 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTU to perform the dot product.
!
                  IF( UPPER ) THEN
                     CSUMJ = ZDOTU( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTU( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     DO 130 I = 1, J - 1
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 140 I = J + 1, N
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  140                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
      &                  GO TO 160
                  END IF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**T *x = 0.
!
                     DO 150 I = 1, N
                        X( I ) = ZERO
  150                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  160             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  170       CONTINUE
!
         ELSE
!
!           Solve A**H * x = b
!
            DO 220 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTC to perform the dot product.
!
                  IF( UPPER ) THEN
                     CSUMJ = ZDOTC( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTC( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     DO 180 I = 1, J - 1
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* &
      &                          X( I )
  180                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 190 I = J + 1, N
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* &
      &                          X( I )
  190                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
      &                  GO TO 210
                  END IF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**H *x = 0.
!
                     DO 200 I = 1, N
                        X( I ) = ZERO
  200                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  210             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  220       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
!
      RETURN
!
!     End of ZLATRS
!
      END SUBROUTINE ZLATRS
!> \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZROT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       DOUBLE PRECISION   C
!       COMPLEX*16         S
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         CX( * ), CY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZROT   applies a plane rotation, where the cos (C) is real and the
!> sin (S) is complex, and the vectors CX and CY are complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vectors CX and CY.
!> \endverbatim
!>
!> \param[in,out] CX
!> \verbatim
!>          CX is COMPLEX*16 array, dimension (N)
!>          On input, the vector X.
!>          On output, CX is overwritten with C*X + S*Y.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of CX.  INCX <> 0.
!> \endverbatim
!>
!> \param[in,out] CY
!> \verbatim
!>          CY is COMPLEX*16 array, dimension (N)
!>          On input, the vector Y.
!>          On output, CY is overwritten with -CONJG(S)*X + C*Y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between successive values of CY.  INCX <> 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX*16
!>          C and S define a rotation
!>             [  C          S  ]
!>             [ -conjg(S)   C  ]
!>          where C*C + S*CONJG(S) = 1.0.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C
      COMPLEX*16         S
!     ..
!     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) &
      &   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) &
      &   GO TO 20
!
!     Code for unequal increments or equal increments not equal to 1
!
      IX = 1
      IY = 1
      IF( INCX.LT.0 ) &
      &   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 ) &
      &   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         STEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - DCONJG( S )*CX( IX )
         CX( IX ) = STEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
!
!     Code for both increments equal to 1
!
   20 CONTINUE
      DO 30 I = 1, N
         STEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - DCONJG( S )*CX( I )
         CX( I ) = STEMP
   30 CONTINUE
      RETURN
      END SUBROUTINE ZROT

!> \brief \b ZTREVC3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTREVC3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!      $                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTREVC3 computes some or all of the right and/or left eigenvectors of
!> a complex upper triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>              T*x = w*x,     (y**H)*T = w*(y**H)
!>
!> where y**H denotes the conjugate transpose of the vector y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix. If Q is the unitary factor that reduces a matrix A to
!> Schur form T, then Q*X and Q*Y are the matrices of right and left
!> eigenvectors of A.
!>
!> This uses a Level 3 BLAS version of the back transformation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  compute right eigenvectors only;
!>          = 'L':  compute left eigenvectors only;
!>          = 'B':  compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A':  compute all right and/or left eigenvectors;
!>          = 'B':  compute all right and/or left eigenvectors,
!>                  backtransformed using the matrices supplied in
!>                  VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          The eigenvector corresponding to the j-th eigenvalue is
!>          computed if SELECT(j) = .TRUE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          The upper triangular matrix T.  T is modified, but restored
!>          on exit.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is COMPLEX*16 array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by ZHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.
!>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX*16 array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by ZHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.
!>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of columns in the arrays VL and/or VR. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns in the arrays VL and/or VR actually
!>          used to store the eigenvectors.
!>          If HOWMNY = 'A' or 'B', M is set to N.
!>          Each selected eigenvector occupies one column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of array WORK. LWORK >= max(1,2*N).
!>          For optimum performance, LWORK >= N + 2*N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (LRWORK)
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of array RWORK. LRWORK >= max(1,N).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the RWORK array, returns
!>          this value as the first entry of the RWORK array, and no error
!>          message related to LRWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The algorithm used in this program is basically backward (forward)
!>  substitution, with scaling to make the the code robust against
!>  possible overflow.
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x| + |y|.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
      &                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
      IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
      &                   WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
      &                     CONE  = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            NBMIN, NBMAX
      PARAMETER          ( NBMIN = 8, NBMAX = 128 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI, IV, MAXWRK, NB
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           LSAME, ILAENV, IZAMAX, DLAMCH, DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLATRS, &
      &                   ZGEMM, DLABAD, ZLASET, ZLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, CONJG, DIMAG, MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      BOTHV  = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV  = LSAME( SIDE, 'L' ) .OR. BOTHV
!
      ALLV  = LSAME( HOWMNY, 'A' )
      OVER  = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )
!
!     Set M to the number of columns required to store the selected
!     eigenvectors.
!
      IF( SOMEV ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) ) &
      &         M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
!
      INFO = 0
      NB = ILAENV( 1, 'ZTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
      MAXWRK = N + 2*N*NB
      WORK(1) = MAXWRK
      RWORK(1) = N
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) THEN
         INFO = -14
      ELSE IF ( LRWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREVC3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
      &   RETURN
!
!     Use blocked version of back-transformation if sufficient workspace.
!     Zero-out the workspace to avoid potential NaN propagation.
!
      IF( OVER .AND. LWORK .GE. N + 2*N*NBMIN ) THEN
         NB = (LWORK - N) / (2*N)
         NB = MIN( NB, NBMAX )
         CALL ZLASET( 'F', N, 1+2*NB, CZERO, CZERO, WORK, N )
      ELSE
         NB = 1
      END IF
!
!     Set the constants to control overflow.
!
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
!
!     Store the diagonal elements of T in working array WORK.
!
      DO 20 I = 1, N
         WORK( I ) = T( I, I )
   20 CONTINUE
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      RWORK( 1 ) = ZERO
      DO 30 J = 2, N
         RWORK( J ) = DZASUM( J-1, T( 1, J ), 1 )
   30 CONTINUE
!
      IF( RIGHTV ) THEN
!
!        ============================================================
!        Compute right eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=NB=1;
!        blocked     version starts with IV=NB, goes down to 1.
!        (Note the "0-th" column is used to store the original diagonal.)
         IV = NB
         IS = M
         DO 80 KI = N, 1, -1
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
      &            GO TO 80
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
!           --------------------------------------------------------
!           Complex right eigenvector
!
            WORK( KI + IV*N ) = CONE
!
!           Form right-hand side.
!
            DO 40 K = 1, KI - 1
               WORK( K + IV*N ) = -T( K, KI )
   40       CONTINUE
!
!           Solve upper triangular system:
!           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
!
            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) &
      &            T( K, K ) = SMIN
   50       CONTINUE
!
            IF( KI.GT.1 ) THEN
               CALL ZLATRS( 'Upper', 'No transpose', 'Non-unit', 'Y', &
      &                      KI-1, T, LDT, WORK( 1 + IV*N ), SCALE, &
      &                      RWORK, INFO )
               WORK( KI + IV*N ) = SCALE
            END IF
!
!           Copy the vector x or Q*x to VR and normalize.
!
            IF( .NOT.OVER ) THEN
!              ------------------------------
!              no back-transform: copy x to VR and normalize.
               CALL ZCOPY( KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 )
!
               II = IZAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               CALL ZDSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
               DO 60 K = KI + 1, N
                  VR( K, IS ) = CZERO
   60          CONTINUE
!
            ELSE IF( NB.EQ.1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.GT.1 ) &
      &            CALL ZGEMV( 'N', N, KI-1, CONE, VR, LDVR, &
      &                        WORK( 1 + IV*N ), 1, DCMPLX( SCALE ), &
      &                        VR( 1, KI ), 1 )
!
               II = IZAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               CALL ZDSCAL( N, REMAX, VR( 1, KI ), 1 )
!
            ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out below vector
               DO K = KI + 1, N
                  WORK( K + IV*N ) = CZERO
               END DO
!
!              Columns IV:NB of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
               IF( (IV.EQ.1) .OR. (KI.EQ.1) ) THEN
                  CALL ZGEMM( 'N', 'N', N, NB-IV+1, KI+NB-IV, CONE, &
      &                        VR, LDVR, &
      &                        WORK( 1 + (IV)*N    ), N, &
      &                        CZERO, &
      &                        WORK( 1 + (NB+IV)*N ), N )
!                 normalize vectors
                  DO K = IV, NB
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     CALL ZDSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                  END DO
                  CALL ZLACPY( 'F', N, NB-IV+1, &
      &                         WORK( 1 + (NB+IV)*N ), N, &
      &                         VR( 1, KI ), LDVR )
                  IV = NB
               ELSE
                  IV = IV - 1
               END IF
            END IF
!
!           Restore the original diagonal elements of T.
!
            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K )
   70       CONTINUE
!
            IS = IS - 1
   80    CONTINUE
      END IF
!
      IF( LEFTV ) THEN
!
!        ============================================================
!        Compute left eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=1;
!        blocked     version starts with IV=1, goes up to NB.
!        (Note the "0-th" column is used to store the original diagonal.)
         IV = 1
         IS = 1
         DO 130 KI = 1, N
!
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
      &            GO TO 130
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
!           --------------------------------------------------------
!           Complex left eigenvector
!
            WORK( KI + IV*N ) = CONE
!
!           Form right-hand side.
!
            DO 90 K = KI + 1, N
               WORK( K + IV*N ) = -CONJG( T( KI, K ) )
   90       CONTINUE
!
!           Solve conjugate-transposed triangular system:
!           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
!
            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) &
      &            T( K, K ) = SMIN
  100       CONTINUE
!
            IF( KI.LT.N ) THEN
               CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', &
      &                      'Y', N-KI, T( KI+1, KI+1 ), LDT, &
      &                      WORK( KI+1 + IV*N ), SCALE, RWORK, INFO )
               WORK( KI + IV*N ) = SCALE
            END IF
!
!           Copy the vector x or Q*x to VL and normalize.
!
            IF( .NOT.OVER ) THEN
!              ------------------------------
!              no back-transform: copy x to VL and normalize.
               CALL ZCOPY( N-KI+1, WORK( KI + IV*N ), 1, VL(KI,IS), 1 )
!
               II = IZAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               CALL ZDSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
!
               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CZERO
  110          CONTINUE
!
            ELSE IF( NB.EQ.1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.LT.N ) &
      &            CALL ZGEMV( 'N', N, N-KI, CONE, VL( 1, KI+1 ), LDVL, &
      &                        WORK( KI+1 + IV*N ), 1, DCMPLX( SCALE ), &
      &                        VL( 1, KI ), 1 )
!
               II = IZAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               CALL ZDSCAL( N, REMAX, VL( 1, KI ), 1 )
!
            ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out above vector
!              could go from KI-NV+1 to KI-1
               DO K = 1, KI - 1
                  WORK( K + IV*N ) = CZERO
               END DO
!
!              Columns 1:IV of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
               IF( (IV.EQ.NB) .OR. (KI.EQ.N) ) THEN
                  CALL ZGEMM( 'N', 'N', N, IV, N-KI+IV, CONE, &
      &                        VL( 1, KI-IV+1 ), LDVL, &
      &                        WORK( KI-IV+1 + (1)*N ), N, &
      &                        CZERO, &
      &                        WORK( 1 + (NB+1)*N ), N )
!                 normalize vectors
                  DO K = 1, IV
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     CALL ZDSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                  END DO
                  CALL ZLACPY( 'F', N, IV, &
      &                         WORK( 1 + (NB+1)*N ), N, &
      &                         VL( 1, KI-IV+1 ), LDVL )
                  IV = 1
               ELSE
                  IV = IV + 1
               END IF
            END IF
!
!           Restore the original diagonal elements of T.
!
            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K )
  120       CONTINUE
!
            IS = IS + 1
  130    CONTINUE
      END IF
!
      RETURN
!
!     End of ZTREVC3
!
      END SUBROUTINE ZTREVC3
!> \brief \b ZTREXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTREXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ
!       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         Q( LDQ, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTREXC reorders the Schur factorization of a complex matrix
!> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
!> is moved to row ILST.
!>
!> The Schur form T is reordered by a unitary similarity transformation
!> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
!> postmultplying it with Z.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V':  update the matrix Q of Schur vectors;
!>          = 'N':  do not update Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!>          If N == 0 arguments ILST and IFST may be any value.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          On entry, the upper triangular matrix T.
!>          On exit, the reordered upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          unitary transformation matrix Z which reorders T.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1, and if
!>          COMPQ = 'V', LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] IFST
!> \verbatim
!>          IFST is INTEGER
!> \endverbatim
!>
!> \param[in] ILST
!> \verbatim
!>          ILST is INTEGER
!>
!>          Specify the reordering of the diagonal elements of T:
!>          The element with row index IFST is moved to row ILST by a
!>          sequence of transpositions between adjacent elements.
!>          1 <= IFST <= N; 1 <= ILST <= N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         Q( LDQ, * ), T( LDT, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            K, M1, M2, M3
      DOUBLE PRECISION   CS
      COMPLEX*16         SN, T11, T22, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF(( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -7
      ELSE IF(( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREXC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.1 .OR. IFST.EQ.ILST ) &
      &   RETURN
!
      IF( IFST.LT.ILST ) THEN
!
!        Move the IFST-th diagonal element forward down the diagonal.
!
         M1 = 0
         M2 = -1
         M3 = 1
      ELSE
!
!        Move the IFST-th diagonal element backward up the diagonal.
!
         M1 = -1
         M2 = 0
         M3 = -1
      END IF
!
      DO 10 K = IFST + M1, ILST + M2, M3
!
!        Interchange the k-th and (k+1)-th diagonal elements.
!
         T11 = T( K, K )
         T22 = T( K+1, K+1 )
!
!        Determine the transformation to perform the interchange.
!
         CALL ZLARTG( T( K, K+1 ), T22-T11, CS, SN, TEMP )
!
!        Apply transformation to the matrix T.
!
         IF( K+2.LE.N ) &
      &      CALL ZROT( N-K-1, T( K, K+2 ), LDT, T( K+1, K+2 ), LDT, CS, &
      &                 SN )
         CALL ZROT( K-1, T( 1, K ), 1, T( 1, K+1 ), 1, CS, &
      &              DCONJG( SN ) )
!
         T( K, K ) = T22
         T( K+1, K+1 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL ZROT( N, Q( 1, K ), 1, Q( 1, K+1 ), 1, CS, &
      &                 DCONJG( SN ) )
         END IF
!
   10 CONTINUE
!
      RETURN
!
!     End of ZTREXC
!
      END SUBROUTINE ZTREXC
!> \brief \b ZTRMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRMM  performs one of the matrix-matrix operations
!>
!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
!>
!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**H.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, N ).
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B,  and  on exit  is overwritten  by the
!>           transformed matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!     ..
!     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
      &         (.NOT.LSAME(TRANSA,'T')) .AND. &
      &         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*A*B.
!
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
!
              IF (UPPER) THEN
                  DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*A.
!
              IF (UPPER) THEN
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
!
              IF (UPPER) THEN
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of ZTRMM
!
      END SUBROUTINE ZTRMM
!> \brief \b ZTRMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRMV  performs one of the matrix-vector operations
!>
!>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
!>
!> where x is an n element vector and  A is an n by n unit, or non-unit,
!> upper or lower triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   x := A*x.
!>
!>              TRANS = 'T' or 't'   x := A**T*x.
!>
!>              TRANS = 'C' or 'c'   x := A**H*x.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N ).
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x. On exit, X is overwritten with the
!>           transformed vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
      &         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := A*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := A**T*x  or  x := A**H*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of ZTRMV
!
      END SUBROUTINE ZTRMV
!> \brief \b ZTRSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSV  solves one of the systems of equations
!>
!>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
!>
!> where b and x are n element vectors and A is an n by n unit, or
!> non-unit, upper or lower triangular matrix.
!>
!> No test for singularity or near-singularity is included in this
!> routine. Such tests must be performed before calling this routine.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the equations to be solved as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   A*x = b.
!>
!>              TRANS = 'T' or 't'   A**T*x = b.
!>
!>              TRANS = 'C' or 'c'   A**H*x = b.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element right-hand side vector b. On exit, X is overwritten
!>           with the solution vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
      &         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRSV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 120 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 130 I = 1,J - 1
                              TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 180 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 190 I = N,J + 1,-1
                              TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  200             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of ZTRSV
!
      END SUBROUTINE ZTRSV
!> \brief \b ZUNG2R
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNG2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the m by n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
      &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
      &   RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
      &                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
      &      CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of ZUNG2R
!
      END SUBROUTINE ZUNG2R
!> \brief \b ZUNGHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNGHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunghr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunghr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunghr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGHR generates a complex unitary matrix Q which is defined as the
!> product of IHI-ILO elementary reflectors of order N, as returned by
!> ZGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI must have the same values as in the previous call
!>          of ZGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by ZGEHRD.
!>          On exit, the N-by-N unitary matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEHRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= IHI-ILO.
!>          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
      &                   ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQR
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NH = IHI - ILO
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'ZUNGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE
!
      IF( NH.GT.0 ) THEN
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
         CALL ZUNGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), &
      &                WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNGHR
!
      END SUBROUTINE ZUNGHR
!> \brief \b ZUNGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNGQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
      &                   LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'ZUNGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) &
      &   CALL ZUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
      &                TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
      &                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL ZLARFB( 'Left', 'No transpose', 'Forward', &
      &                      'Columnwise', M-I+1, N-I-IB+1, IB, &
      &                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
      &                      LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL ZUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
      &                   IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of ZUNGQR
!
      END SUBROUTINE ZUNGQR
!> \brief \b ZUNM2R multiplies a general matrix by the unitary matrix from a QR factorization determined by cgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNM2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunm2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunm2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunm2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNM2R overwrites the general complex m-by-n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**H if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left
!>          = 'R': apply Q or Q**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      &                   WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
      &   RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) or H(i)**H is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) or H(i)**H is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i) or H(i)**H
!
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC, &
      &               WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of ZUNM2R
!
      END SUBROUTINE ZUNM2R
!> \brief \b ZUNMHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmhr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmhr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmhr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMHR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> IHI-ILO elementary reflectors, as returned by ZGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI must have the same values as in the previous call
!>          of ZGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
!>          ILO = 1 and IHI = 0, if M = 0;
!>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
!>          ILO = 1 and IHI = 0, if N = 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by ZGEHRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEHRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, &
      &                   LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNMQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NH = IHI - ILO
      LEFT = LSAME( SIDE, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) &
      &          THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, NQ ) ) THEN
         INFO = -5
      ELSE IF( IHI.LT.MIN( ILO, NQ ) .OR. IHI.GT.NQ ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( LEFT ) THEN
            NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, NH, N, NH, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, NH, NH, -1 )
         END IF
         LWKOPT = NW*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NH.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( LEFT ) THEN
         MI = NH
         NI = N
         I1 = ILO + 1
         I2 = 1
      ELSE
         MI = M
         NI = NH
         I1 = 1
         I2 = ILO + 1
      END IF
!
      CALL ZUNMQR( SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA, &
      &             TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO )
!
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMHR
!
      END SUBROUTINE ZUNMHR
!> \brief \b ZUNMQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMQR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Conjugate transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQRF in the first k columns of its array argument A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      &                   WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
      &                     TSIZE = LDT*NBMAX )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, &
      &                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N, K, &
      &        -1 ) )
         LWKOPT = NW*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.LWKOPT ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQR', SIDE // TRANS, M, N, K, &
      &              -1 ) )
         END IF
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
      &                IINFO )
      ELSE
!
!        Use blocked code
!
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
      &       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL ZLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
      &                   LDA, TAU( I ), WORK( IWT ), LDT )
            IF( LEFT ) THEN
!
!              H or H**H is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H**H is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H**H
!
            CALL ZLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
      &                   IB, A( I, I ), LDA, WORK( IWT ), LDT, &
      &                   C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMQR
!
      END SUBROUTINE ZUNMQR



!> \brief \b DLAMCHF77 deprecated
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!     CHARACTER          CMACH
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMCHF77 determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          Specifies the value to be returned by DLAMCH:
!>          = 'E' or 'e',   DLAMCH := eps
!>          = 'S' or 's ,   DLAMCH := sfmin
!>          = 'B' or 'b',   DLAMCH := base
!>          = 'P' or 'p',   DLAMCH := eps*base
!>          = 'N' or 'n',   DLAMCH := t
!>          = 'R' or 'r',   DLAMCH := rnd
!>          = 'M' or 'm',   DLAMCH := emin
!>          = 'U' or 'u',   DLAMCH := rmin
!>          = 'L' or 'l',   DLAMCH := emax
!>          = 'O' or 'o',   DLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
 
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION dlamch( CMACH )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          cmach
!     ..
!     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            first, lrnd
      INTEGER            beta, imax, imin, it
      DOUBLE PRECISION   base, emax, emin, eps, prec, rmach, rmax, rmin, &
      &                   rnd, sfmin, small, t
!     ..
!     .. External Functions ..
      LOGICAL            lsame
      EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
      EXTERNAL           dlamc2
!     ..
!     .. Save statement ..
      SAVE               first, eps, sfmin, base, t, rnd, emin, rmin, &
      &                   emax, rmax, prec
!     ..
!     .. Data statements ..
      DATA               first / .true. /
!     ..
!     .. Executable Statements ..
!
      IF( first ) THEN
         CALL dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         IF( lrnd ) THEN
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         ELSE
            rnd = zero
            eps = base**( 1-it )
         END IF
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         IF( small.GE.sfmin ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            sfmin = small*( one+eps )
         END IF
      END IF
!
      IF( lsame( cmach, 'E' ) ) THEN
         rmach = eps
      ELSE IF( lsame( cmach, 'S' ) ) THEN
         rmach = sfmin
      ELSE IF( lsame( cmach, 'B' ) ) THEN
         rmach = base
      ELSE IF( lsame( cmach, 'P' ) ) THEN
         rmach = prec
      ELSE IF( lsame( cmach, 'N' ) ) THEN
         rmach = t
      ELSE IF( lsame( cmach, 'R' ) ) THEN
         rmach = rnd
      ELSE IF( lsame( cmach, 'M' ) ) THEN
         rmach = emin
      ELSE IF( lsame( cmach, 'U' ) ) THEN
         rmach = rmin
      ELSE IF( lsame( cmach, 'L' ) ) THEN
         rmach = emax
      ELSE IF( lsame( cmach, 'O' ) ) THEN
         rmach = rmax
      END IF
!
      dlamch = rmach
      first  = .false.
      RETURN
!
!     End of DLAMCH
!
      END FUNCTION dlamch
!
!***********************************************************************
!
!> \brief \b DLAMC1
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC1 determines the machine parameters given by BETA, T, RND, and
!> IEEE1.
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          The number of ( BETA ) digits in the mantissa.
!> \endverbatim
!>
!> \param[out] RND
!> \verbatim
!>          Specifies whether proper rounding  ( RND = .TRUE. )  or
!>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!>          be a reliable guide to the way in which the machine performs
!>          its arithmetic.
!> \endverbatim
!>
!> \param[out] IEEE1
!> \verbatim
!>          Specifies whether rounding appears to be done in the IEEE
!>          'round to nearest' style.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \ingroup auxOTHERauxiliary
!>
!> \details \b Further \b Details
!> \verbatim
!>
!>  The routine is based on the routine  ENVRON  by Malcolm and
!>  incorporates suggestions by Gentleman and Marovich. See
!>
!>     Malcolm M. A. (1972) Algorithms to reveal properties of
!>        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!>
!>     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!>        that reveal properties of floating point arithmetic units.
!>        Comms. of the ACM, 17, 276-277.
!> \endverbatim
!>
      SUBROUTINE dlamc1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
!     ..
!     .. Save statement ..
      SAVE               first, lieee1, lbeta, lrnd, lt
!     ..
!     .. Data statements ..
      DATA               first / .true. /
!     ..
!     .. Executable Statements ..
!
      IF( first ) THEN
         one = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         a = 1
         c = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( c.EQ.one ) THEN
            a = 2*a
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         b = 1
         c = dlamc3( a, b )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( c.EQ.a ) THEN
            b = 2*b
            c = dlamc3( a, b )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         qtr = one / 4
         savec = c
         c = dlamc3( c, -a )
         lbeta = c + qtr
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         b = lbeta
         f = dlamc3( b / 2, -b / 100 )
         c = dlamc3( f, a )
         IF( c.EQ.a ) THEN
            lrnd = .true.
         ELSE
            lrnd = .false.
         END IF
         f = dlamc3( b / 2, b / 100 )
         c = dlamc3( f, a )
         IF( ( lrnd ) .AND. ( c.EQ.a ) ) &
      &      lrnd = .false.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         t1 = dlamc3( b / 2, a )
         t2 = dlamc3( b / 2, savec )
         lieee1 = ( t1.EQ.a ) .AND. ( t2.GT.savec ) .AND. lrnd
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         lt = 0
         a = 1
         c = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( c.EQ.one ) THEN
            lt = lt + 1
            a = a*lbeta
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      first = .false.
      RETURN
!
!     End of DLAMC1
!
      END SUBROUTINE dlamc1
!
!***********************************************************************
!
!> \brief \b DLAMC2
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC2 determines the machine parameters specified in its argument
!> list.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \ingroup auxOTHERauxiliary
!>
!> \param[out] BETA
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          The number of ( BETA ) digits in the mantissa.
!> \endverbatim
!>
!> \param[out] RND
!> \verbatim
!>          Specifies whether proper rounding  ( RND = .TRUE. )  or
!>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!>          be a reliable guide to the way in which the machine performs
!>          its arithmetic.
!> \endverbatim
!>
!> \param[out] EPS
!> \verbatim
!>          The smallest positive number such that
!>             fl( 1.0 - EPS ) .LT. 1.0,
!>          where fl denotes the computed value.
!> \endverbatim
!>
!> \param[out] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow occurs.
!> \endverbatim
!>
!> \param[out] RMIN
!> \verbatim
!>          The smallest normalized number for the machine, given by
!>          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!>          of BETA.
!> \endverbatim
!>
!> \param[out] EMAX
!> \verbatim
!>          The maximum exponent before overflow occurs.
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          The largest positive number for the machine, given by
!>          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!>          value of BETA.
!> \endverbatim
!>
!> \details \b Further \b Details
!> \verbatim
!>
!>  The computation of  EPS  is based on a routine PARANOIA by
!>  W. Kahan of the University of California at Berkeley.
!> \endverbatim
      SUBROUTINE dlamc2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
      &                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
      &                   SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
!     ..
!     .. External Subroutines ..
      EXTERNAL           dlamc1, dlamc4, dlamc5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
!     ..
!     .. Save statement ..
      SAVE               first, iwarn, lbeta, lemax, lemin, leps, lrmax, &
      &                   lrmin, lt
!     ..
!     .. Data statements ..
      DATA               first / .true. / , iwarn / .false. /
!     ..
!     .. Executable Statements ..
!
      IF( first ) THEN
         zero = 0
         one = 1
         two = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL dlamc1( lbeta, lt, lrnd, lieee1 )
!
!        Start to find EPS.
!
         b = lbeta
         a = b**( -lt )
         leps = a
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         b = two / 3
         half = one / 2
         sixth = dlamc3( b, -half )
         third = dlamc3( sixth, sixth )
         b = dlamc3( third, -half )
         b = dlamc3( b, sixth )
         b = abs( b )
         IF( b.LT.leps ) &
      &      b = leps
!
         leps = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( leps.GT.b ) .AND. ( b.GT.zero ) ) THEN
            leps = b
            c = dlamc3( half*leps, ( two**5 )*( leps**2 ) )
            c = dlamc3( half, -c )
            b = dlamc3( half, c )
            c = dlamc3( half, -b )
            b = dlamc3( half, c )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( a.LT.leps ) &
      &      leps = a
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         rbase = one / lbeta
         small = one
         DO 20 i = 1, 3
            small = dlamc3( small*rbase, zero )
   20    CONTINUE
         a = dlamc3( one, small )
         CALL dlamc4( ngpmin, one, lbeta )
         CALL dlamc4( ngnmin, -one, lbeta )
         CALL dlamc4( gpmin, a, lbeta )
         CALL dlamc4( gnmin, -a, lbeta )
         ieee = .false.
!
         IF( ( ngpmin.EQ.ngnmin ) .AND. ( gpmin.EQ.gnmin ) ) THEN
            IF( ngpmin.EQ.gpmin ) THEN
               lemin = ngpmin
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( gpmin-ngpmin ).EQ.3 ) THEN
               lemin = ngpmin - 1 + lt
               ieee = .true.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               lemin = min( ngpmin, gpmin )
!            ( A guess; no known machine )
               iwarn = .true.
            END IF
!
         ELSE IF( ( ngpmin.EQ.gpmin ) .AND. ( ngnmin.EQ.gnmin ) ) THEN
            IF( abs( ngpmin-ngnmin ).EQ.1 ) THEN
               lemin = max( ngpmin, ngnmin )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               lemin = min( ngpmin, ngnmin )
!            ( A guess; no known machine )
               iwarn = .true.
            END IF
!
         ELSE IF( ( abs( ngpmin-ngnmin ).EQ.1 ) .AND. &
      &            ( gpmin.EQ.gnmin ) ) THEN
            IF( ( gpmin-min( ngpmin, ngnmin ) ).EQ.3 ) THEN
               lemin = max( ngpmin, ngnmin ) - 1 + lt
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               lemin = min( ngpmin, ngnmin )
!            ( A guess; no known machine )
               iwarn = .true.
            END IF
!
         ELSE
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
!         ( A guess; no known machine )
            iwarn = .true.
         END IF
         first = .false.
!**
! Comment out this if block if EMIN is ok
         IF( iwarn ) THEN
            first = .true.
            WRITE( 6, fmt = 9999 )lemin
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         ieee = ieee .OR. lieee1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         lrmin = 1
         DO 30 i = 1, 1 - lemin
            lrmin = dlamc3( lrmin*rbase, zero )
   30    CONTINUE
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax )
      END IF
!
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
      &      '  EMIN = ', i8, / &
      &      ' If, after inspection, the value EMIN looks', &
      &      ' acceptable please comment out ', &
      &      / ' the IF block as marked within the code of routine', &
      &      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END SUBROUTINE dlamc2
!
!***********************************************************************
!
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!>
!> \param[in] A
!>
!> \param[in] B
!> \verbatim
!>          The values A and B.
!> \endverbatim
 
      DOUBLE PRECISION FUNCTION dlamc3( A, B )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   a, b
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      dlamc3 = a + b
!
      RETURN
!
!     End of DLAMC3
!
      END FUNCTION dlamc3
!
!***********************************************************************
!
!> \brief \b DLAMC4
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC4 is a service routine for DLAMC2.
!> \endverbatim
!>
!> \param[out] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow, computed by
!>          setting A = START and dividing by BASE until the previous A
!>          can not be recovered.
!> \endverbatim
!>
!> \param[in] START
!> \verbatim
!>          The starting point for determining EMIN.
!> \endverbatim
!>
!> \param[in] BASE
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
      SUBROUTINE dlamc4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
!     ..
!     .. Executable Statements ..
!
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = dlamc3( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( c1.EQ.a ) .AND. ( c2.EQ.a ) .AND. ( d1.EQ.a ) .AND. &
      &    ( d2.EQ.a ) ) THEN
         emin = emin - 1
         a = b1
         b1 = dlamc3( a / base, zero )
         c1 = dlamc3( b1*base, zero )
         d1 = zero
         DO 20 i = 1, base
            d1 = d1 + b1
   20    CONTINUE
         b2 = dlamc3( a*rbase, zero )
         c2 = dlamc3( b2 / rbase, zero )
         d2 = zero
         DO 30 i = 1, base
            d2 = d2 + b2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END SUBROUTINE dlamc4
!
!***********************************************************************
!
!> \brief \b DLAMC5
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC5 attempts to compute RMAX, the largest machine floating-point
!> number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!> approximately to a power of 2.  It will fail on machines where this
!> assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!> EMAX = 28718).  It will also fail if the value supplied for EMIN is
!> too large (i.e. too close to zero), probably with overflow.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          The base of floating-point arithmetic.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          The number of base BETA digits in the mantissa of a
!>          floating-point value.
!> \endverbatim
!>
!> \param[in] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow.
!> \endverbatim
!>
!> \param[in] IEEE
!> \verbatim
!>          A logical flag specifying whether or not the arithmetic
!>          system is thought to comply with the IEEE standard.
!> \endverbatim
!>
!> \param[out] EMAX
!> \verbatim
!>          The largest exponent before overflow
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          The largest machine floating-point number.
!> \endverbatim
!>
      SUBROUTINE dlamc5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          mod
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      lexp = 1
      exbits = 1
   10 CONTINUE
      try = lexp*2
      IF( try.LE.( -emin ) ) THEN
         lexp = try
         exbits = exbits + 1
         GO TO 10
      END IF
      IF( lexp.EQ.-emin ) THEN
         uexp = lexp
      ELSE
         uexp = try
         exbits = exbits + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( uexp+emin ).GT.( -lexp-emin ) ) THEN
         expsum = 2*lexp
      ELSE
         expsum = 2*uexp
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      emax = expsum + emin - 1
      nbits = 1 + exbits + p
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( mod( nbits, 2 ).EQ.1 ) .AND. ( beta.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         emax = emax - 1
      END IF
!
      IF( ieee ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         emax = emax - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      recbas = one / beta
      z = beta - one
      y = zero
      DO 20 i = 1, p
         z = z*recbas
         IF( y.LT.one ) &
      &      oldy = y
         y = dlamc3( y, z )
   20 CONTINUE
      IF( y.GE.one ) &
      &   y = oldy
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 i = 1, emax
         y = dlamc3( y*beta, zero )
   30 CONTINUE
!
      rmax = y
      RETURN
!
!     End of DLAMC5
!
      END SUBROUTINE dlamc5


