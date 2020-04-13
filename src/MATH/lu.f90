      subroutine cgetrf( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The
!                factorization has been completed, but the factor U
!                is exactly singular, and division by zero will occur
!                if it is used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
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
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
             CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,   &
     &                   N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),   &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1,   &
     &                       N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,     &
     &                       A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),   &
     &                       LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of CGETRF
!
      END

      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is
!               used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            ICAMAX
      EXTERNAL           ICAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGERU, CSCAL, CSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL CSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL CGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of CGETF2
!
      END
      SUBROUTINE cgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by CGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by CGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from CGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLASWP, CTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of CGETRS
!
      END
!      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!      CHARACTER*6        SRNAME
!      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! ==================================================================
!
!     .. Executable Statements ..
!
!      WRITE( *, FMT = 9999 )SRNAME, INFO
!
!      STOP
!
! 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ',I2,' had ',  &
!     &      'an illegal value' )
!
!     End of XERBLA
!
!      END
      subroutine  cswap (n,cx,incx,cy,incy)
!
!     interchanges two vectors.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
      subroutine  cscal(n,ca,cx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, linpack,  3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex ca,cx(*)
      integer i,incx,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
      SUBROUTINE CGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n
!  element vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the
!           matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program. LDA must be
!           at least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of CGERU .
!
      END
      SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through
!  K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in
!          positions K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! ==================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1
!     through K2.
!
      IF( INCX.EQ.0 )                                                   &
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of CLASWP
!
      END
      SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,A,LDA,   &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit,
!  or non-unit,  upper or lower triangular matrix  and  op( A )  is
!  one of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the
!           left or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper
!           or lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used
!           in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must
!           be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N
!           must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry,  ALPHA specifies the scalar  alpha. When alpha
!           is zero then  A is not referenced and  B need not be set
!           before entry.
!           Unchanged on exit.
!

!  A - COMPLEX array of DIMENSION ( LDA, k ), where k is m
!           when SIDE = 'L' or 'l' and is n when SIDE = 'R' or 'r'.
!           Before entry with UPLO = 'U' or 'u', the leading k by k
!           upper triangular part of the array A must contain the
!           upper triangular matrix and the strictly lower triangular
!           part of A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k
!           by k lower triangular part of the array  A must contain
!           the lower triangular matrix  and the strictly upper
!           triangular part of A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements
!           of A  are not referenced either,  but are assumed to be
!           unity.  Unchanged on exit.

!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program.  When
!           SIDE = 'L' or 'l'  then LDA  must be at least
!           max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX          array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array
!           B must contain  the  right-hand  side  matrix  B,  and
!           on exit  is overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as
!           declared in  the  calling  (sub)  program.   LDB  must
!           be  at  least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRSM .
!
      END
      SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,B,LDB,  &
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices,
!  with op( A ) an m by k matrix,  op( B )  a  k by n matrix and
!  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the
!           matrix op( A )  and of the  matrix  C.  M  must  be at
!           least  zero. Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the
!           matrix op( B ) and the number of columns of the matrix C.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the
!           matrix op( A ) and the number of rows of the matrix
!           op( B ). K must be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A   - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
!        k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!        part of the array  A  must contain the matrix  A,  otherwise
!        the leading  k by m  part of the array  A  must contain  the
!        matrix A.
!        Unchanged on exit.
!
!  LDA - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!        LDA must be at least  max( 1, m ), otherwise  LDA must be at
!        least  max( 1, k ).
!        Unchanged on exit.
!
!  B   - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
!        n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!        part of the array  B  must contain the matrix  B,  otherwise
!        the leading  n by k  part of the array  B  must contain  the
!        matrix B.
!        Unchanged on exit.
!
!  LDB - INTEGER.
!        On entry, LDB specifies the first dimension of B as declared
!        in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!        LDB must be at least  max( 1, k ), otherwise  LDB must be at
!        least  max( 1, n ).
!        Unchanged on exit.
!
!  BETA   - COMPLEX         .
!        On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!        supplied as zero then C need not be set on input.
!        Unchanged on exit.
!
!  C      - COMPLEX          array of DIMENSION ( LDC, n ).
!        Before entry, the leading  m by n  part of the array  C must
!        contain the matrix  C,  except when  beta  is zero, in which
!        case C need not be set on entry.
!        On exit, the array  C  is overwritten by the  m by n  matrix
!        ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!        On entry, LDC specifies the first dimension of C as declared
!        in  the  calling  (sub)  program.   LDC  must  be  at  least
!        max( 1, m ).
!        Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!  Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!  conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!  B  respectively are to be  transposed but  not conjugated  and set
!  NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!  and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE )))   &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                  TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMM .
!
      END

!     INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,  &
!    &                 N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!     CHARACTER*( * )    NAME, OPTS
!     INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-
!  dependent parameters for the local environment.  See ISPEC for
!  a description of the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!       Specifies the parameter to be returned as the value of
!       ILAENV.
!       = 1: the optimal blocksize; if this value is 1, an unblocked
!            algorithm will give the best performance.
!       = 2: the minimum block size for which the block routine
!            should be used; if the usable block size is less than
!            this value, an unblocked routine should be used.
!       = 3: the crossover point (in a block routine, for N less
!            than this value, an unblocked routine should be used)
!       = 4: the number of shifts, used in the nonsymmetric
!            eigenvalue routines
!       = 5: the minimum column dimension for blocking to be used;
!            rectangular blocks must have dimension at least k by m,
!            where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!       = 6: the crossover point for the SVD (when reducing an m by n
!            matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!            this value, a QR factorization is used first to reduce
!            the matrix to a triangular form.)
!       = 7: the number of processors
!       = 8: the crossover point for the multishift QR and QZ methods
!            for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not
!          all be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal
!                value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from
!  the LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in
!      determining the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the
!      order that they appear in the argument list for NAME.  N1 is
!      used first, N2 second, and so on, and unused problem dimensions
!      are passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity
!      in the calling subroutine.  For example, ILAENV is used to
!      retrieve the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  ==================================================================
!
!     .. Local Scalars ..
!     LOGICAL            CNAME, SNAME
!     CHARACTER*1        C1
!     CHARACTER*2        C2, C4
!     CHARACTER*3        C3
!     CHARACTER*6        SUBNAM
!     INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
!     ILAENV = -1
!     RETURN
!
! 100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
!     ILAENV = 1
!     SUBNAM = NAME
!     IC = ICHAR( SUBNAM( 1:1 ) )
!     IZ = ICHAR( 'Z' )
!     IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
!        IF( IC.GE.97 .AND. IC.LE.122 ) THEN
!           SUBNAM( 1:1 ) = CHAR( IC-32 )
!           DO 10 I = 2, 6
!              IC = ICHAR( SUBNAM( I:I ) )
!              IF( IC.GE.97 .AND. IC.LE.122 )                           &
!    &            SUBNAM( I:I ) = CHAR( IC-32 )
!  10       CONTINUE
!        END IF
!
!     ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
!        IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                         &
!    &       ( IC.GE.145 .AND. IC.LE.153 ) .OR.                         &
!    &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
!           SUBNAM( 1:1 ) = CHAR( IC+64 )
!           DO 20 I = 2, 6
!              IC = ICHAR( SUBNAM( I:I ) )
!              IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                   &
!    &             ( IC.GE.145 .AND. IC.LE.153 ) .OR.                   &
!    &             ( IC.GE.162 .AND. IC.LE.169 ) )                      &
!    &            SUBNAM( I:I ) = CHAR( IC+64 )
!  20       CONTINUE
!        END IF
!
!     ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
!        IF( IC.GE.225 .AND. IC.LE.250 ) THEN
!           SUBNAM( 1:1 ) = CHAR( IC-32 )
!           DO 30 I = 2, 6
!              IC = ICHAR( SUBNAM( I:I ) )
!              IF( IC.GE.225 .AND. IC.LE.250 )                          &
!    &            SUBNAM( I:I ) = CHAR( IC-32 )
!  30       CONTINUE
!        END IF
!     END IF
!
!     C1 = SUBNAM( 1:1 )
!     SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
!     CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
!     IF( .NOT.( CNAME .OR. SNAME ) )                                   &
!    &   RETURN
!     C2 = SUBNAM( 2:3 )
!     C3 = SUBNAM( 4:6 )
!     C4 = C3( 2:3 )
!
!     GO TO ( 110, 200, 300 ) ISPEC
!
! 110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
!     NB = 1
!
!     IF( C2.EQ.'GE' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.    &
!    &            C3.EQ.'QLF' ) THEN
!           IF( SNAME ) THEN
!              NB = 32
!           ELSE
!              NB = 32
!           END IF
!        ELSE IF( C3.EQ.'HRD' ) THEN
!           IF( SNAME ) THEN
!              NB = 32
!           ELSE
!              NB = 32
!           END IF
!        ELSE IF( C3.EQ.'BRD' ) THEN
!           IF( SNAME ) THEN
!              NB = 32
!           ELSE
!              NB = 32
!           END IF
!        ELSE IF( C3.EQ.'TRI' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'PO' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'SY' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
!           NB = 1
!        ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
!           NB = 64
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           NB = 64
!        ELSE IF( C3.EQ.'TRD' ) THEN
!           NB = 1
!        ELSE IF( C3.EQ.'GST' ) THEN
!           NB = 64
!        END IF
!     ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NB = 32
!           END IF
!        ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NB = 32
!           END IF
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NB = 32
!           END IF
!        ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NB = 32
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'GB' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              IF( N4.LE.64 ) THEN
!                 NB = 1
!              ELSE
!                 NB = 32
!              END IF
!           ELSE
!              IF( N4.LE.64 ) THEN
!                 NB = 1
!              ELSE
!                 NB = 32
!              END IF
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'PB' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              IF( N2.LE.64 ) THEN
!                 NB = 1
!              ELSE
!                 NB = 32
!              END IF
!           ELSE
!              IF( N2.LE.64 ) THEN
!                 NB = 1
!              ELSE
!                 NB = 32
!              END IF
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'TR' ) THEN
!        IF( C3.EQ.'TRI' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'LA' ) THEN
!        IF( C3.EQ.'UUM' ) THEN
!           IF( SNAME ) THEN
!              NB = 64
!           ELSE
!              NB = 64
!           END IF
!        END IF
!     ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
!        IF( C3.EQ.'EBZ' ) THEN
!           NB = 1
!        END IF
!     END IF
!     ILAENV = NB
!     RETURN
!
! 200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
!     NBMIN = 2
!     IF( C2.EQ.'GE' ) THEN
!        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
!    &       C3.EQ.'QLF' ) THEN
!           IF( SNAME ) THEN
!              NBMIN = 2
!           ELSE
!              NBMIN = 2
!           END IF
!        ELSE IF( C3.EQ.'HRD' ) THEN
!           IF( SNAME ) THEN
!              NBMIN = 2
!           ELSE
!              NBMIN = 2
!           END IF
!        ELSE IF( C3.EQ.'BRD' ) THEN
!           IF( SNAME ) THEN
!              NBMIN = 2
!           ELSE
!              NBMIN = 2
!           END IF
!        ELSE IF( C3.EQ.'TRI' ) THEN
!           IF( SNAME ) THEN
!              NBMIN = 2
!           ELSE
!              NBMIN = 2
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'SY' ) THEN
!        IF( C3.EQ.'TRF' ) THEN
!           IF( SNAME ) THEN
!              NBMIN = 8
!           ELSE
!              NBMIN = 8
!           END IF
!        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
!           NBMIN = 2
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
!        IF( C3.EQ.'TRD' ) THEN
!           NBMIN = 2
!        END IF
!     ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NBMIN = 2
!           END IF
!        ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NBMIN = 2
!           END IF
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NBMIN = 2
!           END IF
!        ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NBMIN = 2
!           END IF
!        END IF
!     END IF
!     ILAENV = NBMIN
!     RETURN
!
! 300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
!     NX = 0
!     IF( C2.EQ.'GE' ) THEN
!        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
!    &       C3.EQ.'QLF' ) THEN
!           IF( SNAME ) THEN
!              NX = 128
!           ELSE
!              NX = 128
!           END IF
!        ELSE IF( C3.EQ.'HRD' ) THEN
!           IF( SNAME ) THEN
!              NX = 128
!           ELSE
!              NX = 128
!           END IF
!        ELSE IF( C3.EQ.'BRD' ) THEN
!           IF( SNAME ) THEN
!              NX = 128
!           ELSE
!              NX = 128
!           END IF
!        END IF
!     ELSE IF( C2.EQ.'SY' ) THEN
!        IF( SNAME .AND. C3.EQ.'TRD' ) THEN
!           NX = 1
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
!        IF( C3.EQ.'TRD' ) THEN
!           NX = 1
!        END IF
!     ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NX = 128
!           END IF
!        END IF
!     ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
!        IF( C3( 1:1 ).EQ.'G' ) THEN
!           IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
!    &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
!    &          C4.EQ.'BR' ) THEN
!              NX = 128
!           END IF
!        END IF
!     END IF
!     ILAENV = NX
!     RETURN
!
! 400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
!     ILAENV = 6
!     RETURN
!
! 500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
!     ILAENV = 2
!     RETURN
!
! 600 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
!     ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
!     RETURN
!
! 700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
!     ILAENV = 1
!     RETURN
!
! 800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
!     ILAENV = 50
!     RETURN
!
!     End of ILAENV
!
!     END

      integer function icamax(n,cx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
!
      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end

!     LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     January 31, 1994
!
!     .. Scalar Arguments ..
!     CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! ==================================================================
!
!     .. Intrinsic Functions ..
!     INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
!     INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
!     LSAME = CA.EQ.CB
!     IF( LSAME )                                                       &
!    &   RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
!     ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
!     INTA = ICHAR( CA )
!     INTB = ICHAR( CB )
!
!     IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
!        IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
!        IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
!     ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
!        or upper case 'Z'.
!
!        IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
!    &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
!        IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
!    &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
!    &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
!     ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
!        IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
!        IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
!     END IF
!     LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
!     END
!       SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
! *
! *  -- LAPACK routine (version 2.0) --
! *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
! *     Courant Institute, Argonne National Lab, and Rice University
! *     September 30, 1994
! *
! *     .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
! *     ..
! *     .. Array Arguments ..
!       COMPLEX            A( LDA, * )
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  CTRTRI computes the inverse of a complex upper or lower
! *  triangular matrix A.
! *
! *  This is the Level 3 BLAS version of the algorithm.
! *
! *  Arguments
! *  =========
! *
! *  UPLO    (input) CHARACTER*1
! *          = 'U':  A is upper triangular;
! *          = 'L':  A is lower triangular.
! *
! *  DIAG    (input) CHARACTER*1
! *          = 'N':  A is non-unit triangular;
! *          = 'U':  A is unit triangular.
! *
! *  N       (input) INTEGER
! *          The order of the matrix A.  N >= 0.
! *
! *  A       (input/output) COMPLEX array, dimension (LDA,N)
! *       On entry, the triangular matrix A.  If UPLO = 'U', the
! *       leading N-by-N upper triangular part of the array A
! *       contains the upper triangular matrix, and the strictly lower
! *       triangular part of A is not referenced.  If UPLO = 'L', the
! *       leading N-by-N lower triangular part of the array A contains
! *       the lower triangular matrix, and the strictly upper
! *       triangular part of A is not referenced.  If DIAG = 'U', the
! *       diagonal elements of A are also not referenced and are
! *       assumed to be 1.
! *       On exit, the (triangular) inverse of the original matrix, in
! *       the same storage format.
! *
! *  LDA     (input) INTEGER
! *          The leading dimension of the array A.  LDA >= max(1,N).
! *
! *  INFO    (output) INTEGER
! *       = 0: successful exit
! *       < 0: if INFO = -i, the i-th argument had an illegal value
! *       > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
! *            matrix is singular and its inverse can not be computed.
! *
! *  ================================================================
! *
! *     .. Parameters ..
!       COMPLEX            ONE, ZERO
!       PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
!      $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
! *     ..
! *     .. Local Scalars ..
!       LOGICAL            NOUNIT, UPPER
!       INTEGER            J, JB, NB, NN
! *     ..
! *     .. External Functions ..
!       LOGICAL            LSAME
!       INTEGER            ILAENV
!       EXTERNAL           LSAME, ILAENV
! *     ..
! *     .. External Subroutines ..
!       EXTERNAL           CTRMM, CTRSM, CTRTI2, XERBLA
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC          MAX, MIN
! *     ..
! *     .. Executable Statements ..
! *
! *     Test the input parameters.
! *
!       INFO = 0
!       UPPER = LSAME( UPLO, 'U' )
!       NOUNIT = LSAME( DIAG, 'N' )
!       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!          INFO = -1
!       ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!          INFO = -2
!       ELSE IF( N.LT.0 ) THEN
!          INFO = -3
!       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!          INFO = -5
!       END IF
!       IF( INFO.NE.0 ) THEN
!          CALL XERBLA( 'CTRTRI', -INFO )
!          RETURN
!       END IF
! *
! *     Quick return if possible
! *
!       IF( N.EQ.0 )
!      $   RETURN
! *
! *     Check for singularity if non-unit.
! *
!       IF( NOUNIT ) THEN
!          DO 10 INFO = 1, N
!             IF( A( INFO, INFO ).EQ.ZERO )
!      $         RETURN
!    10    CONTINUE
!          INFO = 0
!       END IF
! *
! *     Determine the block size for this environment.
! *
!       NB = ILAENV( 1, 'CTRTRI', UPLO // DIAG, N, -1, -1, -1 )
!       IF( NB.LE.1 .OR. NB.GE.N ) THEN
! *
! *        Use unblocked code
! *
!          CALL CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!       ELSE
! *
! *        Use blocked code
! *
!          IF( UPPER ) THEN
! *
! *           Compute inverse of upper triangular matrix
! *
!             DO 20 J = 1, N, NB
!                JB = MIN( NB, N-J+1 )
! *
! *              Compute rows 1:j-1 of current block column
! *
!             CALL CTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, ONE, A, LDA, A( 1, J ), LDA )
!             CALL CTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
!    20       CONTINUE
!          ELSE
! *
! *           Compute inverse of lower triangular matrix
! *
!             NN = ( ( N-1 ) / NB )*NB + 1
!             DO 30 J = NN, 1, -NB
!                JB = MIN( NB, N-J+1 )
!                IF( J+JB.LE.N ) THEN
! *
! *                 Compute rows j+jb:n of current block column
! *
!                CALL CTRMM( 'Left', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
!      $                     A( J+JB, J ), LDA )
!                CALL CTRSM( 'Right', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, -ONE, A( J, J ), LDA,
!      $                     A( J+JB, J ), LDA )
!                END IF
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
!    30       CONTINUE
!          END IF
!       END IF
! *
!       RETURN
! *
! *     End of CTRTRI
! *
!       END
      SUBROUTINE CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRTI2 computes the inverse of a complex upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower
!          triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!       On entry, the triangular matrix A.  If UPLO = 'U', the
!       leading n by n upper triangular part of the array A contains
!       the upper triangular matrix, and the strictly lower
!       triangular part of A is not referenced.  If UPLO = 'L', the
!       leading n by n lower triangular part of the array A contains
!       the lower triangular matrix, and the strictly upper
!       triangular part of A is not referenced.  If DIAG = 'U', the
!       diagonal elements of A are also not referenced and are
!       assumed to be 1.
!
!       On exit, the (triangular) inverse of the original matrix, in
!       the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX            AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSCAL, CTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRTI2', -INFO )
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Compute inverse of upper triangular matrix.
!
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
!
!           Compute elements 1:j-1 of j-th column.
!
            CALL CTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,     &
     &                  A( 1, J ), 1 )
            CALL CSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
!
!        Compute inverse of lower triangular matrix.
!
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
               CALL CTRMV( 'Lower', 'No transpose', DIAG, N-J,          &
     &                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL CSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of CTRTI2
!
      END
      SUBROUTINE CTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  CTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!        Before entry with  UPLO = 'U' or 'u', the leading n by n
!        upper triangular part of the array A must contain the upper
!        triangular matrix and the strictly lower triangular part of
!        A is not referenced.
!        Before entry with UPLO = 'L' or 'l', the leading n by n
!        lower triangular part of the array A must contain the lower
!        triangular matrix and the strictly upper triangular part of
!        A is not referenced.
!        Note that when  DIAG = 'U' or 'u', the diagonal elements of
!        A are not referenced either, but are assumed to be unity.
!        Unchanged on exit.
!
!  LDA    - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. LDA must be at least
!        max( 1, n ).
!        Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.                            &
     &         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.                            &
     &         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = LSAME( DIAG , 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := A*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := A'*x  or  x := conjg( A' )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 90, I = J - 1, 1, -1
                        TEMP = TEMP + A( I, J )*X( I )
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 100, I = J - 1, 1, -1
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100                CONTINUE
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 120, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 130, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  130                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 150, I = J + 1, N
                        TEMP = TEMP + A( I, J )*X( I )
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 160, I = J + 1, N
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  160                CONTINUE
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 180, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 190, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  190                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRMV .
!
      END
      subroutine caxpy(n,ca,cx,incx,cy,incy)
!
!     constant times a vector plus a vector.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ca
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
      subroutine  ccopy(n,cx,incx,cy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end

      SUBROUTINE CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,           &
     &                   BETA, Y, INCY )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - COMPLEX         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )                   &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )                                               &
     &   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMV .
!
      END
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END
      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dmax
      integer i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      integer function izamax(n,zx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
!
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine  zcopy(n,zx,incx,zy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 4/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
      SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - COMPLEX*16      .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!     and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*DCONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +                                      &
     &                      DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of ZGEMM .
!
      END
      SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX*16       array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX*16       array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGERU .
!
      END
      subroutine  zscal(n,za,zx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex za,zx(*)
      integer i,incx,ix,n
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      subroutine  zswap (n,zx,incx,zy,incy)
!
!     interchanges two vectors.
!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*16         ALPHA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of ZTRSM .
!
      END
      SUBROUTINE CLACON( N, V, X, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
!     ..
!     .. Array Arguments ..
      COMPLEX            V( N ), X( N )
!     ..
!
!  Purpose
!  =======
!
!  CLACON estimates the 1-norm of a square, complex matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) COMPLEX array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) COMPLEX array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         where A' is the conjugate transpose of A, and CLACON must be
!         re-called with all the other parameters unchanged.
!
!  EST    (output) REAL
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to CLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from CLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named CONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  Last modified:  April, 1999
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E0, TWO = 2.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),                    &
     &                   CONE = ( 1.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
      INTEGER            ICMAX1
      REAL               SCSUM1, SLAMCH
      EXTERNAL           ICMAX1, SCSUM1, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      SAFMIN = SLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 90, 120 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 130
      END IF
      EST = SCSUM1( N, X, 1 )
!
      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   40 CONTINUE
      J = ICMAX1( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL CCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SCSUM1( N, V, 1 )
!
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )                                               &
     &   GO TO 100
!
      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   90 CONTINUE
      JLAST = J
      J = ICMAX1( N, X, 1 )
      IF( ( ABS( X( JLAST ) ).NE.ABS( X( J ) ) ) .AND.                  &
     &    ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  120 CONTINUE
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL CCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  130 CONTINUE
      KASE = 0
      RETURN
!
!     End of CLACON
!
      END
      SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTF2 computes an LU factorization of a complex m-by-n band matrix
!  A using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IZAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 )                                               &
     &         CALL ZSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,              &
     &                     AB( KV+1, J ), LDAB-1 )
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL ZSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J )                                            &
     &            CALL ZGERU( KM, JU-J, -ONE, AB( KV+2, J ), 1,         &
     &                        AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),   &
     &                        LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 )                                             &
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of ZGBTF2
!
      END
      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,  &
     &                   JU, K2, KM, KV, NB, NW
      COMPLEX*16         TEMP
!     ..
!     .. Local Arrays ..
      COMPLEX*16         WORK13( LDWORK, NBMAX ),                       &
     &                   WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, IZAMAX
      EXTERNAL           ILAENV, IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBTF2, ZGEMM, ZGERU, ZLASWP,   &
     &                   ZSCAL, ZSWAP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'ZGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IZAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL ZSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,     &
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,   &
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL ZSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,    &
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL ZSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
     &                        1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )                                        &
     &               CALL ZGERU( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,    &
     &                           AB( KV, JJ+1 ), LDAB-1,                &
     &                           AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 )                                       &
     &               INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,            &
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use ZLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL ZLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,     &
     &                      IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J2,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J2,&
     &                           JB, -ONE, WORK31, LDWORK,              &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J3,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
     &                           LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J3,&
     &                           JB, -ONE, WORK31, LDWORK, WORK13,      &
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, WORK31( 1, JJ-J+1 ), 1,               &
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGBTRF
!
      END
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of ZGETF2
!
      END
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
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
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZGETRF
!
      END
      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by ZGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by ZGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of ZGETRS
!
      END
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASWP
!
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,&
     &                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
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
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
     &               ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,&
     &               A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END
      REAL             FUNCTION SCSUM1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  SCSUM1 takes the sum of the absolute values of a complex
!  vector and returns a single precision result.
!
!  Based on SCASUM from the Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX > 0.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL               STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      SCSUM1 = 0.0E0
      STEMP = 0.0E0
      IF( N.LE.0 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 20
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   20 CONTINUE
      DO 30 I = 1, N
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     End of SCSUM1
!
      END
      INTEGER          FUNCTION ICMAX1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  ICMAX1 finds the index of the element whose real part has maximum
!  absolute value.
!
!  Based on ICAMAX from Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX >= 1.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX
      REAL               SMAX
      COMPLEX            ZDUM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
!
!     NEXT LINE IS THE ONLY MODIFICATION.
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
      ICMAX1 = 0
      IF( N.LT.1 )                                                      &
     &   RETURN
      ICMAX1 = 1
      IF( N.EQ.1 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 30
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      SMAX = CABS1( CX( 1 ) )
      IX = IX + INCX
      DO 20 I = 2, N
         IF( CABS1( CX( IX ) ).LE.SMAX )                                &
     &      GO TO 10
         ICMAX1 = I
         SMAX = CABS1( CX( IX ) )
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   30 CONTINUE
      SMAX = CABS1( CX( 1 ) )
      DO 40 I = 2, N
         IF( CABS1( CX( I ) ).LE.SMAX )                                 &
     &      GO TO 40
         ICMAX1 = I
         SMAX = CABS1( CX( I ) )
   40 CONTINUE
      RETURN
!
!     End of ICMAX1
!
      END
!     INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
!     INTEGER            ISPEC
!     REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
!     REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,    &
!    &                   NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
!     IEEECK = 1
!
!     POSINF = ONE / ZERO
!     IF( POSINF.LE.ONE ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     NEGINF = -ONE / ZERO
!     IF( NEGINF.GE.ZERO ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     NEGZRO = ONE / ( NEGINF+ONE )
!     IF( NEGZRO.NE.ZERO ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     NEGINF = ONE / NEGZRO
!     IF( NEGINF.GE.ZERO ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     NEWZRO = NEGZRO + ZERO
!     IF( NEWZRO.NE.ZERO ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     POSINF = ONE / NEWZRO
!     IF( POSINF.LE.ONE ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     NEGINF = NEGINF*POSINF
!     IF( NEGINF.GE.ZERO ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     POSINF = POSINF*POSINF
!     IF( POSINF.LE.ONE ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
!     IF( ISPEC.EQ.0 )                                                  &
!    &   RETURN
!
!     NAN1 = POSINF + NEGINF
!
!     NAN2 = POSINF / NEGINF
!
!     NAN3 = POSINF / POSINF
!
!     NAN4 = POSINF*ZERO
!
!     NAN5 = NEGINF*NEGZRO
!
!     NAN6 = NAN5*0.0
!
!     IF( NAN1.EQ.NAN1 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     IF( NAN2.EQ.NAN2 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     IF( NAN3.EQ.NAN3 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     IF( NAN4.EQ.NAN4 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     IF( NAN5.EQ.NAN5 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     IF( NAN6.EQ.NAN6 ) THEN
!        IEEECK = 0
!        RETURN
!     END IF
!
!     RETURN
!     END
!     REAL             FUNCTION SLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992 
!
!     .. Scalar Arguments ..
!      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
!     REAL               ONE, ZERO
!     PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST, LRND
!     INTEGER            BETA, IMAX, IMIN, IT
!     REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,&
!    &                   RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!     EXTERNAL           SLAMC2
!     ..
!     .. Save statement ..
!     SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,   &
!    &                   EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..

!     IF( FIRST ) THEN
!        FIRST = .FALSE.
!        CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
!        BASE = BETA
!        T = IT
!        IF( LRND ) THEN
!           RND = ONE
!           EPS = ( BASE**( 1-IT ) ) / 2
!        ELSE
!           RND = ZERO
!           EPS = BASE**( 1-IT )
!       END IF
!        PREC = EPS*BASE
!        EMIN = IMIN
!        EMAX = IMAX
!        SFMIN = RMIN
!        SMALL = ONE / RMAX
!        IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
!           SFMIN = SMALL*( ONE+EPS )
!        END IF
!     END IF
!
!     IF( LSAME( CMACH, 'E' ) ) THEN
!        RMACH = EPS
!     ELSE IF( LSAME( CMACH, 'S' ) ) THEN
!        RMACH = SFMIN
!     ELSE IF( LSAME( CMACH, 'B' ) ) THEN
!        RMACH = BASE
!     ELSE IF( LSAME( CMACH, 'P' ) ) THEN
!        RMACH = PREC
!     ELSE IF( LSAME( CMACH, 'N' ) ) THEN
!        RMACH = T
!     ELSE IF( LSAME( CMACH, 'R' ) ) THEN
!        RMACH = RND
!     ELSE IF( LSAME( CMACH, 'M' ) ) THEN
!        RMACH = EMIN
!     ELSE IF( LSAME( CMACH, 'U' ) ) THEN
!        RMACH = RMIN
!     ELSE IF( LSAME( CMACH, 'L' ) ) THEN
!        RMACH = EMAX
!     ELSE IF( LSAME( CMACH, 'O' ) ) THEN
!        RMACH = RMAX
!     END IF
!
!     SLAMCH = RMACH
!     RETURN
!
!     End of SLAMCH
!
!     END

!***********************************************************************
!
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = SLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )                                &
     &      LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of SLAMC1
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,      &
     &                   NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
     &                   SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,&
     &                   LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )                                                &
     &      B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS )                                                &
     &      LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.                   &
     &            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',        &
     &      '  EMIN = ', I8, /                                          &
     &      ' If, after inspection, the value EMIN looks',              &
     &      ' acceptable please comment out ',                          &
     &      / ' the IF block as marked within the code of routine',     &
     &      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of SLAMC2
!
      END
!
!***********************************************************************
!
!     REAL             FUNCTION SLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     REAL               A, B
!     ..
!
!  Purpose
!  =======
!
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
!     SLAMC3 = A + B
!
!     RETURN
!
!     End of SLAMC3
!
!     END
!
!***********************************************************************
!
      SUBROUTINE SLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
!     ..
!
!  Purpose
!  =======
!
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.         &
     &    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of SLAMC4
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
!     ..
!
!  Purpose
!  =======
!
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
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
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )                                                 &
     &      OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )                                                    &
     &   Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of SLAMC5
!
      END
