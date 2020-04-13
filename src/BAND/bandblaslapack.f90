      SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,    &
     &                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
!     ..
!     .. Array Arguments ..
      LOGICAL            BWORK( * )
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
!     ..
!     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
!     ..
!
!  Purpose
!  =======
!
!  CGEES computes for an N-by-N complex nonsymmetric matrix A, the
!  eigenvalues, the Schur form T, and, optionally, the matrix of Schur
!  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
!
!  Optionally, it also orders the eigenvalues on the diagonal of the
!  Schur form so that selected eigenvalues are at the top left.
!  The leading columns of Z then form an orthonormal basis for the
!  invariant subspace corresponding to the selected eigenvalues.

!  A complex matrix is in Schur form if it is upper triangular.
!
!  Arguments
!  =========
!
!  JOBVS   (input) CHARACTER*1
!          = 'N': Schur vectors are not computed;
!          = 'V': Schur vectors are computed.
!
!  SORT    (input) CHARACTER*1
!          Specifies whether or not to order the eigenvalues on the
!          diagonal of the Schur form.
!          = 'N': Eigenvalues are not ordered:
!          = 'S': Eigenvalues are ordered (see SELECT).
!
!  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX argument
!          SELECT must be declared EXTERNAL in the calling subroutine.
!          If SORT = 'S', SELECT is used to select eigenvalues to order
!          to the top left of the Schur form.
!          IF SORT = 'N', SELECT is not referenced.
!          The eigenvalue W(j) is selected if SELECT(W(j)) is true.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten by its Schur form T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  SDIM    (output) INTEGER
!          If SORT = 'N', SDIM = 0.
!          If SORT = 'S', SDIM = number of eigenvalues for which
!                         SELECT is true.
!
!  W       (output) COMPLEX array, dimension (N)
!          W contains the computed eigenvalues, in the same order that
!          they appear on the diagonal of the output Schur form T.
!
!  VS      (output) COMPLEX array, dimension (LDVS,N)
!          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
!          vectors.
!          If JOBVS = 'N', VS is not referenced.
!
!  LDVS    (input) INTEGER
!          The leading dimension of the array VS.  LDVS >= 1; if
!          JOBVS = 'V', LDVS >= N.
!
!  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,2*N).
!          For good performance, LWORK must generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  RWORK   (workspace) REAL array, dimension (N)
!
!  BWORK   (workspace) LOGICAL array, dimension (N)
!          Not referenced if SORT = 'N'.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, and i is
!               <= N:  the QR algorithm failed to compute all the
!                      eigenvalues; elements 1:ILO-1 and i+1:N of W
!                      contain those eigenvalues which have converged;
!                      if JOBVS = 'V', VS contains the matrix which
!                      reduces A to its partially converged Schur form.
!               = N+1: the eigenvalues could not be reordered because
!                      some eigenvalues were too close to separate (the
!                      problem is very ill-conditioned);
!               = N+2: after reordering, roundoff changed values of
!                      some complex eigenvalues so that leading
!                      eigenvalues in the Schur form no longer satisfy
!                      SELECT = .TRUE..  This could also be caused by
!                      underflow due to scaling.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTST, WANTVS
      INTEGER            HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO, &
     &                   ITAU, IWRK, MAXWRK, MINWRK
      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
!     ..
!     .. Local Arrays ..
      REAL               DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY, &
     &                   CLASCL, CTRSEN, CUNGHR, SLABAD, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
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
!       HSWORK refers to the workspace preferred by CHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'CGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N
!
            CALL CHSEQR( 'S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS,      &
     &             WORK, -1, IEVAL )
            HSWORK = WORK( 1 )
!
            IF( .NOT.WANTVS ) THEN
               MAXWRK = MAX( MAXWRK, HSWORK )
            ELSE
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'CUNGHR', &
     &                       ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, HSWORK )
            END IF
         END IF
         WORK( 1 ) = MAXWRK
!
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
!
!     Get machine constants
!
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = CLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )                                                      &
     &   CALL CLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Permute the matrix to make it more nearly triangular
!     (CWorkspace: none)
!     (RWorkspace: need N)
!
      IBAL = 1
      CALL CGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
      ITAU = 1
      IWRK = N + ITAU
      CALL CGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),     &
     &             LWORK-IWRK+1, IERR )
!
      IF( WANTVS ) THEN
!
!        Copy Householder vectors to VS
!
         CALL CLACPY( 'L', N, N, A, LDA, VS, LDVS )
!
!        Generate unitary matrix in VS
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL CUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ),&
     &                LWORK-IWRK+1, IERR )
      END IF
!
      SDIM = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (CWorkspace: need 1, prefer HSWORK (see comments) )
!     (RWorkspace: none)
!
      IWRK = ITAU
      CALL CHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS,        &
     &             WORK( IWRK ), LWORK-IWRK+1, IEVAL )
      IF( IEVAL.GT.0 )                                                  &
     &   INFO = IEVAL
!
!     Sort eigenvalues if desired
!
      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA )                                                   &
     &      CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
         DO 10 I = 1, N
            BWORK( I ) = SELECT( W( I ) )
   10    CONTINUE
!
!        Reorder eigenvalues and transform Schur vectors
!        (CWorkspace: none)
!        (RWorkspace: none)
!
         CALL CTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM,  &
     &                S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND )
      END IF
!
      IF( WANTVS ) THEN
!
!        Undo balancing
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
         CALL CGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS,&
     &                IERR )
      END IF
!
      IF( SCALEA ) THEN
!
!        Undo scaling for the Schur form of A
!
         CALL CLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL CCOPY( N, A, LDA+1, W, 1 )
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of CGEES
!
      END





!> \brief \b CGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CGEBAK + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebak.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebak.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebak.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       REAL               SCALE( * )
!       COMPLEX            V( LDV, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBAK forms the right or left eigenvectors of a complex general
!> matrix by backward transformation on the computed eigenvectors of the
!> balanced matrix output by CGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N', do nothing, return immediately;
!>          = 'P', do backward transformation for permutation only;
!>          = 'S', do backward transformation for scaling only;
!>          = 'B', do backward transformations for both permutation and
!>                 scaling.
!>          JOB must be the same as the argument JOB supplied to CGEBAL.
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
!>          The integers ILO and IHI determined by CGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by CGEBAL.
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
!>          V is COMPLEX array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by CHSEIN or CTREVC.
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
!> \date November 2011
!
!> \ingroup complexGEcomputational
!
!  =====================================================================
      SUBROUTINE CGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,      &
     &                   INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      REAL               SCALE( * )
      COMPLEX            V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      REAL               S
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSSCAL, CSWAP, XERBLA
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
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.     &
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
         CALL XERBLA( 'CGEBAK', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
      IF( M.EQ.0 )                                                      &
     &   RETURN
      IF( LSAME( JOB, 'N' ) )                                           &
     &   RETURN
!
      IF( ILO.EQ.IHI )                                                  &
     &   GO TO 30
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL CSSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL CSSCAL( M, S, V( I, 1 ), LDV )
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
               IF( I.GE.ILO .AND. I.LE.IHI )                            &
     &            GO TO 40
               IF( I.LT.ILO )                                           &
     &            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I )                                             &
     &            GO TO 40
               CALL CSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI )                            &
     &            GO TO 50
               IF( I.LT.ILO )                                           &
     &            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I )                                             &
     &            GO TO 50
               CALL CSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEBAK
!
      END
!> \brief \b CGEBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CGEBAL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebal.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebal.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebal.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               SCALE( * )
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEBAL balances a general complex matrix A.  This involves, first,
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
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
!> \date November 2011
!
!> \ingroup complexGEcomputational
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
      SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      REAL               SCALE( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               SCLFAC
      PARAMETER          ( SCLFAC = 2.0E+0 )
      REAL               FACTOR
      PARAMETER          ( FACTOR = 0.95E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      REAL               C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
     &                   SFMIN2
      COMPLEX            CDUM
!     ..
!     .. External Functions ..
      LOGICAL            SISNAN, LSAME
      INTEGER            ICAMAX
      REAL               SLAMCH
      EXTERNAL           SISNAN, LSAME, ICAMAX, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSSCAL, CSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, MIN, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.     &
     &    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEBAL', -INFO )
         RETURN
      END IF
!
      K = 1
      L = N
!
      IF( N.EQ.0 )                                                      &
     &   GO TO 210
!
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
!
      IF( LSAME( JOB, 'S' ) )                                           &
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
      IF( J.EQ.M )                                                      &
     &   GO TO 30
!
      CALL CSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL CSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
!
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
!
!     Search for rows isolating an eigenvalue and push them down.
!
   40 CONTINUE
      IF( L.EQ.1 )                                                      &
     &   GO TO 210
      L = L - 1
!
   50 CONTINUE
      DO 70 J = L, 1, -1
!
         DO 60 I = 1, L
            IF( I.EQ.J )                                                &
     &         GO TO 60
            IF( REAL( A( J, I ) ).NE.ZERO .OR. AIMAG( A( J, I ) ).NE.   &
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
            IF( I.EQ.J )                                                &
     &         GO TO 100
            IF( REAL( A( I, J ) ).NE.ZERO .OR. AIMAG( A( I, J ) ).NE.   &
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
      IF( LSAME( JOB, 'P' ) )                                           &
     &   GO TO 210
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
      SFMIN1 = SLAMCH( 'S' ) / SLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
!
      DO 200 I = K, L
         C = ZERO
         R = ZERO
!
         DO 150 J = K, L
            IF( J.EQ.I )                                                &
     &         GO TO 150
            C = C + CABS1( A( J, I ) )
            R = R + CABS1( A( I, J ) )
  150    CONTINUE
         ICA = ICAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = ICAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!        Guard against zero C or R due to underflow.
!
         IF( C.EQ.ZERO .OR. R.EQ.ZERO )                                 &
     &      GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR.                 &
     &       MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
            IF( SISNAN( C+F+CA+R+G+RA ) ) THEN
!
!           Exit if NaN to avoid infinite loop
!
            INFO = -3
            CALL XERBLA( 'CGEBAL', -INFO )
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
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR.                    &
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
         IF( ( C+R ).GE.FACTOR*S )                                      &
     &      GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 )                                &
     &         GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F )                              &
     &         GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL CSSCAL( N-K+1, G, A( I, K ), LDA )
         CALL CSSCAL( L, F, A( 1, I ), 1 )
!
  200 CONTINUE
!
      IF( NOCONV )                                                      &
     &   GO TO 140
!
  210 CONTINUE
      ILO = K
      IHI = L
!
      RETURN
!
!     End of CGEBAL
!
      END
!> \brief \b CGEHRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CGEHRD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgehrd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgehrd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgehrd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEHRD reduces a complex general matrix A to upper Hessenberg form H by
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
!>          set by a previous call to CGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          TAU is COMPLEX array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!>          zero.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,N).
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
!> \date November 2011
!
!> \ingroup complexGEcomputational
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
!>  This file is a slight modification of LAPACK-3.0's DGEHRD
!>  subroutine incorporating improvements proposed by Quintana-Orti and
!>  Van de Geijn (2006). (See DLAHR2.)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),                   &
     &                     ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, LDWORK, LWKOPT, NB,      &
     &                   NBMIN, NH, NX
      COMPLEX            EI
!     ..
!     .. Local Arrays ..
      COMPLEX            T( LDT, NBMAX )
!     ..
!     .. External Subroutines ..
      EXTERNAL           CAXPY, CGEHD2, CGEMM, CLAHR2, CLARFB, CTRMM,   &
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
      NB = MIN( NBMAX, ILAENV( 1, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
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
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEHRD', -INFO )
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
      NB = MIN( NBMAX, ILAENV( 1, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code)
!
         NX = MAX( NB, ILAENV( 3, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
!
!           Determine if workspace is large enough for blocked code
!
            IWS = N*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code
!
               NBMIN = MAX( 2, ILAENV( 2, 'CGEHRD', ' ', N, ILO, IHI,   &
     &                 -1 ) )
               IF( LWORK.GE.N*NBMIN ) THEN
                  NB = LWORK / N
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
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V**H
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL CLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT,  &
     &                   WORK, LDWORK )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
!           to 1
!
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL CGEMM( 'No transpose', 'Conjugate transpose',          &
     &                  IHI, IHI-I-IB+1,                                &
     &                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, &
     &                  A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
!
!           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
!           right
!
            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',        &
     &                  'Unit', I, IB-1,                                &
     &                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL CAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1,              &
     &                     A( 1, I+J+1 ), 1 )
   30       CONTINUE
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL CLARFB( 'Left', 'Conjugate transpose', 'Forward',      &
     &                   'Columnwise',                                  &
     &                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, T, LDT, &
     &                   A( I+1, I+IB ), LDA, WORK, LDWORK )
   40    CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL CGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = IWS
!
      RETURN
!
!     End of CGEHRD
!
      END
!> \brief \b CHSEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CHSEQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chseqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chseqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chseqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
!                          WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!       CHARACTER          COMPZ, JOB
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CHSEQR computes the eigenvalues of a Hessenberg matrix H
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
!>           The order of the matrix H.  N .GE. 0.
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
!>           set by a previous call to CGEBAL, and then passed to ZGEHRD
!>           when the matrix output by CGEBAL is reduced to Hessenberg
!>           form. Otherwise ILO and IHI should be set to 1 and N
!>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and JOB = 'S', H contains the upper
!>           triangular matrix T from the Schur decomposition (the
!>           Schur form). If INFO = 0 and JOB = 'E', the contents of
!>           H are unspecified on exit.  (The output value of H when
!>           INFO.GT.0 is given under the description of INFO below.)
!>
!>           Unlike earlier versions of CHSEQR, this subroutine may
!>           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
!>           or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH .GE. max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>           The computed eigenvalues. If JOB = 'S', the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>           If COMPZ = 'N', Z is not referenced.
!>           If COMPZ = 'I', on entry Z need not be set and on exit,
!>           if INFO = 0, Z contains the unitary matrix Z of the Schur
!>           vectors of H.  If COMPZ = 'V', on entry Z must contain an
!>           N-by-N matrix Q, which is assumed to be equal to the unit
!>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
!>           if INFO = 0, Z contains Q*Z.
!>           Normally Q is the unitary matrix generated by CUNGHR
!>           after the call to CGEHRD which formed the Hessenberg matrix
!>           H. (The output value of Z when INFO.GT.0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if COMPZ = 'I' or
!>           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!>           On exit, if INFO = 0, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK .GE. max(1,N)
!>           is sufficient and delivers very good and sometimes
!>           optimal performance.  However, LWORK as large as 11*N
!>           may be required for optimal performance.  A workspace
!>           query is recommended to determine the optimal workspace
!>           size.
!>
!>           If LWORK = -1, then CHSEQR does a workspace query.
!>           In this case, CHSEQR checks the input parameters and
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
!>           .LT. 0:  if INFO = -i, the i-th argument had an illegal
!>                    value
!>           .GT. 0:  if INFO = i, CHSEQR failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO .GT. 0 and JOB = 'E', then on exit, the
!>                remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO .GT. 0 and JOB   = 'S', then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO .GT. 0 and COMPZ = 'V', then on exit
!>
!>                  (final value of Z)  =  (initial value of Z)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO .GT. 0 and COMPZ = 'I', then on exit
!>                      (final value of Z)  = U
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO .GT. 0 and COMPZ = 'N', then Z is not
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
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
!>             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
!>             It is suggested that these defaults be adjusted in order
!>             to attain best performance in each particular
!>             computational environment.
!>
!>            ISPEC=12: The CLAHQR vs CLAQR0 crossover point.
!>                      Default: 75. (Must be at least 11.)
!>
!>            ISPEC=13: Recommended deflation window size.
!>                      This depends on ILO, IHI and NS.  NS is the
!>                      number of simultaneous shifts returned
!>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
!>                      The default for (IHI-ILO+1).LE.500 is NS.
!>                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
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
!>                       CLAHQR and this parameter is ignored.  See
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
      SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,    &
     &                   WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    CLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
!
!     ==== NL allocates some local workspace to help small matrices
!     .    through a rare CLAHQR failure.  NL .GT. NTINY = 11 is
!     .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
!     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
!     .    allows up to six simultaneous shifts and a 16-by-16
!     .    deflation window.  ====
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0e0 )
!     ..
!     .. Local Arrays ..
      COMPLEX            HL( NL, NL ), WORKL( NL )
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
      EXTERNAL           CCOPY, CLACPY, CLAHQR, CLAQR0, CLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     ==== Decode and check the input parameters. ====
!
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = CMPLX( REAL( MAX( 1, N ) ), RZERO )
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
         CALL XERBLA( 'CHSEQR', -INFO )
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
         CALL CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z,&
     &                LDZ, WORK, LWORK, INFO )
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
         WORK( 1 ) = CMPLX( MAX( REAL( WORK( 1 ) ), REAL( MAX( 1,       &
     &               N ) ) ), RZERO )
         RETURN
!
      ELSE
!
!        ==== copy eigenvalues isolated by CGEBAL ====
!
         IF( ILO.GT.1 )                                                 &
     &      CALL CCOPY( ILO-1, H, LDH+1, W, 1 )
         IF( IHI.LT.N )                                                 &
     &      CALL CCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )
!
!        ==== Initialize Z, if requested ====
!
         IF( INITZ )                                                    &
     &      CALL CLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
!
!        ==== Quick return if possible ====
!
         IF( ILO.EQ.IHI ) THEN
            W( ILO ) = H( ILO, ILO )
            RETURN
         END IF
!
!        ==== CLAHQR/CLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'CHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N,    &
     &          ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== CLAQR0 for big matrices; CLAHQR for small ones ====
!
         IF( N.GT.NMIN ) THEN
            CALL CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,&
     &                   Z, LDZ, WORK, LWORK, INFO )
         ELSE
!
!           ==== Small matrix ====
!
            CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,&
     &                   Z, LDZ, INFO )
!
            IF( INFO.GT.0 ) THEN
!
!              ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds
!              .    when CLAHQR fails. ====
!
               KBOT = INFO
!
               IF( N.GE.NL ) THEN
!
!                 ==== Larger matrices have enough subdiagonal scratch
!                 .    space to call CLAQR0 directly. ====
!
                  CALL CLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, W,   &
     &                         ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
!
               ELSE
!
!                 ==== Tiny matrices don't have enough subdiagonal
!                 .    scratch space to benefit from CLAQR0.  Hence,
!                 .    tiny matrices must be copied into a larger
!                 .    array before calling CLAQR0. ====
!
                  CALL CLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL CLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), &
     &                         NL )
                  CALL CLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W,  &
     &                         ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                  IF( WANTT .OR. INFO.NE.0 )                            &
     &               CALL CLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
!
!        ==== Clear out the trash, if necessary. ====
!
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )                    &
     &      CALL CLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
!
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
!
         WORK( 1 ) = CMPLX( MAX( REAL( MAX( 1, N ) ),                   &
     &               REAL( WORK( 1 ) ) ), RZERO )
      END IF
!
!     ==== End of CHSEQR ====
!
      END
!> \brief \b CLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLACPY + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacpy.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacpy.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacpy.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACPY copies all or part of a two-dimensional matrix A to another
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          B is COMPLEX array, dimension (LDB,N)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
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
!     End of CLACPY
!
      END
!> \brief \b CLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLANGE + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clange.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clange.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clange.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               WORK( * )
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex matrix A.
!> \endverbatim
!>
!> \return CLANGE
!> \verbatim
!>
!>    CLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in CLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          CLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          CLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
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
!> \date September 2012
!
!> \ingroup complexGEauxiliary
!
!  =====================================================================
      REAL             FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL               SCALE, SUM, VALUE, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      EXTERNAL           LSAME, SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLASSQ
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
               IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
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
            IF( VALUE.LT.SUM .OR. SISNAN( SUM ) ) VALUE = SUM
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
            IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL CLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      CLANGE = VALUE
      RETURN
!
!     End of CLANGE
!
      END
!> \brief \b CLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLASCL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clascl.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clascl.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clascl.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       REAL               CFROM, CTO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASCL multiplies the M by N complex matrix A by the real scalar
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
!>                  bandwidth KU. See CGBTRF for storage details.
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
!>          CFROM is REAL
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is REAL
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL               CFROM, CTO
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH, SISNAN
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
      ELSE IF( CFROM.EQ.ZERO .OR. SISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( SISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.             &
     &         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.                 &
     &            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )   &
     &             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.                 &
     &            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.                 &
     &            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 )                                          &
     &   RETURN
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' )
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
      IF( .NOT.DONE )                                                   &
     &   GO TO 10
!
      RETURN
!
!     End of CLASCL
!
      END
!> \brief \b CTRSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CTRSEN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsen.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsen.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsen.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,
!                          SEP, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, JOB
!       INTEGER            INFO, LDQ, LDT, LWORK, M, N
!       REAL               S, SEP
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRSEN reorders the Schur factorization of a complex matrix
!> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
!> the leading positions on the diagonal of the upper triangular matrix
!> T, and the leading columns of Q form an orthonormal basis of the
!> corresponding right invariant subspace.
!>
!> Optionally the routine computes the reciprocal condition numbers of
!> the cluster of eigenvalues and/or the invariant subspace.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies whether condition numbers are required for the
!>          cluster of eigenvalues (S) or the invariant subspace (SEP):
!>          = 'N': none;
!>          = 'E': for eigenvalues only (S);
!>          = 'V': for invariant subspace only (SEP);
!>          = 'B': for both eigenvalues and invariant subspace (S and
!>                 SEP).
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V': update the matrix Q of Schur vectors;
!>          = 'N': do not update Q.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          SELECT specifies the eigenvalues in the selected cluster. To
!>          select the j-th eigenvalue, SELECT(j) must be set to .TRUE..
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
!>          T is COMPLEX array, dimension (LDT,N)
!>          On entry, the upper triangular matrix T.
!>          On exit, T is overwritten by the reordered matrix T, with the
!>          selected eigenvalues as the leading diagonal elements.
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
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          unitary transformation matrix which reorders T; the leading M
!>          columns of Q form an orthonormal basis for the specified
!>          invariant subspace.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>          The reordered eigenvalues of T, in the same order as they
!>          appear on the diagonal of T.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the specified invariant subspace.
!>          0 <= M <= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
!>          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
!>          condition number for the selected cluster of eigenvalues.
!>          S cannot underestimate the true reciprocal condition number
!>          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
!>          If JOB = 'N' or 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is REAL
!>          If JOB = 'V' or 'B', SEP is the estimated reciprocal
!>          condition number of the specified invariant subspace. If
!>          M = 0 or N, SEP = norm(T).
!>          If JOB = 'N' or 'E', SEP is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If JOB = 'N', LWORK >= 1;
!>          if JOB = 'E', LWORK = max(1,M*(N-M));
!>          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  CTRSEN first collects the selected eigenvalues by computing a unitary
!>  transformation Z to move them to the top left corner of T. In other
!>  words, the selected eigenvalues are the eigenvalues of T11 in:
!>
!>          Z**H * T * Z = ( T11 T12 ) n1
!>                         (  0  T22 ) n2
!>                            n1  n2
!>
!>  where N = n1+n2. The first
!>  n1 columns of Z span the specified invariant subspace of T.
!>
!>  If T has been obtained from the Schur factorization of a matrix
!>  A = Q*T*Q**H, then the reordered Schur factorization of A is given by
!>  A = (Q*Z)*(Z**H*T*Z)*(Q*Z)**H, and the first n1 columns of Q*Z span the
!>  corresponding invariant subspace of A.
!>
!>  The reciprocal condition number of the average of the eigenvalues of
!>  T11 may be returned in S. S lies between 0 (very badly conditioned)
!>  and 1 (very well conditioned). It is computed as follows. First we
!>  compute R so that
!>
!>                         P = ( I  R ) n1
!>                             ( 0  0 ) n2
!>                               n1 n2
!>
!>  is the projector on the invariant subspace associated with T11.
!>  R is the solution of the Sylvester equation:
!>
!>                        T11*R - R*T22 = T12.
!>
!>  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
!>  the two-norm of M. Then S is computed as the lower bound
!>
!>                      (1 + F-norm(R)**2)**(-1/2)
!>
!>  on the reciprocal of 2-norm(P), the true reciprocal condition number.
!>  S cannot underestimate 1 / 2-norm(P) by more than a factor of
!>  sqrt(N).
!>
!>  An approximate error bound for the computed average of the
!>  eigenvalues of T11 is
!>
!>                         EPS * norm(T) / S
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal condition number of the right invariant subspace
!>  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
!>  SEP is defined as the separation of T11 and T22:
!>
!>                     sep( T11, T22 ) = sigma-min( C )
!>
!>  where sigma-min(C) is the smallest singular value of the
!>  n1*n2-by-n1*n2 matrix
!>
!>     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
!>
!>  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
!>  product. We estimate sigma-min(C) by the reciprocal of an estimate of
!>  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
!>  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
!>
!>  When SEP is small, small changes in T can cause large changes in
!>  the invariant subspace. An approximate bound on the maximum angular
!>  error in the computed right invariant subspace is
!>
!>                      EPS * norm(T) / SEP
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,&
     &                   SEP, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LDQ, LDT, LWORK, M, N
      REAL               S, SEP
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, WANTBH, WANTQ, WANTS, WANTSP
      INTEGER            IERR, K, KASE, KS, LWMIN, N1, N2, NN
      REAL               EST, RNORM, SCALE
!     ..
!     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
      REAL               RWORK( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE
      EXTERNAL           LSAME, CLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLACN2, CLACPY, CTREXC, CTRSYL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )
!
!     Set M to the number of selected eigenvalues.
!
      M = 0
      DO 10 K = 1, N
         IF( SELECT( K ) )                                              &
     &      M = M + 1
   10 CONTINUE
!
      N1 = M
      N2 = N - M
      NN = N1*N2
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
!
      IF( WANTSP ) THEN
         LWMIN = MAX( 1, 2*NN )
      ELSE IF( LSAME( JOB, 'N' ) ) THEN
         LWMIN = 1
      ELSE IF( LSAME( JOB, 'E' ) ) THEN
         LWMIN = MAX( 1, NN )
      END IF
!
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP )   &
     &     THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -14
      END IF
!
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS )                                                    &
     &      S = ONE
         IF( WANTSP )                                                   &
     &      SEP = CLANGE( '1', N, N, T, LDT, RWORK )
         GO TO 40
      END IF
!
!     Collect the selected eigenvalues at the top left corner of T.
!
      KS = 0
      DO 20 K = 1, N
         IF( SELECT( K ) ) THEN
            KS = KS + 1
!
!           Swap the K-th eigenvalue to position KS.
!
            IF( K.NE.KS )                                               &
     &         CALL CTREXC( COMPQ, N, T, LDT, Q, LDQ, K, KS, IERR )
         END IF
   20 CONTINUE
!
      IF( WANTS ) THEN
!
!        Solve the Sylvester equation for R:
!
!           T11*R - R*T22 = scale*T12
!
         CALL CLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL CTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ),    &
     &                LDT, WORK, N1, SCALE, IERR )
!
!        Estimate the reciprocal of the condition number of the cluster
!        of eigenvalues.
!
         RNORM = CLANGE( 'F', N1, N2, WORK, N1, RWORK )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )*            &
     &          SQRT( RNORM ) )
         END IF
      END IF
!
      IF( WANTSP ) THEN
!
!        Estimate sep(T11,T22).
!
         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL CLACN2( NN, WORK( NN+1 ), WORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
!
!              Solve T11*R - R*T22 = scale*X.
!
               CALL CTRSYL( 'N', 'N', -1, N1, N2, T, LDT,               &
     &                      T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,      &
     &                      IERR )
            ELSE
!
!              Solve T11**H*R - R*T22**H = scale*X.
!
               CALL CTRSYL( 'C', 'C', -1, N1, N2, T, LDT,               &
     &                      T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,      &
     &                      IERR )
            END IF
            GO TO 30
         END IF
!
         SEP = SCALE / EST
      END IF
!
   40 CONTINUE
!
!     Copy reordered eigenvalues to W.
!
      DO 50 K = 1, N
         W( K ) = T( K, K )
   50 CONTINUE
!
      WORK( 1 ) = LWMIN
!
      RETURN
!
!     End of CTRSEN
!
      END
!> \brief \b CUNGHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNGHR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunghr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunghr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunghr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNGHR generates a complex unitary matrix Q which is defined as the
!> product of IHI-ILO elementary reflectors of order N, as returned by
!> CGEHRD:
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
!>          of CGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by CGEHRD.
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
!>          TAU is COMPLEX array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEHRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),                   &
     &                   ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CUNGQR, XERBLA
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
         NB = ILAENV( 1, 'CUNGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNGHR', -INFO )
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
         CALL CUNGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ),   &
     &                WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of CUNGHR
!
      END
!> \brief \b SLABAD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SLABAD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slabad.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slabad.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slabad.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLABAD( SMALL, LARGE )
! 
!       .. Scalar Arguments ..
!       REAL               LARGE, SMALL
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLABAD takes as input the values computed by SLAMCH for underflow and
!> overflow, and returns the square root of each of these values if the
!> log of LARGE is sufficiently large.  This subroutine is intended to
!> identify machines with a large exponent range, such as the Crays, and
!> redefine the underflow and overflow limits to be the square roots of
!> the values computed by SLAMCH.  This subroutine is needed because
!> SLAMCH does not compensate for poor arithmetic in the upper half of
!> the exponent range, as is found on a Cray.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SMALL
!> \verbatim
!>          SMALL is REAL
!>          On entry, the underflow threshold as computed by SLAMCH.
!>          On exit, if LOG10(LARGE) is sufficiently large, the square
!>          root of SMALL, otherwise unchanged.
!> \endverbatim
!>
!> \param[in,out] LARGE
!> \verbatim
!>          LARGE is REAL
!>          On entry, the overflow threshold as computed by SLAMCH.
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
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      REAL               LARGE, SMALL
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
      IF( LOG10( LARGE ).GT.2000. ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of SLABAD
!
      END
!> \brief \b CGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CGEHD2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgehd2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgehd2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgehd2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEHD2 reduces a complex general matrix A to upper Hessenberg form H
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
!>          set by a previous call to CGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= max(1,N).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          TAU is COMPLEX array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \date September 2012
!
!> \ingroup complexGEcomputational
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
      SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX            ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARF, CLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
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
         CALL XERBLA( 'CGEHD2', -INFO )
         RETURN
      END IF
!
      DO 10 I = ILO, IHI - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         ALPHA = A( I+1, I )
         CALL CLARFG( IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
         A( I+1, I ) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL CLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),     &
     &               A( 1, I+1 ), LDA, WORK )
!
!        Apply H(i)**H to A(i+1:ihi,i+1:n) from the left
!
         CALL CLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1,                &
     &               CONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK )
!
         A( I+1, I ) = ALPHA
   10 CONTINUE
!
      RETURN
!
!     End of CGEHD2
!
      END
!> \brief \b CLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLACN2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacn2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacn2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacn2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )
! 
!       .. Scalar Arguments ..
!       INTEGER            KASE, N
!       REAL               EST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISAVE( 3 )
!       COMPLEX            V( * ), X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACN2 estimates the 1-norm of a square, complex matrix A.
!> Reverse communication is used for evaluating matrix-vector products.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The order of the matrix.  N >= 1.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (N)
!>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!>         (W is not returned).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
!>         On an intermediate return, X should be overwritten by
!>               A * X,   if KASE=1,
!>               A**H * X,  if KASE=2,
!>         where A**H is the conjugate transpose of A, and CLACN2 must be
!>         re-called with all the other parameters unchanged.
!> \endverbatim
!>
!> \param[in,out] EST
!> \verbatim
!>          EST is REAL
!>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
!>         unchanged from the previous call to CLACN2.
!>         On exit, EST is an estimate (a lower bound) for norm(A). 
!> \endverbatim
!>
!> \param[in,out] KASE
!> \verbatim
!>          KASE is INTEGER
!>         On the initial call to CLACN2, KASE should be 0.
!>         On an intermediate return, KASE will be 1 or 2, indicating
!>         whether X should be overwritten by A * X  or A**H * X.
!>         On the final return from CLACN2, KASE will again be 0.
!> \endverbatim
!>
!> \param[in,out] ISAVE
!> \verbatim
!>          ISAVE is INTEGER array, dimension (3)
!>         ISAVE is used to save variables between calls to SLACN2
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Originally named CONEST, dated March 16, 1988.
!>
!>  Last modified:  April, 1999
!>
!>  This is a thread safe version of CLACON, which uses the array ISAVE
!>  in place of a SAVE statement, as follows:
!>
!>     CLACON     CLACN2
!>      JUMP     ISAVE(1)
!>      J        ISAVE(2)
!>      ITER     ISAVE(3)
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Nick Higham, University of Manchester
!
!> \par References:
!  ================
!>
!>  N.J. Higham, "FORTRAN codes for estimating the one-norm of
!>  a real or complex matrix, with applications to condition estimation",
!>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!>
!  =====================================================================
      SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISAVE( 3 )
      COMPLEX            V( * ), X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER              ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL                 ONE,         TWO
      PARAMETER          ( ONE = 1.0E0, TWO = 2.0E0 )
      COMPLEX              CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),                    &
     &                            CONE = ( 1.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, JLAST
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
!     .. Executable Statements ..
!
      SAFMIN = SLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 90, 120 )ISAVE( 1 )
!
!     ................ ENTRY   (ISAVE( 1 ) = 1)
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
      ISAVE( 1 ) = 2
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   40 CONTINUE
      ISAVE( 2 ) = ICMAX1( N, X, 1 )
      ISAVE( 3 ) = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = CONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 3)
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
      ISAVE( 1 ) = 4
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   90 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = ICMAX1( N, X, 1 )
      IF( ( ABS( X( JLAST ) ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND.         &
     &    ( ISAVE( 3 ).LT.ITMAX ) ) THEN
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE + REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 5)
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
!     End of CLACN2
!
      END
!> \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAHQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAHQR is an auxiliary routine called by CHSEQR to update the
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
!>          CLAHQR works primarily with the Hessenberg submatrix in rows
!>          and columns ILO to IHI, but applies transformations to all of
!>          H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
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
!>          W is COMPLEX array, dimension (N)
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
!>          Z is COMPLEX array, dimension (LDZ,N)
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
!>           =   0: successful exit
!>          .GT. 0: if INFO = i, CLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of W contain
!>                  those eigenvalues which have been successfully
!>                  computed.
!>
!>                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix
!>                  rows and columns ILO thorugh INFO of the final,
!>                  output value of H.
!>
!>                  If INFO .GT. 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthognal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
!>     This is a modified version of CLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,    &
     &                   IHIZ, Z, LDZ, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
!     ..
!
!  =========================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 30 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE, HALF
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0, HALF = 0.5e0 )
      REAL               DAT1
      PARAMETER          ( DAT1 = 3.0e0 / 4.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U, &
     &                   V2, X, Y
      REAL               AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,    &
     &                   SAFMIN, SMLNUM, SX, T2, TST, ULP
      INTEGER            I, I1, I2, ITS, J, JHI, JLO, K, L, M, NH, NZ
!     ..
!     .. Local Arrays ..
      COMPLEX            V( 2 )
!     ..
!     .. External Functions ..
      COMPLEX            CLADIV
      REAL               SLAMCH
      EXTERNAL           CLADIV, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY, CLARFG, CSCAL, SLABAD
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
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
      IF( ILO.LE.IHI-2 )                                                &
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
         IF( AIMAG( H( I, I-1 ) ).NE.RZERO ) THEN
!           ==== The following redundant normalization
!           .    avoids problems with both gradual and
!           .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
            SC = CONJG( SC ) / ABS( SC )
            H( I, I-1 ) = ABS( H( I, I-1 ) )
            CALL CSCAL( JHI-I+1, SC, H( I, I ), LDH )
            CALL CSCAL( MIN( JHI, I+1 )-JLO+1, CONJG( SC ), H( JLO, I ),&
     &                  1 )
            IF( WANTZ )                                                 &
     &         CALL CSCAL( IHIZ-ILOZ+1, CONJG( SC ), Z( ILOZ, I ), 1 )
         END IF
   20 CONTINUE
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( NH ) / ULP )
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
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   30 CONTINUE
      IF( I.LT.ILO )                                                    &
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
            IF( CABS1( H( K, K-1 ) ).LE.SMLNUM )                        &
     &         GO TO 50
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )                                         &
     &            TST = TST + ABS( REAL( H( K-1, K-2 ) ) )
               IF( K+1.LE.IHI )                                         &
     &            TST = TST + ABS( REAL( H( K+1, K ) ) )
            END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some examples.  ====
            IF( ABS( REAL( H( K, K-1 ) ) ).LE.ULP*TST ) THEN
               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               AA = MAX( CABS1( H( K, K ) ),                            &
     &              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( CABS1( H( K, K ) ),                            &
     &              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,                        &
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
         IF( L.GE.I )                                                   &
     &      GO TO 140
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
         IF( ITS.EQ.10 ) THEN
!
!           Exceptional shift.
!
            S = DAT1*ABS( REAL( H( L+1, L ) ) )
            T = S + H( L, L )
         ELSE IF( ITS.EQ.20 ) THEN
!
!           Exceptional shift.
!
            S = DAT1*ABS( REAL( H( I, I-1 ) ) )
            T = S + H( I, I )
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
                  IF( REAL( X / SX )*REAL( Y )+AIMAG( X / SX )*         &
     &                AIMAG( Y ).LT.RZERO )Y = -Y
               END IF
               T = T - U*CLADIV( U, ( X+Y ) )
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
            H21 = REAL( H( M+1, M ) )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = REAL( H( M, M-1 ) )
            IF( ABS( H10 )*ABS( H21 ).LE.ULP*                           &
     &          ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) )       &
     &          GO TO 70
   60    CONTINUE
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = REAL( H( L+1, L ) )
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
!           V(2) is always real before the call to CLARFG, and hence
!           after the call T2 ( = T1*V(2) ) is also real.
!
            IF( K.GT.M )                                                &
     &         CALL CCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL CLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = REAL( T1*V2 )
!
!           Apply G from the left to transform the rows of the matrix
!           in columns K to I2.
!
            DO 80 J = K, I2
               SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J )
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
               H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 )
   90       CONTINUE
!
            IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
               DO 100 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 )
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
               H( M+1, M ) = H( M+1, M )*CONJG( TEMP )
               IF( M+2.LE.I )                                           &
     &            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 110 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J )                                      &
     &                  CALL CSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL CSCAL( J-I1, CONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL CSCAL( NZ, CONJG( TEMP ), Z( ILOZ, J ), 1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
!
!        Ensure that H(I,I-1) is real.
!
         TEMP = H( I, I-1 )
         IF( AIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I )                                               &
     &         CALL CSCAL( I2-I, CONJG( TEMP ), H( I, I+1 ), LDH )
            CALL CSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL CSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
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
!
!     return to start of the main loop with new value of I.
!
      I = L - 1
      GO TO 30
!
  150 CONTINUE
      RETURN
!
!     End of CLAHQR
!
      END
!> \brief \b CLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAHR2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahr2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahr2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahr2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
! 
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LDT, LDY, N, NB
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ),
!      $                   Y( LDY, NB )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)
!> matrix A so that elements below the k-th subdiagonal are zero. The
!> reduction is performed by an unitary similarity transformation
!> Q**H * A * Q. The routine returns the matrices V and T which determine
!> Q as a block reflector I - V*T*v**H, and also the matrix Y = A * V * T.
!>
!> This is an auxiliary routine called by CGEHRD.
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
!>          A is COMPLEX array, dimension (LDA,N-K+1)
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
!>          TAU is COMPLEX array, dimension (NB)
!>          The scalar factors of the elementary reflectors. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,NB)
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
!>          Y is COMPLEX array, dimension (LDY,NB)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
!>  This subroutine is a slight modification of LAPACK-3.0's DLAHRD
!>  incorporating improvements proposed by Quintana-Orti and Van de
!>  Gejin. Note that the entries of A(1:K,2:NB) differ from those
!>  returned by the original LAPACK-3.0's DLAHRD routine. (This
!>  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)
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
      SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ),          &
     &                   Y( LDY, NB )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),                   &
     &                     ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX            EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CGEMM, CGEMV, CLACPY,            &
     &                   CLARFG, CSCAL, CTRMM, CTRMV, CLACGV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 )                                                      &
     &   RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(K+1:N,I)
!
!           Update I-th column of A - Y * V**H
!
            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA ) 
            CALL CGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY,  &
     &                  A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA ) 
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
            CALL CCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL CTRMV( 'Lower', 'Conjugate transpose', 'UNIT',         &
     &                  I-1, A( K+1, 1 ),                               &
     &                  LDA, T( 1, NB ), 1 )
!
!           w := w + V2**H * b2
!
            CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1,            &
     &                  ONE, A( K+I, 1 ),                               &
     &                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
!
!           w := T**H * w
!
            CALL CTRMV( 'Upper', 'Conjugate transpose', 'NON-UNIT',     &
     &                  I-1, T, LDT,                                    &
     &                  T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL CGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE,             &
     &                  A( K+I, 1 ),                                    &
     &                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL CTRMV( 'Lower', 'NO TRANSPOSE',                        &
     &                  'UNIT', I-1,                                    &
     &                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL CAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(I) to annihilate
!        A(K+I+1:N,I)
!
         CALL CLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
     &                TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
!
!        Compute  Y(K+1:N,I)
!
         CALL CGEMV( 'NO TRANSPOSE', N-K, N-K-I+1,                      &
     &               ONE, A( K+1, I+1 ),                                &
     &               LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
         CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1,               &
     &               ONE, A( K+I, 1 ), LDA,                             &
     &               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL CGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE,                    &
     &               Y( K+1, 1 ), LDY,                                  &
     &               T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
         CALL CSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
!
!        Compute T(1:I,I)
!
         CALL CSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL CTRMV( 'Upper', 'No Transpose', 'NON-UNIT',               &
     &               I-1, T, LDT,                                       &
     &               T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
!     Compute Y(1:K,1:NB)
!
      CALL CLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
      CALL CTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE',                     &
     &            'UNIT', K, NB,                                        &
     &            ONE, A( K+1, 1 ), LDA, Y, LDY )
      IF( N.GT.K+NB )                                                   &
     &   CALL CGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K,                 &
     &               NB, N-K-NB, ONE,                                   &
     &               A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y,    &
     &               LDY )
      CALL CTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE',                     &
     &            'NON-UNIT', K, NB,                                    &
     &            ONE, T, LDT, Y, LDY )
!
      RETURN
!
!     End of CLAHR2
!
      END
!> \brief \b CLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR0 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr0.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr0.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr0.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAQR0 computes the eigenvalues of a Hessenberg matrix H
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
!>           The order of the matrix H.  N .GE. 0.
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
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to CGEBAL, and then passed to CGEHRD when the
!>           matrix output by CGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H
!>           contains the upper triangular matrix T from the Schur
!>           decomposition (the Schur form). If INFO = 0 and WANT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO.GT.0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH .GE. max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
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
!>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO.GT.0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK .GE. max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then CLAQR0 does a workspace query.
!>           In this case, CLAQR0 checks the input parameters and
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
!>           .GT. 0:  if INFO = i, CLAQR0 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO .GT. 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO .GT. 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
      SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,    &
     &                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    CLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
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
      REAL               WILK1
      PARAMETER          ( WILK1 = 0.75e0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               TWO
      PARAMETER          ( TWO = 2.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      REAL               S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,   &
     &                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,     &
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
      COMPLEX            ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLACPY, CLAHQR, CLAQR3, CLAQR4, CLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, INT, MAX, MIN, MOD, REAL,   &
     &                   SQRT
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
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
!        ==== Tiny matrices must use CLAHQR. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 )                                              &
     &      CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,    &
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
!        .    point,  N .GT. NTINY = 11, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.3.) ====
!
         NWR = ILAENV( 13, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 11, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to CLAQR3 ====
!
         CALL CLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,   &
     &                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, &
     &                LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(CLAQR5, CLAQR3) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = CMPLX( LWKOPT, 0 )
            RETURN
         END IF
!
!        ==== CLAHQR/CLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
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
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
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
            IF( KBOT.LT.ILO )                                           &
     &         GO TO 80
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO )                                &
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
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT.                  &
     &                CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 )                                       &
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
            CALL CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
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
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-&
     &          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if CLAQR3
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    CLAQR3 above or from the eigenvalues
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
!                 ==== Got NS/2 or fewer shifts? Use CLAQR4 or
!                 .    CLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL CLACPY( 'A', NS, NS, H( KS, KS ), LDH,        &
     &                            H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL CLAQR4( .false., .false., NS, 1, NS,       &
     &                               H( KT, 1 ), LDH, W( KS ), 1, 1,    &
     &                               ZDUM, 1, WORK, LWORK, INF )
                     ELSE
                        CALL CLAHQR( .false., .false., NS, 1, NS,       &
     &                               H( KT, 1 ), LDH, W( KS ), 1, 1,    &
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
                        S = CABS1( H( KBOT-1, KBOT-1 ) ) +              &
     &                      CABS1( H( KBOT, KBOT-1 ) ) +                &
     &                      CABS1( H( KBOT-1, KBOT ) ) +                &
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
                        IF( SORTED )                                    &
     &                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) )   &
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
                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT.            &
     &                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                     W( KBOT-1 ) = W( KBOT )
                  ELSE
                     W( KBOT ) = W( KBOT-1 )
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnatiude
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
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,    &
     &                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK,  &
     &                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH,  &
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
      WORK( 1 ) = CMPLX( LWKOPT, 0 )
!
!     ==== End of CLAQR0 ====
!
      END
!> \brief \b CLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLARFB + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfb.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfb.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfb.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFB applies a complex block reflector H or its transpose H**H to a
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
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The matrix V. See Further Details.
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
!>          T is COMPLEX array, dimension (LDT,K)
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
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          WORK is COMPLEX array, dimension (LDWORK,K)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
      SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,  &
     &                   T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),         &
     &                   WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J, LASTV, LASTC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILACLR, ILACLC
      EXTERNAL           LSAME, ILACLR, ILACLC
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEMM, CLACGV, CTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 )                                          &
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
               LASTV = MAX( K, ILACLR( M, K, V, LDV ) )
               LASTC = ILACLC( LASTV, N, C, LDC )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C1**H
!
               DO 10 J = 1, K
                  CALL CCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( LASTC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**H *V2
!
                  CALL CGEMM( 'Conjugate transpose', 'No transpose',    &
     &                 LASTC, K, LASTV-K, ONE, C( K+1, 1 ), LDC,        &
     &                 V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit',        &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 LASTV-K, LASTC, K, -ONE, V( K+1, 1 ), LDV,       &
     &                 WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 30 J = 1, K
                  DO 20 I = 1, LASTC
                     C( J, I ) = C( J, I ) - CONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILACLR( N, K, V, LDV ) )
               LASTC = ILACLR( M, LASTV, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL CCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL CGEMM( 'No transpose', 'No transpose',           &
     &                 LASTC, K, LASTV-K,                               &
     &                 ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,         &
     &                 ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit',         &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 LASTC, LASTV-K, K,                               &
     &                 -ONE, WORK, LDWORK, V( K+1, 1 ), LDV,            &
     &                 ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, LASTC
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
               LASTC = ILACLC( M, N, C, LDC )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C2**H
!
               DO 70 J = 1, K
                  CALL CCOPY( LASTC, C( M-K+J, 1 ), LDC,                &
     &                 WORK( 1, J ), 1 )
                  CALL CLACGV( LASTC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V( M-K+1, 1 ), LDV,                  &
     &              WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**H*V1
!
                  CALL CGEMM( 'Conjugate transpose', 'No transpose',    &
     &                 LASTC, K, M-K, ONE, C, LDC, V, LDV,              &
     &                 ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit',        &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 M-K, LASTC, K, -ONE, V, LDV, WORK, LDWORK,       &
     &                 ONE, C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V( M-K+1, 1 ), LDV,          &
     &              WORK, LDWORK )
!
!              C2 := C2 - W**H
!
               DO 90 J = 1, K
                  DO 80 I = 1, LASTC
                     C( M-K+J, I ) = C( M-K+J, I ) -                    &
     &                               CONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTC = ILACLR( M, N, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL CCOPY( LASTC, C( 1, N-K+J ), 1,                  &
     &                 WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V( N-K+1, 1 ), LDV,                  &
     &              WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL CGEMM( 'No transpose', 'No transpose',           &
     &                 LASTC, K, N-K,                                   &
     &                 ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit',         &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 LASTC, N-K, K, -ONE, WORK, LDWORK, V, LDV,       &
     &                 ONE, C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V( N-K+1, 1 ), LDV,          &
     &              WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, LASTC
                     C( I, N-K+J ) = C( I, N-K+J )                      &
     &                    - WORK( I, J )
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
               LASTV = MAX( K, ILACLC( K, M, V, LDV ) )
               LASTC = ILACLC( LASTV, N, C, LDC )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C1**H
!
               DO 130 J = 1, K
                  CALL CCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL CLACGV( LASTC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1**H
!
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',     &
     &                     'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**H*V2**H
!
                  CALL CGEMM( 'Conjugate transpose',                    &
     &                 'Conjugate transpose', LASTC, K, LASTV-K,        &
     &                 ONE, C( K+1, 1 ), LDC, V( 1, K+1 ), LDV,         &
     &                 ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL CTRMM( 'Right', 'Upper', TRANST, 'Non-unit',        &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - V2**H * W**H
!
                  CALL CGEMM( 'Conjugate transpose',                    &
     &                 'Conjugate transpose', LASTV-K, LASTC, K,        &
     &                 -ONE, V( 1, K+1 ), LDV, WORK, LDWORK,            &
     &                 ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 150 J = 1, K
                  DO 140 I = 1, LASTC
                     C( J, I ) = C( J, I ) - CONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILACLC( K, N, V, LDV ) )
               LASTC = ILACLR( M, LASTV, C, LDC )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL CCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1**H
!
               CALL CTRMM( 'Right', 'Upper', 'Conjugate transpose',     &
     &                     'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 LASTC, K, LASTV-K, ONE, C( 1, K+1 ), LDC,        &
     &                 V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL CTRMM( 'Right', 'Upper', TRANS, 'Non-unit',         &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL CGEMM( 'No transpose', 'No transpose',           &
     &                 LASTC, LASTV-K, K,                               &
     &                 -ONE, WORK, LDWORK, V( 1, K+1 ), LDV,            &
     &                 ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL CTRMM( 'Right', 'Upper', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, LASTC
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
               LASTC = ILACLC( M, N, C, LDC )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C2**H
!
               DO 190 J = 1, K
                  CALL CCOPY( LASTC, C( M-K+J, 1 ), LDC,                &
     &                 WORK( 1, J ), 1 )
                  CALL CLACGV( LASTC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2**H
!
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V( 1, M-K+1 ), LDV,          &
     &              WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**H * V1**H
!
                  CALL CGEMM( 'Conjugate transpose',                    &
     &                 'Conjugate transpose', LASTC, K, M-K,            &
     &                 ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit',        &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1**H * W**H
!
                  CALL CGEMM( 'Conjugate transpose',                    &
     &                 'Conjugate transpose', M-K, LASTC, K,            &
     &                 -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V( 1, M-K+1 ), LDV,                  &
     &              WORK, LDWORK )
!
!              C2 := C2 - W**H
!
               DO 210 J = 1, K
                  DO 200 I = 1, LASTC
                     C( M-K+J, I ) = C( M-K+J, I ) -                    &
     &                               CONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTC = ILACLR( M, N, C, LDC )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL CCOPY( LASTC, C( 1, N-K+J ), 1,                  &
     &                 WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2**H
!
               CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',     &
     &              'Unit', LASTC, K, ONE, V( 1, N-K+1 ), LDV,          &
     &              WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1**H
!
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',    &
     &                 LASTC, K, N-K, ONE, C, LDC, V, LDV, ONE,         &
     &                 WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit',         &
     &              LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL CGEMM( 'No transpose', 'No transpose',           &
     &                 LASTC, N-K, K, -ONE, WORK, LDWORK, V, LDV,       &
     &                 ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL CTRMM( 'Right', 'Lower', 'No transpose', 'Unit',    &
     &              LASTC, K, ONE, V( 1, N-K+1 ), LDV,                  &
     &              WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, LASTC
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
!     End of CLARFB
!
      END
!> \brief \b CLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLASET + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claset.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claset.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claset.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       COMPLEX            ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASET initializes a 2-D array A to BETA on the diagonal and
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
!>          ALPHA is COMPLEX
!>          All the offdiagonal array elements are set to ALPHA.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>          All the diagonal array elements are set to BETA.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX            ALPHA, BETA
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
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
!     End of CLASET
!
      END
!> \brief \b CLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLASSQ + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/classq.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/classq.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/classq.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASSQ returns the values scl and ssq such that
!>
!>    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
!> assumed to be at least unity and the value of ssq will then satisfy
!>
!>    1.0 .le. ssq .le. ( sumsq + 2*n ).
!>
!> scale is assumed to be non-negative and scl returns the value
!>
!>    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
!>           i
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ respectively.
!> SCALE and SUMSQ are overwritten by scl and ssq respectively.
!>
!> The routine makes only one pass through the vector X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector X.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
!>          The vector x as described above.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector X.
!>          INCX > 0.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with the value  scl .
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is REAL
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with the value  ssq .
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      COMPLEX            X( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      REAL               TEMP1
!     ..
!     .. External Functions ..
      LOGICAL            SISNAN
      EXTERNAL           SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, REAL
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            TEMP1 = ABS( REAL( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.SISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            TEMP1 = ABS( AIMAG( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.SISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 .OR. SISNAN( TEMP1 ) ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
!
      RETURN
!
!     End of CLASSQ
!
      END
      SUBROUTINE CSSCAL(N,SA,CX,INCX)
!     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX CX(*)
!     ..
!
!  Purpose
!  =======
!
!     CSSCAL scales a complex vector by a real constant.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG,CMPLX,REAL
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DO I = 1,N
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      END IF
      RETURN
      END
!> \brief \b CTREXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CTREXC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrexc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrexc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrexc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ
!       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            Q( LDQ, * ), T( LDT, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTREXC reorders the Schur factorization of a complex matrix
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
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,N)
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
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          unitary transformation matrix Z which reorders T.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      COMPLEX            Q( LDQ, * ), T( LDT, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            K, M1, M2, M3
      REAL               CS
      COMPLEX            SN, T11, T22, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARTG, CROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
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
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -7
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTREXC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.1 .OR. IFST.EQ.ILST )                                    &
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
         CALL CLARTG( T( K, K+1 ), T22-T11, CS, SN, TEMP )
!
!        Apply transformation to the matrix T.
!
         IF( K+2.LE.N )                                                 &
     &      CALL CROT( N-K-1, T( K, K+2 ), LDT, T( K+1, K+2 ), LDT, CS, &
     &                 SN )
         CALL CROT( K-1, T( 1, K ), 1, T( 1, K+1 ), 1, CS, CONJG( SN ) )
!
         T( K, K ) = T22
         T( K+1, K+1 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL CROT( N, Q( 1, K ), 1, Q( 1, K+1 ), 1, CS,             &
     &                 CONJG( SN ) )
         END IF
!
   10 CONTINUE
!
      RETURN
!
!     End of CTREXC
!
      END
      SUBROUTINE CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!     .. Scalar Arguments ..
      COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  CTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
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
!              TRANSA = 'T' or 't'   op( A ) = A**T.
!
!              TRANSA = 'C' or 'c'   op( A ) = A**H.
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
!  ALPHA  - COMPLEX         .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
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
!  B      - COMPLEX          array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
      INTRINSIC CONJG,MAX
!     ..
!     .. Local Scalars ..
      COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE= (1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
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
               (.NOT.LSAME(TRANSA,'T')) .AND. &
               (.NOT.LSAME(TRANSA,'C'))) THEN
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
          CALL XERBLA('CTRMM ',INFO)
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
                              IF (NOUNIT) TEMP = TEMP*CONJG(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + CONJG(A(K,I))*B(K,J)
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
                              IF (NOUNIT) TEMP = TEMP*CONJG(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + CONJG(A(K,I))*B(K,J)
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
                                  TEMP = ALPHA*CONJG(A(J,K))
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
                              TEMP = TEMP*CONJG(A(K,K))
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
                                  TEMP = ALPHA*CONJG(A(J,K))
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
                              TEMP = TEMP*CONJG(A(K,K))
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
!     End of CTRMM .
!
      END
!> \brief \b CTRSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CTRSYL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsyl.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsyl.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsyl.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
!                          LDC, SCALE, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          TRANA, TRANB
!       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRSYL solves the complex Sylvester matrix equation:
!>
!>    op(A)*X + X*op(B) = scale*C or
!>    op(A)*X - X*op(B) = scale*C,
!>
!> where op(A) = A or A**H, and A and B are both upper triangular. A is
!> M-by-M and B is N-by-N; the right hand side C and the solution X are
!> M-by-N; and scale is an output scale factor, set <= 1 to avoid
!> overflow in X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANA
!> \verbatim
!>          TRANA is CHARACTER*1
!>          Specifies the option op(A):
!>          = 'N': op(A) = A    (No transpose)
!>          = 'C': op(A) = A**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] TRANB
!> \verbatim
!>          TRANB is CHARACTER*1
!>          Specifies the option op(B):
!>          = 'N': op(B) = B    (No transpose)
!>          = 'C': op(B) = B**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] ISGN
!> \verbatim
!>          ISGN is INTEGER
!>          Specifies the sign in the equation:
!>          = +1: solve op(A)*X + X*op(B) = scale*C
!>          = -1: solve op(A)*X - X*op(B) = scale*C
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A, and the number of rows in the
!>          matrices X and C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B, and the number of columns in the
!>          matrices X and C. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M)
!>          The upper triangular matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The upper triangular matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the M-by-N right hand side matrix C.
!>          On exit, C is overwritten by the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M)
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          The scale factor, scale, set <= 1 to avoid overflow in X.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1: A and B have common or very close eigenvalues; perturbed
!>               values were used to solve the equation (but the matrices
!>               A and B are unchanged).
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
!> \date November 2011
!
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,   &
     &                   LDC, SCALE, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      REAL               SCALE
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            J, K, L
      REAL               BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,      &
     &                   SMLNUM
      COMPLEX            A11, SUML, SUMR, VEC, X11
!     ..
!     .. Local Arrays ..
      REAL               DUM( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE, SLAMCH
      COMPLEX            CDOTC, CDOTU, CLADIV
      EXTERNAL           LSAME, CLANGE, SLAMCH, CDOTC, CDOTU, CLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSSCAL, SLABAD, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     Decode and Test input parameters
!
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
!
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRSYL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Set constants to control overflow
!
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*REAL( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*CLANGE( 'M', M, M, A, LDA, DUM ),         &
     &       EPS*CLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN
!
      IF( NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                        L-1
!          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
!                  I=K+1                      J=1
!
         DO 30 L = 1, N
            DO 20 K = M, 1, -1
!
               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,           &
     &                C( MIN( K+1, M ), L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
!
               SCALOC = ONE
               A11 = A( K, K ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )                               &
     &               SCALOC = ONE / DB
               END IF
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
               IF( SCALOC.NE.ONE ) THEN
                  DO 10 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   10             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
!
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A**H *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1                           L-1
!          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                           J=1
!
         DO 60 L = 1, N
            DO 50 K = 1, M
!
               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
!
               SCALOC = ONE
               A11 = CONJG( A( K, K ) ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )                               &
     &               SCALOC = ONE / DB
               END IF
!
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
!
   50       CONTINUE
   60    CONTINUE
!
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A**H*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        upper-right corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    K-1
!           R(K,L) = SUM [A**H(I,K)*X(I,L)] +
!                    I=1
!                           N
!                     ISGN*SUM [X(K,J)*B**H(L,J)].
!                          J=L+1
!
         DO 90 L = N, 1, -1
            DO 80 K = 1, M
!
               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,           &
     &                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
!
               SCALOC = ONE
               A11 = CONJG( A( K, K )+SGN*B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )                               &
     &               SCALOC = ONE / DB
               END IF
!
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
               IF( SCALOC.NE.ONE ) THEN
                  DO 70 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
!
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                          N
!          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
!                  I=K+1                      J=L+1
!
         DO 120 L = N, 1, -1
            DO 110 K = M, 1, -1
!
               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,           &
     &                C( MIN( K+1, M ), L ), 1 )
               SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,           &
     &                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
!
               SCALOC = ONE
               A11 = A( K, K ) + SGN*CONJG( B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )                               &
     &               SCALOC = ONE / DB
               END IF
!
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
!
               IF( SCALOC.NE.ONE ) THEN
                  DO 100 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
  100             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
!
  110       CONTINUE
  120    CONTINUE
!
      END IF
!
      RETURN
!
!     End of CTRSYL
!
      END
!> \brief \b CUNGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNGQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by CGEQRF.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by CGEQRF in the first k columns of its array
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,       &
     &                   LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARFB, CLARFT, CUNG2R, XERBLA
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
      NB = ILAENV( 1, 'CUNGQR', ' ', M, N, K, -1 )
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
         CALL XERBLA( 'CUNGQR', -INFO )
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
         NX = MAX( 0, ILAENV( 3, 'CUNGQR', ' ', M, N, K, -1 ) )
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
               NBMIN = MAX( 2, ILAENV( 2, 'CUNGQR', ' ', M, N, K, -1 ) )
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
      IF( KK.LT.N )                                                     &
     &   CALL CUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,           &
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
               CALL CLARFT( 'Forward', 'Columnwise', M-I+1, IB,         &
     &                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL CLARFB( 'Left', 'No transpose', 'Forward',          &
     &                      'Columnwise', M-I+1, N-I-IB+1, IB,          &
     &                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
     &                      LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL CUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
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
!     End of CUNGQR
!
      END
!> \brief \b SISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sisnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sisnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sisnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION SISNAN( SIN )
! 
!       .. Scalar Arguments ..
!       REAL               SIN
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIN
!> \verbatim
!>          SIN is REAL
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
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!!  =====================================================================
!      LOGICAL FUNCTION SISNAN( SIN )
!!
!!  -- LAPACK auxiliary routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      REAL               SIN
!!     ..
!!
!!  =====================================================================
!!
!!  .. External Functions ..
!      LOGICAL SLAISNAN
!      EXTERNAL SLAISNAN
!!  ..
!!  .. Executable Statements ..
!      SISNAN = SLAISNAN(SIN,SIN)
!      RETURN
!      END
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
!     ..
!
!  Purpose
!  =======
!
!     forms the dot product of two vectors, conjugating the first
!     vector.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack,  3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX CTEMP
      INTEGER I,IX,IY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG
!     ..
      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(I))*CY(I)
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
            CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTC = CTEMP
      RETURN
      END
      COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
!     ..
!
!  Purpose
!  =======
!
!     CDOTU forms the dot product of two vectors.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX CTEMP
      INTEGER I,IX,IY
!     ..
      CTEMP = (0.0,0.0)
      CDOTU = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
         DO I = 1,N
            CTEMP = CTEMP + CX(I)*CY(I)
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
            CTEMP = CTEMP + CX(IX)*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTU = CTEMP
      RETURN
      END
!> \brief \b CLACGV conjugates a complex vector.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLACGV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacgv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacgv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacgv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLACGV( N, X, INCX )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLACGV conjugates a complex vector of length N.
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
!>          X is COMPLEX array, dimension
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLACGV( N, X, INCX )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            X( * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IOFF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG
!     ..
!     .. Executable Statements ..
!
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )                                                &
     &      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
!
!     End of CLACGV
!
      END
!> \brief \b CLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLADIV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cladiv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cladiv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cladiv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       COMPLEX FUNCTION CLADIV( X, Y )
! 
!       .. Scalar Arguments ..
!       COMPLEX            X, Y
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!> will not overflow on an intermediary step unless the results
!> overflows.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is COMPLEX
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      COMPLEX FUNCTION CLADIV( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      COMPLEX            X, Y
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               ZI, ZR
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
!     ..
!     .. Executable Statements ..
!
      CALL SLADIV( REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR,    &
     &             ZI )
      CLADIV = CMPLX( ZR, ZI )
!
      RETURN
!
!     End of CLADIV
!
      END
!> \brief \b CLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR3 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr3.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr3.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr3.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
!                          NV, WV, LDWV, WORK, LWORK )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
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
!>    CLAQR3 accepts as input an upper Hessenberg matrix
!>    H and performs an unitary similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an unitary similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
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
!>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
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
!>          LDH is integer
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N .LE. LDH
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
!>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the unitary
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is integer
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 .LE. LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is integer
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is integer
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SH
!> \verbatim
!>          SH is COMPLEX array, dimension KBOT
!>          On output, approximate eigenvalues that may
!>          be used for shifts are stored in SH(KBOT-ND-NS+1)
!>          through SR(KBOT-ND).  Converged eigenvalues are
!>          stored in SH(KBOT-ND+1) through SH(KBOT).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is integer scalar
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW .LE. LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is integer scalar
!>          The number of columns of T.  NH.GE.NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is integer
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW .LE. LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is integer
!>          The number of rows of work array WV available for
!>          workspace.  NV.GE.NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is integer
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW .LE. LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension LWORK.
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; CLAQR3
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
     &                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,  &
     &                   NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,   &
     &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),&
     &                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            BETA, CDUM, S, TAU
      REAL               FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,  &
     &                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,      &
     &                   LWKOPT, NMIN
!     ..
!     .. External Functions ..
      REAL               SLAMCH
      INTEGER            ILAENV
      EXTERNAL           SLAMCH, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEHRD, CGEMM, CLACPY, CLAHQR, CLAQR4,  &
     &                   CLARF, CLARFG, CLASET, CTREXC, CUNMHR, SLABAD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, INT, MAX, MIN, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
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
!        ==== Workspace query call to CGEHRD ====
!
         CALL CGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to CUNMHR ====
!
         CALL CUNMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV,  &
     &                WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to CLAQR4 ====
!
         CALL CLAQR4( .true., .true., JW, 1, JW, T, LDT, SH, 1, JW, V,  &
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
         WORK( 1 ) = CMPLX( LWKOPT, 0 )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT )                                                &
     &   RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 )                                                     &
     &   RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N ) / ULP )
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
         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP,            &
     &       KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )                                         &
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
      CALL CLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL CCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL CLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      NMIN = ILAENV( 12, 'CLAQR3', 'SV', JW, 1, JW, LWORK )
      IF( JW.GT.NMIN ) THEN
         CALL CLAQR4( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,&
     &                JW, V, LDV, WORK, LWORK, INFQR )
      ELSE
         CALL CLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,&
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
         IF( FOO.EQ.RZERO )                                             &
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
!           .    way.   (CTREXC can not fail in this case.) ====
!
            IFST = NS
            CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
            ILST = ILST + 1
         END IF
   10 CONTINUE
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 )                                                     &
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
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) )     &
     &            IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST )                                          &
     &         CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
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
            CALL CCOPY( NS, V, LDV, WORK, 1 )
            DO 50 I = 1, NS
               WORK( I ) = CONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            CALL CLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL CLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL CLARF( 'L', NS, JW, WORK, 1, CONJG( TAU ), T, LDT,     &
     &                  WORK( JW+1 ) )
            CALL CLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,              &
     &                  WORK( JW+1 ) )
            CALL CLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,              &
     &                  WORK( JW+1 ) )
!
            CALL CGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),         &
     &                   LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 )                                               &
     &      H( KWTOP, KWTOP-1 ) = S*CONJG( V( 1, 1 ) )
         CALL CLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL CCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),       &
     &               LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO )                                  &
     &      CALL CUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, &
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
            CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),   &
     &                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL CLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   60    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL CGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,          &
     &                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL CLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),     &
     &                      LDH )
   70       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),&
     &                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL CLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),   &
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
      WORK( 1 ) = CMPLX( LWKOPT, 0 )
!
!     ==== End of CLAQR3 ====
!
      END
!> \brief \b CLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR4 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr4.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr4.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr4.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAQR4 implements one level of recursion for CLAQR0.
!>    It is a complete implementation of the small bulge multi-shift
!>    QR algorithm.  It may be called by CLAQR0 and, for large enough
!>    deflation window size, it may be called by CLAQR3.  This
!>    subroutine is identical to CLAQR0 except that it calls CLAQR2
!>    instead of CLAQR3.
!>
!>    CLAQR4 computes the eigenvalues of a Hessenberg matrix H
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
!>           The order of the matrix H.  N .GE. 0.
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
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to CGEBAL, and then passed to CGEHRD when the
!>           matrix output by CGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H
!>           contains the upper triangular matrix T from the Schur
!>           decomposition (the Schur form). If INFO = 0 and WANT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO.GT.0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH .GE. max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
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
!>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO.GT.0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK .GE. max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then CLAQR4 does a workspace query.
!>           In this case, CLAQR4 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!> \verbatim
!>          INFO is INTEGER
!>             =  0:  successful exit
!>           .GT. 0:  if INFO = i, CLAQR4 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO .GT. 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO .GT. 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
      SUBROUTINE CLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,    &
     &                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!
!  ================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    CLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
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
      REAL               WILK1
      PARAMETER          ( WILK1 = 0.75e0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               TWO
      PARAMETER          ( TWO = 2.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      REAL               S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,   &
     &                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,     &
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
      COMPLEX            ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLACPY, CLAHQR, CLAQR2, CLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, INT, MAX, MIN, MOD, REAL,   &
     &                   SQRT
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
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
!        ==== Tiny matrices must use CLAHQR. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 )                                              &
     &      CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,    &
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
!        .    point,  N .GT. NTINY = 11, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.3.) ====
!
         NWR = ILAENV( 13, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 11, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to CLAQR2 ====
!
         CALL CLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,   &
     &                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, &
     &                LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(CLAQR5, CLAQR2) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = CMPLX( LWKOPT, 0 )
            RETURN
         END IF
!
!        ==== CLAHQR/CLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
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
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
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
            IF( KBOT.LT.ILO )                                           &
     &         GO TO 80
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO )                                &
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
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT.                  &
     &                CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 )                                       &
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
            CALL CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
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
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-&
     &          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if CLAQR2
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    CLAQR2 above or from the eigenvalues
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
!                 ==== Got NS/2 or fewer shifts? Use CLAHQR
!                 .    on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL CLACPY( 'A', NS, NS, H( KS, KS ), LDH,        &
     &                            H( KT, 1 ), LDH )
                     CALL CLAHQR( .false., .false., NS, 1, NS,          &
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
                        S = CABS1( H( KBOT-1, KBOT-1 ) ) +              &
     &                      CABS1( H( KBOT, KBOT-1 ) ) +                &
     &                      CABS1( H( KBOT-1, KBOT ) ) +                &
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
                        IF( SORTED )                                    &
     &                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) )   &
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
                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT.            &
     &                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                     W( KBOT-1 ) = W( KBOT )
                  ELSE
                     W( KBOT ) = W( KBOT-1 )
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnatiude
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
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,    &
     &                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK,  &
     &                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH,  &
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
      WORK( 1 ) = CMPLX( LWKOPT, 0 )
!
!     ==== End of CLAQR4 ====
!
      END
!> \brief \b CLAQR5 performs a single small-bulge multi-shift QR sweep.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR5 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr5.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr5.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr5.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
!                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
!                          WV, LDWV, NH, WH, LDWH )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
!      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
!      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAQR5 called by CLAQR0 performs a
!>    single small-bulge multi-shift QR sweep.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is logical scalar
!>             WANTT = .true. if the triangular Schur factor
!>             is being computed.  WANTT is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is logical scalar
!>             WANTZ = .true. if the unitary Schur factor is being
!>             computed.  WANTZ is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] KACC22
!> \verbatim
!>          KACC22 is integer with value 0, 1, or 2.
!>             Specifies the computation mode of far-from-diagonal
!>             orthogonal updates.
!>        = 0: CLAQR5 does not accumulate reflections and does not
!>             use matrix-matrix multiply to update far-from-diagonal
!>             matrix entries.
!>        = 1: CLAQR5 accumulates reflections and uses matrix-matrix
!>             multiply to update the far-from-diagonal matrix entries.
!>        = 2: CLAQR5 accumulates reflections, uses matrix-matrix
!>             multiply to update the far-from-diagonal matrix entries,
!>             and takes advantage of 2-by-2 block structure during
!>             matrix multiplies.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer scalar
!>             N is the order of the Hessenberg matrix H upon which this
!>             subroutine operates.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is integer scalar
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is integer scalar
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
!>          NSHFTS is integer scalar
!>             NSHFTS gives the number of simultaneous shifts.  NSHFTS
!>             must be positive and even.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is COMPLEX array of size (NSHFTS)
!>             S contains the shifts of origin that define the multi-
!>             shift QR sweep.  On output S may be reordered.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array of size (LDH,N)
!>             On input H contains a Hessenberg matrix.  On output a
!>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
!>             to the isolated diagonal block in rows and columns KTOP
!>             through KBOT.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is integer scalar
!>             LDH is the leading dimension of H just as declared in the
!>             calling procedure.  LDH.GE.MAX(1,N).
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
!>             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array of size (LDZ,IHI)
!>             If WANTZ = .TRUE., then the QR Sweep unitary
!>             similarity transformation is accumulated into
!>             Z(ILOZ:IHIZ,ILO:IHI) from the right.
!>             If WANTZ = .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is integer scalar
!>             LDA is the leading dimension of Z just as declared in
!>             the calling procedure. LDZ.GE.N.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array of size (LDV,NSHFTS/2)
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is integer scalar
!>             LDV is the leading dimension of V as declared in the
!>             calling procedure.  LDV.GE.3.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array of size
!>             (LDU,3*NSHFTS-3)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is integer scalar
!>             LDU is the leading dimension of U just as declared in the
!>             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is integer scalar
!>             NH is the number of columns in array WH available for
!>             workspace. NH.GE.1.
!> \endverbatim
!>
!> \param[out] WH
!> \verbatim
!>          WH is COMPLEX array of size (LDWH,NH)
!> \endverbatim
!>
!> \param[in] LDWH
!> \verbatim
!>          LDWH is integer scalar
!>             Leading dimension of WH just as declared in the
!>             calling procedure.  LDWH.GE.3*NSHFTS-3.
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is integer scalar
!>             NV is the number of rows in WV agailable for workspace.
!>             NV.GE.1.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX array of size
!>             (LDWV,3*NSHFTS-3)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is integer scalar
!>             LDWV is the leading dimension of WV as declared in the
!>             in the calling subroutine.  LDWV.GE.NV.
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
!>
!  =====================================================================
      SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,&
     &                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,&
     &                   WV, LDWV, NH, WH, LDWH )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
     &                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), &
     &                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            ALPHA, BETA, CDUM, REFSUM
      REAL               H11, H12, H21, H22, SAFMAX, SAFMIN, SCL,       &
     &                   SMLNUM, TST1, TST2, ULP
      INTEGER            I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,    &
     &                   JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS,  &
     &                   M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL,&
     &                   NS, NU
      LOGICAL            ACCUM, BLK22, BMP22
!     ..
!     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, MOD, REAL
!     ..
!     .. Local Arrays ..
      COMPLEX            VT( 3 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLAQR1, CLARFG, CLASET, CTRMM,  &
     &                   SLABAD
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     ==== If there are no shifts, then there is nothing to do. ====
!
      IF( NSHFTS.LT.2 )                                                 &
     &   RETURN
!
!     ==== If the active block is empty or 1-by-1, then there
!     .    is nothing to do. ====
!
      IF( KTOP.GE.KBOT )                                                &
     &   RETURN
!
!     ==== NSHFTS is supposed to be even, but if it is odd,
!     .    then simply reduce it by one.  ====
!
      NS = NSHFTS - MOD( NSHFTS, 2 )
!
!     ==== Machine constants for deflation ====
!
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N ) / ULP )
!
!     ==== Use accumulated reflections to update far-from-diagonal
!     .    entries ? ====
!
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
!
!     ==== If so, exploit the 2-by-2 block structure? ====
!
      BLK22 = ( NS.GT.2 ) .AND. ( KACC22.EQ.2 )
!
!     ==== clear trash ====
!
      IF( KTOP+2.LE.KBOT )                                              &
     &   H( KTOP+2, KTOP ) = ZERO
!
!     ==== NBMPS = number of 2-shift bulges in the chain ====
!
      NBMPS = NS / 2
!
!     ==== KDU = width of slab ====
!
      KDU = 6*NBMPS - 3
!
!     ==== Create and chase chains of NBMPS bulges ====
!
      DO 210 INCOL = 3*( 1-NBMPS ) + KTOP - 1, KBOT - 2, 3*NBMPS - 2
         NDCOL = INCOL + KDU
         IF( ACCUM )                                                    &
     &      CALL CLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
!
!        ==== Near-the-diagonal bulge chase.  The following loop
!        .    performs the near-the-diagonal part of a small bulge
!        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
!        .    chunk extends from column INCOL to column NDCOL
!        .    (including both column INCOL and column NDCOL). The
!        .    following loop chases a 3*NBMPS column long chain of
!        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
!        .    may be less than KTOP and and NDCOL may be greater than
!        .    KBOT indicating phantom columns from which to chase
!        .    bulges before they are actually introduced or to which
!        .    to chase bulges beyond column KBOT.)  ====
!
         DO 140 KRCOL = INCOL, MIN( INCOL+3*NBMPS-3, KBOT-2 )
!
!           ==== Bulges number MTOP to MBOT are active double implicit
!           .    shift bulges.  There may or may not also be small
!           .    2-by-2 bulge, if there is room.  The inactive bulges
!           .    (if any) must wait until the active bulges have moved
!           .    down the diagonal to make room.  The phantom matrix
!           .    paradigm described above helps keep track.  ====
!
            MTOP = MAX( 1, ( ( KTOP-1 )-KRCOL+2 ) / 3+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL ) / 3 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+3*( M22-1 ) ).EQ.   &
     &              ( KBOT-2 )
!
!           ==== Generate reflections to chase the chain right
!           .    one column.  (The minimum value of K is KTOP-1.) ====
!
            DO 10 M = MTOP, MBOT
               K = KRCOL + 3*( M-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL CLAQR1( 3, H( KTOP, KTOP ), LDH, S( 2*M-1 ),     &
     &                         S( 2*M ), V( 1, M ) )
                  ALPHA = V( 1, M )
                  CALL CLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  CALL CLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
!
!                 ==== A Bulge may collapse because of vigilant
!                 .    deflation or destructive underflow.  In the
!                 .    underflow case, try the two-small-subdiagonals
!                 .    trick to try to reinflate the bulge.  ====
!
                  IF( H( K+3, K ).NE.ZERO .OR. H( K+3, K+1 ).NE.        &
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
                     CALL CLAQR1( 3, H( K+1, K+1 ), LDH, S( 2*M-1 ),    &
     &                            S( 2*M ), VT )
                     ALPHA = VT( 1 )
                     CALL CLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
                     REFSUM = CONJG( VT( 1 ) )*                         &
     &                        ( H( K+1, K )+CONJG( VT( 2 ) )*           &
     &                        H( K+2, K ) )
!
                     IF( CABS1( H( K+2, K )-REFSUM*VT( 2 ) )+           &
     &                   CABS1( REFSUM*VT( 3 ) ).GT.ULP*                &
     &                   ( CABS1( H( K, K ) )+CABS1( H( K+1,            &
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
!                       ==== Stating a new bulge here would
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
   10       CONTINUE
!
!           ==== Generate a 2-by-2 reflection, if needed. ====
!
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 ) THEN
               IF( K.EQ.KTOP-1 ) THEN
                  CALL CLAQR1( 2, H( K+1, K+1 ), LDH, S( 2*M22-1 ),     &
     &                         S( 2*M22 ), V( 1, M22 ) )
                  BETA = V( 1, M22 )
                  CALL CLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  CALL CLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               END IF
            END IF
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
            DO 30 J = MAX( KTOP, KRCOL ), JBOT
               MEND = MIN( MBOT, ( J-KRCOL+2 ) / 3 )
               DO 20 M = MTOP, MEND
                  K = KRCOL + 3*( M-1 )
                  REFSUM = CONJG( V( 1, M ) )*                          &
     &                     ( H( K+1, J )+CONJG( V( 2, M ) )*H( K+2, J )+&
     &                     CONJG( V( 3, M ) )*H( K+3, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
   20          CONTINUE
   30       CONTINUE
            IF( BMP22 ) THEN
               K = KRCOL + 3*( M22-1 )
               DO 40 J = MAX( K+1, KTOP ), JBOT
                  REFSUM = CONJG( V( 1, M22 ) )*                        &
     &                     ( H( K+1, J )+CONJG( V( 2, M22 ) )*          &
     &                     H( K+2, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
   40          CONTINUE
            END IF
!
!           ==== Multiply H by reflections from the right.
!           .    Delay filling in the last row until the
!           .    vigilant deflation check is complete. ====
!
            IF( ACCUM ) THEN
               JTOP = MAX( KTOP, INCOL )
            ELSE IF( WANTT ) THEN
               JTOP = 1
            ELSE
               JTOP = KTOP
            END IF
            DO 80 M = MTOP, MBOT
               IF( V( 1, M ).NE.ZERO ) THEN
                  K = KRCOL + 3*( M-1 )
                  DO 50 J = JTOP, MIN( KBOT, K+3 )
                     REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )*        &
     &                        H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
                     H( J, K+2 ) = H( J, K+2 ) -                        &
     &                             REFSUM*CONJG( V( 2, M ) )
                     H( J, K+3 ) = H( J, K+3 ) -                        &
     &                             REFSUM*CONJG( V( 3, M ) )
   50             CONTINUE
!
                  IF( ACCUM ) THEN
!
!                    ==== Accumulate U. (If necessary, update Z later
!                    .    with with an efficient matrix-matrix
!                    .    multiply.) ====
!
                     KMS = K - INCOL
                     DO 60 J = MAX( 1, KTOP-INCOL ), KDU
                        REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )*   &
     &                           U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                        U( J, KMS+2 ) = U( J, KMS+2 ) -                 &
     &                                  REFSUM*CONJG( V( 2, M ) )
                        U( J, KMS+3 ) = U( J, KMS+3 ) -                 &
     &                                  REFSUM*CONJG( V( 3, M ) )
   60                CONTINUE
                  ELSE IF( WANTZ ) THEN
!
!                    ==== U is not accumulated, so update Z
!                    .    now by multiplying by reflections
!                    .    from the right. ====
!
                     DO 70 J = ILOZ, IHIZ
                        REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )*     &
     &                           Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                        Z( J, K+2 ) = Z( J, K+2 ) -                     &
     &                                REFSUM*CONJG( V( 2, M ) )
                        Z( J, K+3 ) = Z( J, K+3 ) -                     &
     &                                REFSUM*CONJG( V( 3, M ) )
   70                CONTINUE
                  END IF
               END IF
   80       CONTINUE
!
!           ==== Special case: 2-by-2 reflection (if needed) ====
!
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 ) THEN
               IF ( V( 1, M22 ).NE.ZERO ) THEN
                  DO 90 J = JTOP, MIN( KBOT, K+3 )
                     REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )*    &
     &                        H( J, K+2 ) )
                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
                     H( J, K+2 ) = H( J, K+2 ) -                        &
     &                             REFSUM*CONJG( V( 2, M22 ) )
   90             CONTINUE
!
                  IF( ACCUM ) THEN
                     KMS = K - INCOL
                     DO 100 J = MAX( 1, KTOP-INCOL ), KDU
                        REFSUM = V( 1, M22 )*( U( J, KMS+1 )+           &
     &                           V( 2, M22 )*U( J, KMS+2 ) )
                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                        U( J, KMS+2 ) = U( J, KMS+2 ) -                 &
     &                                  REFSUM*CONJG( V( 2, M22 ) )
  100                CONTINUE
                  ELSE IF( WANTZ ) THEN
                     DO 110 J = ILOZ, IHIZ
                        REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )* &
     &                           Z( J, K+2 ) )
                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                        Z( J, K+2 ) = Z( J, K+2 ) -                     &
     &                                REFSUM*CONJG( V( 2, M22 ) )
  110                CONTINUE
                  END IF
               END IF
            END IF
!
!           ==== Vigilant deflation check ====
!
            MSTART = MTOP
            IF( KRCOL+3*( MSTART-1 ).LT.KTOP )                          &
     &         MSTART = MSTART + 1
            MEND = MBOT
            IF( BMP22 )                                                 &
     &         MEND = MEND + 1
            IF( KRCOL.EQ.KBOT-2 )                                       &
     &         MEND = MEND + 1
            DO 120 M = MSTART, MEND
               K = MIN( KBOT-1, KRCOL+3*( M-1 ) )
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
               IF( H( K+1, K ).NE.ZERO ) THEN
                  TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
                  IF( TST1.EQ.RZERO ) THEN
                     IF( K.GE.KTOP+1 )                                  &
     &                  TST1 = TST1 + CABS1( H( K, K-1 ) )
                     IF( K.GE.KTOP+2 )                                  &
     &                  TST1 = TST1 + CABS1( H( K, K-2 ) )
                     IF( K.GE.KTOP+3 )                                  &
     &                  TST1 = TST1 + CABS1( H( K, K-3 ) )
                     IF( K.LE.KBOT-2 )                                  &
     &                  TST1 = TST1 + CABS1( H( K+2, K+1 ) )
                     IF( K.LE.KBOT-3 )                                  &
     &                  TST1 = TST1 + CABS1( H( K+3, K+1 ) )
                     IF( K.LE.KBOT-4 )                                  &
     &                  TST1 = TST1 + CABS1( H( K+4, K+1 ) )
                  END IF
                  IF( CABS1( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) ) &
     &                 THEN
                     H12 = MAX( CABS1( H( K+1, K ) ),                   &
     &                     CABS1( H( K, K+1 ) ) )
                     H21 = MIN( CABS1( H( K+1, K ) ),                   &
     &                     CABS1( H( K, K+1 ) ) )
                     H11 = MAX( CABS1( H( K+1, K+1 ) ),                 &
     &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                     H22 = MIN( CABS1( H( K+1, K+1 ) ),                 &
     &                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )
!
                     IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE.       &
     &                   MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  END IF
               END IF
  120       CONTINUE
!
!           ==== Fill in the last row of each bulge. ====
!
            MEND = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 3 )
            DO 130 M = MTOP, MEND
               K = KRCOL + 3*( M-1 )
               REFSUM = V( 1, M )*V( 3, M )*H( K+4, K+3 )
               H( K+4, K+1 ) = -REFSUM
               H( K+4, K+2 ) = -REFSUM*CONJG( V( 2, M ) )
               H( K+4, K+3 ) = H( K+4, K+3 ) - REFSUM*CONJG( V( 3, M ) )
  130       CONTINUE
!
!           ==== End of near-the-diagonal bulge chase. ====
!
  140    CONTINUE
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
            IF( ( .NOT.BLK22 ) .OR. ( INCOL.LT.KTOP ) .OR.              &
     &          ( NDCOL.GT.KBOT ) .OR. ( NS.LE.2 ) ) THEN
!
!              ==== Updates not exploiting the 2-by-2 block
!              .    structure of U.  K1 and NU keep track of
!              .    the location and size of U in the special
!              .    cases of introducing bulges and chasing
!              .    bulges off the bottom.  In these special
!              .    cases and in case the number of shifts
!              .    is NS = 2, there is no 2-by-2 block
!              .    structure to exploit.  ====
!
               K1 = MAX( 1, KTOP-INCOL )
               NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
!
!              ==== Horizontal Multiply ====
!
               DO 150 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
                  CALL CGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ), &
     &                        LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH,  &
     &                        LDWH )
                  CALL CLACPY( 'ALL', NU, JLEN, WH, LDWH,               &
     &                         H( INCOL+K1, JCOL ), LDH )
  150          CONTINUE
!
!              ==== Vertical multiply ====
!
               DO 160 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
                  JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
                  CALL CGEMM( 'N', 'N', JLEN, NU, NU, ONE,              &
     &                        H( JROW, INCOL+K1 ), LDH, U( K1, K1 ),    &
     &                        LDU, ZERO, WV, LDWV )
                  CALL CLACPY( 'ALL', JLEN, NU, WV, LDWV,               &
     &                         H( JROW, INCOL+K1 ), LDH )
  160          CONTINUE
!
!              ==== Z multiply (also vertical) ====
!
               IF( WANTZ ) THEN
                  DO 170 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
                     CALL CGEMM( 'N', 'N', JLEN, NU, NU, ONE,           &
     &                           Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ), &
     &                           LDU, ZERO, WV, LDWV )
                     CALL CLACPY( 'ALL', JLEN, NU, WV, LDWV,            &
     &                            Z( JROW, INCOL+K1 ), LDZ )
  170             CONTINUE
               END IF
            ELSE
!
!              ==== Updates exploiting U's 2-by-2 block structure.
!              .    (I2, I4, J2, J4 are the last rows and columns
!              .    of the blocks.) ====
!
               I2 = ( KDU+1 ) / 2
               I4 = KDU
               J2 = I4 - I2
               J4 = KDU
!
!              ==== KZS and KNZ deal with the band of zeros
!              .    along the diagonal of one of the triangular
!              .    blocks. ====
!
               KZS = ( J4-J2 ) - ( NS+1 )
               KNZ = NS + 1
!
!              ==== Horizontal multiply ====
!
               DO 180 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
!
!                 ==== Copy bottom of H to top+KZS of scratch ====
!                  (The first KZS rows get multiplied by zero.) ====
!
                  CALL CLACPY( 'ALL', KNZ, JLEN, H( INCOL+1+J2, JCOL ), &
     &                         LDH, WH( KZS+1, 1 ), LDWH )
!
!                 ==== Multiply by U21**H ====
!
                  CALL CLASET( 'ALL', KZS, JLEN, ZERO, ZERO, WH, LDWH )
                  CALL CTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN, ONE,       &
     &                        U( J2+1, 1+KZS ), LDU, WH( KZS+1, 1 ),    &
     &                        LDWH )
!
!                 ==== Multiply top of H by U11**H ====
!
                  CALL CGEMM( 'C', 'N', I2, JLEN, J2, ONE, U, LDU,      &
     &                        H( INCOL+1, JCOL ), LDH, ONE, WH, LDWH )
!
!                 ==== Copy top of H to bottom of WH ====
!
                  CALL CLACPY( 'ALL', J2, JLEN, H( INCOL+1, JCOL ), LDH,&
     &                         WH( I2+1, 1 ), LDWH )
!
!                 ==== Multiply by U21**H ====
!
                  CALL CTRMM( 'L', 'L', 'C', 'N', J2, JLEN, ONE,        &
     &                        U( 1, I2+1 ), LDU, WH( I2+1, 1 ), LDWH )
!
!                 ==== Multiply by U22 ====
!
                  CALL CGEMM( 'C', 'N', I4-I2, JLEN, J4-J2, ONE,        &
     &                        U( J2+1, I2+1 ), LDU,                     &
     &                        H( INCOL+1+J2, JCOL ), LDH, ONE,          &
     &                        WH( I2+1, 1 ), LDWH )
!
!                 ==== Copy it back ====
!
                  CALL CLACPY( 'ALL', KDU, JLEN, WH, LDWH,              &
     &                         H( INCOL+1, JCOL ), LDH )
  180          CONTINUE
!
!              ==== Vertical multiply ====
!
               DO 190 JROW = JTOP, MAX( INCOL, KTOP ) - 1, NV
                  JLEN = MIN( NV, MAX( INCOL, KTOP )-JROW )
!
!                 ==== Copy right of H to scratch (the first KZS
!                 .    columns get multiplied by zero) ====
!
                  CALL CLACPY( 'ALL', JLEN, KNZ, H( JROW, INCOL+1+J2 ), &
     &                         LDH, WV( 1, 1+KZS ), LDWV )
!
!                 ==== Multiply by U21 ====
!
                  CALL CLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, LDWV )
                  CALL CTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,       &
     &                        U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),    &
     &                        LDWV )
!
!                 ==== Multiply by U11 ====
!
                  CALL CGEMM( 'N', 'N', JLEN, I2, J2, ONE,              &
     &                        H( JROW, INCOL+1 ), LDH, U, LDU, ONE, WV, &
     &                        LDWV )
!
!                 ==== Copy left of H to right of scratch ====
!
                  CALL CLACPY( 'ALL', JLEN, J2, H( JROW, INCOL+1 ), LDH,&
     &                         WV( 1, 1+I2 ), LDWV )
!
!                 ==== Multiply by U21 ====
!
                  CALL CTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,     &
     &                        U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), LDWV )
!
!                 ==== Multiply by U22 ====
!
                  CALL CGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,        &
     &                        H( JROW, INCOL+1+J2 ), LDH,               &
     &                        U( J2+1, I2+1 ), LDU, ONE, WV( 1, 1+I2 ), &
     &                        LDWV )
!
!                 ==== Copy it back ====
!
                  CALL CLACPY( 'ALL', JLEN, KDU, WV, LDWV,              &
     &                         H( JROW, INCOL+1 ), LDH )
  190          CONTINUE
!
!              ==== Multiply Z (also vertical) ====
!
               IF( WANTZ ) THEN
                  DO 200 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
!
!                    ==== Copy right of Z to left of scratch (first
!                    .     KZS columns get multiplied by zero) ====
!
                     CALL CLACPY( 'ALL', JLEN, KNZ,                     &
     &                            Z( JROW, INCOL+1+J2 ), LDZ,           &
     &                            WV( 1, 1+KZS ), LDWV )
!
!                    ==== Multiply by U12 ====
!
                     CALL CLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV,     &
     &                            LDWV )
                     CALL CTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,    &
     &                           U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ), &
     &                           LDWV )
!
!                    ==== Multiply by U11 ====
!
                     CALL CGEMM( 'N', 'N', JLEN, I2, J2, ONE,           &
     &                           Z( JROW, INCOL+1 ), LDZ, U, LDU, ONE,  &
     &                           WV, LDWV )
!
!                    ==== Copy left of Z to right of scratch ====
!
                     CALL CLACPY( 'ALL', JLEN, J2, Z( JROW, INCOL+1 ),  &
     &                            LDZ, WV( 1, 1+I2 ), LDWV )
!
!                    ==== Multiply by U21 ====
!
                     CALL CTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,  &
     &                           U( 1, I2+1 ), LDU, WV( 1, 1+I2 ),      &
     &                           LDWV )
!
!                    ==== Multiply by U22 ====
!
                     CALL CGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,     &
     &                           Z( JROW, INCOL+1+J2 ), LDZ,            &
     &                           U( J2+1, I2+1 ), LDU, ONE,             &
     &                           WV( 1, 1+I2 ), LDWV )
!
!                    ==== Copy the result back to Z ====
!
                     CALL CLACPY( 'ALL', JLEN, KDU, WV, LDWV,           &
     &                            Z( JROW, INCOL+1 ), LDZ )
  200             CONTINUE
               END IF
            END IF
         END IF
  210 CONTINUE
!
!     ==== End of CLAQR5 ====
!
      END
!> \brief \b CLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLARF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       COMPLEX            TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), V( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARF applies a complex elementary reflector H to a complex M-by-N
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!>
!> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
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
!>          V is COMPLEX array, dimension
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
!>          TAU is COMPLEX
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          WORK is COMPLEX array, dimension
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX            TAU
!     ..
!     .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMV, CGERC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILACLR, ILACLC
      EXTERNAL           LSAME, ILACLR, ILACLC
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
            LASTC = ILACLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILACLR(M, LASTV, C, LDC)
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
            CALL CGEMV( 'Conjugate transpose', LASTV, LASTC, ONE,       &
     &           C, LDC, V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
!
            CALL CGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            CALL CGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,      &
     &           V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
!
            CALL CGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of CLARF
!
      END
!> \brief \b CLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLARFG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       COMPLEX            ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFG generates a complex elementary reflector H of order n, such
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
!>          ALPHA is COMPLEX
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension
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
!>          TAU is COMPLEX
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            ALPHA, TAU
!     ..
!     .. Array Arguments ..
      COMPLEX            X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      REAL               SCNRM2, SLAMCH, SLAPY3
      COMPLEX            CLADIV
      EXTERNAL           SCNRM2, SLAMCH, SLAPY3, CLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSCAL, CSSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = SCNRM2( N-1, X, INCX )
      ALPHR = REAL( ALPHA )
      ALPHI = AIMAG( ALPHA )
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
         BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
!
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
   10       CONTINUE
            KNT = KNT + 1
            CALL CSSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )                                 &
     &         GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = SCNRM2( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI )
            BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
         CALL CSCAL( N-1, ALPHA, X, INCX )
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
!     End of CLARFG
!
      END
!> \brief \b CLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLARFT + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarft.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarft.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarft.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFT forms the triangular factor T of a complex block reflector H
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
!>          V is COMPLEX array, dimension
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,K)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
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
      SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACGV, CTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
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
                  CALL CGEMV( 'Conjugate transpose', J-I, I-1,          &
     &                        -TAU( I ), V( I+1, 1 ), LDV,              &
     &                        V( I+1, I ), 1,                           &
     &                        ONE, T( 1, I ), 1 )
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
                  CALL CGEMM( 'N', 'C', I-1, 1, J-I, -TAU( I ),         &
     &                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV,       &
     &                        ONE, T( 1, I ), LDT )                  
               END IF
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL CTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
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
                     CALL CGEMV( 'Conjugate transpose', N-K+I-J, K-I,   &
     &                           -TAU( I ), V( J, I+1 ), LDV, V( J, I ),&
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
                     CALL CGEMM( 'N', 'C', K-I, 1, N-K+I-J, -TAU( I ),  &
     &                           V( I+1, J ), LDV, V( I, J ), LDV,      &
     &                           ONE, T( I+1, I ), LDT )                     
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL CTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
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
!     End of CLARFT
!
      END
!> \brief \b CLARTG generates a plane rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLARTG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARTG( F, G, CS, SN, R )
! 
!       .. Scalar Arguments ..
!       REAL               CS
!       COMPLEX            F, G, R, SN
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARTG generates a plane rotation so that
!>
!>    [  CS  SN  ]     [ F ]     [ R ]
!>    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a faster version of the BLAS1 routine CROTG, except for
!> the following differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=1 and SN=0.
!>    If F=0, then CS=0 and SN is chosen so that R is real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is COMPLEX
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is COMPLEX
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is REAL
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is COMPLEX
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX
!>          The nonzero component of the rotated vector.
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      REAL               CS
      COMPLEX            F, G, R, SN
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0E+0, ONE = 1.0E+0, ZERO = 0.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
      INTEGER            COUNT, I
      REAL               D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN,      &
     &                   SAFMN2, SAFMX2, SCALE
      COMPLEX            FF, FS, GS
!     ..
!     .. External Functions ..
      REAL               SLAMCH, SLAPY2
      EXTERNAL           SLAMCH, SLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, INT, LOG, MAX, REAL, &
     &                   SQRT
!     ..
!     .. Statement Functions ..
      REAL               ABS1, ABSSQ
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Statement Function definitions ..
      ABS1( FF ) = MAX( ABS( REAL( FF ) ), ABS( AIMAG( FF ) ) )
      ABSSQ( FF ) = REAL( FF )**2 + AIMAG( FF )**2
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
         SAFMIN = SLAMCH( 'S' )
         EPS = SLAMCH( 'E' )
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /             &
     &            LOG( SLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
!        FIRST = .FALSE.
!     END IF
      SCALE = MAX( ABS1( F ), ABS1( G ) )
      FS = F
      GS = G
      COUNT = 0
      IF( SCALE.GE.SAFMX2 ) THEN
   10    CONTINUE
         COUNT = COUNT + 1
         FS = FS*SAFMN2
         GS = GS*SAFMN2
         SCALE = SCALE*SAFMN2
         IF( SCALE.GE.SAFMX2 )                                          &
     &      GO TO 10
      ELSE IF( SCALE.LE.SAFMN2 ) THEN
         IF( G.EQ.CZERO ) THEN
            CS = ONE
            SN = CZERO
            R = F
            RETURN
         END IF
   20    CONTINUE
         COUNT = COUNT - 1
         FS = FS*SAFMX2
         GS = GS*SAFMX2
         SCALE = SCALE*SAFMX2
         IF( SCALE.LE.SAFMN2 )                                          &
     &      GO TO 20
      END IF
      F2 = ABSSQ( FS )
      G2 = ABSSQ( GS )
      IF( F2.LE.MAX( G2, ONE )*SAFMIN ) THEN
!
!        This is a rare case: F is very small.
!
         IF( F.EQ.CZERO ) THEN
            CS = ZERO
            R = SLAPY2( REAL( G ), AIMAG( G ) )
!           Do complex/real division explicitly with two real divisions
            D = SLAPY2( REAL( GS ), AIMAG( GS ) )
            SN = CMPLX( REAL( GS ) / D, -AIMAG( GS ) / D )
            RETURN
         END IF
         F2S = SLAPY2( REAL( FS ), AIMAG( FS ) )
!        G2 and G2S are accurate
!        G2 is at least SAFMIN, and G2S is at least SAFMN2
         G2S = SQRT( G2 )
!        Error in CS from underflow in F2S is at most
!        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
!        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
!        and so CS .lt. sqrt(SAFMIN)
!        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
!        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
!        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
         CS = F2S / G2S
!        Make sure abs(FF) = 1
!        Do complex/real division explicitly with 2 real divisions
         IF( ABS1( F ).GT.ONE ) THEN
            D = SLAPY2( REAL( F ), AIMAG( F ) )
            FF = CMPLX( REAL( F ) / D, AIMAG( F ) / D )
         ELSE
            DR = SAFMX2*REAL( F )
            DI = SAFMX2*AIMAG( F )
            D = SLAPY2( DR, DI )
            FF = CMPLX( DR / D, DI / D )
         END IF
         SN = FF*CMPLX( REAL( GS ) / G2S, -AIMAG( GS ) / G2S )
         R = CS*F + SN*G
      ELSE
!
!        This is the most common case.
!        Neither F2 nor F2/G2 are less than SAFMIN
!        F2S cannot overflow, and it is accurate
!
         F2S = SQRT( ONE+G2 / F2 )
!        Do the F2S(real)*FS(complex) multiply with two real multiplies
         R = CMPLX( F2S*REAL( FS ), F2S*AIMAG( FS ) )
         CS = ONE / F2S
         D = F2 + G2
!        Do complex/real division explicitly with two real divisions
         SN = CMPLX( REAL( R ) / D, AIMAG( R ) / D )
         SN = SN*CONJG( GS )
         IF( COUNT.NE.0 ) THEN
            IF( COUNT.GT.0 ) THEN
               DO 30 I = 1, COUNT
                  R = R*SAFMX2
   30          CONTINUE
            ELSE
               DO 40 I = 1, -COUNT
                  R = R*SAFMN2
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
!
!     End of CLARTG
!
      END
!> \brief \b CROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CROT + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crot.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crot.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crot.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )
! 
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       REAL               C
!       COMPLEX            S
!       ..
!       .. Array Arguments ..
!       COMPLEX            CX( * ), CY( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CROT   applies a plane rotation, where the cos (C) is real and the
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
!>          CX is COMPLEX array, dimension (N)
!>          On input, the vector X.
!>          On output, CX is overwritten with C*X + S*Y.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of CY.  INCX <> 0.
!> \endverbatim
!>
!> \param[in,out] CY
!> \verbatim
!>          CY is COMPLEX array, dimension (N)
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
!>          C is REAL
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      REAL               C
      COMPLEX            S
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * ), CY( * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX            STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )                                   &
     &   GO TO 20
!
!     Code for unequal increments or equal increments not equal to 1
!
      IX = 1
      IY = 1
      IF( INCX.LT.0 )                                                   &
     &   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )                                                   &
     &   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         STEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - CONJG( S )*CX( IX )
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
         CY( I ) = C*CY( I ) - CONJG( S )*CX( I )
         CX( I ) = STEMP
   30 CONTINUE
      RETURN
      END
!> \brief \b CUNG2R
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNG2R + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cung2r.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cung2r.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cung2r.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNG2R generates an m by n complex matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by CGEQRF.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by CGEQRF in the first k columns of its array
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARF, CSCAL, XERBLA
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
         CALL XERBLA( 'CUNG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 )                                                      &
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
            CALL CLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),     &
     &                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )                                                   &
     &      CALL CSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
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
!     End of CUNG2R
!
      END
!> \brief \b ILACLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILACLC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILACLC( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILACLC scans A for its last non-zero column.
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILACLC( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX          ZERO
      PARAMETER ( ZERO = (0.0E+0, 0.0E+0) )
!     ..
!     .. Local Scalars ..
      INTEGER I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILACLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILACLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILACLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILACLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
!> \brief \b ILACLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILACLR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILACLR( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILACLR scans A for its last non-zero row.
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
!>          A is array, dimension (LDA,N)
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILACLR( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX          ZERO
      PARAMETER ( ZERO = (0.0E+0, 0.0E+0) )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILACLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILACLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILACLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILACLR = MAX( ILACLR, I )
         END DO
      END IF
      RETURN
      END
!> \brief \b SLAISNAN tests input for NaN by comparing two arguments for inequality.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SLAISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaisnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaisnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaisnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )
! 
!       .. Scalar Arguments ..
!       REAL               SIN1, SIN2
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in SISNAN.
!>
!> SLAISNAN checks for NaNs by comparing its two arguments for
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
!> \param[in] SIN1
!> \verbatim
!>          SIN1 is REAL
!> \endverbatim
!>
!> \param[in] SIN2
!> \verbatim
!>          SIN2 is REAL
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
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
!      LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )
!!
!!  -- LAPACK auxiliary routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      REAL               SIN1, SIN2
!!     ..
!!
!!  =====================================================================
!!
!!  .. Executable Statements ..
!      SLAISNAN = (SIN1.NE.SIN2)
!      RETURN
!      END
      SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!     .. Scalar Arguments ..
      COMPLEX ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      COMPLEX A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  CGERC  performs the rank 1 operation
!
!     A := alpha*x*y**H + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
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
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      COMPLEX TEMP
      INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
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
          CALL XERBLA('CGERC ',INFO)
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
                  TEMP = ALPHA*CONJG(Y(JY))
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
                  TEMP = ALPHA*CONJG(Y(JY))
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
!     End of CGERC .
!
      END
!> \brief \b CLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR1 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr1.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr1.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr1.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )
! 
!       .. Scalar Arguments ..
!       COMPLEX            S1, S2
!       INTEGER            LDH, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), V( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      Given a 2-by-2 or 3-by-3 matrix H, CLAQR1 sets v to a
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
!>          N is integer
!>              Order of the matrix H. N must be either 2 or 3.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX array of dimension (LDH,N)
!>              The 2-by-2 or 3-by-3 matrix H in (*).
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is integer
!>              The leading dimension of H as declared in
!>              the calling procedure.  LDH.GE.N
!> \endverbatim
!>
!> \param[in] S1
!> \verbatim
!>          S1 is COMPLEX
!> \endverbatim
!>
!> \param[in] S2
!> \verbatim
!>          S2 is COMPLEX
!>
!>          S1 and S2 are the shifts defining K in (*) above.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array of dimension N
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      COMPLEX            S1, S2
      INTEGER            LDH, N
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), V( * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            CDUM, H21S, H31S
      REAL               S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
      IF( N.EQ.2 ) THEN
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) )
         IF( S.EQ.RZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-S1 )*                 &
     &               ( ( H( 1, 1 )-S2 ) / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 )
         END IF
      ELSE
         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) ) +               &
     &       CABS1( H( 3, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-S1 )*( ( H( 1, 1 )-S2 ) / S ) +        &
     &               H( 1, 2 )*H21S + H( 1, 3 )*H31S
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 ) + H( 2, 3 )*H31S
            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-S1-S2 ) + H21S*H( 3, 2 )
         END IF
      END IF
      END
!> \brief \b CLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CLAQR2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
!                          NV, WV, LDWV, WORK, LWORK )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAQR2 is identical to CLAQR3 except that it avoids
!>    recursion by calling CLAHQR instead of CLAQR4.
!>
!>    Aggressive early deflation:
!>
!>    This subroutine accepts as input an upper Hessenberg matrix
!>    H and performs an unitary similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an unitary similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
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
!>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
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
!>          LDH is integer
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N .LE. LDH
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
!>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the unitary
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is integer
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 .LE. LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is integer
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is integer
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SH
!> \verbatim
!>          SH is COMPLEX array, dimension KBOT
!>          On output, approximate eigenvalues that may
!>          be used for shifts are stored in SH(KBOT-ND-NS+1)
!>          through SR(KBOT-ND).  Converged eigenvalues are
!>          stored in SH(KBOT-ND+1) through SH(KBOT).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is integer scalar
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW .LE. LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is integer scalar
!>          The number of columns of T.  NH.GE.NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is integer
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW .LE. LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is integer
!>          The number of rows of work array WV available for
!>          workspace.  NV.GE.NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is integer
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW .LE. LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension LWORK.
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; CLAQR2
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
!> \date September 2012
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
     &                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,  &
     &                   NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,   &
     &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),&
     &                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),                     &
     &                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
!     ..
!     .. Local Scalars ..
      COMPLEX            BETA, CDUM, S, TAU
      REAL               FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,  &
     &                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT
!     ..
!     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEHRD, CGEMM, CLACPY, CLAHQR, CLARF,   &
     &                   CLARFG, CLASET, CTREXC, CUNMHR, SLABAD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, INT, MAX, MIN, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
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
!        ==== Workspace query call to CGEHRD ====
!
         CALL CGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to CUNMHR ====
!
         CALL CUNMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV,  &
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
         WORK( 1 ) = CMPLX( LWKOPT, 0 )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT )                                                &
     &   RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 )                                                     &
     &   RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N ) / ULP )
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
         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP,            &
     &       KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )                                         &
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
      CALL CLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL CCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL CLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      CALL CLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,   &
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
         IF( FOO.EQ.RZERO )                                             &
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
!           .    way.   (CTREXC can not fail in this case.) ====
!
            IFST = NS
            CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
            ILST = ILST + 1
         END IF
   10 CONTINUE
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 )                                                     &
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
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) )     &
     &            IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST )                                          &
     &         CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
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
            CALL CCOPY( NS, V, LDV, WORK, 1 )
            DO 50 I = 1, NS
               WORK( I ) = CONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            CALL CLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL CLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL CLARF( 'L', NS, JW, WORK, 1, CONJG( TAU ), T, LDT,     &
     &                  WORK( JW+1 ) )
            CALL CLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,              &
     &                  WORK( JW+1 ) )
            CALL CLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,              &
     &                  WORK( JW+1 ) )
!
            CALL CGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),         &
     &                   LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 )                                               &
     &      H( KWTOP, KWTOP-1 ) = S*CONJG( V( 1, 1 ) )
         CALL CLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL CCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),       &
     &               LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO )                                  &
     &      CALL CUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, &
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
            CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),   &
     &                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL CLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   60    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL CGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,          &
     &                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL CLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),     &
     &                      LDH )
   70       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),&
     &                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL CLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),   &
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
      WORK( 1 ) = CMPLX( LWKOPT, 0 )
!
!     ==== End of CLAQR2 ====
!
      END
!> \brief \b CUNMHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNMHR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmhr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmhr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmhr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
!      $                   WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNMHR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> IHI-ILO elementary reflectors, as returned by CGEHRD:
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
!>          of CGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
!>          ILO = 1 and IHI = 0, if M = 0;
!>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
!>          ILO = 1 and IHI = 0, if N = 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by CGEHRD.
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
!>          TAU is COMPLEX array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEHRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,   &
     &                   LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),            &
     &                   WORK( * )
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
      EXTERNAL           ILAENV, LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CUNMQR, XERBLA
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
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) )&
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
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( LEFT ) THEN
            NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, NH, N, NH, -1 )
         ELSE
            NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, NH, NH, -1 )
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMHR', -INFO )
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
      CALL CUNMQR( SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA,       &
     &             TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO )
!
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of CUNMHR
!
      END
      REAL FUNCTION SCNRM2(N,X,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      COMPLEX X(*)
!     ..
!
!  Purpose
!  =======
!
!  SCNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     SCNRM2 := sqrt( x**H*x )
!
!  Further Details
!  ===============
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to CLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL NORM,SCALE,SSQ,TEMP
      INTEGER IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL,SQRT
!     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
!
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (REAL(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(REAL(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (AIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(AIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
!
      SCNRM2 = NORM
      RETURN
!
!     End of SCNRM2.
!
      END
!> \brief \b SLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SLADIV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sladiv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sladiv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sladiv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLADIV( A, B, C, D, P, Q )
! 
!       .. Scalar Arguments ..
!       REAL               A, B, C, D, P, Q
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLADIV performs complex division in  real arithmetic
!>
!>                       a + i*b
!>            p + i*q = ---------
!>                       c + i*d
!>
!> The algorithm is due to Robert L. Smith and can be found
!> in D. Knuth, The art of Computer Programming, Vol.2, p.195
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is REAL
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is REAL
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL
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
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               E, F
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of SLADIV
!
      END
!> \brief \b SLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SLAPY2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLAPY2( X, Y )
! 
!       .. Scalar Arguments ..
!       REAL               X, Y
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is REAL
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL
!>          X and Y specify the values x and y.
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
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
!      REAL             FUNCTION SLAPY2( X, Y )
!!
!!  -- LAPACK auxiliary routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      REAL               X, Y
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      REAL               ZERO
!      PARAMETER          ( ZERO = 0.0E0 )
!      REAL               ONE
!      PARAMETER          ( ONE = 1.0E0 )
!!     ..
!!     .. Local Scalars ..
!      REAL               W, XABS, YABS, Z
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          ABS, MAX, MIN, SQRT
!!     ..
!!     .. Executable Statements ..
!!
!      XABS = ABS( X )
!      YABS = ABS( Y )
!      W = MAX( XABS, YABS )
!      Z = MIN( XABS, YABS )
!      IF( Z.EQ.ZERO ) THEN
!         SLAPY2 = W
!      ELSE
!         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
!      END IF
!      RETURN
!!
!!     End of SLAPY2
!!
!      END
!> \brief \b SLAPY3 returns sqrt(x2+y2+z2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download SLAPY3 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy3.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy3.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy3.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLAPY3( X, Y, Z )
! 
!       .. Scalar Arguments ..
!       REAL               X, Y, Z
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is REAL
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL
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
!> \date September 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      REAL             FUNCTION SLAPY3( X, Y, Z )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      REAL               X, Y, Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
!     ..
!     .. Local Scalars ..
      REAL               W, XABS, YABS, ZABS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
!     W can be zero for max(0,nan,0)
!     adding all three entries together will make sure
!     NaN will not disappear.
         SLAPY3 =  XABS + YABS + ZABS
      ELSE
         SLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+              &
     &            ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of SLAPY3
!
      END


!> \brief \b CUNMQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNMQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
!      $                   WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNMQR overwrites the general complex M-by-N matrix C with
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
!> as returned by CGEQRF. Q is of order M if SIDE = 'L' and of order N
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
!>          A is COMPLEX array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CGEQRF in the first k columns of its array argument A.
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \date November 2011
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),            &
     &                   WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
     &                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      COMPLEX            T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARFB, CLARFT, CUNM2R, XERBLA
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
         NW = N
      ELSE
         NQ = N
         NW = M
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
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, N, K,  &
     &        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CUNMQR', -INFO )
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
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CUNMQR', SIDE // TRANS, M, N, K,&
     &              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,  &
     &                IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.                            &
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
            CALL CLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),&
     &                   LDA, TAU( I ), T, LDT )
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
            CALL CLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,  &
     &                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,  &
     &                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of CUNMQR
!
      END
!> \brief \b CUNM2R multiplies a general matrix by the unitary matrix from a QR factorization determined by cgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download CUNM2R + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunm2r.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunm2r.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunm2r.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNM2R overwrites the general complex m-by-n matrix C with
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
!> as returned by CGEQRF. Q is of order m if SIDE = 'L' and of order n
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
!>          A is COMPLEX array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CGEQRF in the first k columns of its array argument A.
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          WORK is COMPLEX array, dimension
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
!> \date September 2012
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX            AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
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
         CALL XERBLA( 'CUNM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )                              &
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
            TAUI = CONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL CLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC,&
     &               WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of CUNM2R
!
      END
