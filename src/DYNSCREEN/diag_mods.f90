!> \brief \b LA_CONSTANTS is a module for the scaling constants for the compiled Fortran single and double precisions
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date May 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!> Nick Papior, Technical University of Denmark, DK
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
module LA_CONSTANTS
!  -- LAPACK auxiliary module --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!  Standard constants for 
   integer, parameter :: sp = kind(1.e0)

   real(sp), parameter :: szero = 0.0_sp
   real(sp), parameter :: shalf = 0.5_sp
   real(sp), parameter :: sone = 1.0_sp
   real(sp), parameter :: stwo = 2.0_sp
   real(sp), parameter :: sthree = 3.0_sp
   real(sp), parameter :: sfour = 4.0_sp
   real(sp), parameter :: seight = 8.0_sp
   real(sp), parameter :: sten = 10.0_sp
   complex(sp), parameter :: czero = ( 0.0_sp, 0.0_sp )
   complex(sp), parameter :: chalf = ( 0.5_sp, 0.0_sp )
   complex(sp), parameter :: cone = ( 1.0_sp, 0.0_sp )
   character*1, parameter :: sprefix = 'S'
   character*1, parameter :: cprefix = 'C'

!  Scaling constants
   real(sp), parameter :: sulp = epsilon(0._sp)
   real(sp), parameter :: seps = sulp * 0.5_sp
   real(sp), parameter :: ssafmin = real(radix(0._sp),sp)**max( &
      minexponent(0._sp)-1, &
      1-maxexponent(0._sp) &
   )
   real(sp), parameter :: ssafmax = sone / ssafmin
   real(sp), parameter :: ssmlnum = ssafmin / sulp
   real(sp), parameter :: sbignum = ssafmax * sulp
   real(sp), parameter :: srtmin = sqrt(ssmlnum)
   real(sp), parameter :: srtmax = sqrt(sbignum)

!  Blue's scaling constants
   real(sp), parameter :: stsml = real(radix(0._sp), sp)**ceiling( &
       (minexponent(0._sp) - 1) * 0.5_sp)
   real(sp), parameter :: stbig = real(radix(0._sp), sp)**floor( &
       (maxexponent(0._sp) - digits(0._sp) + 1) * 0.5_sp)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly 
   real(sp), parameter :: sssml = real(radix(0._sp), sp)**( - floor( &
       (minexponent(0._sp) - digits(0._sp)) * 0.5_sp))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
   real(sp), parameter :: ssbig = real(radix(0._sp), sp)**( - ceiling( &
       (maxexponent(0._sp) + digits(0._sp) - 1) * 0.5_sp))

!  Standard constants for 
   integer, parameter :: dp = kind(1.d0)

   real(dp), parameter :: dzero = 0.0_dp
   real(dp), parameter :: dhalf = 0.5_dp
   real(dp), parameter :: done = 1.0_dp
   real(dp), parameter :: dtwo = 2.0_dp
   real(dp), parameter :: dthree = 3.0_dp
   real(dp), parameter :: dfour = 4.0_dp
   real(dp), parameter :: deight = 8.0_dp
   real(dp), parameter :: dten = 10.0_dp
   complex(dp), parameter :: zzero = ( 0.0_dp, 0.0_dp )
   complex(dp), parameter :: zhalf = ( 0.5_dp, 0.0_dp )
   complex(dp), parameter :: zone = ( 1.0_dp, 0.0_dp )
   character*1, parameter :: dprefix = 'D'
   character*1, parameter :: zprefix = 'Z'

!  Scaling constants
   real(dp), parameter :: dulp = epsilon(0._dp)
   real(dp), parameter :: deps = dulp * 0.5_dp
   real(dp), parameter :: dsafmin = real(radix(0._dp),dp)**max( &
      minexponent(0._dp)-1, &
      1-maxexponent(0._dp) &
   )
   real(dp), parameter :: dsafmax = done / dsafmin
   real(dp), parameter :: dsmlnum = dsafmin / dulp
   real(dp), parameter :: dbignum = dsafmax * dulp
   real(dp), parameter :: drtmin = sqrt(dsmlnum)
   real(dp), parameter :: drtmax = sqrt(dbignum)

!  Blue's scaling constants
   real(dp), parameter :: dtsml = real(radix(0._dp), dp)**ceiling( &
       (minexponent(0._dp) - 1) * 0.5_dp)
   real(dp), parameter :: dtbig = real(radix(0._dp), dp)**floor( &
       (maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly 
   real(dp), parameter :: dssml = real(radix(0._dp), dp)**( - floor( &
       (minexponent(0._dp) - digits(0._dp)) * 0.5_dp))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
   real(dp), parameter :: dsbig = real(radix(0._dp), dp)**( - ceiling( &
       (maxexponent(0._dp) + digits(0._dp) - 1) * 0.5_dp))

end module LA_CONSTANTS
module LA_XISNAN
   interface LA_ISNAN

   module procedure SISNAN
   module procedure DISNAN

   end interface

contains
   
   logical function SISNAN( x )
   use LA_CONSTANTS, only: wp=>sp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   sisnan = ieee_is_nan(x)
#elif USE_ISNAN
   sisnan = isnan(x)
#else
   sisnan = SLAISNAN(x,x)

   contains
   logical function SLAISNAN( x, y )
   use LA_CONSTANTS, only: wp=>sp
   real(wp) :: x, y
   SLAISNAN = ( x.ne.y )
   end function SLAISNAN
#endif
   end function SISNAN

   logical function DISNAN( x )
   use LA_CONSTANTS, only: wp=>dp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   DISNAN = ieee_is_nan(x)
#elif USE_ISNAN
   DISNAN = isnan(x)
#else
   DISNAN = DLAISNAN(x,x)

   contains
   logical function DLAISNAN( x, y )
   use LA_CONSTANTS, only: wp=>dp
   real(wp) :: x, y
   DLAISNAN = ( x.ne.y )
   end function DLAISNAN
#endif
   end function DISNAN

end module LA_XISNAN
!> \brief \b DZNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE COMPLEX X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DZNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DZNRM2 := sqrt( x**H*x )
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
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>         complex vector with N elements
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup single_blas_level1
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
function DZNRM2( n, x, incx ) 
   integer, parameter :: wp = kind(1.d0)
   real(wp) :: DZNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxN = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   DZNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(real(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ax = abs(aimag(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   DZNRM2 = scl*sqrt( sumsq )
   return
end function
!> \brief \b ZLARTG generates a plane rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARTG( F, G, C, S, R )
!
!       .. Scalar Arguments ..
!       REAL(wp)           C
!       COMPLEX(wp)        F, G, R, S
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARTG generates a plane rotation so that
!>
!>    [  C         S  ] . [ F ]  =  [ R ]
!>    [ -conjg(S)  C  ]   [ G ]     [ 0 ]
!>
!> where C is real and C**2 + |S|**2 = 1.
!>
!> The mathematical formulas used for C and S are
!>
!>    sgn(x) = {  x / |x|,   x != 0
!>             {  1,         x = 0
!>
!>    R = sgn(F) * sqrt(|F|**2 + |G|**2)
!>
!>    C = |F| / sqrt(|F|**2 + |G|**2)
!>
!>    S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
!>
!> When F and G are real, the formulas simplify to C = F/R and
!> S = G/R, and the returned values of C, S, and R should be
!> identical to those returned by DLARTG.
!>
!> The algorithm used to compute these quantities incorporates scaling
!> to avoid overflow or underflow in computing the square root of the
!> sum of squares.
!>
!> This is a faster version of the BLAS1 routine ZROTG, except for
!> the following differences:
!>    F and G are unchanged on return.
!>    If G=0, then C=1 and S=0.
!>    If F=0, then C=0 and S is chosen so that R is real.
!>
!> Below, wp=>dp stands for double precision from LA_CONSTANTS module.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is COMPLEX(wp)
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is COMPLEX(wp)
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL(wp)
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX(wp)
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX(wp)
!>          The nonzero component of the rotated vector.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
subroutine ZLARTG( f, g, c, s, r )
   use LA_CONSTANTS, &
   only: wp=>dp, zero=>dzero, one=>done, two=>dtwo, czero=>zzero, &
         rtmin=>drtmin, rtmax=>drtmax, safmin=>dsafmin, safmax=>dsafmax
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp)           c
   complex(wp)        f, g, r, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
   complex(wp) :: fs, gs, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
!  ..
!  .. Executable Statements ..
!
   if( g == czero ) then
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         g2 = ABSSQ( g )
         d = sqrt( g2 )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         d = sqrt( g2 )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         f2 = ABSSQ( f )
         g2 = ABSSQ( g )
         h2 = f2 + g2
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = f2*p
         s = conjg( g )*( f*p )
         r = f*( h2*p )
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         if( f1*uu < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            vv = one / v
            w = v * uu
            fs = f*vv
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            fs = f*uu
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = ( f2*p )*w
         s = conjg( gs )*( fs*p )
         r = ( fs*( h2*p ) )*u
      end if
   end if
   return
end subroutine
!> \brief \b ZLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlassq.f90">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlassq.f90">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlassq.f90">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       DOUBLE COMPLEX     X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASSQ  returns the values  scl  and  smsq  such that
!>
!>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!> assumed to be non-negative.
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ and
!> scl and smsq are overwritten on SCALE and SUMSQ respectively.
!>
!> If scale * sqrt( sumsq ) > tbig then
!>    we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
!> and if 0 < scale * sqrt( sumsq ) < tsml then
!>    we require:   scale <= sqrt( HUGE ) / ssml       on entry,
!> where
!>    tbig -- upper threshold for values whose square is representable;
!>    sbig -- scaling constant for big numbers; \see la_constants.f90
!>    tsml -- lower threshold for values whose square is representable;
!>    ssml -- scaling constant for small numbers; \see la_constants.f90
!> and
!>    TINY*EPS -- tiniest representable number;
!>    HUGE     -- biggest representable number.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector x.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE COMPLEX array, dimension (1+(N-1)*abs(INCX))
!>          The vector for which a scaled sum of squares is computed.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector x.
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with  scl , the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is DOUBLE PRECISION
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!>          squares from which  scl  has been factored out.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!> Nick Papior, Technical University of Denmark, DK
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
subroutine ZLASSQ( n, x, incx, scl, sumsq )
   use LA_CONSTANTS, &
      only: wp=>dp, zero=>dzero, one=>done, &
            sbig=>dsbig, ssml=>dssml, tbig=>dtbig, tsml=>dtsml
   use LA_XISNAN
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scl, sumsq
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, ymax, ymin
!  ..
!
!  Quick return if possible
!
   if( LA_ISNAN(scl) .or. LA_ISNAN(sumsq) ) return
   if( sumsq == zero ) scl = one
   if( scl == zero ) then
      scl = one
      sumsq = zero
   end if
   if (n <= 0) then
      return
   end if
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(real(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ax = abs(aimag(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Put the existing sum of squares into one of the accumulators
!
   if( sumsq > zero ) then
      ax = scl*sqrt( sumsq )
      if (ax > tbig) then
!        We assume scl >= sqrt( TINY*EPS ) / sbig
         abig = abig + (scl*sbig)**2 * sumsq
      else if (ax < tsml) then
!        We assume scl <= sqrt( HUGE ) / ssml
         if (notbig) asml = asml + (scl*ssml)**2 * sumsq
      else
         amed = amed + scl**2 * sumsq
      end if
   end if
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if (amed > zero .or. LA_ISNAN(amed)) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if (amed > zero .or. LA_ISNAN(amed)) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range or zero
!
      scl = one
      sumsq = amed
   end if
   return
end subroutine
