      REAL*8 FUNCTION BRENT(AX,BX,CX,TOL,XMIN) !KJ removed "f" from argument list
	use fitting, xx=>fitx,yy=>fity,n=>fitmeshsize,m=>fitorder  !KJ
	implicit none !KJ
      integer itmax
	real*8 ax,bx,cx,tol,xmin,cgold,zeps
!Given a function F, and a bracketing triplet of abscissas
!AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
!than both F(AX) and F(CX)), this routine isolates the minimum 
!to a fractional precision of about TOL using Brent's method.
!The abscissa of the minimum is returned in XMIN, and the minimum
!function value is returned as BRENT, the returned function value.
      PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
!Maximum allowed number of iterations; golden ration; and a small
!number which protects against trying to achieve fractional accuracy
!for a minimum that happens to be exactly zero.
      integer iter
	  real*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      call terp (xx,yy,n,m,x,fx) !KJ fx=f(x)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX	                                !main loop
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3  !Test for done here
        IF(ABS(E).GT.TOL1) THEN     !Construct a trial parabolic fit
          R=(X-W)*(FX-FV)
      	Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
      	Q=2.D0*(Q-R)         !bug corrected 07/24/2006 (0.2 instead of 2.D0)
      	IF(Q.GT.0)  P=-P
      	Q=ABS(Q)
      	ETEMP=E
          E=D
      	IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  P.GE.Q*(B-X))  GOTO 1
!   The above conditions determine the acceptability of the 
!   parabolic fit. Here it is o.k.:
          D=P/Q
      	U=X+D
      	IF(U-A.LT.TOL2.OR.B-U.LT.TOL2)  D=SIGN(TOL1,XM-X)
      	GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        call terp (xx,yy,n,m,u,fu) !KJ fu=f(u)  !This the one function evaluation per iteration
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN      
      	  A=X
          ELSE
	      B=X
	    ENDIF
	    V=W
	    FV=FW
	    W=X
	    FW=FX
      	X =U      
      	FX=FU
        ELSE
          IF(U.LT.X) THEN
	       A=U      
      	ELSE
	       B=U
      	ENDIF
	    IF(FU.LE.FW.OR.W.EQ.X) THEN
	       V=W
	       FV=FW
	       W=U
	       FW=FU
          ELSE IF(FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
	       V=U
	       FV=FU
          ENDIF
        ENDIF
11       CONTINUE
      Pause ' Brent exceed maximum iterations.'
3       XMIN=X   !exit section
        BRENT=FX
        RETURN
        END
  
!end of file Brent.f90
