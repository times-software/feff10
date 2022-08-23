SUBROUTINE sor(a,b,c,d,e,f,u,jmax,rjac)
INTEGER jmax,MAXITS
DOUBLE PRECISION rjac,a(jmax,jmax),b(jmax,jmax),
* c(jmax,jmax),d(jmax,jmax),e(jmax,jmax),
* f(jmax,jmax),u(jmax,jmax),EPS
PARAMETER (MAXITS=1000,EPS=1.d-5)
Successive overrelaxation solution of equation (19.5.25) with Chebyshev acceleration. a,
b, c, d, e, and f are input as the coefficients of the equation, each dimensioned to the
grid size JMAX × JMAX. u is input as the initial guess to the solution, usually zero, and
returns with the final value. rjac is input as the spectral radius of the Jacobi iteration,
or an estimate of it
INTEGER ipass,j,jsw,l,lsw,n
DOUBLE PRECISION anorm,anormf,
* omega,resid Double precision is a good idea for JMAX bigger than about 25.
anormf=0.d0 Compute initial norm of residual and terminate iteration when
norm has been reduced by a factor EPS.do 12 j=2,jmax-1
do 11 l=2,jmax-1
anormf=anormf+abs(f(j,l)) Assumes initial u is zero.
enddo 11
enddo 12
omega=1.d0
do 16 n=1,MAXITS
anorm=0.d0
jsw=1
do 15 ipass=1,2 Odd-even ordering.
lsw=jsw
do 14 j=2,jmax-1
do 13 l=lsw+1,jmax-1,2
resid=a(j,l)*u(j+1,l)+b(j,l)*u(j-1,l)+
* c(j,l)*u(j,l+1)+d(j,l)*u(j,l-1)+
* e(j,l)*u(j,l)-f(j,l)
anorm=anorm+abs(resid)
u(j,l)=u(j,l)-omega*resid/e(j,l)
enddo 13
lsw=3-lsw
enddo 14
jsw=3-jsw
if(n.eq.1.and.ipass.eq.1) then
omega=1.d0/(1.d0-.5d0*rjac**2)
else
omega=1.d0/(1.d0-.25d0*rjac**2*omega)
endif
enddo 15
if(anorm.lt.EPS*anormf)return
enddo 16
pause ’MAXITS exceeded in sor’
END
