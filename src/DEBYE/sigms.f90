!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sigms.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!     program sigms.f
!
!     calculates debye-waller factors for each multiple
!     scattering path using Debye-Model correlations
!
!     files:  input  pathd_all.dat  multiple scattering path data
!             output fort.3  sig**2 vs path
!                    fort.2  long output
!
!     version 1  (29 july 91)
!
!     coded by j. rehr
!     path data from s. zabinsky
!
!     modified to use pdata.inp, Dec 1991, siz
!     Subroutine version, Dec 1991, siz
!
!---------------------------------------------------------------------

      subroutine sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
!               tk temperature in degrees K
!               thetad debye temp in degrees K
!               rs=wigner seitz or norman radius in bohr, averaged
!                  over entire problem
!                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
!                  (sum is over all atoms in the problem)
!               nlegx used in dimensions of rat and iz
!               nleg nlegs in path
!               rat positions of each atom in path
!               iz atomic number of each atom in path
!               NB Units of distance in this routine
!                  are angstroms, including sig**2.  rs is in bohr.
!               sig2 is output debye waller factor
      use constants

      implicit double precision (a-h,o-z)

!     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)
!#mn
       external dist

!      parameters
!               x = k_d*R   (distance parameter)
!               R distance in angstroms
!               y = hbar omegad/kT = thetad/t
!               thetad debye temp in degrees K
!               tk temperature in degrees K
!               k_d = (6*pi**2 N/V) = debye wave number
!               N/V=1/(4pi/3rs**3)
!               rs=wigner seitz or norman radius in bohr
!               ami, amj masses at sites i and j in amu
!               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

!     Note:  There are nleg atoms including the central atom
!            index 0 and index nleg both refer to central atom,
!            which makes special code unnecessary later.

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

!        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+        &
     &          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+        &
     &          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

!        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

!     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn(rij,cij,thetad,tk,iz1,iz2,rsavg)
!     subroutine calculates correlation function
!     c(ri,rj)=<xi xj> in the Debye approximation
!
!             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
!              (hbar/2w_k)*coth(beta hbar w_k/2)
!             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
!
!      parameters
!               x = k_d*R   (distance parameter)
!               R distance in angstroms
!               y = hbar omegad/kT = thetad/t
!               thetad debye temp in degrees K
!               tk temperature in degrees K
!               k_d = (6*pi**2 N/V) = debye wave number
!               N/V=1/(4pi/3rs**3)
!               rs=wigner seitz or norman radius in bohr
!               ami, amj masses at sites i and j in amu
!               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
!
!      solution by numerical integration
!
      use constants
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

!     con=hbar**2/kB*amu)*10**20   in ang**2 units
!     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
!     kB = 1.380 6581 d-23
      parameter (con = 48.508459393094)
!#mn
       external atwtd

!     external fn
!     rij=2.55
!     tk=295
!     thetad=315
!     ami=amj=63.55 at wt for Cu
!     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
!     thetad in degrees K, t temperature in degrees K
!     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
!     call numerical integration
      call bingrt (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
!     fn=(sin(wx)/x)*coth(wy/2)
!     change code to allow t=0 without bombing
!     fn=2/y
      fn=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
      emwy=0.
      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
      emwy=exp(-w/yinv)
      fn=fac*(1+emwy)/(1-emwy)
      return
      end
!-----------------------------------------------
      subroutine bingrt (b, eps, n)
!     subroutine calculates integrals between [0,1]
!      b = int_0^1 f(z) dz
!     by trapezoidal rule and binary refinement
!     (romberg integration)
!     coded by j rehr (10 Feb 92)
!     see, e.g., numerical recipes for discussion
!     and a much fancier version
!-----------------------------------------------
!     del=dz  itn=2**n tol=1.e-5
!     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
      character*512 slog
!     external fn
!     error is approximately 2**(-2n) ~ 10**(-.6n)
!     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn(zero)+fn(one))/2
      bo=bn
 10   continue
!     nth iteration
!     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn(zi)
!     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
!     cancel leading error terms b=[4b-bn]/3
!     note: this is the first term in the
!     neville table - remaining errors were
!     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(slog,50) n,itn, b,eps
      call wlog(slog)
 50   format(' not converged, n,itn,b,eps=',                            &
     &  2i4,2e14.6)
      return
 60   continue
      return
      end
