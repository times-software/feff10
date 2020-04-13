!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mksat.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function xmkesat(w)
! Find the extrinsic satellite spectral function.  This is the extrinsic 
! spectral function with the quasiparticle pole subtracted off and the 
! quasiparticle broadening removed.
! input: w - energy (omega)
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       omp - plasma frequency omega_p
!       se - real part of the on shell self energy
!       width - absolute value of the imaginargy part of the 
!            on shell self energy plus the core hole broadening
!       se2 - real part of the self energy at energy w
!       xise - imaginary part of the self energy at energy w
!       z1 - real part of the renormalization constant
!       z1i - imaginary part of the renormalization constant
      implicit none
      integer it1,it2
      double precision w,etot,emain,etothi,etotlo,z1m
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      z1m=sqrt(z1**2+z1i**2)
!      etot=-z1*(width-xise)+z1i*(w+se-se2)
      etot=-(width-xise)
      etot=etot/((w+se-se2)**2+(width-xise)**2)
!      it1=int(xmkesat)
!      it2=int(2*xmkesat)
!      if(it1.eq.it2.and.it1.gt.5) then
!        etot=-sqrt(z1**2+z1i**2)*(width-xise)
!        etothi=etot/((w+omp*1.d-3+se-se2)**2+(width-xise)**2)
!        etotlo=etot/((w-omp*1.d-3+se-se2)**2+(width-xise)**2)
!        etot=(etothi+etotlo)/2.d0
!      endif
!      emain=0.d0
      emain=-z1i/(w*pi*z1m)
      emain=emain*exp(-(w/(2*omp))**2)
      xmkesat=etot/(pi*z1m)-emain
!      it1=int(xmkesat)
!      it2=int(2*xmkesat)
!      if(it1.eq.it2.and.it1.gt.5) then
!        write(6,*) 'nan error: xmkesat'
!        write(6,*) z1,z1i,width,xise,w,se,se2,etot,xmkesat
!      endif
      return
500   format(1x,5(e12.5,1x))
      end

      double precision function xmkgwext(w)
! Find the extrinsic satellite, with full broadening and all quasiparticle
! contributions.
! input: w - energy (omega)
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       se - real part of the on shell self energy
!       se2 - real part of the self energy at energy w
!       xise - imaginary part of the self energy at energy w
      implicit none
      integer it1,it2
      double precision w,etot,emain
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      etot=xise/(pi*((w+se-se2)**2+xise**2))
      emain=0.d0
      xmkgwext=etot-emain
      it1=int(xmkgwext)
      it2=int(2*xmkgwext)
      if(it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'nan error: xmkgwext'
        write(6,*) z1,z1i,width,xise,w,se,se2,etot,xmkgwext
      endif
      return
500   format(1x,5(e12.5,1x))
      end

      double precision function xintxsat(q)
! the integrand of the interference spectral function
! input: q - momentum to be integrated over
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       omp - plasma frequency omega_p
!       ek - bare photoelectron kinetic energy = pk**2/2
!       brd - global broadening parameter to stabilize logarithms
!       ac2 - additional accuracy parameter
!       wp2 - omega prime, an additional energy variable
      implicit none
      integer numcal,maxns,nsing
      double precision wwq,q,xk,vpp2,wdisp,xfact,xloren,                &
     &       eps,ac2,wp2,dw1,abr,rlr,xsing,error,tol
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision se,ce,width,z1,z1i,se2,xise
      common /energies/ se,ce,width,z1,z1i,se2,xise
      common /ff/ ac2,wp2
      external vpp2,wdisp
      wwq=wdisp(q)
      tol=2.d-1*omp
      if (ek-wp2.ge.0.d0) then
        xk=sqrt(2.d0*(ek-wp2))
        xfact=log(((wwq-q**2/2+xk*q)**2+tol**2)/                        &
     &            ((wwq-q**2/2-xk*q)**2+tol**2))/2.d0
        xloren=ac2/(pi*((wp2-wwq)**2+ac2**2))
        xintxsat=q*vpp2(q)*xloren*xfact/(wwq*xk)
      else
        xk=sqrt(-2.d0*(ek-wp2))
        xfact=datan(xk*q/(wwq-q**2/2))
        xloren=ac2/(pi*((wp2-wwq)**2+ac2**2))
        xintxsat=q*vpp2(q)*xloren*xfact/(wwq*xk)
      endif
      return
      end

      double precision function xintisat(q)
! the integrand of the intrinsic spectral function
! input: q - momentum to be integrated over
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ac2 - additional accuracy parameter
!       wp2 - omega prime, an additional energy variable
      implicit none
      integer numcal,maxns,nsing
      double precision wwq,q,xk,vpp2,wdisp,xfact,xloren,                &
     &       ac2,wp2,dw1,abr,rlr,xsing,error
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external vpp2,wdisp
      wwq=wdisp(q)
      xloren=ac2/(((wp2-wwq)**2+ac2**2)*pi)
      xintisat=q**2*vpp2(q)*xloren/wwq**2
      return
      end

      double precision function xmkxsat(w,width)
! generates the interference spectral function
! input: w - energy (omega)
!        width - additional broadening
! input from common blocks
!       omp - plasma frequency omega_p
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp2 - omega prime, an additional energy variable
      implicit none
      integer nsing,numcal,maxns,i
      double precision xintxsat,grater,abr,rlr,xsing(20),error,             &   !KJ bugfix 11-2011 xsing->xsing(20)
     &       qk,wp,acc,pi,ef,xmu,qf,omp,ompl,wt,w,pk,width,             &
     &       ekp,ek,ac2,wp2,q2,qwidth,qmin,qmax,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external xintxsat,grater
      wp2=w
      ac2=width
      rlr=acc
      abr=omp*acc
      nsing=0
      q2=sqrt(max(2*(w-ompl),ac2))
      qwidth=10.d0*ac2/q2
      qmin=max(0.d0,q2-qwidth)
      qmax=q2+qwidth
!      do i=1,1000
!      enddo
!      write(6,*) 'xmkxsat'
      xmkxsat=grater(xintxsat,qmin,q2,                                  &
     &             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkxsat=xmkxsat+grater(xintxsat,q2,qmax,                          &
     &             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkxsat=xmkxsat/(2.d0*pi)**2
      return
      end

      double precision function xmkisat(w,width)
! generates the intrinsic spectral function
! input: w - energy (omega)
!        width - additional broadening
! input from common blocks
!       omp - plasma frequency omega_p
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp2 - omega prime, an additional energy variable
      implicit none
      integer nsing,numcal,maxns
      double precision xintisat,grater,abr,rlr,xsing(20),error,             &    !KJ bugfix 11-2011 xsing->xsing(20)
     &       qk,wp,acc,pi,ef,xmu,qf,omp,ompl,wt,z,w,pk,width,           &
     &       ekp,ek,ac2,wp2,q2,qmin,qmax,qwidth,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /ff/ ac2,wp2
      external xintisat,grater
      wp2=w
      ac2=width
      rlr=acc
      abr=omp*acc
      nsing=0
      if ((w-ompl).gt.ac2) then
        q2=dsqrt(2.d0*(w-ompl))
      else
        q2=dsqrt(2.d0*ac2)
      endif
      qwidth=10.d0*min(q2,ac2/q2)
      qmin=max(0.d0,q2-qwidth)
      qmax=q2+qwidth
      xmkisat=grater(xintisat,0.d0,q2,                                  &
     &             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkisat=xmkisat+grater(xintisat,q2,qmax,                          &
     &             abr,rlr,nsing,xsing,error,numcal,maxns)
      xmkisat=xmkisat/(2.d0*pi**2)
      return
      end

!      double precision function xmkisat(w,width)
!*     generates the intrinsic spectral function
!      implicit none
!      double precision w,wdisp,qdisp,dwdq,q,vpp2,width
!      double precision pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
!      common /convsf/ pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
!      external wdisp,dwdq,qdisp,vpp2
!      xmkisat=0.d0
!      if (w.gt.omp) then
!        q=qdisp(w)
!        xmkisat=q**2*vpp2(q)/(2.d0*pi**2*w**2*dwdq(q))
!      endif
!      return
!      end
!
!      double precision function xmkxsat(w)
!*     generates the interference spectral function
!      implicit none
!      double precision w,q,wdisp,qdisp,dwdq,vpp2,xk,xfact,
!     2                 eta
!      complex*16 cfac,coni
!      parameter (coni=(0,1))
!      double precision pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
!      common /convsf/ pi,ef,xmu,qf,omp,ekp,ek,pk,acc,brd,adisp
!      external wdisp,dwdq,qdisp,vpp2
!      xmkxsat=0.d0
!      eta=1.d-2*omp
!      if (w.gt.omp.and.ek.gt.w) then
!        q=qdisp(w)
!        xk=dsqrt(2.d0*(ek-w))
!        xfact=log(((w-q**2/2+xk*q)**2+(omp*acc)**2)/
!     2            ((w-q**2/2-xk*q)**2+(omp*acc)**2))/2.d0
!*        cfac=log((w-q**2/2+xk*q-coni*eta)/(w-q**2/2-xk*q-coni*eta))
!*        xfact=dble(cfac)
!        xmkxsat=q*vpp2(q)*xfact/((2.d0*pi)**2*xk*w*dwdq(q))
!      endif
!      return
!      end

      double precision function xmkak(w)
! This function returns the interference contribution to the 
! quasiparticle.
! input: w - energy (omega)
! input from common blocks
!       omp - plasma frequency omega_p
!       ek - bare photoelectron kinetic energy = pk**2/2
!       acc - global accuracy parameter 
!       wmax - highest allowed energy
! common block control of subprograms
!       pkq - momentum variable used in the integrand but kept 
!             fixed throughout the integration
      implicit none
      integer i,j,jj,nsing,numcal,maxns,it1,it2
      double precision w,qmax,xintak,                                   &
     &                 abr,rlr,error,xsing(20),grater   !KJ bugfix 11-2011 xsing->xsing(20)
      integer nnpts
      double precision pkq
      common /funct2/ pkq
      double precision wmin,wmax,wmin1,wmax1
      common /limits/ wmin,wmax,wmin1,wmax1
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external xintak
      if (w.gt.0.d0) then
        rlr=acc
        abr=dsqrt(omp)*acc
        qmax=dsqrt(2*wmax)
        pkq=dsqrt(2*ek)
        nsing=0
!        write(6,*) 'xmkak'
        xmkak=grater(xintak,abr,qmax,abr,rlr,nsing,xsing,               &
     &                 error,numcal,maxns)
      else
        xmkak=0.d0
      endif
      return
      end

      double precision function xintak(q)
! Integrand for the function xmkak.
! input: q - momentum to be integrated over.
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       omp - plasma frequency omega_p
!       pkq - momentum variable used in the integrand but kept 
!             fixed throughout the integration
      implicit none
      double precision wq,q,pkq,xlog,wdisp,vpp2,eps
      common /funct2/ pkq
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external vpp2,wdisp
      wq=wdisp(q)
      eps=1.d-1
      xlog=(wq+q**2/2.d0+pkq*q)**2+(ompl*eps)**2
      xlog=xlog/((wq+q**2/2.d0-pkq*q)**2+(ompl*eps)**2)
      xlog=log(xlog)/2.d0
      xintak=q*vpp2(q)*xlog/(wq*pkq*4.d0*pi**2)
      return
      end
