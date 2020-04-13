!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: senergies.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine renergies(w,rbeta)
! Calculates the real part of the photelectron self energy due to 
! pole ipl in the inverse dielectric function epsilon^{-1}.
! input: w - energy (omega)
! output: rbeta - the real part of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - width of pole in epsilon^(-1)
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,rbeta,wmax,wmin,           &  !KJ ,beta,grater,exchange
     &                 abr,rlr,xsing(20),error,qmax  !KJ bugfix 11-2011 xsing->xsing(20)
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      !KJ double precision rseint1,rseint2,rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      double precision,external :: beta,grater,exchange,rseint1,rseint2,rseint3
      integer ijkwrite
      common /morewrite/ ijkwrite
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
!     rlr for plasmon pole should be 10^-7 to eliminate numerical 
!     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        rbeta=grater(rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        rbeta=rbeta+grater(rseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,      &
     &             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,pk-qf,pk+qf,abr,rlr,nsing,xsing,     &
     &             error,numcal,maxns)
      elseif (pk.lt.qf) then
        rbeta=grater(rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,qf-pk,pk+qf,abr,rlr,nsing,xsing,     &
     &             error,numcal,maxns)
        if (lowq.ne.0) rbeta=rbeta+grater(rseint3,0.d0,                 &
     &             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        rbeta=grater(rseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,          &
     &             error,numcal,maxns)
        rbeta=rbeta+grater(rseint2,0.d0,2.d0*qf,abr,rlr,nsing,xsing,    &
     &             error,numcal,maxns)
      endif
      rbeta=-rbeta*omp**2/(2.d0*pi*pk)
!      rbeta=rbeta+exchange(pk)
      return
      end

      subroutine brsigma(w,rbeta,xibeta)
! Calculates the broadened photelectron self energy due to 
! pole ipl in the inverse dielectric function epsilon^{-1}.
! input: w - energy (omega)
! output: rbeta - the real part of the self energy
!         xibeta - the imaginary part of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - linewidth of pole in epsilon^{-1}
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns,nq,nqq
      parameter (nqq=5)
      double precision w,rbeta,xibeta,grater,qdisp,wmax,wmin,           &
     &                 abr,rlr,xsing(20),error
      double precision qmax,qlimh,qliml,qh,q0,q1,q2,q3,qsing(nqq)
      double precision sig1r,sig2r,sig3r,sig4r,sig5r,                   &
     &                 sig6r,sig7r,sig8r,sig9r,sig10r,                  &
     &                 sig1i,sig2i,sig3i,sig4i,sig5i,                   &
     &                 sig6i,sig7i,sig8i,sig9i,sig10i
      double precision fqlogr1,fqlogr2,fqlogr3,fqlogr4,                 &
     &                 fqlogi1,fqlogi2,fqlogi3,fqlogi4,                 &
     &                 fqatnr1,fqatnr2,fqatnr3,fqatnr4,                 &
     &                 fqatni1,fqatni2,fqatni3,fqatni4
      double complex cfact,csigma
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer iwrite,jj
      common /flag/ iwrite,jj
      integer lowq
      common /belowqf/ lowq
      external grater,qdisp,                                            &
     &         fqlogr1,fqlogr2,fqlogr3,fqlogr4,                         &
     &         fqlogi1,fqlogi2,fqlogi3,fqlogi4,                         &
     &         fqatnr1,fqatnr2,fqatnr3,fqatnr4,                         &
     &         fqatni1,fqatni2,fqatni3,fqatni4
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
!     rlr for plasmon pole should be 10^-7 to eliminate numerical 
!     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      qlimh=pk+qf
      qliml=abs(pk-qf)
      qh=qdisp(max(wp-ef,ompl))
      q0=qdisp(max(ef-wp,ompl))
      call qlimits(wp,pk,ompl,adisp,qmax,nq,q1,q2,q3)
      qsing(1)=q0
      qsing(2)=q1
      qsing(3)=q2
      qsing(4)=q3
      qsing(5)=qh
      sig1r=0.d0
      sig1i=0.d0
      sig2r=0.d0
      sig2i=0.d0
      sig3r=0.d0
      sig3i=0.d0
      sig4r=0.d0
      sig4i=0.d0
      sig5r=0.d0
      sig5i=0.d0
      sig6r=0.d0
      sig6i=0.d0
      sig7r=0.d0
      sig7i=0.d0
      sig8r=0.d0
      sig8i=0.d0
      sig9r=0.d0
      sig9i=0.d0
      sig10r=0.d0
      sig10i=0.d0
      call findsing(qlimh,qmax,nqq,qsing,nsing,xsing)
      sig1r=grater(fqlogr1,qlimh,qmax,abr,rlr,nsing,xsing,              &
     &           error,numcal,maxns)
      sig1i=grater(fqlogi1,qlimh,qmax,abr,rlr,nsing,xsing,              &
     &           error,numcal,maxns)
      sig2r=grater(fqatnr1,qlimh,qmax,abr,rlr,nsing,xsing,              &
     &           error,numcal,maxns)
      sig2i=grater(fqatni1,qlimh,qmax,abr,rlr,nsing,xsing,              &
     &           error,numcal,maxns)
      call findsing(qliml,qlimh,nqq,qsing,nsing,xsing)
      sig3r=grater(fqlogr2,qliml,qlimh,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      sig3i=grater(fqlogi2,qliml,qlimh,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      sig4r=grater(fqatnr2,qliml,qlimh,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      sig4i=grater(fqatni2,qliml,qlimh,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      if (lowq.ne.0) then
        sig7r=grater(fqlogr3,qliml,qlimh,abr,rlr,nsing,xsing,           &
     &           error,numcal,maxns)
        sig7i=grater(fqlogi3,qliml,qlimh,abr,rlr,nsing,xsing,           &
     &           error,numcal,maxns)
        sig8r=grater(fqatnr3,qliml,qlimh,abr,rlr,nsing,xsing,           &
     &           error,numcal,maxns)
        sig8i=grater(fqatni3,qliml,qlimh,abr,rlr,nsing,xsing,           &
     &           error,numcal,maxns)
      endif
      call findsing(0.d0,qliml,nqq,qsing,nsing,xsing)
      if (pk.gt.qf) then
        sig5r=grater(fqlogr1,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        sig5i=grater(fqlogi1,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        sig6r=grater(fqatnr1,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        sig6i=grater(fqatni1,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
      endif
      if (pk.lt.qf.and.lowq.ne.0) then
        sig9r=grater(fqlogr4,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        sig9i=grater(fqlogi4,0.d0,qliml,abr,rlr,nsing,xsing,            &
     &             error,numcal,maxns)
        sig10r=grater(fqatnr4,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        sig10i=grater(fqatni4,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
      endif
      rbeta=(sig1r+sig3r+sig5r+sig7r+sig9r)*omp**2/(4.d0*pi*pk)         &
     &     +(sig2r+sig4r+sig6r+sig8r+sig10r)*omp**2/(2.d0*pi*pk)
      xibeta=(sig1i+sig3i+sig5i+sig7i+sig9i)*omp**2/(4.d0*pi*pk)        &
     &     -(sig2i+sig4i+sig6i+sig8i+sig10i)*omp**2/(2.d0*pi*pk)
      cfact=(1,0)-(brd/ompl)*(0,1)
      csigma=((1,0)*rbeta+(0,1)*xibeta)*cfact
      rbeta=dble(csigma)
      xibeta=dimag(csigma)
!      rbeta=rbeta+exchange(pk)
      return
      end

      double precision function fqlogr1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr1=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi1=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=(ef-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr2=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=(ef-wp+wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp+wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi2=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/((ef-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr3=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/((ef-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi3=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogr4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogr4=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqlogi4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xlog=((pk-q)**2/2.d0-wp-wq)**2+(brd)**2
      xlog=xlog/(((pk+q)**2/2.d0-wp-wq)**2+(brd)**2)
      xlog=log(xlog)
      fqlogi4=xpole*xlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function fqatni1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp-wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatni1=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp-wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatnr1=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp-wq-ef)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatni2=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp-wq-ef)/brd)
      xatan=xatan-atan((wp-wq-(pk+q)**2/2.d0)/brd)
      fqatnr2=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-ef)/brd)
      fqatni3=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-ef)/brd)
      fqatnr3=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatni4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-(pk+q)**2/2.d0)/brd)
      fqatni4=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function fqatnr4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,xatan
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xatan=atan((wp+wq-(pk-q)**2/2.d0)/brd)
      xatan=xatan-atan((wp+wq-(pk+q)**2/2.d0)/brd)
      fqatnr4=xpole*xatan/dsqrt(q**2+omp*acc)
      return
      end

      subroutine dbrsigma(w,drbeta,dibeta)
! Calculates the derivative of the broadened photelectron self energy 
! with respect to frequency w due to 
! pole ipl in the inverse dielectric function epsilon^{-1}.
! input: w - energy (omega)
! output: drbeta - the real part of the self energy derivative
!         dibeta - the imaginary part of the self energy derivative
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - linewidth of pole in epsilon^{-1}
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns,nq,nqq
      parameter (nqq=5)
      double precision w,drbeta,dibeta,wmax,wmin,          &
     &                 abr,rlr,xsing(20),error !KJ ,grater,qdisp
      double precision qmax,qlimh,qliml,qh,q0,q1,q2,q3,qsing(nqq)
      double precision dsig1r,dsig2r,dsig3r,dsig4r,dsig5r,              &
     &                 dsig6r,dsig7r,dsig8r,dsig9r,dsig10r,             &
     &                 dsig1i,dsig2i,dsig3i,dsig4i,dsig5i,              &
     &                 dsig6i,dsig7i,dsig8i,dsig9i,dsig10i
      double precision fqlogr1,fqlogr2,fqlogr3,fqlogr4,                 &
     &                 fqlogi1,fqlogi2,fqlogi3,fqlogi4,                 &
     &                 fqatnr1,fqatnr2,fqatnr3,fqatnr4,                 &
     &                 fqatni1,fqatni2,fqatni3,fqatni4
      double complex cfact,csigma
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer iwrite,jj
      common /flag/ iwrite,jj
      integer lowq
      common /belowqf/ lowq
      double precision,external :: grater,qdisp,                                            & !KJ added "double precision"
     &         dqlogr1,dqlogr2,dqlogr3,dqlogr4,                         &
     &         dqlogi1,dqlogi2,dqlogi3,dqlogi4,                         &
     &         dqatnr1,dqatnr2,dqatnr3,dqatnr4,                         &
     &         dqatni1,dqatni2,dqatni3,dqatni4
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)+pk+qf
!     rlr for plasmon pole should be 10^-7 to eliminate numerical 
!     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      qlimh=pk+qf
      qliml=abs(pk-qf)
      qh=qdisp(max(wp-ef,ompl))
      q0=qdisp(max(ef-wp,ompl))
      call qlimits(wp,pk,ompl,adisp,qmax,nq,q1,q2,q3)
      qsing(1)=q0
      qsing(2)=q1
      qsing(3)=q2
      qsing(4)=q3
      qsing(5)=qh
      dsig1r=0.d0
      dsig1i=0.d0
      dsig2r=0.d0
      dsig2i=0.d0
      dsig3r=0.d0
      dsig3i=0.d0
      dsig4r=0.d0
      dsig4i=0.d0
      dsig5r=0.d0
      dsig5i=0.d0
      dsig6r=0.d0
      dsig6i=0.d0
      dsig7r=0.d0
      dsig7i=0.d0
      dsig8r=0.d0
      dsig8i=0.d0
      dsig9r=0.d0
      dsig9i=0.d0
      dsig10r=0.d0
      dsig10i=0.d0
      call findsing(qlimh,qmax,nqq,qsing,nsing,xsing)
      dsig1r=grater(dqlogr1,qlimh,qmax,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      dsig1i=grater(dqlogi1,qlimh,qmax,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      dsig2r=grater(dqatnr1,qlimh,qmax,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      dsig2i=grater(dqatni1,qlimh,qmax,abr,rlr,nsing,xsing,             &
     &           error,numcal,maxns)
      call findsing(qliml,qlimh,nqq,qsing,nsing,xsing)
      dsig3r=grater(dqlogr2,qliml,qlimh,abr,rlr,nsing,xsing,            &
     &           error,numcal,maxns)
      dsig3i=grater(dqlogi2,qliml,qlimh,abr,rlr,nsing,xsing,            &
     &           error,numcal,maxns)
      dsig4r=grater(dqatnr2,qliml,qlimh,abr,rlr,nsing,xsing,            &
     &           error,numcal,maxns)
      dsig4i=grater(dqatni2,qliml,qlimh,abr,rlr,nsing,xsing,            &
     &           error,numcal,maxns)
      if (lowq.ne.0) then
        dsig7r=grater(dqlogr3,qliml,qlimh,abr,rlr,nsing,xsing,          &
     &           error,numcal,maxns)
        dsig7i=grater(dqlogi3,qliml,qlimh,abr,rlr,nsing,xsing,          &
     &           error,numcal,maxns)
        dsig8r=grater(dqatnr3,qliml,qlimh,abr,rlr,nsing,xsing,          &
     &           error,numcal,maxns)
        dsig8i=grater(dqatni3,qliml,qlimh,abr,rlr,nsing,xsing,          &
     &           error,numcal,maxns)
      endif
      call findsing(0.d0,qliml,nqq,qsing,nsing,xsing)
      if (pk.gt.qf) then
        dsig5r=grater(dqlogr1,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        dsig5i=grater(dqlogi1,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        dsig6r=grater(dqatnr1,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        dsig6i=grater(dqatni1,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
      endif
      if (pk.lt.qf.and.lowq.ne.0) then
        dsig9r=grater(dqlogr4,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        dsig9i=grater(dqlogi4,0.d0,qliml,abr,rlr,nsing,xsing,           &
     &             error,numcal,maxns)
        dsig10r=grater(dqatnr4,0.d0,qliml,abr,rlr,nsing,xsing,          &
     &             error,numcal,maxns)
        dsig10i=grater(dqatni4,0.d0,qliml,abr,rlr,nsing,xsing,          &
     &             error,numcal,maxns)
      endif
      drbeta=(dsig1r+dsig3r+dsig5r+dsig7r+dsig9r)*omp**2/(2.d0*pi*pk)   &
     &     +(dsig2r+dsig4r+dsig6r+dsig8r+dsig10r)*omp**2/(2.d0*pi*pk)
      dibeta=(dsig1i+dsig3i+dsig5i+dsig7i+dsig9i)*omp**2/(2.d0*pi*pk)   &
     &     -(dsig2i+dsig4i+dsig6i+dsig8i+dsig10i)*omp**2/(2.d0*pi*pk)
      cfact=(1,0)-(brd/ompl)*(0,1)
      csigma=((1,0)*drbeta+(0,1)*dibeta)*cfact
      drbeta=dble(csigma)
      dibeta=dimag(csigma)
!      rbeta=rbeta+exchange(pk)
      return
      end

      double precision function dqlogr1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr1=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi1=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogr2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr2=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi2=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogr3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr3=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi3=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end


      double precision function dqlogr4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogr4=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqlogi4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxlog,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxlog=xw1/((xw1)**2+(brd)**2)-xw2/((xw2)**2+(brd)**2)
      dqlogi4=xpole*dxlog/dsqrt(q**2+ompl*acc)
      return
      end

      double precision function dqatni1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni1=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr1(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-(pk-q)**2/2.d0
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr1=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni2=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr2(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp-wq-ef
      xw2=wp-wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr2=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni3=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr3(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-ef
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr3=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatni4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=wq/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatni4=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      double precision function dqatnr4(q)
! Integrand for one of the intergals in calculating the 
! self energy in subroutine brsigma.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
!      integer it1,it2,i
      double precision q,wq,xpole,dxatan,xw1,xw2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xpole=brd/(wq**2+brd**2)
      xw1=wp+wq-(pk-q)**2/2.d0
      xw2=wp+wq-(pk+q)**2/2.d0
      dxatan=brd/((xw1)**2+(brd)**2)-brd/((xw2)**2+(brd)**2)
      dqatnr4=xpole*dxatan/dsqrt(q**2+omp*acc)
      return
      end

      subroutine findsing(ql,qh,nqq,qsing,nsing,xsing)
! finds which values in the array qsing are between ql and qh, and puts
! them in the array xsing.  Returns number of elements in xsing as nsing
      implicit none
      integer i,j,k,nqq,nsing,it
      double precision ql,qh,qsing(nqq),xsing(20),store
      integer iwrite,jj
      common /flag/ iwrite,jj
      nsing=0
      do i=1,nqq
        if (qsing(i).gt.ql.and.qsing(i).lt.qh) then
          nsing=nsing+1
          xsing(nsing)=qsing(i)
        elseif (qsing(i).lt.ql.and.qsing(i).gt.qh) then
          nsing=nsing+1
          xsing(nsing)=qsing(i)
        endif
      enddo
      if (nsing.le.1) goto 40
      j=2
 20   k=j
 30   continue
      if (xsing(k-1).gt.xsing(k)) then
        store=xsing(k-1)
        xsing(k-1)=xsing(k)
        xsing(k)=store
        k=k-1
        if (k.gt.1) goto 30
      endif
      j=j+1
      if (j.le.nsing) goto 20
 40   return
      end

      double precision function exchange(qk)
! compute the HF exchange potential for the free electron gas
! input: qk - photoelectron momentum
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       qf - Fermi momentum
      implicit none
      integer i,ii,j,jj
      double precision qk
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      if (qk.eq.qf) then
        exchange=-qf/pi
      else
        exchange=-(1.d0/pi)*(qf+                                        &
     &  ((qf**2-qk**2)/(2.d0*qk))*log(dabs((qk+qf)/(qk-qf))))
      endif
      return
      end

      subroutine xienergies(w,xibeta)
! Calculates the imaginary part of the photelectron self energy.
! input: w - energy (omega)
! output: xibeta - the imaginary part of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,xibeta,beta,grater,wmax,wmin,         &
     &                 abr,rlr,xsing,error
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,grater
      xibeta=-beta(w)*pi
      return
      end

      subroutine drenergies(w,drbeta)
! Calculates the real part of the derivative of the photelectron 
! self energy with respect to omega.
! input: w - energy (omega)
! output: drbeta - the imaginary part of the derivative of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,drbeta,wmax,wmin,                  & !KJ ,beta,grater
     &                 abr,rlr,xsing(20),error,qmax  !KJ bugfix 11-2011 xsing->xsing(20)
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      !KJdouble precision rseint1,rseint2,rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      double precision,external :: beta,grater,drseint1,drseint2,drseint3
      wp=w+ekp
      qmax=1.d2*dsqrt(ompl)
!     rlr for plasmon pole should be 10^-7 to eliminate numerical 
!     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        drbeta=grater(drseint1,pk+qf,qmax,abr,rlr,nsing,xsing,          &
     &             error,numcal,maxns)
        drbeta=drbeta+grater(drseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,   &
     &             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,pk-qf,pk+qf,abr,rlr,nsing,xsing,  &
     &             error,numcal,maxns)
      elseif (pk.lt.qf) then
        drbeta=grater(drseint1,pk+qf,qmax,abr,rlr,nsing,xsing,          &
     &             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,qf-pk,pk+qf,abr,rlr,nsing,xsing,  &
     &             error,numcal,maxns)
        if (lowq.ne.0) drbeta=drbeta+grater(drseint3,0.d0,              &
     &             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        drbeta=grater(drseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,        &
     &             error,numcal,maxns)
        drbeta=drbeta+grater(drseint2,0.d0,2.d0*qf,abr,rlr,nsing,       &
     &             xsing,error,numcal,maxns)
      endif
      drbeta=drbeta*omp**2/(2.d0*pi*pk)
      return
      end

      subroutine dienergies(w,dibeta)
! Calculates the imaginary part of the derivative of the photelectron 
! self energy with respect to omega.
! input: w - energy (omega)
! output: dibeta - the imaginary part of the derivative of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
      implicit none
      integer i,npts,nnpts
      double precision w,dibeta
      integer it1,it2,nq
      double precision qh,q1,q2,q3,q0,wh,wq1,wq2,wq3,wq0,qmax
! the q's are limiting values of momenta, the w's are energies
! corresponding to these limiting momenta
      double precision qdisp,A,wdisp,test1,test2
      double precision dq0dw,dq1dw,dq2dw,dq3dw,dqhdw,xfact
! the dqdw's are derivatives of the momentum limits with
! changing energy
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      A=adisp
      q1=0.d0
      q2=0.d0
      q3=0.d0
!      xfact=0.d0
      dibeta=0.d0
! Must find the momenta which limit the 
! final integration (done analytically).
! Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
! Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
! and omega(q)+omega-(k-q)^2/2=0 for q>0.
      qmax=1.d6*qf
      call qlimits(w+ekp,pk,ompl,A,qmax,nq,q1,q2,q3)
! find derivatives of q's
      wq1=wdisp(q1)
      if (w+ekp-ef.gt.ompl) then
        dqhdw=(w+ekp-ef)/(qh*dsqrt(A**2+(w+ekp-ef)**2-ompl**2))
        dq0dw=0.d0
      elseif(ef-w-ekp.gt.ompl) then
        dq0dw=-(ef-w-ekp)/(q0*dsqrt(A**2+(ef-w-ekp)**2-ompl**2))
        dqhdw=0.d0
      else
        dqhdw=0.d0
        dq0dw=0.d0
      endif
      test1=(pk+q1)**2/2.d0-(w+ekp)+wq1
      test2=(pk-q1)**2/2.d0-(w+ekp)+wq1
      if (q1.ge.qh) then
        q1=qh
        dq1dw=dqhdw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq1dw=wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)
      else
        dq1dw=wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)
      endif
      wq2=wdisp(q2)
      test1=(pk+q2)**2/2.d0-(w+ekp)+wq2
      test2=(pk-q2)**2/2.d0-(w+ekp)+wq2
      if (q2.ge.qh) then
        q2=qh
        dq2dw=dqhdw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq2dw=wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)
      else
        dq2dw=wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)
      endif
      wq3=wdisp(q3)
      test1=(pk+q3)**2/2.d0-(w+ekp)-wq3
      test2=(pk-q3)**2/2.d0-(w+ekp)-wq3
      if (q3.ge.q0) then
        q3=q0
        dq3dw=dq0dw
      elseif (dabs(test1).lt.dabs(test2)) then
        dq3dw=wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)
      else
        dq3dw=wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)
      endif
! find derivative of imaginary part of self energy
! find contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        xfact=A*q1*(1.d0/wq1+1.d0/ompl)+q1**3/(2.d0*wq1)
        xfact=2.d0/q1-xfact/(ompl+wq1+A*q1**2/(2.d0*ompl))
        dibeta=dibeta+omp**2/(4.d0*pk*ompl)*dq1dw*xfact
        xfact=A*q2*(1.d0/wq2+1.d0/ompl)+q2**3/(2.d0*wq2)
        xfact=2.d0/q2-xfact/(ompl+wq2+A*q2**2/(2.d0*ompl))
        dibeta=dibeta-omp**2/(4.d0*pk*ompl)*dq2dw*xfact
      endif
! find contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        xfact=A*q0*(1.d0/wq0+1.d0/ompl)+q0**3/(2.d0*wq0)
        xfact=2.d0/q0-xfact/(ompl+wq0+A*q0**2/(2.d0*ompl))
        dibeta=dibeta+omp**2/(4.d0*pk*ompl)*dq0dw*xfact
        xfact=A*q3*(1.d0/wq3+1.d0/ompl)+q3**3/(2.d0*wq3)
        xfact=2.d0/q3-xfact/(ompl+wq3+A*q3**2/(2.d0*ompl))
        dibeta=dibeta-omp**2/(4.d0*pk*ompl)*dq3dw*xfact
      endif
      dibeta=dibeta
! Test that dibeta is a number (stop for nan results).
      it1=int(dibeta)
      it2=int(2.d0*dibeta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'dienergies ',dibeta
        stop
      endif
      dibeta=dibeta
      return
      end

      subroutine d2renergies(w,d2rbeta)
! Calculates the real part of the second derivative of the photelectron 
! self energy with respect to omega.
! input: w - energy (omega)
! output: d2rbeta - the real part of the second derivative of the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable
      implicit none
      integer i,npts,nnpts,nsing,numcal,maxns
      double precision w,d2rbeta,beta,grater,wmax,wmin,                 &
     &                 abr,rlr,xsing(20),error,qmax  !KJ bugfix 11-2011 xsing->xsing(20)
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision d2rseint1,d2rseint2,d2rseint3
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      integer lowq
      common /belowqf/ lowq
      external beta,grater,d2rseint1,d2rseint2,d2rseint3
      wp=w+ekp
      wmin=-100*(pi*beta(0.d0)+ompl)+w
      wmax=100*(pi*beta(0.d0)+ompl)+w
      qmax=1.d2*dsqrt(ompl)
!     rlr for plasmon pole should be 10^-7 to eliminate numerical 
!     artifacts
      rlr=1.d-7
      abr=rlr/1.d3
      nsing=0
      if (pk.gt.qf) then
        d2rbeta=grater(d2rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,        &
     &             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint1,0.d0,pk-qf,abr,rlr,nsing,xsing,&
     &             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,pk-qf,pk+qf,abr,rlr,nsing,     &
     &             xsing,error,numcal,maxns)
      elseif (pk.lt.qf) then
        d2rbeta=grater(d2rseint1,pk+qf,qmax,abr,rlr,nsing,xsing,        &
     &             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,qf-pk,pk+qf,abr,rlr,nsing,     &
     &             xsing,error,numcal,maxns)
        if (lowq.ne.0) d2rbeta=d2rbeta+grater(d2rseint3,0.d0,           &
     &             qf-pk,abr,rlr,nsing,xsing,error,numcal,maxns)
      else
        d2rbeta=grater(d2rseint1,2.d0*qf,qmax,abr,rlr,nsing,xsing,      &
     &             error,numcal,maxns)
        d2rbeta=d2rbeta+grater(d2rseint2,0.d0,2.d0*qf,abr,rlr,nsing,    &
     &             xsing,error,numcal,maxns)
      endif
      d2rbeta=d2rbeta*omp**2/(2.d0*pi*pk)
      return
      end

      subroutine d2ienergies(w,d2ibeta)
! Calculates the imaginary part of the second derivative of the photelectron 
! self energy with respect to omega.
! input: w - energy (omega)
! output: d2ibeta - the imaginary part of the second derivative of 
!                   the self energy
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      integer i,npts,nnpts
      double precision w,d2ibeta
      integer it1,it2,nq
      double precision qh,q1,q2,q3,q0,wh,wq1,wq2,wq3,wq0,qmax
! the q's are limiting values of momenta, the w's are energies
! corresponding to these limiting momenta
      double precision qdisp,A,wdisp,test1,test2,www,d2fact,xfact
      double precision dq0dw,dq1dw,dq2dw,dq3dw,dqhdw,dwqdw
! the dqdw's are derivatives of the momentum limits with
! changing energy
      double precision d2q0dw2,d2q1dw2,d2q2dw2,d2q3dw2,d2qhdw2
! the d2qdw2's are second derivatives of the momentum limits with
! changing energy (but you've already figured that out by now, right?)
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp,d2fact
      A=adisp
      q1=0.d0
      q2=0.d0
      q3=0.d0
      d2ibeta=0.d0
! Must find the momenta which limit the 
! final integration (done analytically).
! Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
! Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
! and omega(q)+omega-(k-q)^2/2=0 for q>0.
      qmax=1.d6*qf
      call qlimits(w+ekp,pk,ompl,A,qmax,nq,q1,q2,q3)
! find derivatives of q's
      wq1=wdisp(q1)
      if (w+ekp-ef.gt.ompl) then
        www=A**2+(w+ekp-ef)**2-ompl**2
        dqhdw=(w+ekp-ef)/(qh*dsqrt(www))
        dq0dw=0.d0
        d2qhdw2=1/(qh*dsqrt(www))-(w+ekp-ef)**2/(qh*dsqrt(www)**3)      &
     &          -(w+ekp-ef)*dqhdw/(qh**2*dsqrt(www))
        d2q0dw2=0.d0
      elseif(ef-w-ekp.gt.ompl) then
        www=A**2+(w+ekp-ef)**2-ompl**2
        dq0dw=(w+ekp-ef)/(q0*dsqrt(www))
        dqhdw=0.d0
        d2q0dw2=1/(q0*dsqrt(www))-(w+ekp-ef)**2/(q0*dsqrt(www)**3)      &
     &          -(w+ekp-ef)*dq0dw/(q0**2*dsqrt(www))
        d2qhdw2=0.d0
      else
        dqhdw=0.d0
        dq0dw=0.d0
        d2qhdw2=0.d0
        d2q0dw2=0.d0
      endif
      test1=(pk+q1)**2/2.d0-(w+ekp)+wq1
      test2=(pk-q1)**2/2.d0-(w+ekp)+wq1
      if (q1.ge.qh) then
        q1=qh
        dq1dw=dqhdw
        d2q1dw2=d2qhdw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq1dw=wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)
        dwqdw=(A*q1+q1**3/2.d0)/wq1*dq1dw
        d2q1dw2=dwqdw/((q1+pk)*wq1+A*q1+q1**3/2.d0)                     &
     &  -wq1/((q1+pk)*wq1+A*q1+q1**3/2.d0)                              &
     &  *((wq1+A+3*q1**2/2.d0)*dq1dw+(q1+pk)*dwqdw)
      else
        dq1dw=wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)
        dwqdw=(A*q1+q1**3/2.d0)/wq1*dq1dw
        d2q1dw2=dwqdw/((q1-pk)*wq1+A*q1+q1**3/2.d0)                     &
     &  -wq1/((q1-pk)*wq1+A*q1+q1**3/2.d0)                              &
     &  *((wq1+A+3*q1**2/2.d0)*dq1dw+(q1-pk)*dwqdw)
      endif
      wq2=wdisp(q2)
      test1=(pk+q2)**2/2.d0-(w+ekp)+wq2
      test2=(pk-q2)**2/2.d0-(w+ekp)+wq2
      if (q2.ge.qh) then
        q2=qh
        dq2dw=dqhdw
        d2q2dw2=d2qhdw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq2dw=wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)
        dwqdw=(A*q2+q2**3/2.d0)/wq2*dq2dw
        d2q2dw2=dwqdw/((q2+pk)*wq2+A*q2+q2**3/2.d0)                     &
     &  -wq2/((q2+pk)*wq2+A*q2+q2**3/2.d0)                              &
     &  *((wq2+A+3*q2**2/2.d0)*dq2dw+(q2+pk)*dwqdw)
      else
        dq2dw=wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)
        dwqdw=(A*q2+q2**3/2.d0)/wq2*dq2dw
        d2q2dw2=dwqdw/((q2-pk)*wq2+A*q2+q2**3/2.d0)                     &
     &  -wq2/((q2-pk)*wq2+A*q2+q2**3/2.d0)                              &
     &  *((wq2+A+3*q2**2/2.d0)*dq2dw+(q2-pk)*dwqdw)
      endif
      wq3=wdisp(q3)
      test1=(pk+q3)**2/2.d0-(w+ekp)-wq3
      test2=(pk-q3)**2/2.d0-(w+ekp)-wq3
      if (q3.ge.q0) then
        q3=q0
        dq3dw=dq0dw
        d2q3dw2=d2q0dw2
      elseif (dabs(test1).lt.dabs(test2)) then
        dq3dw=wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)
        dwqdw=(A*q3+q3**3/2.d0)/wq3*dq3dw
        d2q3dw2=dwqdw/((q3+pk)*wq3-A*q3-q3**3/2.d0)                     &
     &  -wq3/((q3+pk)*wq3-A*q3-q3**3/2.d0)                              &
     &  *((wq3-A-3*q3**2/2.d0)*dq3dw+(q3+pk)*dwqdw)
      else
        dq3dw=wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)
        dwqdw=(A*q3+q3**3/2.d0)/wq3*dq3dw
        d2q3dw2=dwqdw/((q3-pk)*wq3-A*q3-q3**3/2.d0)                     &
     &  -wq3/((q3-pk)*wq3-A*q3-q3**3/2.d0)                              &
     &  *((wq3-A-3*q3**2/2.d0)*dq3dw+(q3-pk)*dwqdw)
      endif
! find second derivative of imaginary part of self energy
! find contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        d2ibeta=d2ibeta+d2fact(q1,wq1,dq1dw,d2q1dw2)
        d2ibeta=d2ibeta-d2fact(q2,wq2,dq2dw,d2q2dw2)
      endif
! find contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        d2ibeta=d2ibeta+d2fact(q0,wq0,dq0dw,d2q0dw2)
        d2ibeta=d2ibeta-d2fact(q3,wq3,dq3dw,d2q3dw2)
      endif
! Test that d2ibeta is a number (stop for nan results).
      it1=int(d2ibeta)
      it2=int(2.d0*d2ibeta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2ienergies ',d2ibeta
        stop
      endif
      d2ibeta=d2ibeta
      return
      end

      double precision function d2fact(q,wq,dqdw,d2qdw2)
! the terms in a sum used to calculate the imaginary part of
! the second derivative of the self energy in the subroutine
! d2ienergies.
! input: q - momentum limit
!        wq - energy corresponding to q
!        dqdw - derivative of momentum limit with respect to energy
!        d2qdw2 - second derviative of momentum limit
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      double precision q,wq,dqdw,d2qdw2,xfact1,xfact2,dwqdw
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,A
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,A
      dwqdw=(A*q+q**3/2.d0)/wq*dqdw
      xfact1=((A*(1.d0/wq+1.d0/ompl)+3.d0*q**2/(2.d0*wq))*dqdw          &
     &  -dwqdw*(q**3+2.d0*A*q)/(2*wq**2))                               &
     &  /(ompl+wq+A*q**2/(2.d0*ompl))
      xfact1=-xfact1+(dwqdw+A*q*dqdw/ompl)                              &
     &  *(A*q*(1.d0/wq+1.d0/ompl)+q**3/(2.d0*wq))                       &
     &  /(ompl+wq+A*q**2/(2.d0*ompl))**2
      xfact1=xfact1-2.d0*dqdw/q**2
      xfact2=A*q*(1.d0/wq+1.d0/ompl)+q**3/(2.d0*wq)
      xfact2=2.d0/q-xfact2/(ompl+wq+A*q**2/(2.d0*ompl))
      d2fact=omp**2/(4.d0*pk*ompl)*(dqdw*xfact1+d2qdw2*xfact2)
      return
      end

      double precision function rseint1(q)
! Integrand for one of the intergals in calculating the real part
! of the self energy in subroutine renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      integer ijkwrite
      common /morewrite/ ijkwrite
      wq=wdisp(q)
      xlog=((pk+q)**2/2.d0-wp+wq)**2+(acc*ompl)**2
      xlog=xlog/(((pk-q)**2/2.d0-wp+wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint1=xlog/(wq*dsqrt(q**2+ompl*acc))
!* NAN detector
!      it1=int(rseint1)
!      it2=int(2.d0*rseint1)
!      if (it1.eq.it2.and.it1.gt.5) then
!        write(6,*) 'rseint1 ',rseint1,q,wq,xlog,wp,pk,acc,ompl
!        stop
!      endif
      return
      end

      double precision function rseint2(q)
! Integrand for one of the intergals in calculating the real part
! of the self energy in subroutine renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xlog=1.d0
      if (lowq.ne.0) then
        xlog=(ef-wp-wq)**2+(acc*ompl)**2
        xlog=xlog/(((pk-q)**2/2.d0-wp-wq)**2+(acc*ompl)**2)
      endif
      xlog=xlog*(((pk+q)**2/2.d0-wp+wq)**2+(acc*ompl)**2)
      xlog=xlog/((ef-wp+wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint2=xlog/(wq*dsqrt(q**2+ompl*acc))
!* NAN detector
!      it1=int(rseint2)
!      it2=int(2.d0*rseint2)
!      if (it1.eq.it2.and.it1.gt.5) then
!        write(6,*) 'rseint2 ',rseint2,q,wq,xlog,wp,pk,acc,ompl
!        stop
!      endif
      return
      end

      double precision function rseint3(q)
! Integrand for one of the intergals in calculating the real part
! of the self energy in subroutine renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xlog
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xlog=((pk+q)**2/2.d0-wp-wq)**2+(acc*ompl)**2
      xlog=xlog/(((pk-q)**2/2.d0-wp-wq)**2+(acc*ompl)**2)
      xlog=log(xlog)/2.d0
      rseint3=xlog/(wq*dsqrt(q**2+ompl*acc))
!* NAN detector
!      it1=int(rseint3)
!      it2=int(2.d0*rseint3)
!      if (it1.eq.it2.and.it1.gt.5) then
!        write(6,*) 'rseint3 ',rseint3,q,wq,xlog,wp,pk,acc,ompl
!        stop
!      endif
      return
      end

      double precision function drseint1(q)
! Integrand for one of the intergals in calculating the derivative 
! of the real part of the self energy in subroutine drenergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp+wq
      xfact=xnum1/(xnum1**2+brd**2)
      xnum2=(pk-q)**2/2.d0-wp+wq
      xfact=xfact-xnum2/(xnum2**2+brd**2)
      drseint1=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint1)
      it2=int(2.d0*drseint1)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint1 ',drseint1,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function drseint2(q)
! Integrand for one of the intergals in calculating the derivative 
! of the real part of the self energy in subroutine drenergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2,xnum3,xnum4
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xfact=0.d0
      if (lowq.ne.0) then
        xnum1=ef-wp-wq
        xfact=xnum1/(xnum1**2+brd**2)
        xnum2=(pk-q)**2/2.d0-wp-wq
        xfact=xfact-xnum2/(xnum2**2+brd**2)
      endif
      xnum3=(pk+q)**2/2.d0-wp+wq
      xfact=xfact+xnum3/(xnum3**2+brd**2)
      xnum4=ef-wp+wq
      xfact=xfact-xnum4/(xnum4**2+brd**2)
      drseint2=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint2)
      it2=int(2.d0*drseint2)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint2 ',drseint2,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function drseint3(q)
! Integrand for one of the intergals in calculating the derivative 
! of the real part of the self energy in subroutine drenergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp-wq
      xfact=xnum1/(xnum1**2+brd**2)
      xnum2=(pk-q)**2/2.d0-wp-wq
      xfact=xfact-xnum2/sqrt(xnum2**2+brd**2)
      drseint3=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(drseint3)
      it2=int(2.d0*drseint3)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'drseint3 ',drseint3,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function d2rseint1(q)
! Integrand for one of the intergals in calculating the second derivative 
! of the real part of the self energy in subroutine d2renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp+wq
      xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
      xnum2=(pk-q)**2/2.d0-wp+wq
      xfact=xfact-(xnum2**2-brd**2)/(xnum2**2+brd**2)**2
      d2rseint1=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint1)
      it2=int(2.d0*d2rseint1)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint1 ',d2rseint1,q,wq,xfact,wp,pk,acc,ompl
        stop
      endif
      return
      end

      double precision function d2rseint2(q)
! Integrand for one of the intergals in calculating the second derivative 
! of the real part of the self energy in subroutine d2renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2,xnum3,xnum4
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      wq=wdisp(q)
      xfact=0.d0
      if (lowq.ne.0) then
        xnum1=ef-wp-wq
        xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
        xnum2=(pk-q)**2/2.d0-wp-wq
        xfact=xfact-(xnum2**2-brd**2)/(xnum2**2+brd**2)**2
      endif
      xnum3=(pk+q)**2/2.d0-wp+wq
      xfact=xfact+(xnum3**2-brd**2)/(xnum3**2+brd**2)**2
      xnum4=ef-wp+wq
      xfact=xfact-(xnum4**2-brd**2)/(xnum4**2+brd**2)**2
      d2rseint2=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint2)
      it2=int(2.d0*d2rseint2)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint2 ',d2rseint2
        stop
      endif
      return
      end

      double precision function d2rseint3(q)
! Integrand for one of the intergals in calculating the second derivative 
! of the real part of the self energy in subroutine d2renergies.
! input: q - momentum (variable to be integrated over)
! input from common blocks
!       ompl - energy of pole ipl in epsilon^{-1}
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
! common block control of subprograms
!       ac2 - additional accuracy parameter
!       wp - omega prime, an additional energy variable to be held
!            constant durring the integration
      implicit none
      integer it1,it2,i
      double precision q,wq,xfact,xnum1,xnum2
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      double precision ac2,wp
      common /ff/ ac2,wp
      double precision qdisp,wdisp
      external qdisp,wdisp
      wq=wdisp(q)
      xnum1=(pk+q)**2/2.d0-wp-wq
      xfact=(xnum1**2-brd**2)/(xnum1**2+brd**2)**2
      xnum2=(pk-q)**2/2.d0-wp-wq
      xfact=xfact-(xnum2**2-brd**2)/sqrt(xnum2**2+brd**2)**2
      d2rseint3=xfact/(wq*dsqrt(q**2+ompl*acc))
      it1=int(d2rseint3)
      it2=int(2.d0*d2rseint3)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'd2rseint3 ',d2rseint3
        stop
      endif
      return
      end

      double precision function beta(w)
! the extrinsic beta function
! beta(k,w)=(1/pi)*|Im(self energy(k,w+(k^2)/2))|*theta((w+(k^2)/2)-xmu)
! input: w - energy (omega)
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent 
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole ipl in epsilon^{-1}
!       wt - weight of pole ipl in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of 
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       pk - photoelectron momentum
!       acc - global accuracy parameter 
!       brd - global broadening parameter to stabilize logarithms
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
      implicit none
      integer it1,it2,nq,i
      double precision w,q1,q2,q3,wth,wq0,wq1,wq2,wq3
      double precision q0,qh,qdisp,A,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
!      common /temp/ q1,q2
      integer lowq
      common /belowqf/ lowq
      external qdisp,wdisp
      A=adisp
      q1=0.d0
      q2=0.d0
! Must find the momenta which limit the
! final integration (done analytically).
      beta=0.d0
! Find limit due to Fermi level
      qh=qdisp(max(w+ekp-ef,ompl))
      q0=qdisp(max(ef-w-ekp,ompl))
! Find roots of omega(q)-omega+(q-k)^2/2=0, omega(q)-omega+(q+k)^2/2=0,
! and omega(q)+omega-(k-q)^2/2=0.
      call qlimits(w+ekp,pk,ompl,A,qh,nq,q1,q2,q3)
! Calculate beta
! Calculate contributions from above Fermi momentum
      if (nq.eq.3) then
        q1=dsqrt(q1**2+acc*ompl)
        q2=dsqrt(q2**2+acc*ompl)
        wq1=wdisp(q1)
        wq2=wdisp(q2)
        beta=beta+omp**2/(4*pi*pk*ompl)                                 &
     &       *dlog(q2**2/(ompl+wq2+A*q2**2/(2.d0*ompl))                 &
     &       *(ompl+wq1+A*q1**2/(2.d0*ompl))/q1**2)
!        beta=beta+ompl/(4*pi*pk)
!     2       *dlog(q2**2/(ompl+wq2+A*q2**2/(2.d0*ompl)))
      endif
! Calculate contributions from below Fermi momentum
      if (q3.lt.q0.and.lowq.ne.0) then
        q0=dsqrt(q0**2+acc*ompl)
        q3=dsqrt(q3**2+acc*ompl)
        wq0=wdisp(q0)
        wq3=wdisp(q3)
        beta=beta-omp**2/(4*pi*pk*ompl)                                 &
     &       *dlog(q0**2/(ompl+wq0+A*q0**2/(2.d0*ompl))                 &
     &       *(ompl+wq3+A*q3**2/(2.d0*ompl))/q3**2)
      endif
! Test that beta is a number (stop for nan results).
      it1=int(beta)
      it2=int(2.d0*beta)
      if (it1.eq.it2.and.it1.gt.5) then
        write(6,*) 'beta ',beta,q1,q2,q3,q0
        write(6,*) 'beta ',pk,w,ekp,acc,ompl
        write(6,*) nq
        write(6,*) q1,q2
        stop
      endif
      return
      end
      
