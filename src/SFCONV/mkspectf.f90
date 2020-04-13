!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mkspectf.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Computes spectral function window centered on resonance.  
!  Spectral function "deconvoluted" to remove broadening of main peak
      subroutine mkspectf(rs,pk,gammach,xreduc,wpts,spectf,             &
     &                    weights,isattype,npl,nplmax,                  &
     &                    plengy,plwt,plbrd,brpole) 
! input: rs - radius that contains one electron
!        pk - photoelectron momentum
!        gammach - core hole lifetime broadening
!        xreduc - add-hoc interference reduction
!        isattype - approximation to use when computing satellite.
!        npl - number of poles in epsilon^{-1}
!        nplmax - maximum number of poles for dimensioning arrays
!        plengy - energy of each pole in epsilon^{-1}
!        plwt - weight of each pole in epsilon^{-1}
!        plbrd - broadening of each pole in epsilon^{-1}
!        brpole - true if poles are to be calculated with broadening
! output: wpts - array of energy points on which the spectral
!                function is specified.  array is of length npts
!         spectf - array of values of the spectral function.
!                array is of length npts
!         weights - array containing the total spectral weight of 
!                various contributions to the spectral function
!         weights(1) - extrinsic quasiparticle weight
!         weights(2) - asymmetry of quasiparticle
!         weights(3) - interference quasiparticle weight
!         weights(4) - extrinsic satellite weight
!         weights(5) - interference satellite weight
!         weights(6) - intrinsic satellite weight
!         weights(7) - weight of the extrinsic satellite clipped 
!                      to include only that part in the satellite region
!         weights(8) - weight of the extrinsic satellite clipped
!                      to include that part in the vicinity of the
!                      quasiparticle
      implicit none
      logical brpole
      integer npts,i,ii,j,jj,k,iemax,iixmax,ishift,isetold,             &
     &        iset,iset2,iswitch,isattype,iqph,iqpl,                    &
     &        isetd,isetd2,iswd,iswd2,jwidth,ipl
!             npts - number of grid points for spectral function
!             i,ii,j,jj,k - counters and dummy variables
!             iixmax - grid point of highest value of the
!                   interference spectral function
!             ishift - number of grid points to shift the 
!                   extrinsic spectral function
!             iset,iset2,isetd,isetd2,isetold - keeping track
!                   of trigger conditions to separate the 
!                   extrinsic satellite into "satellite" and
!                   "main peak" regions
!             iswd,iswd2,iswitch - keeping track of grid points
!                   at which the above triggering conditions are met
!             iqph,iqpl - array points that bound the quasiparticle
!             jwidth - extra region tacked onto integration for 
!                   lorentzian broadening of the self energy
!             ipl - counter for loops over poles
!      parameter(npts=80)
      parameter(npts=112)
      integer npl,nplmax
      double precision rs,pk,gammach,xreduc,                            &
     &                 spectf(8,npts),wpts(npts)
      double precision plengy(nplmax),plwt(nplmax),plbrd(nplmax)
      double precision conc,xa,expa,ak,aangstrom,eV,dsat,d2sat,         &
     &                 specttmp(npts),wswitch,swidth,wd,wd2,            &
     &                 dsatold,d2satold,dsath,wsearch,esfhi,esflo
!             conc - electron concentration
!             xa - dimensionless coupling constant of the electron gas
!             expa - the base of the natural logarithm raised to the
!                  power xa (e**xa)
!             ak - interference contribution to the quasiparticle
!             aangstrom - convert from angstrom to Bohr units of length
!             eV - convert from electron volts to Hartree units of energy
!             dsat,d2sat,dsatold,d2satold,dsath - finite difference 
!                  derivatives of the extrinsic satellite, to find
!                  the triggers to separate the "main peak" from the
!                  "satellite" structure of the extrinsic satellite 
!                  spectral function
!             specttmp - interpolated spectral function for shifting
!                  the extrinsic spectral function without numerical
!                  jitter
!             wsearch - energy used for finding the shift in the
!                  extrinsic spectral function
!             wswitch,wd,wdold - keep track of energy at which the trigger
!                  for separating the extrinsic spectral function
!                  is met
!             esfhi,esflo - used to bracket the high point of the
!                  extrinsic spectral function for its shift
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      double precision xfact1,xfact2,wemax,wshift,wshift2,              &
     &                 whi,wlo,whi2,wlo2,ehi,elo,emax,ehi1,ehi2,        &
     &                 elo1,elo2,delta
!             xfact1,xfact2 - dummy variables for keeping track of 
!                  lengthly computations
!             wemax - the energy at the maximum value of the 
!                  extrinsic spectral function
!             wshift,wshift2 - amounts to shift the extrinsic 
!                  spectral function
!             whi,wlo,whi2,wlo2 - used to bracket the search region for the
!                  maximum of the extrinsic spectral function
!             ehi,elo,ehi1,ehi2,elo1,elo2 - values of the extrinsic 
!                  spectral function at the above bracket points
!             emax - the maximum value of the extrinsic spectral function
!             delta - a small value
      double precision sefr(npts),sefi(npts),                           &
     &                 sef2r(npts),sef2i(npts),sumr,sumi,brpl,          &
     &                 wh,wl,w2,wh2,wl2,dw2
!             sefr - array of the real part of the self energy 
!                  values on the energy grid wpts
!             sefi - array of the imaginary part of the self energy 
!                  values on the energy grid wpts
!             sef2r,sef2i - self energy after lorentzian broadening
!             sumr,sumi - running integration counters
!             brpl - the broadening of a given pole in epsilon^{-1}
!             wh,wl,wh2,wl2 - high and low values of energy intervals,
!                  used for integration bounds over that interval
!             w2 - an energy variable
!             dw2 - an energy interval
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       xmu - chemical potential = Fermi energy + self consistent
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of ipl'th pole in epsilon^{-1}
!       wt - weight of ipl'th pole in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       qpk - photoelectron momentum
!       acc - global accuracy parameter
!       brd - width of pole in epsilon^(-1)
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
      double precision wmin,wmax,wmin1,wmax1,dw,dw1,w,wlim(0:npts)
      common /limits/ wmin,wmax,wmin1,wmax1
!       wmin,wmax,wmin1,wmax1 - a bunch of extreme values of the energy
!       dw,dw1 - energy intervals
!       w - energy
!       wlim - energy values separating the energy gridpoints
      double precision se,xibeta,ce,width,xnn,xaa,se2,xise,             &
     &                 z0,z0i,zm,z1,z1i,z1m,xrz,xiz,x2a,x2n,            &
     &                 sef0,se0,yr,yi,zzr,zzi,zzm,qpengy,qpwidth
!       se - real part of the on shell self energy
!       xibeta - imaginary part of the on shell self energy, primarily 
!             used only to find width
!       ce - core energy
!       width - quasiparticle broadening
!       xnn - on shell energy derivative of the real part of the self
!             energy
!       xaa - on shell energy derivative of the imaginary part of the self
!             energy
!       se2 - real part of the self energy calculated at energy w
!       xise - imaginary part of the self energy calculated at energy w
!       z0 - approximation for the real part of the renormalization constant
!       z0i - approximation for the imaginary part of the 
!             renormalization constant
!       zm - magnitude of the above approximation for the 
!             renormalization constant
!       z1,z1i,z1m - a better approximation for the renormalization constant 
!       xrz,xiz - variables for intermediate steps in the calculation of
!             the renormalization constant
!       x2n,x2a - the real and imaginary parts, respectively, of the
!             second energy derivative of the on shell self energy.
!             These variables are not currently used in this version
!             of the code, but if needed, here they are.
!       sef0 - self energy at the fermi level
!       se0 - approximation for the self energy
!       yr - real part of second derivative of self energy
!       yi - imaginary part of second derivative of self energy
!       zzr,zzi,zzm - refined approximation for renormalization constant
!       qpengy - refined approximation for quasiparticle energy
!       qpwidth - refined approximation for quasiparticle width
      double complex zzc,yyc,xxc,discr,epc
!       zzc - complex renormalization constant
!       yyc - second derivative of self energy
!       xxc - 1-derivative of self energy
!       discr - square root of complex discriminant
!       epc - complex energy of quasiparticle pole
      double precision sxnn,sxaa,sse,sxise
!       for sums over the poles in epsilon^{-1} when computing the 
!       self energy and its derivatives
      common /energies/ se,ce,width,z1,z1i,se2,xise
      double precision xsat,xmain,emain,esat,xisat,esat2
!       the values of the interference satellite, interference 
!       quasiparticle, extrinsic quasiparticle, extrinsic satellite,
!       intrinsic satellite, and quasiparticle asymmetry spectral 
!       functions, respectively
      double precision xwidth, xiwidth
!       xwidth - broadening of the interference satellite
!       xiwidth - broadening of the intrinsic satellite
      double precision wtemain,wtesat,wtxmain,wtxsat,wtisat,            &
     &                 wtmesat,wtmemain,wtesat2,weights(8)
!       the total spectral weights of the extrinsic quasiparticle,
!       extrinsic satellite, interference quasiparticle, 
!       interference satellite, intrinsic satellite, clipped 
!       "satellite" part of the extrinsic satellite, clipped
!       "quasiparticle" part of the extrinsic satellite, and quasipartcle
!       asymmetry terms in the spectral function, respectively, in 
!       addition to the array of spectral weights
      double precision swtcorr,satwt,swtfac
!       swtcorr - the weight of the negative regions of the satellite
!       satwt - the satellite weight (neglecting factors of expa, as above)
!       swtfac - renormalization factor to keep satellite weight the 
!          same after negative parts are chopped off
      double precision grater,abr,rlr,xsing(20),error  !KJ bugfix 11-2011 xsing->xsing(20)
      integer nsing,numcal,maxns
!       Grater is an integration function.  See its description for
!       the meaning of these variables.
      double precision xmkesat,beta,xmkak,xmkxsat,xmkisat, &
           & xmkgwext,xmkpsat
!       These are functions called in the execution of this subroutine.
      integer iwrite,jcount
      common /flag/ iwrite,jcount
!       These are used to flag a particular run of this subroutine
!       to write intermediate results to a file.
      integer lowq
!       lowq - If not equal to zero, calculate contributions
!          to self energy from below Fermi level.
      common /belowqf/ lowq
      external beta,grater,xmkesat,xmkxsat,xmkisat,xmkak,xmkgwext
      double precision w3,q3,q4
      common /fqso2/ q3
      common /fw/ w
!     Josh - comment out extraneous externals
!      double precision wintr1,winti1,wintr2,winti2,xdumr,xdumi
!      external wintr1,winti1,wintr2,winti2
!      double precision srrpaw1,sirpaw1,rpadsfr,rpadsfi
!      external srrpaw1,sirpaw1,rpadsfr,rpadsfi
!      double precision rpapolr,rpapoli
!     external rpapolr,rpapoli
!      double precision xtestr,xtesti,qq,wplas,dq
!      external xtestr,xtesti,wplas
!      double precision rpasigr,rpasigi
!      external rpasigr,rpasigi
!     End Josh
      double precision exchange
      external exchange
      integer iw
      integer ijkwrite
      common /morewrite/ ijkwrite
! test functions

! compute initial values
      ijkwrite=0
      acc=1.d-4
      pi=dacos(-1.d0)
      qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
      ef=qf*qf/2.d0
      conc=3.d0/(4.d0*pi*(rs**3))
      omp=dsqrt(4.d0*pi*conc)
      qpk=qf
      ek=ef
      ekp=ef
      adisp=2.d0*ef/3.d0
      sef0=0.d0
!     Compute self energy at the fermi level
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(0.d0,sse,sxise)
        else
          call renergies(0.d0,sse)
        endif
!        call ppset(rs,pi,qf,ef,omp)
        sef0=sef0+sse*wt
      enddo
      sef0=sef0+exchange(qf)
!      call ppset(rs,pi,qf,ef,omp)
      xmu=ef+sef0
      ekp=xmu

! first estimate for on shell energies and self energies
      qpk=pk
      ek=pk*pk/2.d0
      ekp=ek
      se0=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(0.d0,sse,sxise)
          call dbrsigma(0.d0,sxnn,sxaa)
        else
          call renergies(0.d0,sse)
          call drenergies(0.d0,sxnn)
          call dienergies(0.d0,sxaa)
        endif
!        call ppset(rs,pi,qf,ef,omp)
        se0=se0+sse*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se0=se0+exchange(pk)
!      call ppset(rs,pi,qf,ef,omp)
      xrz=1.d0-xnn
      xiz=-xaa
      z0=xrz/(xrz**2+xiz**2)
      z0i=-xiz/(xrz**2+xiz**2)
      zm=dsqrt(z0**2+z0i**2)
      ekp=ek+sef0+z0*(se0-sef0)

! refined estimate for on shell self energies, using the self energy at the 
! Fermi level to increase self consistency
      se=0.d0
      xibeta=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(-sef0,sse,sxise)
          call dbrsigma(-sef0,sxnn,sxaa)
        else
          call renergies(-sef0,sse)
          call xienergies(-sef0,sxise)
          call drenergies(-sef0,sxnn)
          call dienergies(-sef0,sxaa)
        endif
        se=se+sse*wt
        xibeta=xibeta+sxise*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se=se+exchange(pk)
      width=dabs(xibeta)+gammach

! calculate renormalization constant
      xa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        xa=xa+3*wt*(omp/ompl)**2/(8*dsqrt(2.d0*ompl))
      enddo
      expa=exp(-xa)
      xrz=1.d0-xnn
      xiz=-xaa
      z1=xrz/(xrz**2+xiz**2)
      z1i=-xiz/(xrz**2+xiz**2)
      z1m=sqrt(z1**2+z1i**2)

! Find energy gridpoints
      dw=omp/30.d0
      iqph=54
      iqpl=53
      wpts(iqph)=dw*1.d-2
      wpts(iqpl)=-dw*1.d-2
      wpts(iqph+1)=dw*2.d-2
      wpts(iqpl-1)=-dw*2.d-2
      do i=1,30
        wpts(i+1+iqph)=i*dw
        wpts(iqpl-1-i)=-i*dw
      enddo
      do i=1,3
        wpts(i+31+iqph)=wpts(31+iqph)+i*dw
        wpts(iqpl-31-i)=wpts(iqpl-31)-i*dw
      enddo
      do i=1,3
        wpts(i+34+iqph)=wpts(34+iqph)+(2*i)*dw
        wpts(iqpl-34-i)=wpts(iqpl-33)-(2*i)*dw
      enddo
      do i=1,3
        wpts(i+37+iqph)=wpts(37+iqph)+(4*i)*dw
        wpts(iqpl-37-i)=wpts(iqpl-36)-(4*i)*dw
      enddo
      do i=1,12
        wpts(i+40+iqph)=wpts(40+iqph)+(10*i)*dw
        wpts(iqpl-40-i)=wpts(iqpl-39)-(10*i)*dw
      enddo
      do i=1,6
        wpts(i+52+iqph)=wpts(52+iqph)+(30*i)*dw
      enddo
      do i=1,npts-1
        wlim(i)=(wpts(i)+wpts(i+1))/2.d0
      enddo
      wlim(0)=2.d0*wpts(1)-wpts(2)
      wlim(npts)=2.d0*wpts(npts)-wpts(npts-1)
      wmin=wpts(1)
      wmax=wpts(npts)
      wmax1=ekp+wmax
      wmin1=ekp+wmin

!     compute contributoin to off shell self energies from the poles of loss fcn
!     on energy gridpoints wpts.
      do i=1,npts
        sefr(i)=0.d0
        sefi(i)=0.d0
      enddo
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        do i=1,npts
          w=wpts(i)+ekp
          dw1=wlim(i)-wlim(i-1)
!          dw1=wpts(i)-wpts(i-1)
  
!          call renergies(w-ekp-sef0,sse)
!          call xienergies(w-ekp-sef0,sxise)
          if (brpole) then
            call brsigma(w-ekp-se,sse,sxise)
          else
            call renergies(w-ekp-se,sse)
            call xienergies(w-ekp-se,sxise)
          endif
          sefr(i)=sefr(i)+sse*wt
          sefi(i)=sefi(i)+sxise*wt
        enddo
!        call ppset(rs,pi,qf,ef,omp)
      enddo

!     Further refinement of on shell self energies, third level of self consistency.
      se2=0.d0
      xise=0.d0
      xnn=0.d0
      xaa=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        if (brpole) then
          call brsigma(-se,sse,sxise)
          call dbrsigma(-se,sxnn,sxaa)
        else
          call renergies(-se,sse)
          call xienergies(-se,sxise)
          call drenergies(-se,sxnn)
          call dienergies(-se,sxaa)
        endif
        se2=se2+sse*wt
        xise=xise+sxise*wt
        xnn=xnn+sxnn*wt
        xaa=xaa+sxaa*wt
      enddo
      se=se2+exchange(pk)
      width=dabs(xise)+gammach

!     Add Hartree-Fock term to off shell self-energy
      do i=1,npts
        sefr(i)=sefr(i)+exchange(pk)
        sefi(i)=dabs(sefi(i))+gammach
        if (iwrite.eq.jcount) write(24,500) wpts(i),sefr(i),sefi(i)
      enddo
!      se=(sefr(iqph)+sefr(iqpl))/2
!      width=(sefi(iqph)+sefi(iqpl))/2
!      xnn=(sefr(iqph)-sefr(iqpl))/(wpts(iqph)-wpts(iqpl))
!      xaa=(sefi(iqph)-sefi(iqpl))/(wpts(iqph)-wpts(iqpl))
      xrz=1.d0-xnn
      xiz=-xaa
      z1=xrz/(xrz**2+xiz**2)
      z1i=-xiz/(xrz**2+xiz**2)
      z1m=sqrt(z1**2+z1i**2)
      qpengy=ekp+width*z1i
      qpwidth=width*z1
      zzr=z1
      zzi=z1i
      zzm=z1m

!      write(68,*) pk/qf,se,width,z1,z1i

! correct for endpoint effects
      ak=0.d0
      do ipl=1,npl
        call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
        ak=ak+xmkak(ekp)*xreduc*wt
!        call ppset(rs,pi,qf,ef,omp)
      enddo
      wtemain=(datan(wlim(0)/width)+pi/2.d0)/pi                         &
     &        +(pi/2.d0-datan(wlim(npts)/width))/pi
      wtxmain=2.d0*wtemain*zm*z1*ak
      wtemain=wtemain*z1*expa
      wtesat=0.d0
      wtesat2=0.d0
      wtxsat=0.d0
      wtisat=0.d0
      do i=1,npts
        do ii=1,8
          spectf(ii,i)=0.d0
        enddo
      enddo

! compute spectral function
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)

        se2=sefr(i)
        xise=sefi(i)

! This form for the extrinsic quasiparticle is calculated at the
! correct quasiparticle pole, for much better cancelation of the
! quasiparticle structure from the xmkgwext function, below.
        emain=z1*(datan((wlim(i)-qpengy+ekp)/qpwidth)                   &
     &        -datan((wlim(i-1)-qpengy+ekp)/qpwidth))/(pi*dw1)          &
     &        -z1i*dlog((qpwidth**2+(wlim(i)-qpengy+ekp)**2)/           &
     &         (qpwidth**2+(wlim(i-1)-qpengy+ekp)**2))/(2.d0*pi*dw1)    &
     &         *exp(-((w-qpengy)/(2*omp))**2)
        xmain=2.d0*zm*ak*emain
        wtemain=wtemain+emain*expa*dw1
        wtxmain=wtxmain+xmain*expa*dw1
! Depending on the value of isattype, different approximations
! can be used for the satellite.
        if (isattype.eq.1) then
          esat=xmkgwext(w-ekp)-emain
        elseif (isattype.eq.2) then
          esat=(xise-width-(w-ekp)*xaa)/(pi*(w-ekp)**2)
        elseif (isattype.eq.3) then
! generate full extrinsic spectral function - do not use for spectroscopy!
          esat=xmkgwext(w-ekp)
        else
          esat=xmkesat(w-ekp)
        endif
        xsat=0.d0
        xisat=0.d0
        do ipl=1,npl
          call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
          xwidth=max(5.d0*dw,brd)
          xiwidth=max(2.d0*dw,brd)
          if (isattype.eq.3) then
            xwidth=xwidth+width
            xiwidth=xiwidth+width
          endif
          xsat=xsat+xmkxsat(w-ekp,xwidth)*xreduc*wt
          xisat=xisat+xmkisat(w-ekp,xiwidth)*wt
!          call ppset(rs,pi,qf,ef,omp)
        enddo
        wtxsat=wtxsat+xsat*dw1*expa
        wtesat=wtesat+esat*dw1*expa
        wtisat=wtisat+xisat*dw1*expa
        spectf(1,i)=emain
        spectf(2,i)=esat
        spectf(3,i)=xmain
        spectf(4,i)=xsat
        spectf(5,i)=xisat
        spectf(6,i)=(esat+xisat-2.d0*xsat)
        if (isattype.eq.3) then
          spectf(6,i)=spectf(6,i)+xmain
        endif
      enddo
      spectf(2,iqpl)=(spectf(2,iqpl)+spectf(2,iqph))/2
      spectf(2,iqph)=spectf(2,iqpl)

!     separate quasiparticle-like structure from satellite in non-delta
!     function part of extrinsic spectral function.
      iset=0
      iset2=0
      dsatold=0.d0
      d2satold=0.d0
      isetd=0
      isetd2=0
      do ii=2,npts-1
        isetold=iset2
        i=npts+1-ii
        w=wpts(i)+ekp
        dsat=(spectf(2,i)-spectf(2,i-1))/(wpts(i)-wpts(i-1))
        dsath=(spectf(2,i+1)-spectf(2,i))/(wpts(i+1)-wpts(i))
        d2sat=(dsath-dsat)/(wlim(i)-wlim(i-1))
        if (dsat.gt.0.d0.and.spectf(2,i).gt.0.d0                        &
     &      .and.iset.eq.0.d0) iemax=i
!        if (spectf(5,i).lt.spectf(5,i+1)                                &
!     &      .and.spectf(5,i+1).gt.spectf(5,i+2)) iixmax=i
        ! Josh - Changed the if statement so that array spectf is
        ! not out of bounds.
        if (spectf(5,i-1).lt.spectf(5,i)                                &
     &      .and.spectf(5,i).gt.spectf(5,i+1)) iixmax=i
        if (dsat.gt.0.d0.and.spectf(2,i).gt.0.d0) iset=1
        if (beta(0.d0).gt.0.d0) then
          if (dsat.lt.0.d0.and.dsatold.ge.0.d0                          &
     &        .and.iset.eq.1.and.isetd.eq.0) then
             isetd=1
             iswd=i
             wd=w
          endif
          if (d2sat.gt.0.d0.and.d2satold.le.0.d0                        &
     &        .and.iset.eq.1.and.isetd2.eq.0) then
             isetd2=1
             iswd2=i
             wd2=w
          endif
        else
          if (dsat.lt.0.d0.and.dsatold.ge.0.d0                          &
     &        .and.iset.eq.1.and.isetd.eq.0.and.w.gt.0.d0) then
             isetd=1
             iswd=i
             wd=w
          endif
          if (d2sat.gt.0.d0.and.d2satold.le.0.d0                        &
     &        .and.iset.eq.1.and.isetd2.eq.0) then
             isetd2=1
             iswd2=i
             wd2=w
          endif
        endif
      enddo
      if (isetd.eq.1) then
        iswitch=iswd
        wswitch=wd
      else
        iswitch=iswd2
        wswitch=wd2
      endif
!      swidth=(omp/4.d0)*spectf(2,iswitch)/spectf(2,iemax)
      swidth=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        if (swidth.gt.0.d0) then
          spectf(7,i)=spectf(2,i)/(1.d0+exp((w-wswitch)/swidth))
          spectf(8,i)=spectf(2,i)/(1.d0+exp((wswitch-w)/swidth))
        elseif (i.ge.iswitch) then
          spectf(8,i)=spectf(2,i)
        else
          spectf(7,i)=spectf(2,i)
        endif
      enddo

!* find the energy at which the extrinsic satellite reaches its 
!* maximum value.
!      wemax=wpts(iemax)+ekp
!      whi=wpts(iemax+1)+ekp
!      wlo=wpts(iemax-1)+ekp
!      emax=spectf(2,iemax)
!      do i=1,8
!        whi2=(whi+wemax)/2.d0
!        wlo2=(wlo+wemax)/2.d0
!        if (whi2-ekp.ne.0.d0) then
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          ehi=xmkesat(whi2-ekp)
!        else
!          delta=(whi2-wemax)/1.d3
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          ehi1=xmkesat(whi2-ekp+delta)
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          ehi2=xmkesat(whi2-ekp-delta)
!          ehi=(ehi1+ehi2)/2.d0
!        endif
!        if (wlo2-ekp.ne.0.d0) then
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          elo=xmkesat(wlo2-ekp)
!        else
!          delta=(wemax-wlo2)/1.d3
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          elo1=xmkesat(wlo2-ekp+delta)
!          se2=0.d0
!          xise=0.d0
!          do ipl=1,npl
!            call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
!            call brsigma(whi2-ekp-sef0,sse,sxise)
!*            call ppset(rs,pi,qf,ef,omp)
!            se2=se2+sse*wt
!            xise=xise+sxise*wt
!          enddo
!          se2=se2+exchange(pk)
!          xise=dabs(xise)+gammach
!          elo2=xmkesat(wlo2-ekp-delta)
!          elo=(elo1+elo2)/2.d0
!        endif
!        if (ehi.gt.emax.and.ehi.gt.elo) then
!          wlo=wemax
!          wemax=whi2
!          emax=ehi
!        elseif(elo.gt.emax.and.elo.ge.ehi) then
!          whi=wemax
!          wemax=wlo2
!          emax=elo
!        else
!          whi=whi2
!          wlo=wlo2
!        endif
!      enddo
!      wemax=wemax-ekp
!
!* if the extrinsic spectral function peaks beyond the plasma frequency,
!* shift the spectral function so it peaks at the plasma frequency.
!*      if(wemax.lt.omp) then
!      if(.false.) then
!        wshift=omp-wemax
!        ishift=53-iemax
!*       here, wpts(53)=omp
!        wshift2=wpts(iemax)-wemax
!        if (wshift.lt.0.d0) then
!          ishift=ishift-1
!          wshift2=wpts(iemax+1)-wemax
!        endif
!*        ishift=int(wshift/dw)
!*        wshift2=wshift-ishift*dw
!        do i=1,ishift
!          specttmp(i)=0.d0
!        enddo
!        do i=1,npts
!          wsearch=wpts(i)-wshift
!          if (wpts(1).gt.wsearch) then
!            specttmp(i)=0.d0
!          else
!            do j=1,npts
!              if (wpts(j).gt.wsearch) then
!                whi=wpts(j)
!                wlo=wpts(j-1)
!                esfhi=spectf(2,j)
!                esflo=spectf(2,j-1)
!                wshift2=wpts(j)-wsearch
!                specttmp(i)=esflo+(esfhi-esflo)*wshift2/(whi-wlo)
!                goto 100
!              endif
!            enddo
!          endif
!100       continue
!        enddo
!      endif


! Find the weights for the clipped parts of the extrinsic satellite
      wtmesat=0.d0
      wtmemain=0.d0
      do i=1,npts
        wtmesat=wtmesat+spectf(8,i)*expa*dw
        wtmemain=wtmemain+spectf(7,i)*expa*dw
      enddo
!* Weight the satellite terms by the power of the renormalization
!* constant reflecting how many powers of the extrinsic Green's function
!* appear in the original expressions
!      wtxsat=wtxsat*(zm+wtmemain)**2
!      wtisat=wtisat*(zm+wtmemain)
      do i=1,npts
!        spectf(4,i)=spectf(4,i)*(zm+wtmemain)**2
!        spectf(5,i)=spectf(5,i)*(zm+wtmemain)
        spectf(6,i)=spectf(2,i)-2.d0*spectf(4,i)+spectf(5,i)
      enddo
      
! eliminate regions of negative spectral weight, keep track of 
! the total satellite weight before correction and weight of negative
! regions removed for later renormalization.
      swtcorr=0.d0
      satwt=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)
        satwt=satwt+spectf(6,i)*dw1
        if (spectf(6,i).lt.0.d0) then
          swtcorr=swtcorr+spectf(6,i)*dw1
          spectf(6,i)=0.d0
          spectf(4,i)=(spectf(2,i)+spectf(5,i))/2.d0
        endif
      enddo
! renormalize to keep satellite weight the same.
      swtfac=satwt/(satwt-swtcorr)
      if (swtfac.lt.0.d0) swtfac=0.d0
      wtxsat=0.d0
      wtisat=0.d0
      wtesat=0.d0
      wtmesat=0.d0
      wtmemain=0.d0
      do i=1,npts
        w=wpts(i)+ekp
        dw1=wlim(i)-wlim(i-1)
        spectf(4,i)=(spectf(2,i)+spectf(5,i)-spectf(6,i)*swtfac)/2.d0
        wtxsat=wtxsat+spectf(4,i)*dw1*expa
        wtesat=wtesat+spectf(2,i)*dw1*expa
        wtisat=wtisat+spectf(5,i)*dw1*expa
        wtmesat=wtmesat+spectf(8,i)*expa*dw
        wtmemain=wtmemain+spectf(7,i)*expa*dw
      enddo
        
! Write out array of weights
      weights(1)=z1*expa
      write(66,*) z1,expa
      weights(2)=z1i*expa
      weights(3)=2.d0*z1*zm*ak*xreduc*expa
      weights(4)=wtesat
      weights(5)=wtxsat
      weights(6)=wtisat
      weights(7)=wtmesat
      weights(8)=wtmemain

!      write(6,*) xxc,yyc,discr,abs(discr),
!     2           -2.d0*width*(0,1)*yyc,zzc,abs(zzc),z1,z1i
!      z1=zzr
!      z1i=zzi
!      z1m=zzm

      if (iwrite.eq.jcount) then
        do i=1,npts
          write(12,500) wpts(i),spectf(1,i),spectf(3,i),                &
     &                  spectf(7,i),spectf(8,i)
          write(13,500) wpts(i),spectf(2,i),spectf(4,i),                &
     &                  spectf(5,i),                                    &
     &                  spectf(2,i)+spectf(5,i)-2.d0*spectf(4,i)
!     2                  spectf(5,i),spectf(6,i)
        enddo
      endif

      return
 500  format(1x,5(e12.5,1x))
 700  format(1x,7(f10.5,1x))
      end
