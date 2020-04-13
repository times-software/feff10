!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: so2conv.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2012/05/30 00:55:55 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program so2conv, by Luke Campbell (2002)
      subroutine so2conv  !KJ I stuffed all the calling args into modules 7-09

! This program computes many body effects on an XAS spectrum given in
! one or more FEFF output files of the form chi.dat,chipnnnn.dat,
! xmu.dat, or feffnnnn.dat.  The files read are from the most recent FEFF 
! from the directory in which this (sub) program is run, as determined
! by the files feff.inp and list.dat.  The files are overwritten
! with the computed spectra corrected for many body effects.
      use sfconv_inp
	  use errorfile
      implicit none
      integer npts,npts2,j,nsfpts,nqpts,maxfiles,ifiles,ifeffct
!      parameter(npts=400,npts2=401,nsfpts=80,nqpts=66,maxfiles=200)
      parameter(npts=112,npts2=401,nsfpts=112,nqpts=66,maxfiles=200)
!       npts - The number of points of the uniform energy grid used for 
!              the spectral function.
!       npts2 - The maximum number of data points in the data files
!       j - Used to store the number of data points in a data file.
!       nsfpts - Number of points of the minimal energy grid used for 
!              the spectral function.
!       nqpts - Number of points in the minimal momentum grid used for
!              the spectral function.
!       maxfiles - The maximum expected number of FEFF files to be read.
!       ifiles - The actual number of FEFF files specified.
!       ifeffct - Number of data points in a feffnnnn.dat file.
      integer i,ii,jj,ik,ikk,ipl
!       i,ii,jj,ik,ikk,ipl - counters for loops.
      integer last,lfirst,isearch,ifirst(maxfiles),ilast(maxfiles)
!       last - Collum number of last non-blank character in text string.
!       lfirst - Collum number of first non-blank character in text string.
!       isearch - Another counter for loops, used in text string searches.
!       ifirst - Collum number of start of file name.
!       ilast - Collum number of end of file name.
!KJ INPUT FROM MODULE      integer ispec,ipr6
!       ispec - Identifies type of spectroscopy.  0 for EXAFS, 1 for XANES.
!       ipr6 - Identifes which files need to be convoluted.
      integer iwrite,iwrite2,intout,iwrconv,iwrcfile,iout,ipw,ipwc  ! ,ipsk,ipse !KJ INPUT FROM MODULE
      double precision den,den2   !,cen !KJ INPUT FROM MODULE
!KJ INPUT FROM MODULE      character*12 cfname
!       iwrite - Write the spectral function and self energy 
!              to a file at this point on the minimal grid.
!       iwrite2 - Write the spectral function to a file at this point 
!              on the uniform grid.
!       intout - Set to one to write out the running integration of 
!              the convolution with the spectral function.
!       iwrconv - The data point at which the running integration is
!              to be written.
!       iwrcfile - The file for which the running integration is to be 
!              written.
!       ipw - Set greater than zero to write a file of integrated spectral
!              weights.
!       ipwc - Set greater than zero to write a file of integrated spectra
!              weights with the extrinsic spectral function separated into
!              quasiparticle and satellite terms.
!       ipse - Set greater to zero to write a file of the on shell
!              self energy.
!       ipsk - Set greater than zero to write the spectral function 
!              and self energy at a specified momentum.
!       iout - Flags whether the file specfunct.dat needs to
!              be created.
!       cen -  closest energy point to write out running convolution
!       den -  difference in gridpoint energy and cen (above)
!       cfname - name of file shown in running convolution (NULL if none)
!KJ INPUT FROM MODULE      double precision wsigk
!       wsigk - Momentum at which to write self energy file if 
!              ipsk specified.
      integer iedge,ios,itype(maxfiles),iasym,iasymt,isattype,isattypt, &
     &        icut,nleg,npi
!       iedge - Pointer to the Fermi energy in xmu.dat.
!       ios - Greater than 1 upon failure to open a file.
!       itype - Identifies type of file read.
!       iasym - Set to 1 to include quasiparticle phase as an
!              asymmetric 1/omega term to the extrinsic satellite
!              rather than as a complex spectral weight.  This is
!              necessary when convoluting with a real valued
!              function or one whose imaginary part is not known as
!              is the case with the xmu.dat file.
!       iasymt - Check if spectral function read from file uses
!              the same method of calculating the quasiparticle phase.
!       isattype - Indicates what approximations to use for 
!              satellite.
!       isattypt - Check if spectral function read from file uses
!              the same approximation.
!       icut - Set to 0 to avoid cutting off the spectral function at
!              energies where excitations are energetically forbidden
!              in the convolution.
!       nleg - Number of legs in the scattering path, not used for 
!              any calculations.
!       npi - Number of 2 pi jumps in phase.
      integer ipath,nlegs
!       Data from file list.dat.  Only ipath is used.
!       ipath - Path index.
      double precision sig2,ampratio,reff
!       Data from file list.dat.  Not used in this program except 
!       as dummy variables.
      character*80 cfirst,cinp,cblankl
      character*12 cfile(maxfiles)
      character*4 cpath
!       cfirst - Lines from FEFF files.
!       cinp - Lines from so2conv.inp.
!       cblankl - A long line of blanks.
!       cfile - The names of the FEFF files to read.
!       cpath - text containing path index
      character*20 cblanks
!       cblanks - A short line of blanks.
      character*10 vintx,mutx
!       Text containing the interstitial potential and the chemical
!       potential, respectively.
      character*9 gchtx,kftx
!       Text containing the core hole lifetime and the Fermi level,
!       respectively.
      character*5 rstx
!       Text containing the interstitial electron density parameter R_s.
      double precision aangstrom,eV
!       Conversion constants.
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      double precision dw,w,wp,dk,dk2,dp,delta
!       Various energy and momenta and their intervals.
      double precision rs,vint,deg,Rnn,gammach,conc,brpl1,brpl,cmu,     &
     &                 ckf,rst,gcht,brplt
!       rs - Electron density parameter R_s.
!       vint - Interstitial potential.
!       deg - Path degeneracy, number of identical scattering paths.
!       Rnn - Half length of scattering path.
!       gammach - Core hole lifetime.
!       conc - Electron concentration.
!       brpl1 - plasmon broadening in units of plasma frequency
!       brpl - plasmon broadening
!       cmu - Chemical potential from FEFF file.
!       ckf - Fermi level from FEFF file.
!       rst - Test to see if rs matches old value.
!       gcht - Test to see if gammach matches old value.
!       brplt - Test to see if brpl1 matches old value.
      double precision sef0,se0,se,ce,width,z1,z1i,se2,xise,            &
     &                 zkk,dsig,xpkg(npts2),seg(npts2),ekpg(npts2),     &
     &                 sse,sxise
!       sef0 - Self energy at the Fermi level Sigma(E_F,k_F)
!       se0 - Non-self consistent on shell self energy.
!       se - Real part of on shell self energy.
!       ce - Electron gas core hole self energy.
!       width - Quasiparticle lifetime.
!       z1 - Real part of renormalization constant.
!       z1i - Imaginary part of renormalization constant.
!       se2 - Real part of self energy (not neccessarily on shell).
!       xise - Imaginary part of self energy (not neccessarily on shell).
!       zkk - Momentum derivative renormalization constant.
!       dsig - Difference in self energies, to estimate zkk.
!       xpkg - Array of momenta, to estimate zkk.
!       seg - Array of self energies, to estimate zkk.
!       ekpg - Array of energies, to estimate zkk.
!       sse - self energy due to given pole in epsilon^{-1}
      common /energies/ se,ce,width,z1,z1i,se2,xise
      double precision xktest,xktold
!       Used for checking the location of the fermi level in xmu.dat.
      double precision xreduc,xfact1
!       xreduc - Ad hoc overall uniform amplitude reduction.  Obsolete.
!       xfact1 - Dummy variable for storing intermediate calculations.
      double precision pk(npts2),pgrid(nqpts),epts(nsfpts),wpts(npts),  &
     &                 pthresh
!       pk - Array of momentum values at data points from FEFF files.
!       pgrid - Array of momentum values on minimal momentum grid.
!       epts - Array of energy values on minimal energy grid.
!       wpts - Array of energy variables on uniform energy grid.
!       pthresh - Estimate of momentum threshold for plasmon creation.
      double precision spectf(8,nsfpts),cspec(npts),weights(8)
!       spectf - Spectral function computed on minimal grid.
!       cspec - Spectral function interpolated onto uniform grid.
!       Weights - Spectral weights of the components of the spectral function.
      double precision dphsum,s02sum,xnorm,ww
!       Used in finite element integration for averaging over
!            data points on uniform momentum grid to get output 
!            for momentum points in feffnnnn.dat file.
!       dphsum - Phase shift integration.
!       s02sum - Amplitude reduction integration.
!       xnorm - Normalization.
!       ww - Weighting or importance function.
      double precision chirr,chiii,so2mag,phaseshft,                    &
     &                 phshftold,phrmu,phrmu0,                          &
     &                 xchir,xchii,phchir,phchii,xmu2,xmu02,            &
     &                 phmu,phmu0,xchi2,phchi,rmu2,rmu02,ximu2
!       chirr - Real part of chi.
!       chiii - Imaginary part of chi.
!       so2mag - Amplitude reduction.
!       phaseshift - Phaseshift.
!       phaseshiftold - Used to check for jumps of 2 pi in phase.
!       xchir - Intermediate real part of chi.
!       xchii - Intermediate imaginary part of chi.
!       phchir - Phase of convolution of real part of chi.
!       phchii - Phase of convolution of imaginary part of chi.
!       xmu2 - Convolution of XANES spectrum with spectral function.
!       xmu02 - Convolution of background with spectral function.
!       xchi2 - Convolution of XAFS signal from xmu.dat.
!       phmu - 'Phase' of XANES spectrum, should be zero, but the 
!           subroutine needs a variable to put it in.
!       phmu0 - 'Phase' of background, should be zero, but the 
!           subroutine needs a variable to put it in.
!       phchi - 'Phase' of XAFS signal, should be zero, but the 
!           subroutine needs a variable to put it in.
      double precision emsf(nqpts,nsfpts),essf(nqpts,nsfpts),           &
     &                 xmsf(nqpts,nsfpts),xssf(nqpts,nsfpts),           &
     &                 xissf(nqpts,nsfpts),escsf(nqpts,nsfpts),         &
     &                 engrid(nqpts,nsfpts),wgts(nqpts,8),              &
     &                 sfinfo(nqpts,8)
!       emsf - Array of extrinsic quasiparticle spectral function
!            values for file specfunct.dat.
!       essf - Array of extrinsic satellite spectral function
!            values for file specfunct.dat.
!       xmsf - Array of interference quasiparticle spectral function
!            values for file specfunct.dat.
!       xssf - Array of interference satellite spectral function
!            values for file specfunct.dat.
!       xissf - Array of intrinsic satellite spectral function
!            values for file specfunct.dat.
!       essf - Array of clipped extrinsic satellite spectral function
!            values for file specfunct.dat.
!       engrid - Energy grid for spectral function in file specfunct.dat.
!       wgts - Array of spectral weights in file specfunct.dat.
!       sfinfo - Array of miscellaneous information in file specfunct.dat.
      double precision xk(npts2),chi(npts2),xmag(npts2),phase(npts2),   &
     &                 phm2kr(npts2),epts2(npts2),chir(npts2),          &
     &                 chii(npts2),e1(npts2),xmu(npts2),xmu0(npts2),    &
     &                 xk2(npts2),caph(npts2),xmfeff(npts2),            &
     &                 phfeff(npts2),redfac(npts2),xlam(npts2),         &
     &                 realck2(npts2),caph2(npts2),xmfeff2(npts2),      &
     &                 phfeff2(npts2),redfac2(npts2),xlam2(npts2)
      double precision rmu(npts2),rmu0(npts2),ximu(npts2)
!       xk - EXAFS wavenumber from chi.dat or chipnnnn.dat.
!       chi - EXAFS signal from chi.dat or chipnnnn.dat.
!       xmag - Magnitude of chi from chi.dat or chipnnnn.dat.
!       phase - Phase of chi from chi.dat or chipnnnn.dat.
!       phm2kr - PHase Minus 2 K R, phase of chi with dominant 
!           k dependant oscillation removed from chi.dat or chipnnnn.dat.
!       epts2 - Energy from files, measured from edge.
!       chir - Real part of EXAFS signal.
!       chii - Imaginary part of EXAFS signal.
!       e1 - Energy from xmu.dat.
!       xmu - XANES signal from xmu.dat.
!       xmu0 - Embedded atom background signal.
!       xk2 - EXAFS wavenumber from feffnnnn.dat
!       caph - Central atom phase shift interpolated onto uniform 
!           momentum grid.
!       xmfeff - Magnitude of f_(eff) interpolated onto uniform
!           momentum grid.
!       phfeff - Phase of f_(eff) interpolated onto uniform 
!           momentum grid.
!       redfac - Reduction factor on uniform momentum grid.
!       xlam - Mean free path of photoelectron, interpolated
!           onto uniform momentum grid.
!       realck2 - Real part of complex momentum from feffnnnn.dat.
!       caph2 - Central atom phase shift on feffnnnn.dat grid.
!       xmfeff2 - Magnitude of f_(eff) from feffnnnn.dat.
!       phfeff2 - Phase of f_(eff) from feffnnnn.dat.
!       redfac2 - Reduction factor feffnnnn.dat grid
!       xlam2 - Mean free path of photoelectron from feffnnnn.dat.
      double precision s02list(npts2),phlist(npts2)
!       s02list - Array of amplitude reduction values on uniform 
!           momentum grid.
!       phlist - Array of phase shifts on uniform momentum grid.
      double precision emain,esat,xmain,xsat,xisat,sat,eclip,esatr
!       emain - Extrinsic quasiparticle spectral function.
!       esat - Extrinsic satellite spectral function.
!       xmain - Interference quasiparticle spectral function.
!       xsat - Interference satellite spectral function.
!       xisat - Intrinsic spectral function.
!       sat - Total satellite spectral function.
!       eclip - Anomalous extrinsic satellite in region of quasiparticle.
!       esatr - Everything left in the extrinsic satellite after eclip is
!         subtracted off.
      integer npl,nplmax,nplt
      parameter (nplmax=5000)
!       npl - number of poles in epsilon^{-1}
!       nplmax - maximum number of poles in epsilon^{-1} for array dimensioning
      double precision plengy(nplmax),plwt(nplmax),oscstr(nplmax),      &
     &                 plbrd(nplmax),epswt,                             &
     &                 plengyt(nplmax),plwtt(nplmax),plbrdt(nplmax)
!       plengy - energy of poles in epsilon^{-1}
!       plwt - weight of poles in epsilon^{-1}, such that sum(plwt)=1
!       oscstr - oscilator strength f of resonances epsilon^{-1}
!       plbrd - broadening of poles in epsilon^{-1}
!       epswt - sum of weight of poles, to enforce sum(plwt)=1 {a consequence of
!               int(0,oo) d_omega omega epsilon^{-1}(q,omega)=-(pi/2)*omp}
      logical brpole
!       brpole - true if pole broadening is to be calculated
      double precision pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ef - Fermi energy
!       fmu - chemical potential = Fermi energy + self consistent
!             on shell self energy at the Fermi level
!       qf - Fermi momentum
!       omp - plasma frequency omega_p
!       ompl - energy of pole in epsilon^{-1}
!       wt - weight of pole in epsilon^{-1}
!       ekp - photoelectron energy = bare kinetic energy + real part of
!             on shell self energy
!       ek - bare photoelectron kinetic energy = pk**2/2
!       qpk - photoelectron momentum
!       acc - global accuracy parameter
!       brd - width of pole in epsilon^(-1)
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=omp+adisp*q**2+q**4/4
      common /convsf/ pi,ef,fmu,qf,omp,ompl,wt,ekp,ek,qpk,acc,brd,adisp
!       Common block one is used in most functions and subroutines.
      common /flag/ iwrite,jj
!       Common block flag tells subroutine  mkspectf when to print out
!       interesting information.
      double precision qthresh,beta,exchange
!       Functions which may be called.
      integer lowq,lowqt
!       lowq - Set equal to 0 to avoid calculating self energy from 
!              contributions below Fermi level.  When considering 
!              interference effects, lowq should be set equal to zero.
!              When finding quasiparticle properties, lowq should not
!              be zero.
!       lowqt - Value of lowq from previous run.
      common /belowqf/ lowq
      external qthresh,beta,exchange


      integer elnes,ipmin,ipmax,ipstep  !KJ added these variables 2-06
        character*12 f1,f2  !KJ dummy filenames  2-06
        integer kilast(10) !KJ stupid fix for stupid setup of this program 2-06
        integer ip !KJ local index 2-06
!     elnes : we're in EELS mode if elnes=1
!     ipmin,ipmax,ipstep : range of polarization components to calculate

!  !KJ Next section added to read ELNES variables 2-06     
!     read eels.inp
      elnes=0
      open(file='eels.inp',unit=3,status='old',err=900)
        read(3,*,err=900,end=900) 
        read(3,20,err=900,end=900) elnes
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,20,err=900,end=900) ipmin,ipstep,ipmax
  20  format (20i4)
      close(3)
      goto 901
900   continue
      elnes=0
901   continue
      if(elnes.eq.0) then
        ipstep=1
        ipmax=1
        ipmin=1
      endif
      do i=1,10
      kilast(i)=9
      enddo
      kilast(1)=7
!  !KJ end my changes



!  !KJ big loop for all polarization types 2-06 :
      do ip=ipmin,ipmax,ipstep
      if(ip.eq.1) then
          f1(1:12)='chi.dat     '
          f2(1:12)='xmu.dat     '
        elseif(ip.eq.10) then
          f1(1:12)='chi10.dat   '
          f2(1:12)='xmu10.dat   '
        elseif(ip.gt.1.and.ip.lt.10) then
          f1(1:4)='chi0'
          f1(5:5)= char(48+ip)
          f1(6:12)='.dat   '
          f2(1:4)='xmu0'
          f2(5:5)= char(48+ip)
          f2(6:12)='.dat   '
        else
          stop 'crazy ip in ff2xmu'
        endif
!  !KJ end my changes


      isattype=0
      iwrite=0
      iwrite2=0
      iwrconv=0
      iwrcfile=0
      ipw=0
      ipwc=0
!      ipse=0
!      ipsk=0
      lowq=0
      icut=0
      brpole=.true.
!      brpole=.false.
!      wsigk=1.d0
      if (abs(ispec).eq.1) then
        cfile(1)= f2  !KJ 2-06    'xmu.dat     '
        ilast(1)= kilast(ip)  !KJ 2-06     7
        itype(1)=2
        cfile(2)= f1  !KJ 2-06    'chi.dat     '
        ilast(2)=7
        itype(2)=1
        ifiles=2
      elseif (abs(ispec).eq.0) then
!       Josh - Added convolution of xmu.dat for path expansion
        cfile(1)= f2  !KJ 2-06    'xmu.dat     '
        ilast(1)= kilast(ip)  !KJ 2-06   7
        itype(1)=2
!       Josh END
        cfile(2)= f1  !KJ 2-06    'chi.dat     '
        ilast(2)= kilast(ip)  !KJ 2-06   7
        itype(2)=1
        ik=2
        if (ipr6.ge.2) then
          open(unit=1,file='list.dat',status='old',iostat=ios)
          if (ios.gt.0) goto 10
 2        read(1,'(A)') cinp
          if (cinp(6:14).ne.'---------') goto 2
          read(1,'(A)') cinp
 3          read(1,'(A)',end=10) cinp
            ik=ik+1
            lfirst=1
            last=1
 4          if (cinp(lfirst:lfirst).eq.' ') then
              lfirst=lfirst+1
              last=lfirst
              goto 4
            endif
 5          if (cinp(last:last).ne.' ') then
              last=last+1
              goto 5
            endif
            if (last-lfirst.eq.3) then
              cpath='0'//cinp(lfirst:last-1)
            elseif (last-lfirst.eq.2) then
              cpath='00'//cinp(lfirst:last-1)
            elseif (last-lfirst.eq.1) then
              cpath='000'//cinp(lfirst:last-1)
            else
              cpath=cinp(lfirst:last-1)
            endif
!            open(unit=23,file='temp',status='unknown')
!            close(unit=23)
!            open(unit=23,file='temp',status='old')
!            read(23,*) cpath
!            close(unit=23,status='delete')
!            last=4
! 4          if (cpath(last:last).eq.' ') then
!              last=last-1
!              goto 4
!            endif
!            if (last.eq.3) then
!              cpath='0'//cpath(:last)
!            elseif (last.eq.2) then
!              cpath='00'//cpath(:last)
!            elseif (last.eq.1) then
!              cpath='000'//cpath(:last)
!            endif
            if (ipr6.ge.2) then
              cfile(ik+1)='chip'//cpath//'.dat'
              itype(ik+1)=1
            endif
            if (ipr6.ge.3) then
              cfile(ik+1)='feff'//cpath//'.dat'
              itype(ik+1)=3
            endif
            ilast(ik+1)=12
          goto 3
 10       close(unit=1)
        endif
        ifiles=ik+1
      endif
      do i=1,ifiles
        if (cfile(i).eq.cfname) iwrcfile=i
      enddo

      if (ipsk.ne.0) then
        iwrite=0
        iwrite2=0
      elseif (iwrite.ne.0) then
        iwrite2=0
      endif
        
! Open diagnostic and auxilliary output files.
      if (iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0) then
        open(unit=12,status='unknown',file='qpsf.dat')
        open(unit=13,status='unknown',file='satsf.dat')
      endif
      if (ipsk.ne.0.or.iwrite.ne.0)                                     &
     &  open(unit=24,status='unknown',file='sigma.dat')
      if (ipwc.ne.0) open(unit=14,status='unknown',file='weightscl.dat')
      if (ipw.ne.0) open(unit=15,status='unknown',file='weights.dat')
      if (ipse.ne.0)                                                    &
     & open(unit=18,status='unknown',file='selfenergy.dat')

      do ikk=1,ifiles
! Open input and output files.
        open(unit=21,status='old',                                      &
     &       file=cfile(ikk)(:ilast(ikk)),iostat=ios)
        open(unit=16,status='unknown',                                  &
     &       file='mbheader.temp')
!        if (itype(ikk).eq.1) then
!          open(unit=17,status='unknown',
!     2       file=cfile(ikk)(:ilast(ikk)-4)//'so2.dat')
!        endif
        if (ios.gt.0) goto 640

 12     read(21,'(A)') cfirst
        last=80
! Trim off excess white space.
 13     if (cfirst(last:last).eq.' ') then
          last=last-1
          goto 13
        endif
! Josh - Check that convolution has not been done on this file.
        if(cfirst(1:last).eq.'# Convoluted with A(omega).') then
          print '(A)', 'WARNING: '// cfile(ikk) 
          print '(A)', 'has already been convoluted in a previous calculation.'
          print '(A)', 'Rerun module 6 if you still wish to proceed.'
		  call WipeErrorfileAtFinish !KJ this is a legitimate stop
          stop
        end if
! Find physical parameters of the material needed for computation.
        do isearch=1,last-7
          if(cfirst(isearch:isearch+6).eq.'Gam_ch=') then
            gchtx=cfirst(isearch+7:isearch+15)
          elseif(cfirst(isearch:isearch+6).eq.'Rs_int=') then
            rstx=cfirst(isearch+8:isearch+12)
          elseif(cfirst(isearch:isearch+4).eq.'Vint=') then
            vintx=cfirst(isearch+5:isearch+14)
          elseif(cfirst(isearch:isearch+2).eq.'Mu=') then
            mutx=cfirst(isearch+3:isearch+12)
          elseif(cfirst(isearch:isearch+2).eq.'kf=') then
            kftx=cfirst(isearch+3:isearch+11)
          endif
        enddo
! Write header lines.
        if ((iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0).and.ikk.eq.1)    &
     &  then
          write(12,'(A)') cfirst(1:last)
          write(13,'(A)') cfirst(1:last)
        endif
        if (ikk.eq.1.and.(ipsk.ne.0.or.iwrite.ne.0))                    &
     &    write(24,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipwc.ne.0) write(14,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipw.ne.0) write(15,'(A)') cfirst(1:last)
        write(16,'(A)') cfirst(1:last)
        if (ikk.eq.1.and.ipse.ne.0) write(18,'(A)') cfirst(1:last)
        if (cfirst(6:14).ne.'---------') goto 12
! feffnnnn.dat files have more information lines after the line with
! dashes than the other file types.
        if (itype(ikk).eq.3) then
 14       read(21,'(A)') cfirst
          last=80
 15        if (cfirst(last:last).eq.' ') then
            last=last-1
            goto 15
          endif
! If first character is a number, keep it.  Otherwise it might be a
! comment character, so dump it.
          if (cfirst(1:1).eq.'1'.or.cfirst(1:1).eq.'2'                  &
     &        .or.cfirst(1:1).eq.'3'.or.cfirst(1:1).eq.'4'              &
     &        .or.cfirst(1:1).eq.'5'.or.cfirst(1:1).eq.'6'              &
     &        .or.cfirst(1:1).eq.'7'.or.cfirst(1:1).eq.'8'              &
     &        .or.cfirst(1:1).eq.'9'.or.cfirst(1:1).eq.'0') then
            lfirst=1
          else
            lfirst=2
          endif
! Search for scattering half length, read into variable Rnn
          do isearch=1,last
            if(cfirst(isearch:isearch+3).eq.'reff') then
              open(unit=23,file='temp',status='unknown')
              write(23,*) cfirst(lfirst:last)
              close(unit=23)
              open(unit=23,file='temp',status='old')
              read(23,*) nleg,deg,Rnn
              close(unit=23,status='delete')
              Rnn=Rnn*aangstrom
! @# indicates end of comment lines.
            elseif(cfirst(isearch:isearch+1).eq.'@#') then
              goto 16
            endif
          enddo
          write(16,'(A)') cfirst(1:last)
          goto 14
        else 
          read(21,'(A)') cfirst
        endif
        last=80
 16     if (cfirst(last:last).eq.' ') then
          last=last-1
          goto 16
        endif
! Write column labels.
        if (ikk.eq.1) then
          if (iwrite.ne.0.or.iwrite2.ne.0.or.ipsk.ne.0) then
            write(12,'(A)') '#  omega       delta ext    interf.'       &
     &       // '      main ext     sat ext'
            write(13,'(A)') '#  omega       extrinsic    interf.'       &
     &       // '      intrinsic    total'
          endif
          if (iwrite.ne.0.or.ipsk.ne.0)                                 &
     &    write(24,'(A)') '#   omega      Re(Sigma)    -Im(Sigma)'
          if (ipwc.ne.0)                                                &
     &    write(14,'(A)') '#   energy     ext q.p.   inter q.p. ext sat'&
     &     // '    inter sat  intrin sat'
          if (ipw.ne.0)                                                 &
     &    write(15,'(A)') '#   energy     ext q.p.   inter q.p. ext sat'&
     &     // '    inter sat  intrin sat'
          if (ipse.ne.0)                                                &
     &    write(18,'(A)') '#   energy      Re(Sigma)   -Im(Sigma)'      &
     &     //'        Re(Z)        Im(Z)'
        endif
        write(16,'(A)') cfirst(1:last)
!        if (itype(ikk).eq.1) then
!          write(17,'(A)') '#     k          So2        phase shift'
!        endif
! Find numerical values of core hole lifetime, Wigner-Seitz radius, 
! interstitial potential, chemical potential, and Fermi momentum.
        open(unit=23,file='temp',status='unknown')
        write(23,*) gchtx
        write(23,*) rstx
        write(23,*) vintx
        write(23,*) mutx
        write(23,*) kftx
        close(unit=23)
        open(unit=23,file='temp',status='old')
        read(23,*) gammach
        read(23,*) rs
        read(23,*) vint
        read(23,*) cmu
        read(23,*) ckf
        close(unit=23,status='delete')
        gammach=(gammach/2.d0)*eV
        cmu=(cmu-vint)*eV
        vint=vint*eV
        ckf=ckf/aangstrom
        xreduc=1.d0
! Compute important material properties and constants.
        pi=dacos(-1.d0)
        qf=((9.d0*pi/4.d0)**(1.d0/3.d0))/rs
        ef=qf*qf/2.d0
        conc=3.d0/(4.d0*pi*(rs**3))
        omp=dsqrt(4.d0*pi*conc)
        adisp=2.d0*ef/3.d0
        ekp=ef
        qpk=qf
        acc=1.d-4

! Read pole expansion for epsilon^{-1}
        call rdeps(omp,nplmax,npl,plengy,oscstr,plbrd)
        epswt=0.d0
        do ipl=1,npl
          plwt(ipl)=dabs(oscstr(ipl)*plengy(ipl)**2/(omp**2))
          epswt=epswt+plwt(ipl)
        enddo
        open (unit=88,file='apl.dat', status='unknown')
        do ipl=1,npl
          write(88,'(5f10.5)') plengy(ipl)/eV,oscstr(ipl)*plengy(ipl)
        enddo
        close (88)

        sef0=0.d0
        do ipl=1,npl
          call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
          if (brpole) then
            call brsigma(0.d0,sse,sxise)
          else 
            call renergies(0.d0,sse)
          endif
          sef0=sef0+sse*wt
        enddo
        sef0=sef0+exchange(qf)

! Estimate threshold for plasmon creation.
        pthresh=qthresh(omp,adisp,ef,qf)
! Find minimal momentum grid.
        do i=1,10
          dp=(pthresh-qf)/10.d0
          pgrid(i)=qf+i*dp
        enddo
        do i=1,30
          dp=0.25d0*pthresh/30.d0
          pgrid(i+10)=pgrid(10)+i*dp
        enddo
        do i=1,10
          dp=0.75d0*pthresh/10.d0
          pgrid(i+40)=pgrid(40)+i*dp
        enddo
        do i=1,10
          dp=2.d0*pthresh/10.d0
          pgrid(i+50)=pgrid(50)+i*dp
        enddo
        pgrid(61)=5.d0*pthresh
        pgrid(62)=7.d0*pthresh
        pgrid(63)=1.d1*pthresh
        pgrid(64)=3.d1*pthresh
        pgrid(65)=1.d2*pthresh
        pgrid(66)=3.d2*pthresh
  
! Read information from files.
        if (itype(ikk).eq.1) then
          iasym=0
          do i=1,npts2
            epts2(i)=0.d0
            xk(i)=0.d0
            pk(i)=0.d0
            chi(i)=0.d0
            xmag(i)=0.d0
            phase(i)=0.d0
            phm2kr(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=25) xk(i),chi(i),xmag(i),phase(i)
            xk(i)=xk(i)/aangstrom
            chir(i)=xmag(i)*cos(phase(i))
            chii(i)=xmag(i)*sin(phase(i))
            j=i
          enddo
        elseif (itype(ikk).eq.2) then
          iasym=1
          do i=1,npts2
            e1(i)=0.d0
            epts2(i)=0.d0
            xk(i)=0.d0
            pk(i)=0.d0
            xmu(i)=0.d0
            xmu0(i)=0.d0
            chi(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=25) e1(i),epts2(i),xk(i),xmu(i),xmu0(i),chi(i)
            e1(i)=e1(i)*eV
            epts2(i)=epts2(i)*eV
            xktest=dabs(xk(i))
            if (xktest.lt.xktold) iedge=i
            xktold=xktest
            xk(i)=xk(i)/aangstrom
            j=i
          enddo
          dw=epts2(j)-epts2(j-1)
        elseif (itype(ikk).eq.3) then
          iasym=0
          do i=1,npts2
            epts2(i)=0.d0
            xk(i)=0.05d0*(i-1)/aangstrom
            pk(i)=0.d0
            chi(i)=0.d0
            xmag(i)=0.d0
            phase(i)=0.d0
            phm2kr(i)=0.d0
          enddo
          xktold=0.d0
          do i=1,npts2
            read(21,*,end=17) xk2(i),caph2(i),xmfeff2(i),phfeff2(i),    &
     &          redfac2(i),xlam2(i),realck2(i)
            xk2(i)=xk2(i)/aangstrom
            xmfeff2(i)=xmfeff2(i)*aangstrom
            xlam2(i)=xlam2(i)*aangstrom
            realck2(i)=realck2(i)/aangstrom
            j=i
          enddo
 17       ifeffct=j
! Interpolate data from feffnnnn.dat grid to the uniform grid used
! in chi.dat.
          do i=1,npts2
            do jj=1,j-1
              if (xk(i).ge.xk2(jj).and.xk(i).lt.xk2(jj+1)) then
                caph(i)=caph2(jj)+(caph2(jj+1)-caph2(jj))               &
     &                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                xmfeff(i)=xmfeff2(jj)+(xmfeff2(jj+1)-xmfeff2(jj))       &
     &                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                phfeff(i)=phfeff2(jj)+(phfeff2(jj+1)-phfeff2(jj))       &
     &                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                redfac(i)=redfac2(jj)+(redfac2(jj+1)-redfac2(jj))       &
     &                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                xlam(i)=xlam2(jj)+(xlam2(jj+1)-xlam2(jj))               &
     &                  *(xk(i)-xk2(jj))/(xk2(jj+1)-xk2(jj))
                goto 18
              endif
            enddo
            if (xk(i).eq.xk2(jj)) then
              caph(i)=caph2(jj)
              xmfeff(i)=xmfeff2(jj)
              phfeff(i)=phfeff2(jj)
              redfac(i)=redfac2(jj)
              xlam(i)=xlam2(jj)
            endif
! Compute exafs signal Chi from feffnnnn.dat info.
 18         xmag(i)=(deg*xmfeff(i)*redfac(i)*exp(-2.d0*Rnn/xlam(i)))    &
     &              /(xk(i)*Rnn**2)
            phm2kr(i)=phfeff(i)+caph(i)
            phase(i)=phm2kr(i)+2.d0*xk(i)*Rnn
            chir(i)=xmag(i)*cos(phase(i))
            chii(i)=xmag(i)*sin(phase(i))
          enddo
          xmag(1)=xmag(2)+(xk(1)-xk(2))*(xmag(3)-xmag(2))/(xk(3)-xk(2))
          chir(1)=xmag(1)*cos(phase(1))
          chii(1)=xmag(1)*sin(phase(1))
          j=npts2
        endif
 25     continue

        close(unit=21)
        close(unit=16)
        open(unit=21,status='old',file='mbheader.temp',iostat=ios)
        if (ios.gt.0) goto 640
        open(unit=16,status='unknown',                                  &
     &       file=cfile(ikk)(:ilast(ikk)))
!	Josh - added to header so that we can check if a previous convolution
!	has been performed.
        write(16,'(A)') '# Convoluted with A(omega).'
 27       read(21,'(A)',end=35) cfirst
          last=80
 30       if (cfirst(last:last).eq.' ') then
            last=last-1
            goto 30
          endif
          write(16,'(A)') cfirst(1:last)
          goto 27
 35     close(unit=21,status='delete')
          
! Need to find photoelectron momentum.
! First, find kinetic energy.
        do i=1,j
          if (xk(i).ge.0.d0) then
            ekp=xk(i)**2/2.d0+cmu
          else
            ekp=-xk(i)**2/2.d0+cmu
          endif
          ekpg(i)=ekp
          if (itype(ikk).eq.1.or.itype(ikk).eq.3) epts2(i)=ekp
! Make arrays of 0 order self energies and 0 order estimates for momentum.
          if (ekp.ge.0.d0) then
            qpk=sqrt(qf**2+2.d0*(ekp-fmu))
            xpkg(i)=qpk
            se0=0.d0
            do ipl=1,npl
              call plset(ipl,nplmax,plengy,plwt,plbrd,ompl,wt,brd)
              if (brpole) then
                call brsigma(0.d0,sse,sxise)
              else
                call renergies(0.d0,sse)
              endif
              se0=se0+sse*wt
            enddo
            seg(i)=se0+exchange(qpk)
          endif
        enddo
! Estimate momentum derivative renormalization constant.
        do i=1,j
          if (ekpg(i).ge.0) then
            if (i.ne.1.and.i.ne.j) then
              dsig=seg(i+1)-seg(i-1)
              dk2=xpkg(i+1)**2/2.d0-xpkg(i-1)**2/2.d0
            elseif (i.eq.1) then
              dsig=seg(i+1)-seg(i)
              dk2=xpkg(i+1)**2/2.d0-xpkg(i)**2/2.d0
            else
              dsig=seg(i)-seg(i-1)
              dk2=xpkg(i)**2/2.d0-xpkg(i-1)**2/2.d0
            endif
            zkk=1.d0/(1.d0+dsig/dk2)
! Find array of photoelectrom momenta.
            pk(i)=sqrt(xpkg(i)**2-2.d0*zkk*(seg(i)-sef0))
          endif
        enddo
! if running convolution is written, choose momentum point closest to desired energy
        if (iwrcfile.ne.0) then
          den=abs(ekpg(1)-cen)
          iwrconv=1
          do i=2,j
            den2=abs(ekpg(i)-cen)
            if (den2.lt.den) then
              den=den2
              iwrconv=i
            endif
          enddo
        endif
  
! Pad out the ends of the arrays so the convolution does not screw up the
! end data points.
        if (itype(ikk).eq.1.or.itype(ikk).eq.3) then
          dw=epts2(j)-epts2(j-1)
          do i=j+1,npts2
            epts2(i)=epts2(i-1)+dw
          enddo
        elseif (itype(ikk).eq.2) then
          dw=epts2(j)-epts2(j-1)
          do i=j,npts2
            e1(i)=e1(i-1)+dw
            epts2(i)=epts2(i-1)+dw
            xmu0(i)=xmu0(j)
            xmu(i)=xmu0(j)
            chi(i)=0.d0
          enddo
          call mkrmu(xmu,xmu0,rmu,epts2,npts2)
!          call mkrmu(xmu0,xmu0,rmu0,epts2,npts2)
          do i=1,npts2
            ximu(i)=xmu(i)-xmu0(i)
          enddo
        endif

! Find the spectral function and associated data.
        iout=0
! If specfunct.dat does not exist, or if it is for a material with 
! different electron gas properties, set flag 'iout' to recompute 
! spectral function.  Otherwise, read in spectral function.
        open(unit=23,file='specfunct.dat',status='old',                 &
     &       access='sequential',form='unformatted',iostat=ios)
        if (ios.gt.0) then
          iout=1
        else
          read(23) rst,gcht,iasymt,isattypt,lowqt,nplt
          read(23) (plengyt(ii), ii=1,nplmax)
          read(23) (plbrdt(ii), ii=1,nplmax)
          read(23) (plwtt(ii), ii=1,nplmax)
          read(23) ((sfinfo(jj,ii), jj=1,nqpts), ii=1,8)
          read(23) ((wgts(jj,ii), jj=1,nqpts), ii=1,8)
          read(23) ((emsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((essf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xmsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xssf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((xissf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((escsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
          read(23) ((engrid(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        endif
        close(unit=23)
        if (rst.ne.rs) iout=1
        if (gcht.ne.gammach) iout=1
        if (iasymt.ne.iasym) iout=1
        if (lowqt.ne.lowq) iout=1
        if (isattypt.ne.isattype) iout=1
        if (nplt.ne.npl) iout=1
        do jj=1,npl
          if (plengy(jj).ne.plengyt(jj)) iout=1
          if (plbrd(jj).ne.plbrdt(jj)) iout=1
          if (plwt(jj).ne.plwtt(jj)) iout=1
        enddo
        do jj=1,nqpts
          if (real(sfinfo(jj,1)).ne.real(pgrid(jj))) iout=1
        enddo
        if (iwrite.ne.0.and.ikk.eq.1) iout=1

        npi=0
        if (iout.eq.1) then
          write(6,*) 'computing spectral function'
        endif
        if(ipsk.ne.0.and.ikk.eq.1) then
          jj=0
          !PRINT*, "Line 949."
          call mkspectf(rs,wsigk*qf,gammach,xreduc,                     &
     &                  epts,spectf,weights,isattype,                   &
     &                  npl,nplmax,plengy,plwt,plbrd,brpole)
          !PRINT*, "Line 953."

!          do ipl=1,npl
!            write(6,*) ipl,plengy(ipl),plbrd(ipl),plwt(ipl),
!     2        plengy(ipl)+sef0-se,sef0,se
!          enddo
  
        endif
        do jj=1,nqpts
          if (iout.eq.1) then
! Spectral function has not yet been computed, or was computed for the
! wrong material.  Compute spectral function.
            call mkspectf(rs,pgrid(jj),gammach,xreduc,                  &
     &                    epts,spectf,weights,isattype,                 &
     &                  npl,nplmax,plengy,plwt,plbrd,brpole)
            !write(6,'(I3,A10)') int(real(jj)/real(nqpts)*100.0),        &
   !  &                '% computed'
            sfinfo(jj,1)=pgrid(jj)
            sfinfo(jj,2)=ekp
            sfinfo(jj,3)=ek
            sfinfo(jj,4)=se
            sfinfo(jj,5)=ce
            sfinfo(jj,6)=width
            sfinfo(jj,7)=z1
            sfinfo(jj,8)=z1i
            do ii=1,8
              wgts(jj,ii)=weights(ii)
            enddo
            do ii=1,nsfpts
              emsf(jj,ii)=spectf(1,ii)
              essf(jj,ii)=spectf(2,ii)
              xmsf(jj,ii)=spectf(3,ii)
              xssf(jj,ii)=spectf(4,ii)
              xissf(jj,ii)=spectf(5,ii)
              escsf(jj,ii)=spectf(8,ii)
              engrid(jj,ii)=epts(ii)
            enddo
          else
! Correct spectral function has already been computed.
            ekp=sfinfo(jj,2)
            ek=sfinfo(jj,3)
            se=sfinfo(jj,4)
            ce=sfinfo(jj,5)
            width=sfinfo(jj,6)
            z1=sfinfo(jj,7)
            z1i=sfinfo(jj,8)
            do ii=1,8
              weights(ii)=wgts(jj,ii)
            enddo
          endif
        enddo
        write(6,*) 'convoluting file '                                  &
     &           //cfile(ikk)(:ilast(ikk))
! Interpolate spectral function onto uniform momentum grid.
        do jj=1,j
          do i=1,nqpts-1
            if (pk(jj).ge.pgrid(i).and.pk(jj).lt.pgrid(i+1)) then
              delta=(pk(jj)-pgrid(i))/(pgrid(i+1)-pgrid(i))
              do ii=1,nsfpts
                epts(ii)=engrid(i,ii)                                   &
     &               +delta*(engrid(i+1,ii)-engrid(i,ii))
                spectf(1,ii)=emsf(i,ii)+delta*(emsf(i+1,ii)-emsf(i,ii))
                spectf(2,ii)=essf(i,ii)+delta*(essf(i+1,ii)-essf(i,ii))
                spectf(3,ii)=xmsf(i,ii)+delta*(xmsf(i+1,ii)-xmsf(i,ii))
                spectf(4,ii)=xssf(i,ii)+delta*(xssf(i+1,ii)-xssf(i,ii)) 
                spectf(5,ii)=xissf(i,ii)                                &
     &               +delta*(xissf(i+1,ii)-xissf(i,ii))
                spectf(6,ii)=essf(i,ii)+xissf(i,ii)-2.d0*xssf(i,ii)     &
     &               +delta*(essf(i+1,ii)+xissf(i+1,ii)                 &
     &                     -2.d0*xssf(i+1,ii)                           &
     &                     -(essf(i,ii)+xissf(i,ii)-2.d0*xssf(i,ii)))
                spectf(7,ii)=essf(i,ii)-escsf(i,ii)                     &
     &               +delta*(essf(i+1,ii)-escsf(i+1,ii)                 &
     &                     -(essf(i,ii)-escsf(i,ii)))
                spectf(8,ii)=escsf(i,ii)                                &
     &               +delta*(escsf(i+1,ii)-escsf(i,ii))
              enddo
              do ii=1,8
                weights(ii)=wgts(i,ii)+delta*(wgts(i+1,ii)-wgts(i,ii))
              enddo
              se=sfinfo(i,4)+delta*(sfinfo(i+1,4)-sfinfo(i,4))
              ce=sfinfo(i,5)+delta*(sfinfo(i+1,5)-sfinfo(i,5))
              width=sfinfo(i,6)+delta*(sfinfo(i+1,6)-sfinfo(i,6))
              z1=sfinfo(i,7)+delta*(sfinfo(i+1,7)-sfinfo(i,7))
              z1i=sfinfo(i,8)+delta*(sfinfo(i+1,8)-sfinfo(i,8))
            endif
          enddo
! Take special care with the endpoints.
          if (pk(jj).ge.pgrid(nqpts)) then
            do ii=1,nsfpts
              epts(ii)=engrid(nqpts,ii)
              spectf(1,ii)=emsf(nqpts,ii)
              spectf(2,ii)=essf(nqpts,ii)
              spectf(3,ii)=xmsf(nqpts,ii)
              spectf(4,ii)=xssf(nqpts,ii)
              spectf(5,ii)=xissf(nqpts,ii)
              spectf(6,ii)=essf(nqpts,ii)+xissf(nqpts,ii)               &
     &                     -2.d0*xssf(nqpts,ii)
              spectf(7,ii)=essf(nqpts,ii)-escsf(nqpts,ii)
              spectf(8,ii)=escsf(nqpts,ii)
            enddo
            do ii=1,8
              weights(ii)=wgts(nqpts,ii)
            enddo
            se=sfinfo(nqpts,4)
            ce=sfinfo(nqpts,5)
            width=sfinfo(nqpts,6)
            z1=sfinfo(nqpts,7)
            z1i=sfinfo(nqpts,8)
          elseif (pk(jj).lt.pgrid(1)) then
            do ii=1,nsfpts
              epts(ii)=engrid(nqpts,ii)
              spectf(1,ii)=emsf(1,ii)
              spectf(2,ii)=essf(1,ii)
              spectf(3,ii)=xmsf(1,ii)
              spectf(4,ii)=xssf(1,ii)
              spectf(5,ii)=xissf(1,ii)
              spectf(6,ii)=essf(1,ii)+xissf(1,ii)-2.d0*xssf(1,ii)
              spectf(7,ii)=essf(1,ii)-escsf(1,ii)
              spectf(8,ii)=escsf(1,ii)
            enddo
            do ii=1,8
              weights(ii)=wgts(1,ii)
            enddo
            se=sfinfo(1,4)
            ce=sfinfo(1,5)
            width=sfinfo(1,6)
            z1=sfinfo(1,7)
            z1i=sfinfo(1,8)
          endif
! Interpolate spectral function onto uniform energy grid.
          call interpsf(npts,epts,wpts,spectf,cspec)
! Write information to diagnostic and supplementary files.
          if (ikk.eq.1.and.ipwc.ne.0) then
            write(14,701) epts2(jj),weights(1)+weights(8),              &
     &                  weights(3),weights(7),weights(5),weights(6)
          endif
          if (ikk.eq.1.and.ipw.ne.0) then
            write(15,701) epts2(jj),weights(1),weights(3)/2.d0,         &
     &                  weights(4),weights(5),weights(6)
          endif
          if (ikk.eq.1.and.ipse.ne.0) then
            write(18,500) epts2(jj),se,width,z1,z1i
          endif
          do i=1,nsfpts
            w=epts(i)
            emain=spectf(1,i)
            esat =spectf(2,i)
            xmain=spectf(3,i)
            xsat =spectf(4,i)
            xisat =spectf(5,i)
            sat  =spectf(6,i)
            eclip=spectf(7,i)
            esatr=spectf(8,i)
            if (jj.eq.iwrite2.and.ikk.eq.1) then
              write(12,500) w,emain,xmain,eclip,esatr
              write(13,500) w,esat,-2.d0*xsat,xisat,                    &
     &                      (esat+xisat-2.d0*xsat)
            endif
          enddo
          wp=epts2(jj)
! Begin convolution.
          if (itype(ikk).eq.1.or.itype(ikk).eq.3) then
            if (jj.eq.iwrconv.and.ikk.eq.iwrcfile) then
              intout=1
              open (unit=28,file='realint.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu,gammach,npts2,epts2,chir,npts,           &
     &                wpts,cspec,weights,xchir,phchir,iasym,            &
     &                1,intout,omp)
            if (jj.eq.iwrconv.and.ikk.eq.iwrcfile) then
              intout=1
              open (unit=28,file='imagint.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu,gammach,npts2,epts2,chii,npts,           &
     &                wpts,cspec,weights,xchii,phchii,iasym,            &
     &                1,intout,omp)
! Find real and imaginary many body EXAFS signal chi.
            chirr=xchir*cos(phchir)-xchii*sin(phchii)
            chiii=xchii*cos(phchii)+xchir*sin(phchir)
! Compute phase shift, remove jumps of 2 pi.
            phaseshft=datan2(chiii,chirr)
            if(phaseshft-phshftold.gt.5.d0) then
              npi=npi+2
            elseif(phaseshft-phshftold.lt.-5.d0) then
              npi=npi-2
            endif
            phshftold=phaseshft
            phaseshft=phaseshft-pi*npi
            if (itype(ikk).eq.1) then
! Write output file of many body EXAFS signal.
              write(16,630) xk(jj)*aangstrom,                           &
     &                chiii,dsqrt(chirr**2+chiii**2),                   &
     &                phaseshft,                                        &
     &                phaseshft+phm2kr(jj)-phase(jj)
            endif
 630        format (1x, f10.4, 3x, 4(1pe13.6,1x))
! Find amplitude reduction and phase shift.
            so2mag=dsqrt(chirr**2+chiii**2)/xmag(jj)
            phaseshft=datan2(chiii,chirr)
            if(phaseshft-phshftold.gt.5.d0) then
              npi=npi+2
            elseif(phaseshft-phshftold.lt.-5.d0) then
              npi=npi-2
            endif
            phshftold=phaseshft
            phaseshft=phaseshft-pi*npi-phase(jj)
            s02list(jj)=so2mag
            phlist(jj)=phaseshft
!            if (itype(ikk).eq.1) then
!* Write output file of amplitude reduction and phase shift.
!              write(17,500) xk(jj)*aangstrom,
!     2                so2mag,phaseshft
!            endif
          elseif (itype(ikk).eq.2) then
! Find many body convolution on absorption signal, embedded atom background,
! and fine structure.
!            call sfconv(wp,cmu+vint,gammach,npts2,epts2,chi,npts,
!     2                wpts,cspec,weights,xchi2,phchi,iasym,1,0,omp)
            if (iasym.eq.0) then
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='imagint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,ximu,npts,    &
     &            wpts,cspec,weights,ximu2,phmu,iasym,icut,intout,omp)
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='realint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,rmu,npts,     &
     &            wpts,cspec,weights,rmu2,phrmu,iasym,icut,intout,omp)
            else
              if (jj.eq.iwrconv) then
                intout=1
                open (unit=28,file='realint.dat',status='unknown')
              else
                intout=0
              endif
              call sfconv(wp,cmu+vint,gammach,npts2,epts2,xmu,npts,     &
     &            wpts,cspec,weights,xmu2,phrmu,iasym,icut,intout,omp)
            endif
            if (jj.eq.iwrconv) then
              intout=1
              open (unit=28,file='xmu0int.dat',status='unknown')
            else
              intout=0
            endif
            call sfconv(wp,cmu+vint,gammach,npts2,epts2,xmu0,npts,      &
     &            wpts,cspec,weights,xmu02,phmu0,iasym,icut,intout,omp)
            if (iasym.eq.0) then
              xmu2=ximu2*cos(phmu)+rmu2*sin(phrmu)+xmu02
            endif
            write(16,800)                                               &
     &                e1(jj)/eV,epts2(jj)/eV,xk(jj)*aangstrom,          &
!     2                xmu2,xmu02,(xmu2-xmu02)/xmu02    !KJ this is the original instruction  3-06
     &                xmu2,xmu02,(xmu2-xmu02)   !KJ removed normalization of chi   3-06  
          endif
        enddo
        if (itype(ikk).eq.3) then
          dk=0.05d0
! Average over nearby points on uniform momentum grid to find
! values for coarser grid in feffnnnn.dat files.  Use a trianguar
! weighting function to ensure each point contributes its full weight
! to the average.
          do jj=1,ifeffct
            s02sum=0.d0
            dphsum=0.d0
            xnorm=0.d0
            do i=1,j
              if (xk2(jj).eq.xk(i)) then
                ww=1.d0
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              elseif (xk(i).gt.xk2(jj-1).and.xk(i).le.xk2(jj)           &
     &                .and.xk2(jj-1).ne.xk2(jj)) then
                ww=(xk(i)-xk2(jj-1))/(xk2(jj)-xk2(jj-1))
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              elseif(xk(i).gt.xk2(jj).and.xk(i).lt.xk2(jj+1)            &
     &                .and.xk2(jj+1).ne.xk2(jj)) then
                ww=(xk2(jj+1)-xk(i))/(xk2(jj+1)-xk2(jj))
                s02sum=s02sum+s02list(i)*ww*dk
                dphsum=dphsum+phlist(i)*ww*dk
                xnorm=xnorm+ww*dk
              else
                ww=0.d0
              endif
            enddo
! Find many body values of phase and reduction.
            redfac2(jj)=redfac2(jj)*s02sum/xnorm
            caph2(jj)=caph2(jj)+dphsum/xnorm
! Write convoluted feffnnnnc.dat.
            write(16,400) xk2(jj)*aangstrom,caph2(jj),                  &
     &          xmfeff2(jj)/aangstrom,phfeff2(jj),redfac2(jj),          &
     &          xlam2(jj)/aangstrom,realck2(jj)*aangstrom
  400       format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,             &
     &                            2(1pe11.4,1x))
          enddo
        endif
! Write out the spectral function to a data file.
        open(unit=23,file='specfunct.dat',status='unknown',             &
     &       access='sequential',form='unformatted',iostat=ios)
        write(23) rs,gammach,iasym,isattype,lowq,npl
        write(23) (plengy(ii), ii=1,nplmax)
        write(23) (plbrd(ii), ii=1,nplmax)
        write(23) (plwt(ii), ii=1,nplmax)
        write(23) ((sfinfo(jj,ii), jj=1,nqpts), ii=1,8)
        write(23) ((wgts(jj,ii), jj=1,nqpts), ii=1,8)
        write(23) ((emsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((essf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xmsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xssf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((xissf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((escsf(jj,ii), jj=1,nqpts), ii=1,nsfpts)
        write(23) ((engrid(jj,ii), jj=1,nqpts), ii=1,nsfpts)
! Do a little housecleaning.
        close(unit=23)
        close(unit=21)
        close(unit=16)
        close(unit=17)
 640    continue
      enddo
 500  format(1x,5(e12.5,1x))
 700  format(1x,7(f10.5,1x))
 701  format(1x,e10.5,1x,6(f10.5,1x))
 800  format(1x,2f11.3,f8.3,1p,3e13.5)
! End of program.


      enddo  !KJ end of loop do ip=ipmin,ipmax,ipstep  2-06


!     sub-program so2conv
!     stop
      return

      end
