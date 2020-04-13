!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: csigma.f90,v $:
! $Revision: 1.8 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CSigma(Energy, Mu, Rs, ReSig, ImSig, WpScl, Gamma, AmpFac)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Written by Josh Kas
!     This subroutine calculates the self energy Sigma(k(Energy),Energy)
!     based on an electron gas model of epsilon^-1.
!      
!     Solve: k0**2 = 2*Energy - 2*Mu -
!                    2*(Sigma(k0,Energy)-Sigma(kFermi,EFermi))
!            
!     Steps:
!
!            1. k0**2  = 2*(Energy-Mu) + SigmaF (SigmaF is self energy at fermi level).
!            2. Sigma0 = Sigma(k0,Energy)
!            3. k1**2  = 
!                  k0**2 - 2*(Sigma0-SigmaF)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
	  use dimsmod,only: MxPole !KJ 7-09 replaced parameter statement
!     Parameters:
!     MxPole - Maximum number of poles      
!KJ      INTEGER MxPole
!KJ      PARAMETER(MxPole=1000)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Input:
!     Energy - Energy at which to evaluate Sigma
!     Mu     - Fermi energy as calculated by FEFF
!     Rs     - R sub s (sphere of radius Rs holds charge e)
!     WpScl  - Scale Wp in interstitial by WpScl
!     Gamma  - Use broadening Gamma when calculating Sigma
!     AmpFac - Use amplitude AmpFac for plasmon pole.
!     UseBP  - If true, use broadened pole SE.
!     Note: Atomic units are used.
      DOUBLE PRECISION Rs, WpScl(MxPole), Gamma(MxPole), AmpFac(MxPole),&
     &     Mu
      COMPLEX*16 Energy
      LOGICAL UseBP
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output:
!     ReSig  - Re[Sigma(Energy,k(Energy))]
!     ImSig  - Im[Sigma(Energy,k(Energy))]
      DOUBLE PRECISION ReSig, ImSig
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables:
!     kFermi - Fermi momentum calculated from Rs via electron gas
!              approximation.
!     EFermi - Fermi energy = kFermi^2/2
!     Wp     - Electron gas plasmon frequency.
!     Gam    - broadening for broadened pole.
!     ckF    - complex variable to store kFermi
!     ck0    - complex momentum
!     SigmaF - Self energy at the fermi energy and fermi momentum
!              Does not include the Hartree Fock part.
!     Sigma0 - Single pole self energy evaluated at
!              Wp = Wp(ipole)/Wp(R_Interstitial)*Wp(Rs)
!     dSgdE  - Derivative of sigma w.r.t. Energy
!     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
!     RelEn  - Energy relative to the fermi energy from FEFF      
!              Relen = Energy - Mu + EFermi
!     SigTot - Total many pole deltaSigma = Sigma(E,k(E))-Sigma(EF,kF)
!     DelHF  - Sigma_HartreeFock(k) - Sigma_HartreeFock(kF)      
      DOUBLE PRECISION kFermi, EFermi, Wp, Gam
      COMPLEX*16 ckF, ck0, Sig, SigmaF, Sigma0,                         &
     &     RelEn, SigTot, DelHF
      INTEGER i1, i2

!     Parameters:
      DOUBLE PRECISION DPZero, h
      PARAMETER(DPZero = 0.d0, h = 1.d-3)
      INTEGER MxIter
      LOGICAL MPole
      PARAMETER(MxIter = 1)

!     Externals:
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
!     Initialization
      UseBP = .TRUE.
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      SigmaF = 0.d0 
      Gam = 0.d0

!     Loop1: Start self consistency loop. Disabled for now.
      DO i2 = 1, MxIter
!        Loop2: Loop over poles to find SigmaF
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 5
            
!           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)
            Gam = Gamma(i1)
!           find Sigma_Fermi (SigmaF)
            ckF = kFermi*1.00001d0
            RelEn = EFermi
            SigmaF = SigmaF + Sigma1(ckF,RelEn,Wp,Gam,AmpFac(i1),       &
     &           kFermi,EFermi,.TRUE.,.FALSE.)
         END DO
 5       CONTINUE
         
!        Loop3: Loop over poles
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 10
!           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)
            Gam = Gamma(i1)
!           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
!           Find Sigma0
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,         &
     &           EFermi,.TRUE.,.FALSE.)
            
            SigTot = SigTot + Sigma0
            
!        End loop over poles.            
         END DO
 10      CONTINUE
         
!     End self-consistency loop
      END DO

!     Form delta sigma and retur.n
      SigTot = SigTot - SigmaF
      DelHF = HFExc(ck0,EFermi,kFermi) - HFExc(ckF,EFermi,kFermi)
      SigTot = SigTot + DelHF
      
      ReSig = DBLE(SigTot)
      ImSig = DIMAG(SigTot)
      
      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION Sigma1(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi,OnShll,      &
     &   UseBP)
!     Written by Josh Kas
!     Function Sigma calculates the energy dependent part
!     of Sigma(ck,Energy) from H.L. electron gas model.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
!     Input:
!     ck     - complex momentum
!     Energy - Energy
!     Wi     - Plasmon pole energy
!     Gamma  - Broadening of plasmon pole.
!     Amp    - Amplitude of plasmon pole.
!     kFermi - Fermi momentum
!     EFermi - Fermi energy
!              This is used when calculating dSigma/dE
!     OnShll - Logical flag for on shell or off shell calculation
!     UseBP  - If true, use broadened pole SE.
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
      LOGICAL OnShll, UseBP
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output:
!     Sigma(ck,Energy)
      COMPLEX*16 Sigma1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables:
!     NSing  - Number of singularities of integrand (used in cgratr)
!     NCalls - Number of fcn calls used in cgratr
!     MaxNR  - Max number of regions used in cgratr
!     DPPar  - Array of double precision parameters passed to cgratr
!              to be used in the functions r1, r2, and r3
!     CPar   - Array of complex parameters passed to cgratr
!              to be used in the functions r1, r2, and r3
!     Limit1 - Lower limit of integration
!     Limit2 - Upper limit of integration
!     HLInt1 - Integral of r2 (first integral in eq. 13 of H.L.)
!     HLInt2 - Integral of r1 (second integral in eq. 13 of H.L.)
!     HLInt3 - Integral of r1 or r3 (3rd or 4th integral)
!     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Wp, Beta, WPlus
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3
      DOUBLE PRECISION     XSing(20)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Loop variables:
      INTEGER i1, i2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Parameters:
!     ZeroPl - lower bound Limit1
!     Inf    - Upper bound of Limit2
!     AbsErr - absolute error used by cgratr
!     RelErr - Relative error used by cgratr
!     Error  - used for error codes by cgratr      
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Externals:
!     Externals:
!     cgratr      - integration routine
!     dr1,dr2,dr3 - functions to integrate
!     HFExc       - Calculates Hartree Fock exchange      
      COMPLEX*16 cgratr, r1, r2, r3, HFExc, Intgrl, bpr1, bpr2, bpr3
      EXTERNAL cgratr, r1, r2, r3, bpr1, bpr2, bpr3, HFExc, Intgrl
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0

!     DPPar is array of dp parameters to evaluate functions in cgratr.
!     Everything is in dimensionless units.
!     1. xwg
      DPPar(1) = (Wi/EFermi)
!     2. xgam
      DPPar(2) = gamma/EFermi
!     3. DPPar(3) =  1 -> On shell calculation
!        DPPar(3) = -1 -> Off shell calculation
      IF(OnShll) THEN
         DPPar(3) = 1.d0
      ELSE
         DPPar(3) = -1.d0
      END IF

!     4. xeg (gap energy)
      DPPar(4) = 0.d0
!     CPar is array of complex parameters to evaluate functions in cgratr.
!     ck in dimensionless units.
      CPar(1) = ck/kFermi
!     2. ce (complex energy)
      CPar(2) = Energy/EFermi
      
!     Josh - This is a possible fix for functions that overlap zero by a large
!     amount so that Wp does not equal Wi. 
!      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
!      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
!      Wp =  2.d0*Gamma*LOG(Gamma/Beta) + 2.d0*ATAN2(Wi,Gamma) -
!     &     Gamma*LOG(Gamma**2 + Wi**2)
!      Wp = SQRT(Wi*Wp)
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     Calculate integrals in eq. 13 of H.L.
!     1)
!     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
!     from ck + kFermi to Inf.
      Limit1 = DBLE(ck)/kFermi+1.d0
      Limit2 = Limit1*Inf
      
!     Find singularities in r2
      iFcn = 2
!      CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
!      DO i2 = NSing, 1
!         XSing(i2+1) = XSing(i2)
!      END DO
!      XSing(1) = Limit1
!      NSing = NSing+1
!     Calculate integral
      IF(UseBP) THEN
         HLInt1 = cgratr(bpr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,   &
     &        NSing, XSing,Error,NCalls,MaxNR)
      ELSE
         HLInt1 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,     &
     &        NSing, XSing,Error,NCalls,MaxNR)
      END IF
!      DO i2 = 1, NSing
!         XSing(i2) = (1.d0,0)*XSing(i2)
!      END DO

!     2)
!     Integral { ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
!                     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
!     From ck - kFermi to ck + kFermi
      Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
      Limit2 = DBLE(ck)/kFermi+1.d0

!     Find singularities in r1
!      IFcn = 1
!      CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)

!     Calculate integral
      IF(UseBP) THEN
         HLInt2 = cgratr(bpr1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,   &
     &        NSing, XSing,Error,NCalls,MaxNR)
      ELSE
         HLInt2 = cgratr(r1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,     &
     &        NSing, XSing,Error,NCalls,MaxNR)
      END IF
!      DO i2 = 1, NSing
!         XSing(i2) = (1.d0,0)*XSing(i2)
!      END DO      
         
!     3)
!     Theta(kFermi-Re(ck)) *
!     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] } +
!     Theta(Re(ck)-kFermi) *
!     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
!     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
      Limit1 = ZeroPl
      Limit2 = ABS(DBLE(ck)/kfermi-1.d0)

!     If ck = kFermi, HLInt3 = 0
      IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.                           &
     &     (DBLE(Limit2).le.DBLE(Limit1))) THEN
         HLInt3 = 0.d0
!     If ck < kFermi, HLint3 = Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }
      ELSEIF(DBLE(ck).lt.kFermi) THEN
         Limit2 = 1.d0 - DBLE(ck)/kFermi
            
!     Find singularities in r3
         iFcn = 3         
!         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
!     Calculate integral
         IF(UseBP) THEN
            HLInt3 = cgratr(bpr3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,&
     &           NSing,XSing,Error,NCalls,MaxNR)
         ELSE
            HLInt3 = cgratr(r3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,  &
     &           NSing,XSing,Error,NCalls,MaxNR)
         END IF         
!         DO i2 = 1, NSing
!            XSing(i2) = (1.d0,0)*XSing(i2)
!         END DO
         
!     Else ck > kFermi, HLint3 = Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
      ELSE
         Limit2 = DBLE(ck)/kFermi - 1.d0
            
!     Find singularities in r2
!         iFcn = 2
!         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
!     Calculate integral
         IF(UseBP) THEN
            HLInt3 = cgratr(bpr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,&
     &           NSing,XSing,Error,NCalls,MaxNR)
         ELSE
            HLInt3 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,  &
     &           NSing,XSing,Error,NCalls,MaxNR)
         END IF
!         DO i2 = 1, NSing
!            XSing(i2) = (1.d0,0)*XSing(i2)
!         END DO
      END IF
      
      Sigma1 = - Amp*Wi**2/(2.d0*pi*EFermi*ck)*                         &
     &        (HLInt1 + HLInt2 + HLInt3)

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION dSigma(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi)
!     Written by Josh Kas
!     Function dSigma calculates dSigma(ck,Energy)/dE from H.L.
!     electron gas model.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
!     Input:
!     ck     - complex momentum
!     Energy - Energy
!     Wi     - Plasmon pole energy
!     Gamma  - Broadening of plasmon pole.
!     Amp    - Amplitude of plasmon pole.
!     kFermi - Fermi momentum
!     EFermi - Fermi energy
!              This is used when calculating dSigma/dE
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output:
!     Sigma(ck,Energy)
      COMPLEX*16 dSigma
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables:
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Local Variables:
!     Wp     - Plasmon frequency, for broadened poles, Wp is not the
!              same as Wi.
!     Beta   - g_i for the negative weight pole at zero. 
!              Adding the negative weight pole at zero corrects the
!              diverging sum rule for epsilon^-1. This is irrelevant
!              for unbroadened poles.
!     NSing  - Number of singularities of integrand (used in cgratr)
!     NCalls - Number of fcn calls used in cgratr
!     MaxNR  - Max number of regions used in cgratr
!     DPPar  - Array of double precision parameters passed to cgratr
!              to be used in the functions r1, r2, and r3
!     CPar   - Array of complex parameters passed to cgratr
!              to be used in the functions r1, r2, and r3
!     Limit1 - Lower limit of integration
!     Limit2 - Upper limit of integration
!     HLInt1 - Integral of dr2 (derivative of first integral in eq. 13 of H.L.)
!     HLInt2 - Integral of dr1 (derivative second integral in eq. 13 of H.L.)
!     HLInt3 - Integral of dr1 or dr3 (derivative of 3rd or 4th integral)
!     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Wp, Beta
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3
      DOUBLE PRECISION  XSing(20)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     Loop Variables:
      INTEGER i1, i2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!     Parameters:
!     ZeroPl - lower bound Limit1
!     Inf    - Upper bound of Limit2
!     AbsErr - absolute error used by cgratr
!     RelErr - Relative error used by cgratr
!     Error  - used for error codes by cgratr
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Externals:
!     cgratr      - integration routine
!     dr1,dr2,dr3 - functions to integrate
!     HFExc       - Calculates Hartree Fock exchange
      COMPLEX*16 cgratr, dr1, dr2, dr3, HFExc
      EXTERNAL cgratr, dr1, dr2, dr3, HFExc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
!     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0
!     DPPar is array of dp parameters to evaluate functions in cgratr.
!     Everything is in dimensionless units.
!     1. xwg
      DPPar(1) = Wi/EFermi
!     2. xgam
      DPPar(2) = Gamma/EFermi
!     3. xe
      DPPar(3) = Energy/EFermi
!     4. xeg
      DPPar(4) = 0.d0
!     CPar is array of complex parameters to evaluate functions in cgratr.
!     ck in dimensionless units.
      CPar(1) = ck/kFermi
!     2. ce (complex energy)
      CPar(2) = Energy/EFermi + coni*DPPar(2)
      
!     Josh - This is a possible fix for functions that overlap zero by a large
!     amount so that Wp does not equal Wi. 
!     Wp= pi*Wi/2 + Wi*ArcTan[Wi/Gamma] - Gamma*Log[Beta] + Gamma*Log[Gamma] - 
!               1/2*Gamma*Log[Wi**2 + Gamma**2]
!      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
!      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
!      Wp =  2*Gamma*LOG(Gamma/Beta) + 2*ATAN2(Wi,Gamma) -
!     &     Gamma*LOG(Gamma**2 + Wi**2)
!      Wp = SQRT(Wi*Wp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
!     Calculate derivatives of integrals in eq. 13 of H.L.      
!     1)
!     Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
!     from ck + kFermi to Inf.      
      Limit1 = ck/kFermi+1.d0
      Limit2 = Limit1*Inf
!     Find singularities in dr2
      iFcn = 2
      CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         
!     Calculate integral
      HLInt1 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,       &
     &     NSing, XSing,Error,NCalls,MaxNR)
      DO i2 = 1, NSing
         XSing(i2) = (1.d0,0)*XSing(i2)
      END DO
      
!     2)
!     Integral d/dE{ ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
!     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
!     From ck - kFermi to ck + kFermi
      Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
      Limit2 = ck/kFermi+1.d0
         
!     Find singularities in dr1
      iFcn = 1
      CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
      HLInt2 = cgratr(dr1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,       &
     &     NSing, XSing,Error,NCalls,MaxNR)
      DO i2 = 1, NSing
         XSing(i2) = (1.d0,0)*XSing(i2)
      END DO
         
!     3)
!     Theta(kFermi-Re(ck)) *
!     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] +
!     Theta(Re(ck)-kFermi) *
!     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)]
!     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
      Limit1 = ZeroPl
      Limit2 = ABS(DBLE(ck)/kfermi-1.d0)
         
!     If ck = kFermi, HLInt3 = 0      
      IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.                           &
     &     (DBLE(Limit2).le.DBLE(Limit1))) THEN         
         HLInt3 = 0.d0
!     If ck < kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }
      ELSEIF(DBLE(ck).lt.kFermi) THEN
         Limit2 = 1.d0 - DBLE(ck)/kFermi
         
!     Find singularities in r3        
         iFcn = 3
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
!     Calculate integral         
         HLInt3 = cgratr(dr3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,    &
     &        NSing,XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO
!     Else ck > kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
      ELSE
         Limit2 = DBLE(ck)/kFermi - 1.d0
         
!     Find singularities in r2         
         iFcn = 2
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         
!     Calculate integral
         HLInt3 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,    &
     &        NSing,XSing,Error,NCalls,MaxNR)
!         DO i2 = 1, NSing
!            XSing(i2) = (1.d0,0)*XSing(i2)
!         END DO
      END IF
      
      dSigma = - Amp*Wi**2/(2.d0*pi*EFermi**2*ck)*                      &
     &     (HLInt1 + HLInt2 + HLInt3)

      RETURN
      END      

      FUNCTION HFExc(ckIn, EFermi, kFermi)
!     returns dirac-hara hartree-fock exchange
!     ck - complex momentum in units of kFermi
      use constants
      COMPLEX*16 ckIn, ck, HFExc, c
      DOUBLE PRECISION EFermi, kFermi
      ck = ckIn/kFermi
      c=-2.d0*EFermi/(pi*kFermi)
      IF(ABS(ck-1.d0).le.0.00001d0) THEN
         HFExc = c
      ELSE
         HFExc = c*(1.d0+(1.d0/ck-ck)* log( (1.d0+ck)/(ck-1.d0) )/2.d0)
      END IF
      RETURN
      END
      
!****************************************************************************
!     the following function routines are used for evaluating integals and
!     their derivatives.
!****************************************************************************
      complex*16 function r1(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)
      
!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,t1,t2,q,ck,xe
      external fq

      ck=CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
      xeg = DPPar(4)
!     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
!     print*,'fqq=', fqq
      fiq=1./(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10
      t1=(a1*a2)
      t2=(a3*a4)
!     print*,'a1,a2,a3,a4,t1,t2',a1,a2,a3,a4,t1,t2
!      t1=t1/t2
      r1=fiq*(log(a1)+log(a2)-log(a3)-log(a4))
!     Test with r=E
!      r1=xe
!     print*,'r1 return to cgratr', r1
      return
      end
!****************************************************************************
      complex*16 function dr1(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)
      
!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,t1,t2,q,ck,xe
      external fq

      ck=CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
      xeg = DPPar(4)
!     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
!     print*,'fqq=', fqq
      fiq=1.d0/(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10
!     print*,'a1,a2,a3,a4,t1,t2',a1,a2,a3,a4,t1,t2
!      t1=t1/t2
      dr1 = -fiq*(1.d0/a1+1.d0/a2-1.d0/a3-1.d0/a4)
!     Test with r=E
!      dr1=1.d0
!      write(51,*) dble(q), dble(dr1)
!     print*,'r1 return to cgratr', r1
      return
      end
!**********************************************************************
      complex*16 function r2(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)
      
!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
!      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      r2=fiq*(log(a1)-log(a2))
!      r2 = fiq*(log(a1/a2))
!     Test with r=E
!      r2=xe

30    return
      end
!**********************************************************************
      complex*16 function dr2(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)
      
!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
!      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      dr2=-fiq*(1.d0/a1-1.d0/a2)
!     Test with r=E
!      dr2=1.d0      
!      write(52,*) dble(q), dble(dr2)
30    return
      end  
!**********************************************************************
      complex*16 function r3(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)

!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
!      xeg = DPPar(4)
!     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      r3=fiq*(log(a1) - log(a2))
!     Test with r=E
!      r3=xe      
30    return
      end
!**********************************************************************
      complex*16 function dr3(q,dppar,cpar)
      use constants
      implicit double precision (a-h,o-z)

!     Input:
      double precision dppar(4)
!     dppar contains:
!     xe, xeg, xwg, xgam
      complex*16 cpar(2)

!     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
!     3. xe
      xe = CPar(2)
!     4. xeg
!      xeg = DPPar(4)
!     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      dr3=-fiq*(1.d0/a1-1.d0/a2)
!     Test with r=E
!      dr3=1.d0
!      write(53,*) dble(q), dble(dr2)
30    return
      end     
!**********************************************************************
      complex*16 function fq(q,dppar)
      use constants
      implicit double precision (a-h,o-z)
      complex*16 q
      double precision dppar(4)

!     1. xwg
      xwg = DPPar(1)
!     2. xgam
      xgam = DPPar(2)
!     Here I am going to change the dispersion relation to 
!     wq = wp + 1/2 * q**2 
!     This makes calculation of broadened poles easier. 
!      fq = xwg-coni*xgam + q**2
      
!     fq(q)=w1(q)=sqrt(w1**2+((omega(q)-omega_p)/omega_f)**2)
!     omega(q)**2=omega_p**2+omega_g**2(q)
!     fq(q)=xwg+a2*q**2+a4*q**4    xwg=(w1/ef)**2
!     electron gas parameters xwg=wp**2 a2=4/3, a4=1
!     uncomment the following 4 lines to use the old dispersion relation.
      a2=4.d0/3.d0
      a4=1.d0
      fq=(xwg-coni*xgam)**2 + a2*q*q + a4*q**4
      fq=sqrt(fq)
      return
      end

      FUNCTION Intgrl(func, a, b, NPts, DpPar, CPar)
!     Function Integ integrates func by trapezoidal rule from a to b
!     using NPts steps size.
      COMPLEX*16 Intgrl, a, b, CPar(10)
      DOUBLE PRECISION DpPar(10)
      INTEGER NPts
      COMPLEX*16 func
      EXTERNAL func

      INTEGER i1
      COMPLEX*16 dx, sum, y1, y2, x1, x2


      sum = 0.d0
      dx = (b-a)/DBLE(NPts-1)
      
      x1 = a
      y1 = func(x1, DpPar, CPar)
      
      DO i1 = 1, NPts
         x2 = a + i1*dx
         y2 = func(x2, DpPar, CPar)
         sum = sum + (y1+y2)
         y1 = y2
         x1 = x2
      END DO

      Intgrl = sum*dx/2.d0

      RETURN
      END
               
!*********************************************************************
!   This is Steve White's rewrite of Mike Teter's integration routine.  
!   Modified by J. Rehr for complex integration.
!   The following is a listing of the arguments in the initial function 
!   statement:
!     fn    -- routine requires external function statement in MAIN
!     xmin  -- lower limit
!     xmax  -- upper limit
!     abr   -- absolute tolerable error
!     rlr   -- relative tolerable error
!     nsing -- number of singularities or regions requiring 
!     special attention
!     xsing -- array of locations of singularities or endpoints
!     of special regions
!     error -- output for routine error messages
!     numcal-- the number of times fn was called
!     maxns -- the maximum number of regions being considered simultaneously
!     function cgratr(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
!     fn declared double precision
!     double precision function cgratr(fn,xmin,xmax,abr,rlr,
!     fn declared complex*16
      
      complex*16 function cgratr(fn,dppar,cpar,xmin,xmax,abr,rlr,       &
     &     nsing,xsing,error,numcal,maxns)
      implicit double precision (a-h,o-z)
      parameter (mx=1500)
      integer nsing
      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1
      complex*16 xleft(mx), cpar(10)
      double precision dppar(10)
      external fn
!     dimension xleft(mx),fval(3,mx),dx(3),wt(3)
      dimension wt9(9),dx(3),wt(3)
      dimension xsing(20)
      logical atsing
      save dx,wt,wt9
      data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
      data wt/0.277777777777777778  ,0.4444444444444444444  ,           &
     &     0.2777777777777777778  /
      data wt9/0.0616938806304841571  ,0.108384229110206161  ,          &
     &     0.0398463603260281088  ,0.175209035316976464  ,              &
     &     0.229732989232610220  ,0.175209035316976464  ,               &
     &     0.0398463603260281088  ,0.108384229110206161  ,              &
     &     0.0616938806304841571  /
!     nstack is the number of different intervals into which the 
!     integration region is currently divided. The number of regions can
!     grow if more accuracy is needed by dividing the right-most region
!     into three regions. The number of regions shrinks when the integral
!     over the right-most region is accurate enough, in which case that
!     integral is added to the total (stored in cgratr) and the region
!     is removed from consideration (and a new region is the right-most).
      nstack=nsing+1
      maxns = nstack
      error=0.  
      cgratr=0.  
!     The array xleft stores the boundary points of the regions.
!     The singular points just govern the initial placement of the regions.
      xleft(1)=xmin
      xleft(nsing+2)=xmax
      if(nsing.gt.0) then
         do 9 j=1,nsing
            xleft(j+1)=xsing(j)
 9       continue
      endif
!     For each region, calculate the function and store at three selected points.
      do 1 jj=1,nstack
         del=xleft(jj+1)-xleft(jj)
!     print*, 'fn call j= ,'
         do 1 j=1,3
!     print*, 'fn call in cgratr j= ',j
            fval(j,jj)=fn(xleft(jj)+del*dx(j),dppar,cpar)
 1    continue
!     print*, 'output of fn call, fval(j,jj)',fval(j,jj)
      numcal = nstack * 3
 6    continue
      if(nstack+3.ge.mx) then
         write(*,*) ' TOO MANY REGIONS'
         stop 0006
      endif
!     Divide the rightmost region into three subregions.  
      del=xleft(nstack+1)-xleft(nstack)
      xleft(nstack+3)=xleft(nstack+1)
      xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
      xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
!     The three data points already found for the region become the 
!     middle data points (number 2 in first index of fval) for each region.
      fval(2,nstack+2)=fval(3,nstack)
      fval(2,nstack+1)=fval(2,nstack)
      fval(2,nstack)=fval(1,nstack)
!     Now do the integral over the right-most region in two different ways-
!     a three point integral (valu) over each of the three subregions 
!     and a more accurate nine-point integral (value) over whole region.
!     valu is used only for the error estimate.
      icount=0
      value=0.  
      valu=0.  
      do 3 j=nstack,nstack+2
         del1=xleft(j+1)-xleft(j)
!     print*, 'fn call 2'
         fval(1,j)=fn(xleft(j)+dx(1)*del1,dppar,cpar)
         fval(3,j)=fn(xleft(j)+dx(3)*del1,dppar,cpar)
!     print*, 'fn call 2'
         numcal = numcal + 2
         do 5 k=1,3
            icount=icount+1
            value=value+wt9(icount)*fval(k,j)*del
            valu=valu+fval(k,j)*wt(k)*del1
 5       continue
 3    continue
      dif=abs(value-valu)
!     If the following condition is true, add in this integral to the total,
!     and reduce the number of regions under consideration.
      frac = del / (xmax - xmin)
      atsing = .false.
      if(frac .le. 1.0e-8) atsing = .true.
      if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.              &
     &     (atsing .and.                                                &
     &     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
!     The following commented out line is Teeter's old error criterion.
!     if(dif.le.abr.or.dif.le.rlr*abs(value))then
         cgratr=cgratr+value
         error=error+abs(dif)
         nstack=nstack-1
!     If no more regions, we are done.
         if(nstack.le.0) return
      else
!     If the integration is insufficiently accurate, make each of the 
!     three subregions of the right-most region into regions.
!     On next pass the right-most of these is the new current region.
         nstack=nstack+2
         maxns = max(maxns,nstack)
      endif
      go to 6
      end
