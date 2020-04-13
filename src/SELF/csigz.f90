!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: csigz.f90,v $:
! $Revision: 1.10 $
! $Author: hebhop $
! $Date: 2010/04/13 18:05:47 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CSigZ(Energy, Mu, Rs, SigTot, ZTot, WpScl, Gamma,      &
     &     AmpFac, EGap, NPoles, OnShll, UseBP)
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
!            3. Find derivative w.r.t. E dSgdE
!            4. k1**2  = 
!                  k0**2 - 2*(Sigma0-SigmaF)/(1-dSgdE)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
	  use dimsmod,only:MxPole !KJ 7-09 replaces parameter statement
!     Parameters:
!     MxPole - Maximum number of poles
!KJ      INTEGER MxPole
!KJ      PARAMETER(MxPole=1000)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Input:
!     Energy - Energy at which to evaluate Sigma
!     Mu     - Fermi energy as calculated by FEFF.
!     Rs     - R sub s (sphere of radius Rs holds charge e)
!     WpScl  - Scale Wp in interstitial by WpScl
!     Gamma  - Use broadening Gamma when calculating Sigma
!     AmpFac - Use amplitude AmpFac for plasmon pole.
!     NPoles - Number of poles.
      DOUBLE PRECISION Rs, WpScl(MxPole), Gamma(MxPole), AmpFac(MxPole), Mu, EGap
      INTEGER NPoles
      COMPLEX*16 Energy
      LOGICAL UseBP
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output:
!     ReSig  - Re[Sigma(k,e)]
!     ImSig  - Im[Sigma(k,e)]
!     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
!     Note: Atomic units are used.
      DOUBLE PRECISION ReSig, ImSig
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      DOUBLE PRECISION kFermi, EFermi, Wp, WpOverEf, Gam, xk0, Ei, EiTot
      COMPLEX*16 ckF, ck0, ckP, SigmaF, Sigma0, SigmaP, dSgdE, ZTot,    &
     &     RelEn, RelEnP, SigTot, DelHF
!     Loop variables
      INTEGER i1, i2

!     Parameters:
      DOUBLE PRECISION DPZero
      PARAMETER(DPZero = 0.d0)
      INTEGER MxIter
      PARAMETER(MxIter = 1)
      LOGICAL OnShll

!     Externals:
!     Sigma1 - calculates the energy dependent part of self energy
!              for a single pole.
!     dSigma - calculates derivative of self energy w.r.t energy.
!     HFExt  - calculates Hartree Fock exchange
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
!     Initialization
      ZTot = 0.d0
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      dSgdE = 0.d0
      SigmaF = 0.d0
      Gam = 0.d0
      EiTot  = 0.d0
      
!     Loop1: Start self consistency loop.
!     This does not seem to work, so MxIter = 1
      DO i2 = 1, MxIter
         
         dsgdE = 0.d0
!        Loop3: Loop over poles         
         DO i1 = 1, NPoles
!           Wp is in Hartrees
            Wp = WpScl(i1)
            IF(UseBP) THEN
               Gam = Gamma(i1)
            ELSE
               Gam = 0.d0
            END IF
!           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi + EGap/2.d0
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
!           Find Sigma0 = Sigma(ck0,E); ck0=SQRT(2*(Energy-Mu))
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,EFermi,.TRUE.,UseBP)
            
!            write(71,*) DBLE(RelEn-EFermi), DBLE(dSgdE), DIMAG(dSgDE)
            
!           dSgdE = dSgdE + 
!    &          dSigma(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,EFermi)
!            RelEnP = 
!     &           RelEn*(1.d0 + MAX(MIN(
!     &           DBLE(wp/(RelEn+0.00001d0)*0.05d0) ,0.01d0),0.001d0))
            IF(.FALSE.) THEN
               ! Use numerical derivative for now.
               RelEnP = RelEn*1.001d0
               SigmaP = Sigma1(ck0,RelEnP,Wp,Gam,AmpFac(i1),kFermi,     &
     &              EFermi,.FALSE.,UseBP)
               dSgdE = dSgdE + (SigmaP - Sigma0)/(RelEnP-RelEn)
            ELSE
               dSgdE = dSgdE + dSigma(ck0,RelEn,Wp,Gam,AmpFac(i1),      &
     &              kFermi,EFermi)
            END IF

!            IF(.NOT.UseBP) THEN
             IF(.FALSE.) THEN
               xk0 = DBLE(ck0)/kFermi
               WpOverEf = Wp/EFermi
               CALL quinn (xk0, Rs, WpOverEf, EFermi, Ei)
               Ei = Ei*AmpFac(i1)
!              Add Quinn correction
               IF(DIMAG(Sigma0).ge.Ei) Sigma0 = DBLE(Sigma0) + coni*Ei
            END IF

!           Uncomment the following line to print derivative            
!            write(72,*) DBLE(RelEn-EFermi), DBLE(dSgdE2), DIMAG(dSgDE2)         

!           SigTot is sum of poles
            SigTot = SigTot + Sigma0
            
!        End Loop3: loop over poles.            
         END DO         
         
!     End Loop1: self-consistency loop      
      END DO

      ! HF part of the self energy.
      SigTot  = SigTot + HFExc(ck0,EFermi,kFermi)
!      SigTot = S
!     Form ZTot and return Re and Im parts of Sigma.
      IF(.TRUE.) THEN
         ZTot = 1.d0/(1.d0-dSgdE)
      ELSE
         ZTot = 1.d0
      END IF

      RETURN
      END

