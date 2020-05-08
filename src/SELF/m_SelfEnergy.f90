!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_SelfEnergy.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SelfEnergyMod
  USE IOMod
  USE ErrorMod
  USE constants

  TYPE SEInput
     ! User Input:
     ! Run self-energy > 0
     ! 1 use many-pole model
     ! 2 use broadened pole model.
     INTEGER RunMode
     ! User specified EpsInv
     DOUBLE PRECISION Eps0
     ! User specified loss function.
     CHARACTER(30) LossFile
     ! User specified pole representation.
     CHARACTER(30) PoleFile
     ! EGridFile - File that contains energy grid for output.
     CHARACTER(30) EGridFile
  END TYPE SEInput

  TYPE SEData
     ! Data used by SelfEnergy subroutine
     ! Fermi energy, Wigner-Sietz radius
     DOUBLE PRECISION Rs
     ! Number of energy points.
     INTEGER NEPts
     ! Energy grid relative to vacuum.
     COMPLEX*16,ALLOCATABLE :: EGrid(:)
     ! Number of poles.
     INTEGER NPoles
     ! Pole data. Position, width, and amplitude.
     DOUBLE PRECISION,ALLOCATABLE :: Omega(:), Width(:), g(:)
     ! Self energy
     COMPLEX*16,ALLOCATABLE :: SigDat(:), Z(:)
  END TYPE SEData
     
CONTAINS

  SUBROUTINE ReadSEInp(InData)
    TYPE(SEInput),INTENT(OUT) :: InData

    CHARACTER(14) SEFile
    CHARACTER(200) ErrorMessage(5)

    SEFile = 'selfenergy.inp'
    ErrorMessage(1) = "File: selfenergy.inp"
    ErrorMessage(2) = "Comment lines begin with #*!cC."
    ErrorMessage(3) = "Data lines are ordered as follows:"
    ErrorMessage(4) = "Integer RunMode: 0 (Many Pole) or 1 (Broadened Pole)."
    ErrorMessage(5) = "Strings: LossFile PoleFile EgridFile"
    
    
    CALL ReadData(SEFile, Int1 = InData%RunMode,ErrorMessage = ErrorMessage)
    CALL ReadData(SEFile, String1 = InData%LossFile, String2 = &
         & InData%PoleFile, String3 = InData%EGridFile)
    CALL ReadData(SEFile, Double1 = InData%Eps0)
  END SUBROUTINE ReadSEInp

  SUBROUTINE ReadSEData(InData,SEData1)
    TYPE(SEData),INTENT(OUT) :: SEData1
    TYPE(SEInput),INTENT(IN) :: InData
    DOUBLE PRECISION Energy
    INTEGER iError

    ! Read NEPts from egrid.dat
    CALL ReadData(InData%EGridFile, Int1 = SEData1%NEPts)    
    !Allocate space for energy grid.
    ALLOCATE(SEData1%EGrid(SEData1%NEPts), STAT = iError)
    CALL CheckAllocation(iError,'Allocation error occured during allocation of energy grid array.')
    ! Read energy points from egrid.dat
    CALL ReadArrayData(InData%EGridFile, DComplex1 = SEData1%EGrid)
    ! Close file.
    CALL CloseFl(InData%EGridFile)

    ! Read Rs from PoleFile.
    CALL ReadData(InData%PoleFile, Double1 = SEData1%Rs)
    ! Read NPoles from PoleFile
    CALL ReadData(InData%PoleFile, Int1 = SEData1%NPoles)
    !Allocate space for pole data.
    ALLOCATE(SEData1%Omega(SEData1%NPoles), SEData1%Width(SEData1%NPoles), SEData1%g(SEData1%NPoles), STAT = iError)    
    CALL CheckAllocation(iError,'Allocation error occured during allocation of pole data arrays.')
    ! Read pole data.
    CALL ReadArrayData(InData%PoleFile, Double1 = SEData1%Omega, Double2 = SEData1%Width, Double3 = SEData1%g)    
    
    ! Close file.
    CALL CloseFl(InData%PoleFile)

    ! Unit conversion from eV to Hartree.
    SEData1%EGrid(:) = SEData1%EGrid(:)/hart
    SEData1%Omega(:) = SEData1%Omega(:)/hart
    SEData1%Width(:) = SEData1%Width(:)/hart

    ! Allocate space
    ALLOCATE(SEData1%SigDat(SEData1%NEPts))
    ALLOCATE(SEData1%Z(SEData1%NEPts))
  END SUBROUTINE ReadSEData
  
  SUBROUTINE SEnergy(SEData1, InData)
    TYPE(SEData),INTENT(INOUT) :: SEData1
    TYPE(SEInput),INTENT(IN) :: InData
    COMPLEX*16 SigmaF, Sigma0, Z, En
    INTEGER iE
    LOGICAL UseBP
    
    
    SELECT CASE(InData%RunMode)
       CASE(0)
          UseBP = .FALSE.
       CASE(1)
          UseBP = .TRUE.
       CASE DEFAULT
          UseBP = .FALSE.
    END SELECT

    ! Find Sigma(EFermi)
    En = 0.001d0
    CALL CSigZ(En, 0.d0, SEData1%Rs, SigmaF, Z, SEData1%Omega, SEData1%Width, SEData1%g, 0.d0, SEData1%NPoles, .TRUE., UseBP)

    ! Now calculate Sigma(E) and Z(E)
    CALL wlog('Beginning calculation of the self energy.')
    DO iE = 1, SEData1%NEPts
       !PRINT '(A9,i4,A4,i5)', '   Point ', iE, ' of ', SEData1%NEPts
       IF(DBLE(SEData1%EGrid(iE)).le.0.d0) THEN
          En = 0.001d0 + (0.d0,1.d0)*SEData1%EGrid(iE)
       ELSE
          En = SEData1%EGrid(iE)
       END IF          
       CALL CSigZ(En, 0.d0, SEData1%Rs, Sigma0, Z, SEData1%Omega, SEData1%Width, SEData1%g, 0.d0,SEData1%NPoles, .TRUE., UseBP) !KJ 12-2011 was missing argument "EGap" ; I inserted "0.d0"
!SUBROUTINE  CSigZ(c16Energy, r8Mu, r8Rs, c16SigTot, c16ZTot, r8(MxPole)WpScl, r8(MxPole)Gamma, r8(MxPole)AmpFac, r8EGap, iNPoles, lOnShll, lUseBP)
       SEData1%SigDat(iE) = (Sigma0 - SigmaF)
       SEData1%Z = Z
    END DO
  END SUBROUTINE SEnergy

  SUBROUTINE MkExc(Energy, Loss, Eps0, Wi, gi, NPoles)
    REAL(8), INTENT(IN)  :: Energy(:), Loss(:), Eps0
    INTEGER, INTENT(INOUT)  :: NPoles
    REAL(8), INTENT(OUT) ::  Wi(NPoles), gi(NPoles)
    INTEGER, PARAMETER :: NFine = 50000
    REAL(8)  M1, MM1, M1Init, MM1Init, frac, Cnt, EMin, EMax, a, Delta(NPoles), &
         & deltaE, EFine(NFine), LFine(NFine)
    INTEGER i1, i2, i3, iPole, NE, iMin, iMax
    
    ! Interpolate loss function onto finely spaced grid.
    NE = SIZE(Loss)
    deltaE = (MIN(Energy(NE),1000.d0) - Energy(1))/NFine
    EFine(1) = 0.d0
    LFine(1) = 0.d0
    DO i1 = 1, NFine -1 !KJ 2/12  NFine
       EFine(i1+1) = deltaE*i1
       IF(EFine(i1+1).GT.Energy(1)) THEN
          CALL terp(Energy,Loss,NE,1,EFine(i1+1),LFine(i1+1))
       ELSE
          LFine(i1+1) = EFine(i1+1)*Loss(1)/Energy(1)
       END IF
    END DO

    NE = NFine
    ! Find the total integral of LFine and LFine(E)/E.
    MM1 = 0.d0
    DO i1 = 1, NE - 1
       IF(EFine(i1).NE.0.d0) THEN
          MM1 = MM1 + 0.5d0*(LFine(i1+1)/EFine(i1+1) &
            & + LFine(i1)/EFine(i1))*(EFine(i1+1) - EFine(i1))
       ELSE
          MM1 = MM1 + LFine(i1+1)
       END IF
    END DO

    ! Put poles evenly in area of MM1
    frac = MM1/(NPoles)
    iMax = 1
    iMin = 1
    EMax = EFine(1)
    M1  = 0.d0
    MM1 = 0.d0
    iPole = 1
    DO i1 = 1, NE - 1
       IF(EFine(i1).NE.0.d0) THEN
          MM1 = MM1 + 0.5d0*(LFine(i1+1)/EFine(i1+1) &
               & + LFine(i1)/EFine(i1))*(EFine(i1+1) - EFine(i1))
       ELSE
          MM1 = MM1 + LFine(i1+1)
       END IF
       M1 = M1 + 0.5d0*(LFine(i1+1)*EFine(i1+1) &
            & + LFine(i1)*EFine(i1))*(EFine(i1+1) - EFine(i1))
       IF((MM1.ge.frac).OR.(i1.EQ.NE-1)) THEN 
          ! We have overshot EMax, OK, just won't use all the poles
          ! Also, use one last pole to finish catch the tail of the integral
          ! even if it doesn't sum to frac.
          gi(iPole) = 2.d0/pi*MM1
          Wi(iPole) = SQRT(M1/MM1)
          Delta(iPole) = EFine(i1+1) - EFine(iMin)
          iPole = iPole + 1
          iMin  = i1 + 1
          MM1   = 0.d0
          M1    = 0.d0
       END IF
      
       IF(iPole.GT.NPoles) EXIT
    END DO
    iPole = iPole - 1
    NPoles = iPole

    ! Now scale poles and amplitudes to obtain correct Eps0
    IF((Eps0.GT.-1.5d0).AND.Eps0.LT.0) THEN
       ! Metallic case. Eps0 = inf
       frac = 1/SUM(gi(1:NPoles))
    ELSEIF(Eps0.GT.1.d0) THEN
       frac = (1 - 1.d0/Eps0)/SUM(gi(1:NPoles))
    ELSE
       frac = 1.d0
    END IF
    gi(1:NPoles) = frac*gi(1:NPoles)
    Wi(1:NPoles) = Wi(1:NPoles)/SQRT(frac)

    DO iPole = 1, NPoles
       CALL WriteData('exc.dat', Double1 = Wi(iPole), Double2 = 0.1d0, &
         Double3 = gi(iPole), Double4 = pi/2.d0*gi(iPole)*Wi(iPole)/Delta(iPole))
    END DO
    CALL CloseFl('exc.dat')
  END SUBROUTINE MkExc

END MODULE SelfEnergyMod
