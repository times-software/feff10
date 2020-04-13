!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: eps2exc.f90,v $:
! $Revision: 1.1 $
! $Author: hebhop $
! $Date: 2010/10/18 19:40:41 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Josh Kas - Standalone version of mkexc called eps2exc
! Written 10/2010
PROGRAM eps2exc
  USE IOMOD
  use DimsMod, only: MxPole
  USE SelfEnergyMod

  implicit double precision (a-h, o-z)

  integer iPl, ipole, NPoles, iErr
  double precision, allocatable :: Energy(:), Loss(:)
  double precision WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole), Eps0

  ! Get input from user.
  PRINT*, 'Enter number of poles:'
  READ*, NPoles
  PRINT*, 'Enter static dielectric constant.'
  PRINT*, 'Enter -1 for a metal (eps0 = -infinity),'
  PRINT*, 'Enter -2 to leave the dielectric constant unchanged.'
  READ*, Eps0

  WpCorr(:) = -1.d30
  CALL OpenFl('loss.dat', FileStatus = 'OLD')
  NData = NumberOfLines('loss.dat')
  
  ALLOCATE(Energy(NData), Loss(NData), STAT = iErr)
  IF(iErr.NE.0) THEN
     PRINT*, "Error in allocation of Loss and Energy arrays:", iErr
     STOP
  END IF
  
  CALL ReadArrayData('loss.dat', Double1 = Energy, Double2 = Loss)
  CALL CloseFl('loss.dat')

  CALL MkExc(Energy, Loss, Eps0, WpCorr, AmpFac, NPoles)
  DEALLOCATE(Energy,Loss)

END PROGRAM eps2exc
