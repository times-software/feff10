!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: calcse.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM CalcSE
      USE SelfEnergyMod

      TYPE(SEInput) SEInput1
      TYPE(SEData) SEData1
      ! Read input file.
      CALL ReadSEInp(SEInput1)

      CALL ReadSEData(SEInput1,SEData1)

      ! Calculate self-energy
      CALL SEnergy(SEData1, SEInput1)

    END PROGRAM CalcSE
         
         
