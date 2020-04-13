! This program takes a list of xy data files and interpolates the data from each onto an x grid which is the 
! union of all of the x grids, then adds all of y data.
! Written by J. J. Kas - 1/25/2010
!PROGRAM AddXY
SUBROUTINE AddEps(Files,NFiles,Weights,print_eps)
  USE IOMod
  USE IOFiles
  IMPLICIT NONE
  CHARACTER*(*) Files(NFiles)
  LOGICAL print_eps
  INTEGER NFiles, iFile, iData, NDataTot, iTot, iSort, iNew, iErr
  REAL(8), ALLOCATABLE :: E(:,:), Eps1(:,:), Eps2(:,:), ETot(:), Eps1Tot(:),Eps2Tot(:)
  REAL(8) xTmp, yTmp, Weights(NFiles)
  INTEGER,ALLOCATABLE :: NData(:), iUnit(:)
  ! Get the filenames.
  !PRINT*, "# Enter the number of files."
  !READ*, NFiles
  ALLOCATE(iUnit(NFiles))

  !DO iFile = 1, NFiles
  !   PRINT*, '# Enter file #', iFile 
  !   READ '(a)', Files(iFile)
  !END DO

  ! Get the number of data points for each file.
  ALLOCATE(NData(NFiles))
  NDataTot = 0
  DO iFile = 1, NFiles
     NData(iFile) = 0
     CALL OpenFl(Files(iFile))
     CALL GetIOFileInfo(Files(iFile), UnitNumber = iUnit(iFile))
     ! Read comments from header. Assumes that comments are ONLY in header and starts with #.
     CALL RdCmt(iUnit(iFile),'#cC!',iErr)
     IF(iErr.EQ.-1) CALL par_stop("Invalid file found when forming epsilon: stopping")
     DO
        READ(iUnit(iFile),*,END = 5)
        NData(iFile) = NData(iFile) + 1
     END DO
5    CONTINUE
     NDataTot = NDataTot + NData(iFile)
     CALL CloseFl(Files(iFile))
  END DO
  ! Make space for x data.
  ALLOCATE(E(MAXVAL(NData),NFiles),Eps1(MAXVAL(NData),NFiles),Eps2(MAXVAL(NData),NFiles))
  ALLOCATE(ETot(NDataTot),Eps1Tot(NDataTot),Eps2Tot(NDataTot))

  ! Read in x and y data.
  iTot = 0
  DO iFile = 1, NFiles
     CALL OpenFl(Files(iFile))
     CALL GetIOFileInfo(Files(iFile), UnitNumber = iUnit(iFile))
     
     CALL RdCmt(iUnit(iFile),'#   ',iErr)
     IF(iErr.EQ.-1) CALL par_stop("Invalid file found when forming epsilon: stopping")
     DO iData = 1, NData(iFile)
        READ(iUnit(iFile),*) E(iData,iFile), Eps1(iData,iFile), Eps2(iData,iFile)
        iTot = iTot + 1
        ! Keep xTot sorted and unique.
        xTmp = E(iData,iFile)
        iNew = iTot
        IF(iTot.EQ.1) ETot(1) = xTmp
        DO iSort = iTot - 1, 1, -1
           IF( xTmp.LT.ETot(iSort) ) THEN
              ! Change iNew.
              iNew = iSort
           ELSEIF (xTmp.EQ.ETot(iSort)) THEN
              ! Don't keep this point.
              iTot = iTot -1
           ELSE
              ! xTot is sorted. Fill it and exit loop.
              IF(iNew.NE.iTot) THEN
                 ! Shift xTot
                 ETot(iNew + 1:iTot) = ETot(iNew:iTot - 1)
                 ETot(iNew) = xTmp
              ELSE
                 ETot(iNew) = xTmp
              END IF
              EXIT
           END IF
        END DO
     END DO
     CLOSE(iUnit(iFile))
  END DO
  NDataTot = iTot
  Eps1Tot(:) = 0.d0
  Eps2Tot(:) = 0.d0
  DO iTot = 1, NDataTot
     DO iFile = 1, NFiles
        ! Interpolate each data set onto xTot, and add to yTot
        CALL terp(E(:,iFile),Eps1(:,iFile),NData(iFile),1,ETot(iTot),yTmp)
        Eps1Tot(iTot) = Eps1Tot(iTot) + yTmp*Weights(iFile)
        CALL terp(E(:,iFile),Eps2(:,iFile),NData(iFile),1,ETot(iTot),yTmp)
        Eps2Tot(iTot) = Eps2Tot(iTot) + yTmp*Weights(iFile)
     END DO
  END DO
  
  ! Open output file and write data.
  CALL OpenFl('loss.dat')
  CALL GetIOFileInfo('loss.dat', UnitNumber = iUnit(1))
  WRITE(iUnit(1),'(A)') '# E(eV)    Loss'
  DO iTot = 1, NDataTot
     WRITE(iUnit(1),*), ETot(iTot), Eps2Tot(iTot)/((Eps1Tot(iTot)+1)**2 + Eps2Tot(iTot)**2)
  END DO
  CALL CloseFl('loss.dat')
  IF(print_eps) THEN
     CALL OpenFl('epsilon.dat')
     CALL GetIOFileInfo('epsilon.dat', UnitNumber = iUnit(1))
     WRITE(iUnit(1),'(A)') '# E(eV)    eps1    eps2'
     DO iTot = 1, NDataTot
        WRITE(iUnit(1),*), ETot(iTot), Eps1Tot(iTot)+1.d0, Eps2Tot(iTot)
     END DO
     CALL CloseFl('epsilon.dat')
  END IF
  DEALLOCATE(iUnit,NData,E,Eps1,Eps2,Eps1Tot,Eps2Tot)
  RETURN
END SUBROUTINE AddEps
