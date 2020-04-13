!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: TestIOFiles.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM TestIOFiles
  USE IOFiles
  
  ! Define a few IOFiles to work with
  TYPE(IOFile) IOF1, IOF2, IOF3
  CHARACTER c

  ! Initialize them.
  PRINT*, 'Initializing files'
  CALL InitIOFile(IOF1,'file1')
  CALL InitIOFile(IOF2,'file2')
  CALL InitIOFile(IOF3,'file3')

  ! Print the IOFile info for each. They should all be the 
  ! same except for filename.
  PRINT*, 'Printing file information.'
  CALL PrintIOFileInfo(IOF1)
  CALL PrintIOFileInfo(IOF2)
  CALL PrintIOFileInfo(IOF3)

  ! Now lets set IOF2=IOF1. This tests CopyIOFile
  IOF2=IOF1
  PRINT*, IOF2%FileName

  ! Reset it.
  IOF2%FileName = 'file2'

  ! Now lets test AddIOFile. This adds the files to the filestack.
  CALL AddIOFile('file1')
  CALL AddIOFile('file2')
  CALL AddIOFile('file3')
  
  ! Print the stack.
  CALL PrintFileStackInfo

  ! Lets test DeleteIOFile. This deletes the file from the filestack.
  CALL DeleteIOFile('file2')
  CALL PrintFileStackInfo
  CALL AddIOFile('file2')
  
  ! Now lets test IndexIOFile.
  PRINT*, 'This should say 2: ', IndexIOFile('file3')

  ! Test adding a file that already exists.
  CALL AddIOFile('file1') ! This should do nothing.
  PRINT*, 'This should say 1: ', IndexIOFile('file1')
  CALL PrintFileStackInfo

  ! Try SetIOFileInfo and GetIOFileInfo
  CALL SetIOFileInfo('file1',UnitNumber=1,NSections=1,FileAction='READWRITE',EOF=.TRUE.)

  CALL DeleteIOFile('file2')
  CALL DeleteIOFile('file3')

  CALL PrintFileStackInfo
END PROGRAM TestIOFiles
  
