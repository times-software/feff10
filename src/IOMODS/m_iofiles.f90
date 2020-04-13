!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_iofiles.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE IOFiles
  USE ErrorMod
  IMPLICIT NONE
  ! MaxLength - Maximum length of filename strings
  ! MaxSections - Maximum number of sections in a file
  INTEGER,PRIVATE :: MaxLength, MaxSections, NError
  PARAMETER( MaxLength = 300, MaxSections = 2000, NError = 2 )
  
  ! Default format strings for Integer, Real, Double, Complex, 
  ! Double Complex, and String
  CHARACTER(300) DefaultIFormat, DefaultRFormat, DefaultDFormat, &
       & DefaultCFormat, DefaultDCFormat, DefaultSFormat
  PARAMETER ( DefaultIFormat  = '(I10)' )
  PARAMETER ( DefaultRFormat  = '(E20.10)' )
  PARAMETER ( DefaultDFormat  = '(E20.10)' )
  PARAMETER ( DefaultCFormat  = '(1E20.10,1X,1E20.10)' )
  PARAMETER ( DefaultDCFormat = '(1E20.10,1X,1E20.10)' )
  PARAMETER ( DefaultSFormat  = '(A)' )
  
  ! Default file format is text
  CHARACTER(3) DefaultFileFormat
  PARAMETER ( DefaultFileFormat = 'TXT' )


  TYPE IOFile
     ! Filename - name of file.
     ! DataTypeLine - type description of the columns of data that will 
     ! be in a section. Section changes when data types change or user 
     ! specifies a new section.
     CHARACTER(300) FileName, DataTypeLine(MaxSections)
     ! Format strings for each type of data (Integer, Real, Double, 
     ! Complex, Complex*16, and String
     CHARACTER(100) IFormat, RFormat, DFormat, CFormat, DCFormat, SFormat
     ! Specify whether you want output to be ASCII ('txt') or packed 
     ! ascii ('pad'). Can add new formats later, i.e. 'xml' ...
     CHARACTER(5) FileFormat
     ! IOFileAction: READ, WRITE, READWRITE
     CHARACTER(9) IOFileAction
     ! NLines     - number of lines in a section
     ! TotLines   - total number of lines in file.
     ! NSections  - Number of sections (Write) or current section (Read).
     ! UnitNumber - unit number of the file.
     INTEGER NLines(MaxSections), TotLines, NSections, UnitNumber
     ! IsInitialized - this IOFile was initialized if true.
     ! ExistedOnOpen - file existed when we opened it if true.
     ! IsOpen        - file is open if true.
     ! EOF           - End of file has been reached (Read) if true.
     LOGICAL IsInitialized, ExistedOnOpen, IsOpen, EOF
     ! ReadError and WriteError define possible errors that could occur, and what action to take
     ! if each error does occur.
     TYPE(ErrorType) IOErrors(NError)
  END TYPE IOFile

  TYPE(IOFile),PRIVATE,SAVE :: FileStack(100)
  INTEGER,PRIVATE,Save :: NFiles = 0

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE CopyIOFile
  END INTERFACE

CONTAINS

  ! Set IOFile1 = IOFile2
  SUBROUTINE CopyIOFile(IOFile1, IOFile2)
    TYPE(IOFile),INTENT(OUT) :: IOFile1
    TYPE(IOFile),INTENT(IN)  :: IOFile2

    IOFile1%FileName = IOFile2%FileName
    IOFile1%DataTypeLine(:) = IOFile2%DataTypeLine(:) 

    IOFile1%IFormat  = IOFile2%IFormat
    IOFile1%RFormat  = IOFile2%RFormat
    IOFile1%DFormat  = IOFile2%DFormat
    IOFile1%CFormat  = IOFile2%CFormat
    IOFile1%DCFormat = IOFile2%DCFormat
    IOFile1%SFormat  = IOFile2%SFormat

    IOFile1%FileFormat = IOFile2%FileFormat

    IOFile1%IOFileAction = IOFile2%IOFileAction

    IOFile1%NLines(:) = IOFile2%NLines(:)
    IOFile1%TotLines = IOFile2%TotLines
    IOFile1%NSections = IOFile2%NSections
    IOFile1%UnitNumber = IOFile2%UnitNumber

    IOFile1%IsInitialized = IOFile2%IsInitialized
    IOFile1%ExistedOnOpen = IOFile2%ExistedOnOpen
    IOFile1%IsOpen = IOFile2%IsOpen
    IOFile1%EOF = IOFile2%EOF


  END SUBROUTINE COPYIOFILE

  ! Initialize IOFile IOF. FileName must be specified.
  ! All other arguments are optional.
  SUBROUTINE InitIOFile(IOF, FileName)
    TYPE(IOFile),INTENT(INOUT) :: IOF
    CHARACTER*(*),INTENT(IN) :: FileName
    ! Loop variables
    INTEGER i1
    ! Set FileName, NSections, NLines
    IOF%FileName  = FileName
    IOF%NSections = 0
    IOF%NLines(:) = 0
    IOF%TotLines  = 0
    IOF%NSections = 0
    IOF%DataTypeLine(:)='Unknown'
    IOF%IsInitialized = .TRUE.
    IOF%ExistedOnOpen = .FALSE.
    IOF%IFormat = DefaultIFormat
    IOF%RFormat = DefaultRFormat
    IOF%DFormat = DefaultDFormat
    IOF%CFormat = DefaultCFormat
    IOF%DCFormat = DefaultDCFormat
    IOF%SFormat = DefaultSFormat
    IOF%FileFormat = DefaultFileFormat
    IOF%EOF = .FALSE.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Error definitions: Not set up yet.
    ! ArrayTooSmall - Error evaluates true if arrays passed to read routines are not big enough
    !                 to hold the data in the file.
    IOF%IOErrors(1)%ErrorName = 'ArrayTooSmall'
    IOF%IOErrors(1)%Message(1)    = 'Array too small to hold data in file ' // FileName // '.'
    IOF%IOErrors(1)%Action        = 'WARN' 
    IOF%IOErrors(1)%ErrorOccured  = .FALSE.
    IOF%IOErrors(1)%NErrLines     = 1

    ! ArrayTooLarge - Error evaluates true if arrays passed to read routines are larger than
    !                 the data in the file.
    IOF%IOErrors(2)%ErrorName = 'ArrayTooLarge'
    IOF%IOErrors(2)%Message(1)    = 'Array too large when compared to data in file ' // FileName // '.'
    IOF%IOErrors(2)%Action        = 'WARN' 
    IOF%IOErrors(2)%ErrorOccured  = .FALSE.
    IOF%IOErrors(2)%NErrLines     = 1

  END SUBROUTINE InitIOFile
    
  ! Add a file to the stack
  SUBROUTINE AddIOFile(FileName)
    CHARACTER*(*),INTENT(IN) :: FileName
    TYPE(IOFile) IOF

    ! If file exists, do nothing
    IF(IndexIOFile(FileName).ne.0) RETURN 

    ! initialize IOFile
    CALL InitIOFile(IOF, FileName)
    
    ! Add IOFile to stack.
    NFiles = NFiles + 1
    FileStack(NFiles) = IOF
    FileStack(NFiles)%IsInitialized = .TRUE.
  END SUBROUTINE AddIOFile

  ! Find the index of a file
  INTEGER FUNCTION IndexIOFile(FileName)
    CHARACTER*(*) FileName
    INTEGER i1
    IndexIOFile = 0
    DO i1 = 1, NFiles
       IF(TRIM(FileStack(i1)%FileName).eq.TRIM(FileName)) THEN
          IndexIOFile = i1
          RETURN
       END IF
    END DO
  END FUNCTION IndexIOFile


  ! Remove a file from the stack
  SUBROUTINE DeleteIOFile(FileName)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER(10) TmpStr

    ! Local variable index is index of file to delete
    INTEGER Index
    INTEGER i1
    
    TmpStr = TRIM(ADJUSTL(FileName))
    CALL Upper(TmpStr)
    
    IF(TmpStr(1:3).eq.'ALL') THEN       
       NFiles = 0
    ELSE
       Index = IndexIOFile(FileName)
       ! If file doesn't exist, do nothing
       IF(Index.eq.0) RETURN
       
       ! Shift remaining files left
       DO i1 = Index, NFiles - 1
          FileStack(i1) = FileStack(i1+1)
       END DO
       
       NFiles = NFiles - 1
    END IF
  END SUBROUTINE DeleteIOFile

  ! Set options for IOFile
  SUBROUTINE SetIOFileInfo(FileName,UnitNumber, DataTypeLine, IFormat, RFormat, &
       & DFormat, CFormat, DCFormat, SFormat, FileFormat, ExistedOnOpen, IsOpen, &
       & FileAction, NSections,EOF)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: DataTypeLine, IFormat, RFormat, DFormat, &
         & CFormat, DCFormat, SFormat, FileFormat, FileAction
    LOGICAL, INTENT(IN), OPTIONAL :: ExistedOnOpen, IsOpen, EOF
    INTEGER, INTENT(IN), OPTIONAL :: UnitNumber, NSections

    CHARACTER(1000) TmpStr
    INTEGER IOFIndex

    IOFIndex = IndexIOFile(FileName)
    IF(IOFIndex.eq.0) THEN
       CALL Error('Error in function: SetIOFileOptions', StopProgram = .FALSE.)
       CALL Error('File: ' // FileName // ' does not exist.')
    END IF

    ! If arguments are present, set them.
    IF(PRESENT(UnitNumber)) FileStack(IOFIndex)%UnitNumber = UnitNumber
    IF(PRESENT(NSections)) FileStack(IOFIndex)%NSections = NSections
    IF(PRESENT(DataTypeLine)) THEN
       IF(FileStack(IOFIndex)%NSections.gt.0) FileStack(IOFIndex)%DataTypeLine(FileStack(IOFIndex)%NSections) = DataTypeLine
    END IF
    IF(PRESENT(IFormat))  FileStack(IOFIndex)%IFormat  = IFormat
    IF(PRESENT(RFormat))  FileStack(IOFIndex)%RFormat  = RFormat
    IF(PRESENT(DFormat))  FileStack(IOFIndex)%DFormat  = DFormat
    IF(PRESENT(CFormat))  FileStack(IOFIndex)%CFormat  = CFormat
    IF(PRESENT(DCFormat)) FileStack(IOFIndex)%DCFormat = DCFormat
    IF(PRESENT(SFormat))  FileStack(IOFIndex)%SFormat  = SFormat
    IF(PRESENT(FileFormat)) THEN
       TmpStr = TRIM(ADJUSTL(FileFormat))
       CALL Upper(TmpStr)
       FileStack(IOFIndex)%FileFormat  = TRIM(TmpStr)
    END IF
    IF(PRESENT(ExistedOnOpen)) FileStack(IOFIndex)%ExistedOnOpen  = ExistedOnOpen
    IF(PRESENT(IsOpen)) FileStack(IOFIndex)%IsOpen  = IsOpen
    IF(PRESENT(FileAction)) THEN
       TmpStr = TRIM(ADJUSTL(FileAction))
       CALL Upper(TmpStr)
       FileStack(IOFIndex)%IOFileAction = TRIM(TmpStr)
    END IF
    IF(PRESENT(EOF)) FileStack(IOFIndex)%EOF = EOF
  END SUBROUTINE SetIOFileInfo

  SUBROUTINE GetIOFileInfo(FileName, UnitNumber, DataTypeLine, IFormat, RFormat, &
       & DFormat, CFormat, DCFormat, SFormat, FileFormat, ExistedOnOpen, IsOpen, &
       & FileAction, NSections, EOF)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: DataTypeLine, IFormat, RFormat, DFormat, &
         & CFormat, DCFormat, SFormat, FileFormat, FileAction
    LOGICAL, INTENT(OUT), OPTIONAL :: ExistedOnOpen, IsOpen, EOF
    INTEGER, INTENT(OUT), OPTIONAL :: UnitNumber, NSections

    INTEGER IOFIndex

    IOFIndex = IndexIOFile(FileName)

    IF(IOFIndex.eq.0) THEN
       CALL Error('Error in function: GetIOFileInfo', StopProgram = .FALSE.)
       CALL Error('File: ' // FileName // ' does not exist.')
    END IF

    ! If arguments are present, set them.
    IF(PRESENT(UnitNumber)) UnitNumber = FileStack(IOFIndex)%UnitNumber
    IF(PRESENT(DataTypeLine)) THEN
       IF(FileStack(IOFIndex)%NSections.gt.0) THEN
          DataTypeLine = FileStack(IOFIndex)%DataTypeLine(FileStack(IOFIndex)%NSections)
       END IF
    END IF
    IF(PRESENT(IFormat))  IFormat  = FileStack(IOFIndex)%IFormat
    IF(PRESENT(RFormat))  RFormat  = FileStack(IOFIndex)%RFormat
    IF(PRESENT(DFormat))  DFormat  = FileStack(IOFIndex)%DFormat
    IF(PRESENT(CFormat))  CFormat  = FileStack(IOFIndex)%CFormat
    IF(PRESENT(DCFormat)) DCFormat = FileStack(IOFIndex)%DCFormat
    IF(PRESENT(SFormat))  SFormat  = FileStack(IOFIndex)%SFormat 

    IF(PRESENT(FileFormat)) FileFormat = FileStack(IOFIndex)%FileFormat
    IF(PRESENT(ExistedOnOpen)) ExistedOnOpen = FileStack(IOFIndex)%ExistedOnOpen
    IF(PRESENT(IsOpen)) IsOpen = FileStack(IOFIndex)%IsOpen
    IF(PRESENT(NSections)) NSections = FileStack(IOFIndex)%NSections
    IF(PRESENT(EOF)) EOF = FileStack(IOFIndex)%EOF
    IF(PRESENT(FileAction)) FileAction = FileStack(IOFIndex)%IOFileAction
  END SUBROUTINE GetIOFileInfo
  
  SUBROUTINE PrintIOFileInfo(IOF)
    TYPE(IOFile) IOF
    INTEGER i1

    PRINT '(A)', 'Information for file ' // TRIM(IOF%FileName)
    PRINT *, 'File is initialized ?: ', IOF%IsInitialized
    PRINT*, 'File is open?: ', IOF%IsOpen
    PRINT*, 'Files existed on open ?: ', IOF%ExistedOnOpen
    PRINT*, 'Unit Number: ', IOF%UnitNumber
    PRINT*, 'End of file?: ', IOF%EOF
    PRINT '(2A)', 'IFormat: ', IOF%IFormat
    PRINT '(2A)', 'RFormat: ', IOF%RFormat
    PRINT '(2A)', 'DFormat: ', IOF%DFormat
    PRINT '(2A)', 'CFormat: ', IOF%CFormat
    PRINT '(2A)', 'DCFormat: ',IOF%DCFormat
    PRINT '(2A)', 'SFormat: ', IOF%SFormat
    PRINT '(2A)', 'FileFormat: ', IOF%FileFormat
    PRINT '(2A)', 'FileAction: ', IOF%IOFileAction
    PRINT '(A,I10)', 'Total number of Lines: ', IOF%TotLines
    PRINT '(A,I10)', 'Number of sections: ', IOF%NSections

    PRINT '(A)', 'Section information:'
    PRINT '(A)', 'Section  NLines  DataTypeLine'
    DO i1 = 1, IOF%NSections
       PRINT '(2I5,1X,A)', i1, IOF%NLines(i1), IOF%DataTypeLine(i1)
    END DO
  END SUBROUTINE PrintIOFileInfo
  
  SUBROUTINE PrintFileStackInfo
    INTEGER i1
    PRINT*, 'Number of files: ', NFiles
    DO i1 = 1, NFiles
       CALL PrintIOFileInfo(FileStack(i1))
    END DO
  END SUBROUTINE PrintFileStackInfo

  LOGICAL FUNCTION ReadingNewSection(FileName,NSections)
    CHARACTER*(*), INTENT(IN) :: FileName
    INTEGER,INTENT(IN) :: NSections
    INTEGER IOFIndex

    IOFIndex = IndexIOFile(FileName)

    IF(IOFIndex.eq.0) THEN
       CALL Error('Error in function IsNewSection. File: ' &
            & // FileName // ' does not exist.')
    ELSE
       
       IF(FileStack(IOFIndex)%NSections.eq.NSections+1) THEN
          ReadingNewSection = .TRUE.
       ELSE
          ReadingNewSection = .FALSE.
       END IF
    END IF
  END FUNCTION ReadingNewSection
    
  LOGICAL FUNCTION WritingNewSection(FileName,DataTypeLine)
    CHARACTER*(*), INTENT(IN) :: FileName, DataTypeLine
    INTEGER IOFIndex, NSections

!    write(*,*) 'filename,datatypeline',filename,datatypeline
!	write(*,*) 'filestack',filestack
    IOFIndex = IndexIOFile(FileName)
    IF(IOFIndex.eq.0) THEN
       CALL Error('Error in function IsNewSection. File: ' &
            & // FileName // ' does not exist.')
    ELSE
       NSections = FileStack(IOFIndex)%NSections
       IF(NSections.le.0) THEN
          WritingNewSection = .FALSE.
          RETURN
       END IF
       IF(TRIM(DataTypeLine).eq. &
            & TRIM(FileStack(IOFIndex)%DataTypeLine(NSections)).or.&
            & LEN_TRIM(DataTypeLine).eq.0) THEN
          WritingNewSection = .FALSE.
       ELSE
          WritingNewSection = .TRUE.
       END IF
    END IF
  END FUNCTION WritingNewSection

  SUBROUTINE GetFileStackInfo(NumberOfFiles, FileNames, UnitNumbers)
    INTEGER,INTENT(OUT) :: NumberOfFiles
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: FileNames(:)
    INTEGER,INTENT(OUT),OPTIONAL :: UnitNumbers(:)

    INTEGER i1

    IF(PRESENT(FileNames)) THEN
       ! Check the size of the array.
       IF(SIZE(FileNames).gt.NFiles) CALL Error('Array to hold file names passed to ' // &
            & 'GetFileStackInfo is not large enough.')
       
       DO i1 = 1, NFiles
          FileNames(i1) = FileStack(i1)%FileName
       END DO
    END IF

    IF(PRESENT(UnitNumbers)) THEN
       ! Check the size of the array.
       IF(SIZE(UnitNumbers).gt.NFiles) CALL Error('Array to hold file names passed to ' // &
            & 'GetFileStackInfo is not large enough.')
       
       DO i1 = 1, NFiles
          UnitNumbers(i1) = FileStack(i1)%UnitNumber
       END DO
    END IF
    NumberOfFiles = NFiles
  END SUBROUTINE GetFileStackInfo

! Find the index of an ioerror.
  INTEGER FUNCTION IndexIOError(Filename,ErrorName)
    CHARACTER*(*) FileName, ErrorName
    INTEGER IOFIndex, i1
    
    IOFIndex = IndexIOFile(FileName)
    ! If file doesn't exist, error
    IF(IOFIndex.eq.0) CALL Error('ERROR: IOFile does not exist for file ' // FileName // '.')

    DO i1 = 1, NError
       IF(TRIM(FileStack(IOFIndex)%IOErrors(i1)%ErrorName).eq.TRIM(ErrorName)) THEN
          IndexIOError = i1
          RETURN
       END IF
    END DO
  END FUNCTION IndexIOError

  SUBROUTINE SetIOErrorInfo(FileName, ErrorName, Message, Action, ErrorOccured)
    CHARACTER*(*),INTENT(IN) :: FileName, ErrorName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Message(:), Action
    LOGICAL,INTENT(IN),OPTIONAL :: ErrorOccured
    INTEGER IOErrorIndex, IOFIndex, i1, NMessage, NErrLines

    IOFIndex = IndexIOFile(FileName)
    IOErrorIndex = IndexIOError(FileName,ErrorName)
    NErrLines = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%NErrLines

    IF(PRESENT(Message)) THEN
       NMessage = MAX(SIZE(Message)+NErrLines,MxErrLines)
       DO i1 = NErrLines, NMessage
          FileStack(IOFIndex)%IOErrors(IOErrorIndex)%Message(i1) = Message(i1-NErrLines+1)
       END DO
    END IF
    FileStack(IOFIndex)%IOErrors(IOErrorIndex)%NErrLines = NMessage
    IF(PRESENT(Action))  FileStack(IOFIndex)%IOErrors(IOErrorIndex)%Action = Action
    IF(PRESENT(ErrorOccured)) FileStack(IOFIndex)%IOErrors(IOErrorIndex)%ErrorOccured = ErrorOccured
  END SUBROUTINE SetIOErrorInfo

SUBROUTINE GetIOErrorInfo(FileName, ErrorName, Message, Action, ErrorOccured, NErrLines)
    CHARACTER*(*),INTENT(IN) :: FileName, ErrorName
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Message(:), Action
    LOGICAL,INTENT(OUT),OPTIONAL :: ErrorOccured
    INTEGER,INTENT(OUT),OPTIONAL :: NErrLines
    INTEGER IOErrorIndex, IOFIndex, i1, NLines

    IOFIndex = IndexIOFile(FileName)
    IOErrorIndex = IndexIOError(FileName,ErrorName)
    NLines = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%NErrLines

    IF(PRESENT(Message)) THEN
       DO i1 = 1, NLines
          Message(i1) = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%Message(i1)
       END DO
    END IF

    IF(PRESENT(Action))  Action = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%Action
    IF(PRESENT(ErrorOccured))  ErrorOccured = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%ErrorOccured
    IF(PRESENT(NErrLines))  NErrLines = FileStack(IOFIndex)%IOErrors(IOErrorIndex)%NErrLines 
  END SUBROUTINE GetIOErrorInfo
END MODULE IOFiles
