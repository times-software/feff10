!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_iomod.f90,v $:
! $Revision: 1.11 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE IOMod
  USE ErrorMod
  USE IOFiles
  USE PADIO

  IMPLICIT NONE

  INTEGER,PRIVATE :: npadxDefault, npadrDefault, MaxStrLen, iFileTypeDefault
  PARAMETER(npadxDefault = 8, npadrDefault = 8, MaxStrLen = 80, iFileTypeDefault = 1)
  LOGICAL,PRIVATE,SAVE :: WriteDataDescription(100)
  CHARACTER(30),PRIVATE :: DefaultFileStatus, DefaultFilePosition, &
       & DefaultFileAction, DefaultCommentCharacters
  PARAMETER(DefaultFileStatus = 'UNKNOWN', DefaultFilePosition = 'REWIND')
  PARAMETER(DefaultFileAction = 'READWRITE', DefaultCommentCharacters = '#!cC')

  INTERFACE Write2D
     MODULE PROCEDURE WriteInt2D
     MODULE PROCEDURE WriteReal2D
     MODULE PROCEDURE WriteDouble2D
     MODULE PROCEDURE WriteComplex2D
     MODULE PROCEDURE WriteDComplex2D
     MODULE PROCEDURE WriteString2D
  END INTERFACE

  INTERFACE Read2D
     MODULE PROCEDURE ReadInt2D
     MODULE PROCEDURE ReadReal2D
     MODULE PROCEDURE ReadDouble2D
     MODULE PROCEDURE ReadComplex2D
     MODULE PROCEDURE ReadDComplex2D
     MODULE PROCEDURE ReadString2D
  END INTERFACE

CONTAINS

  ! OpenFl finds the first available unit and opens the file if it is not already open.
  ! In addition, a new IOFile is created and added to the FileStack. The FileStack holds 
  ! valuable information about the files that are open. See IOFiles.f90 for more info.
  ! If the file is already open, OpenFl does nothing.
  !
  ! FileName.
  ! FileName     - Name of file to open.
  ! FileStatus   - Optional: Can be any legal fortran 90 file status.
  !                Argument is passed directly to OPEN(iUnit, STATUS = FileStatus ...
  ! FilePosition - Optional: Can be any legal file position.
  !                Agument is passed directly to OPEN(iUnit, POSITION = FilePosition ...
  ! FileAction   - This can be 'READ', 'WRITE', or 'WREADWRITE'. Case insensitive.
  SUBROUTINE OpenFl(FileName, FileStatus, FilePosition, FileAction, ErrorMessage)
    CHARACTER*(*) FileName
    CHARACTER*(*),OPTIONAL :: FileStatus, FilePosition, FileAction, ErrorMessage(:)

    LOGICAL FileIsOpen
    INTEGER,SAVE :: iOpenUnit
    INTEGER IOError, iUnit
    LOGICAL :: UnitIsUsed, FileExists
    CHARACTER(300) messg, FStatus, FPosition, FAction

    FStatus   = DefaultFileStatus
    FPosition = DefaultFilePosition
    FAction   = DefaultFileAction

    IF(PRESENT(FileStatus)) FStatus = FileStatus
    IF(PRESENT(FilePosition)) FPosition = FilePosition
    IF(PRESENT(FileAction)) FAction = FileAction

    ! Is file already open?
    INQUIRE(FILE=FileName,OPENED=FileIsOpen,EXIST=FileExists)
    IF(FileIsOpen) THEN
       ! This file is already open. Do nothing.
       RETURN
    ELSE
       ! File was not open. Find first available unit.
       DO iOpenUnit = 5, 100
          INQUIRE(UNIT=iOpenUnit,OPENED=UnitIsUsed)
          IF(.not.UnitIsUsed) THEN
             iUnit = iOpenUnit
             EXIT
          END IF
       END DO

       ! Open file.
       OPEN(UNIT=iUnit,FILE=FileName,STATUS=FStatus,POSITION=FPosition,ACTION=FAction,IOSTAT=IOError)
       IF(IOError.ne.0) THEN
          ! OPEN passed back an error
          CALL Error('Error opening file: ' // FileName, StopProgram = .FALSE.)
          IF(PRESENT(ErrorMessage)) CALL Error(ErrorMessage, StopProgram = .FALSE.)
          WRITE(messg,'(A,i4)') 'OPEN returned error number ', IOError
          CALL Error(messg)
       END IF

       ! Create new IOFile
       CALL AddIOFile(FileName)

       ! Set file options 
       CALL SetIOFileInfo(FileName,ExistedOnOpen = FileExists, UnitNumber = iUnit,IsOpen = .TRUE., &
            & FileAction = FAction)
    END IF
  END SUBROUTINE OpenFl

  ! Subroutine CloseFl closes the file associated with FileName, and deletes the IOFile
  ! from the filestack. If FileName = 'ALL' (case insensitive), all open files are closed
  ! and deleted from the filestack.
  SUBROUTINE CloseFl(FileName)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER(10) FileAction, TmpStr
    INTEGER iUnit, iUnits(100), NFiles, i1

    TmpStr = TRIM(ADJUSTL(FileName))
    CALL Upper(TmpStr)

    IF(TmpStr.ne.'ALL') THEN
       CALL GetIOFileInfo(FileName,UnitNumber = iUnit)
       CLOSE(iUnit)
       CALL DeleteIOFile(FileName)
    ELSE
       CALL GetFileStackInfo(NFiles, UnitNumbers = iUnits)
       DO i1 = 1, NFiles
          CLOSE(iUnits(i1))
       END DO
       CALL DeleteIOFile(FileName)
    END IF
  END SUBROUTINE CloseFl


  ! Subroutine WriteData writes scalar data to the file filename.
  ! Up to 10 pieces of data can be written each specified by using
  ! the optional arguments Int1 - Int10, Real1 - Real10, etc.  
  ! Thus you can write all integer data, all real data, etc. or a mixture 
  ! of various types. The column that each data piece goes in is set by
  ! specifying the number of each type, for example, a call to write a
  ! real number followed by an integer is as follows:
  !
  ! CALL WriteData(FileName, Real1 = RealVar, Int2 = IntVar) 
  !
  ! To write a complex variable followed by a string variable then a real 
  ! variable,
  !
  ! CALL WriteData(FileName, Complex1 = ComplexVar, String2 = StringVar, Real3 = RealVar)
  !
  ! Assigning two variables to the same column will produce an error, i.e.
  !
  ! CALL WriteData(FileName, Complex1 = ComplexVar, Int1 = IntVar)
  !
  ! will produce an error since we are trying to write two variables to the first column.
  !
  ! WriteData will also write a section liine (beginning with #SN#) and a data type line
  ! (beginning with #DT#). New sections are started when the types or number of the data 
  ! arguments are changed, or the optional ForceNewSection argument is set to true.
  ! In addition, the user can specify an array of headers to be written above the data,
  ! as well as an array of column labels to describe the data. 
  !
  ! Arguments
  ! FileName - Name of file to write to. Note that if the file is already open for reading,
  !            you will get an error. Note that this is the only non-optional argument.
  !
  ! The following arguments give the data to be written to the file. Note that all of 
  ! these arguments are optional.
  ! Int1, Real1, ... String1    - Argument to go in the first column.
  ! Int2, Real2, ... String2    - Argument to go in the second column.
  ! .
  ! .
  ! .
  ! Int10, Real10, ... String10 - Argument to go in the tenth column.
  !
  ! Headers - Optional: Array of headers to be written out above the data. Note that
  !           this argument can be given without any data, and the headers will be written
  !           alone. Also, headers will only be written at the beginning of a section, so
  !           that the user can set headers and use WriteData in a loop, and the headers
  !           will only be written in the first call inside the loop.
  !
  ! ColumnLabels - Optional: Array of strings to describe the data being passed to 
  !                WriteData. These will be written directly above the data in a formatted
  !                manner. Note that column labels will not extend past the length of the 
  !                data in a particular column. Thus the length of each column label is 
  !                restricted by the formats which are used to write that data. The default 
  !                formats can be found in IOFiles.f90, and may be set for a particular 
  !                file by using the SetIOFileInfo subroutine.
  !
  ! WriteDataInHeader - When set to .TRUE. data will be written with #HD# inserted at the 
  !                     beginning of the line. This is usefull for files which contain
  !                     various information, including scalar information (which would be
  !                     written in the header) as well a set of arrays written in columns
  !                     which the user would like to be able to plot.
  !
  ! ForceNewSection - In the case where two sets of data have the same types. The user may
  !                   force a new section so that the section line, data type line,any 
  !                   headers, and the column labels will be written.
  SUBROUTINE WriteData(FileName,                                      &
       & Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, Int9, Int10, &
       & Real1, Real2, Real3, Real4, Real5, Real6, Real7, Real8,      &
       & Real9, Real10,                                               &
       & Double1, Double2, Double3, Double4, Double5, Double6,        &
       & Double7, Double8, Double9, Double10,                         &
       & Complex1, Complex2, Complex3, Complex4, Complex5, Complex6,  &
       & Complex7, Complex8, Complex9, Complex10,                     &
       & DComplex1, DComplex2, DComplex3, DComplex4, DComplex5,       &
       & DComplex6, DComplex7, DComplex8, DComplex9, DComplex10,      &
       & String1, String2, String3, String4, String5, String6,        &
       & String7, String8, String9, String10,                         &
       & Headers, ColumnLabels, WriteDataInHeader, ForceNewSection,   &
       & FileType)
    ! Passed variables
    CHARACTER*(*), INTENT(IN) :: FileName
    INTEGER, INTENT(IN), OPTIONAL :: Int1, Int2, Int3, Int4, Int5,    &
         & Int6, Int7, Int8, Int9, Int10
    REAL, INTENT(IN), OPTIONAL :: Real1, Real2, Real3, Real4, Real5,  &
         & Real6, Real7, Real8, Real9, Real10
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: Double1, Double2,       &
         & Double3, Double4, Double5, Double6, Double7, Double8,      &
         & Double9, Double10
    COMPLEX, INTENT(IN), OPTIONAL :: Complex1, Complex2, Complex3,    &
         & Complex4, Complex5,                                        &
         & Complex6, Complex7, Complex8, Complex9, Complex10
    COMPLEX*16, INTENT(IN), OPTIONAL :: DComplex1, DComplex2,         &
         & DComplex3, DComplex4, DComplex5, DComplex6, DComplex7,     &
         & DComplex8, DComplex9, DComplex10
    CHARACTER*(*), INTENT(IN), OPTIONAL :: String1, String2, String3, &
         & String4, String5, String6, String7, String8, String9,      &
         & String10

    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), ColumnLabels(:), &
         & FileType
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader, ForceNewSection

    ! Local variables
    ! NumArgs - Number of data arguments passed
    ! iUnit   - unit number associated with FileName
    ! NColumnLabels - will hold the size of the ColumnLabels array.
    INTEGER NumArgs, iUnit, NColumnLabels

    ! IsPresent  - If IsPresent(i) is true then the ith argument has already been defined.
    ! HeaderData - If true, write data in header format (beginning with #HD#).
    ! txt, pad  - These are here in anticipation of enabling various formats for output,
    !              including text, paded ascii to start, and possibly further formats 
    !              such as xml.
    ! NewSection - If true, force a new section.
    ! FlType     - 
    LOGICAL IsPresent(10), HeaderData, txt, pad, NewSection
    ! DataArray - Array of strings to hold the data. Data is first written to DataArray,
    !             then DataArray(1:NumArgs) is written to the file.
    CHARACTER(MaxStrLen) DataArray(10)
    ! DataTypeLine - The data type line lists the types of data that have been specified
    !                by the user to be written to file. This line will be written above
    !                the data, and above any column label line if present. The line 
    !                begins with #DT#
    CHARACTER(120) DataTypeLine
    ! IntFormat, RealFormat, ... StringFormat - These hold the format strings which will
    !                                           be used to write the data to the 
    !                                           character array DataArray and eventually
    !                                           to the file.
    !
    ! FileAction - String variable that holds the action associated with the IOFile
    !              specified by FileName.
    CHARACTER(20) IntFormat, RealFormat, DoubleFormat, ComplexFormat, &
         & DComplexFormat, StringFormat, FileAction
    ! ColumnLabelFormat - This will be used to hold the format string to write the
    !                     column labels if any.
    CHARACTER(80) ColumnLabelFormat
    ! CLFormats - This holds the separate formats for each column label. They will be
    !             combined to create ColumnLabelFormat
    CHARACTER(10) CLFormats(20)
    ! TmpStr - used when creating column label formats.
    CHARACTER(4) TmpStr, FlType

    ! NSection - number of the current section.
    INTEGER NSections

    ! Loop variables
    INTEGER i1

    ! Initialization
    HeaderData    = .FALSE.
    IsPresent(:)  = .FALSE.
    NumArgs       = 0
    txt           = .TRUE.
    pad           = .FALSE.
    DataTypeLine  = ''
    FlType = 'TXT'

    ! Set the type of file. Default is txt.
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          txt = .TRUE.
          pad = .FALSE.
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          pad = .TRUE.
          txt = .FALSE.
       ELSE
          CALL Error('Error: Illegal file type passed to WriteData')
       END IF
    END IF

!KJ next section was commented out "!NS" but I activated it again 7-09
    IF(PRESENT(ForceNewSection)) THEN
       NewSection = ForceNewSection
    ELSE
       NewSection = .FALSE.
    END IF
!KJ

    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader

    IF(PRESENT(ColumnLabels)) THEN
       NColumnLabels = 0
    ELSE
       NColumnLabels = -1
    END IF

    ! Open the file.
    CALL OpenFl(FileName, FileStatus = 'REPLACE', FileAction = 'WRITE')

    ! Get the formats for each data type.
    CALL GetIOFileInfo(FileName, IFormat = IntFormat, RFormat = RealFormat, &
         & DFormat = DoubleFormat, CFormat = ComplexFormat, &
         & DCFormat = DComplexFormat, SFormat = StringFormat, FileAction = FileAction)

    ! Check that this file was not already opened for READ.
    IF(TRIM(FileAction).eq.'READ') THEN 
       CALL Error('Error: file ' // FileName // &
            & ' was already opened for reading.', StopProgram = .FALSE.)
       CALL Error("Please use CALL CloseFl('" // FileName // "') before writing.")
    END IF
    
    ! The SetArgType(Arg,ColumnNumber) subroutines checks the presence of Arg. If Arg
    ! is present, the data is written to DataArray(ColumnNumber). Also, NumArgs is incremented,
    ! CLFormats(ColumnNumber) is set, and the data type is appended to DataTypeLine.
    ! In addition, various errors are checked.
    ! Find the 1st argument.
    IF(PRESENT(Int1)) CALL SetArgInt(Int1, 1)
    IF(PRESENT(Real1)) CALL SetArgReal(Real1, 1)
    IF(PRESENT(Double1)) CALL SetArgDouble(Double1, 1)
    IF(PRESENT(Complex1)) CALL SetArgComplex(Complex1, 1)
    IF(PRESENT(DComplex1)) CALL SetArgDComplex(DComplex1, 1)
    IF(PRESENT(String1)) CALL SetArgString(String1, 1)

    ! Find the 2nd argument.
    IF(PRESENT(Int2)) CALL SetArgInt(Int2, 2)
    IF(PRESENT(Real2)) CALL SetArgReal(Real2, 2)
    IF(PRESENT(Double2)) CALL SetArgDouble(Double2, 2)
    IF(PRESENT(Complex2)) CALL SetArgComplex(Complex2, 2)
    IF(PRESENT(DComplex2)) CALL SetArgDComplex(DComplex2, 2)
    IF(PRESENT(String2)) CALL SetArgString(String2, 2)

    ! Find the 3rd argument.
    IF(PRESENT(Int3)) CALL SetArgInt(Int3, 3)
    IF(PRESENT(Real3)) CALL SetArgReal(Real3, 3)
    IF(PRESENT(Double3)) CALL SetArgDouble(Double3, 3)
    IF(PRESENT(Complex3)) CALL SetArgComplex(Complex3, 3)
    IF(PRESENT(DComplex3)) CALL SetArgDComplex(DComplex3, 3)
    IF(PRESENT(String3)) CALL SetArgString(String3, 3)

    ! Find the 4th argument.
    IF(PRESENT(Int4)) CALL SetArgInt(Int4, 4)
    IF(PRESENT(Real4)) CALL SetArgReal(Real4, 4)
    IF(PRESENT(Double4)) CALL SetArgDouble(Double4, 4)
    IF(PRESENT(Complex4)) CALL SetArgComplex(Complex4, 4)
    IF(PRESENT(DComplex4)) CALL SetArgDComplex(DComplex4, 4)
    IF(PRESENT(String4)) CALL SetArgString(String4, 4)

    ! Find the 5th argument.
    IF(PRESENT(Int5)) CALL SetArgInt(Int5, 5)
    IF(PRESENT(Real5)) CALL SetArgReal(Real5, 5)
    IF(PRESENT(Double5)) CALL SetArgDouble(Double5, 5)
    IF(PRESENT(Complex5)) CALL SetArgComplex(Complex5, 5)
    IF(PRESENT(DComplex5)) CALL SetArgDComplex(DComplex5, 5)
    IF(PRESENT(String5)) CALL SetArgString(String5, 5)

    ! Find the 6th argument.
    IF(PRESENT(Int6)) CALL SetArgInt(Int6, 6)
    IF(PRESENT(Real6)) CALL SetArgReal(Real6, 6)
    IF(PRESENT(Double6)) CALL SetArgDouble(Double6, 6)
    IF(PRESENT(Complex6)) CALL SetArgComplex(Complex6, 6)
    IF(PRESENT(DComplex6)) CALL SetArgDComplex(DComplex6, 6)
    IF(PRESENT(String6)) CALL SetArgString(String6, 6)

    ! Find the 7th argument.
    IF(PRESENT(Int7)) CALL SetArgInt(Int7, 7)
    IF(PRESENT(Real7)) CALL SetArgReal(Real7, 7)
    IF(PRESENT(Double7)) CALL SetArgDouble(Double7, 7)
    IF(PRESENT(Complex7)) CALL SetArgComplex(Complex7, 7)
    IF(PRESENT(DComplex7)) CALL SetArgDComplex(DComplex7, 7)
    IF(PRESENT(String7)) CALL SetArgString(String7, 7)

    ! Find the 8th argument.
    IF(PRESENT(Int8)) CALL SetArgInt(Int8, 8)
    IF(PRESENT(Real8)) CALL SetArgReal(Real8, 8)
    IF(PRESENT(Double8)) CALL SetArgDouble(Double8, 8)
    IF(PRESENT(Complex8)) CALL SetArgComplex(Complex8, 8)
    IF(PRESENT(DComplex8)) CALL SetArgDComplex(DComplex8, 8)
    IF(PRESENT(String8)) CALL SetArgString(String8, 8)

    ! Find the 9th argument.
    IF(PRESENT(Int9)) CALL SetArgInt(Int9, 9)
    IF(PRESENT(Real9)) CALL SetArgReal(Real9, 9)
    IF(PRESENT(Double9)) CALL SetArgDouble(Double9, 9)
    IF(PRESENT(Complex9)) CALL SetArgComplex(Complex9, 9)
    IF(PRESENT(DComplex9)) CALL SetArgDComplex(DComplex9, 9)
    IF(PRESENT(String9)) CALL SetArgString(String9, 9)

    ! Find the 10th argument.
    IF(PRESENT(Int10)) CALL SetArgInt(Int10, 10)
    IF(PRESENT(Real10)) CALL SetArgReal(Real10, 10)
    IF(PRESENT(Double10)) CALL SetArgDouble(Double10, 10)
    IF(PRESENT(Complex10)) CALL SetArgComplex(Complex10, 10)
    IF(PRESENT(DComplex10)) CALL SetArgDComplex(DComplex10, 10)
    IF(PRESENT(String10)) CALL SetArgString(String10, 10)

    ! Get the unit number and section number.
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit, NSections = NSections)

    ! If data type line is different that last, or user has forced a new section
    ! write a new section header and data type line. Also write any headers and
    ! column labels that have been specified.
    IF(WritingNewSection(FileName,DataTypeLine).or.NewSection) THEN

       NSections = NSections + 1

       WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections
       WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
       WRITE(iUnit,'(A)') '#H#'

       ! Only write a data type line if it is not blank       
       IF(LEN_TRIM(DataTypeLine).ne.0) THEN
          WRITE(iUnit,'(A)') '#H# The following data types are written in this section.'
          WRITE(iUnit,'(2A)') '#DT# ', TRIM(DataTypeLine)
       END IF

       ! If headers are present, write them.
       IF(PRESENT(Headers)) THEN
          DO i1 = 1, SIZE(Headers)
             ! Write only non-blank, non-null headers.
             IF((LEN_TRIM(Headers(i1)).gt.0).and.(ICHAR(Headers(i1)).ne.0)) &
                  & WRITE(iUnit,'(2A)') '#H# ', TRIM(Headers(i1))
          END DO
       END IF

       ! If column labels are present, get the columnlable format and write them.
       IF(PRESENT(ColumnLabels)) THEN
          ColumnLabelFormat = '(A5, '
          DO i1 = 1, NColumnLabels - 1
             ColumnLabelFormat = TRIM(ColumnLabelFormat) // TRIM(CLFormats(i1)) // ', '
          END DO
          ColumnLabelFormat = TRIM(ColumnLabelFormat) // TRIM(CLFormats(NColumnLabels)) // ')'
          WRITE(iUnit,ColumnLabelFormat) '#CL# ', (TRIM(ColumnLabels(i1)) // ' ', i1 = 1, NColumnLabels)
       END IF

       ! Set the section number and data type line.
       CALL SetIOFileInfo(FileName, DataTypeLine = DataTypeLine, NSections = NSections)

    ELSEIF(PRESENT(Headers).and.(NumArgs.eq.0)) THEN
       ! If only headers are present, write the headers without changing the section number.
       DO i1 = 1, SIZE(Headers)
          ! Write only non-blank, non-null headers.
          IF((LEN_TRIM(Headers(i1)).gt.0).and.(ICHAR(Headers(i1)).ne.0)) &
               & WRITE(iUnit,'(2A)') '#H# ', TRIM(Headers(i1))
       END DO
       CALL SetIOFileInfo(FileName, DataTypeLine = DataTypeLine)
    END IF


    ! Write data line.
    IF(NumArgs.gt.0) THEN
       IF(HeaderData) THEN
          WRITE(iUnit, '(11A)') '#HD# ', (TRIM(DataArray(i1)) // ' ', i1 = 1, NumArgs)
       ELSE
          WRITE(iUnit, '(11A)') (TRIM(DataArray(i1)) // ' ', i1 = 1, NumArgs)
       END IF
    END IF
  CONTAINS

    ! The SetArgType(Arg,ColumnNumber) subroutines check the presence of Arg. If Arg
    ! is present, the data is written to DataArray(ColumnNumber). Also, NumArgs is incremented,
    ! CLFormats(ColumnNumber) is set, and the data type is appended to DataTypeLine.
    ! In addition, various errors are checked.

    SUBROUTINE SetArgInt(Arg, iArg)
      INTEGER,INTENT(IN) :: Arg
      INTEGER,INTENT(IN) :: iArg

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Int')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), IntFormat) Arg
      ELSEIF(pad) THEN
         CALL WritePAD(DataArray(iArg),Arg)
      END IF

      ! Set the columnlabel format string if ColumnLabels is passed.
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels) THEN
            NColumnLabels = NColumnLabels + 1
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) - 4
            ELSE
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) + 1
            END IF
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF

    END SUBROUTINE SetArgInt

    SUBROUTINE SetArgReal( Arg, iArg)
      REAL,INTENT(IN) :: Arg
      INTEGER,INTENT(IN) :: iArg
      REAL TmpArg(2)
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Real')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), RealFormat) Arg
      ELSEIF(pad) THEN
         TmpArg(1) = Arg
         CALL WritePAD(DataArray(iArg),npadrDefault,Arg)
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels) THEN
            NColumnLabels = NColumnLabels + 1
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) - 4
            ELSE
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) + 1
            END IF
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF
    END SUBROUTINE SetArgReal

    SUBROUTINE SetArgDouble( Arg, iArg)
      DOUBLE PRECISION,INTENT(IN) :: Arg
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Double')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), DoubleFormat) Arg
      ELSEIF(pad) THEN
         CALL WritePAD(DataArray(iArg),npadrDefault,Arg)
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels) THEN
            NColumnLabels = NColumnLabels + 1
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) - 4
            ELSE
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) + 1
            END IF
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF
    END SUBROUTINE SetArgDouble

    SUBROUTINE SetArgComplex( Arg, iArg)
      COMPLEX,INTENT(IN) :: Arg
      INTEGER,INTENT(IN) :: iArg

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Complex')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), ComplexFormat) REAL(Arg), IMAG(Arg)
      ELSEIF(pad) THEN
         CALL WritePAD(DataArray(iArg),npadrDefault,Arg)
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels+1) THEN
            NColumnLabels = NColumnLabels + 2
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) - 4
               CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) + 1
               CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
            ELSE
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) + 1
               CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
               CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
            END IF
         END IF
      END IF
    END SUBROUTINE SetArgComplex

    SUBROUTINE SetArgDComplex( Arg, iArg)
      COMPLEX*16,INTENT(IN) :: Arg
      INTEGER,INTENT(IN) :: iArg

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'DComplex')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), DComplexFormat) DBLE(Arg), DIMAG(Arg)
      ELSEIF(pad) THEN
         CALL WritePAD(DataArray(iArg),npadrDefault,Arg)
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels+1) THEN
            NColumnLabels = NColumnLabels + 2
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) - 4
               CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) + 1
               CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
            ELSE
               WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(iArg))/2) + 1
               CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
               CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
            END IF
         END IF
      END IF
    END SUBROUTINE SetArgDComplex

    SUBROUTINE SetArgString( Arg, iArg)
      CHARACTER*(*),INTENT(IN) :: Arg
      CHARACTER(MaxStrLen) TmpStr
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'String')

      TmpStr = TRIM(ADJUSTL(Arg))
      IF(LEN_TRIM(Arg).gt.MaxStrLen) CALL Error('Warning: String argument is too long in WriteData. ' // &
           & 'Will be truncated.', StopProgram = .FALSE.)

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         WRITE(DataArray(iArg), StringFormat) TRIM(TmpStr)
      ELSEIF(pad) THEN
         CALL WritePAD(DataArray(iArg),TRIM(TmpStr))
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         IF(SIZE(ColumnLabels).gt.NColumnLabels) THEN
            NColumnLabels = NColumnLabels + 1
            IF((iArg.eq.1).and.txt) THEN
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) - 4
            ELSE
               WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(iArg)) + 1
            END IF
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF
    END SUBROUTINE SetArgString

    ! CheckForErrors checks for the following errors and stops if it finds one:
    ! 
    ! 1. 'Error in WriteData multiple arguments in place' - This error occurs 
    !    when the user has specified that two arguments should be written into
    !    one column, i.e. WriteData(FileName, Int1 = I, Real1 = R, ...)
    !
    ! 2. 'Too many arguments passed to WriteData. Max is 10.' - user has specified
    !    more than 10 aguments. This one is redundant, and error 1. should always 
    !    occur if this happens.
    !
    ! In addition, CheckForErrors appends the type to DataTypeLine
    SUBROUTINE CheckForErrors(iArg,DTL)
      INTEGER,INTENT(IN) :: iArg
      CHARACTER*(*) DTL
      CHARACTER(2) ArgNum
      INTEGER iAllocateError
      ArgNum=''
      WRITE(ArgNum,'(I2)') iArg

      IF(IsPresent(iArg)) THEN
         CALL Error('Error in WriteData multiple arguments in place ' // TRIM(ArgNum) // '.')
      END IF

      ! Set IsPresent(iArg)
      IsPresent(iArg) = .TRUE.

      ! Check that number of arguments is not too large. Redundant, but oh well.
      IF(NumArgs.ge.iArg) CALL Error('Too many arguments passed to WriteData. Max is 10.')
      NumArgs = NumArgs + 1

      ! Set data type line
      DataTypeLine = TRIM(DataTypeLine) // ' ' // TRIM(DTL)
    END SUBROUTINE CheckForErrors

  END SUBROUTINE WriteData


  ! SUBROUTINE ReadData reads scalar data from the file filename.
  ! Up to 10 pieces of data can be read, each specified by using
  ! the optional arguments Int1 - Int10, Real1 - Real10, etc.  
  ! Thus you can read all integer data, all real data, etc. or a mixture 
  ! of various types. The column that each data piece comes from is set by
  ! specifying the number of each type, for example, a call to read a
  ! real number followed by an integer is as follows:
  !
  ! CALL ReadData(FileName, Real1 = RealVar, Int2 = IntVar) 
  !
  ! To read a complex variable followed by a string variable then a real 
  ! variable,
  !
  ! CALL ReadData(FileName, Complex1 = ComplexVar, String2 = StringVar, Real3 = RealVar)
  !
  ! Assigning two variables to the same column will produce an error, i.e.
  !
  ! CALL ReadData(FileName, Complex1 = ComplexVar, Int1 = IntVar)
  !
  ! will produce an error since we are trying to write two variables to the first column.
  !
  ! ReadData will also read section lines (beginning with #SN#) if they exist and a data 
  ! type lines (beginning with #DT#) if they exist. If the arguments passed to the routine
  ! do not match with the data type line for the current section, or (in the case where no
  ! data type line is present) do not match with the types found in the file, an error will
  ! occur. In addition, the user can specify an array of headers to be read from above the data.
  ! These headers will include any column labels that are defined.
  !
  ! Arguments
  ! FileName - Name of file to read from. Note that if the file is already open for writing,
  !            you will get an error. Note that this is the only non-optional argument.
  !
  ! The following arguments give the data to be read from the file. Note that all of 
  ! these arguments are optional.
  ! Int1, Real1, ... String1    - Argument to come from the first column.
  ! Int2, Real2, ... String2    - Argument to come from second column.
  ! .
  ! .
  ! .
  ! Int10, Real10, ... String10 - Argument to come from the tenth column.
  !
  ! SectionNumber - data will be read from the section number supplied if any. If no section
  !                 number is supplied, data is read sequentially from the current section.
  !                 Note: If the section number supplied is equal to the current section 
  !                 number, data continues to be read sequentially from this section.
  !
  ! Headers - Optional: Array of headers to be read from above the data. Note that
  !           this argument can be given without any data, and the headers will be read
  !           alone. 
  !
  ! ReadDataFromHeader - When set to .TRUE., lines beginning with #HD# will be interpretted
  !                     as data. This is usefull for files which contain
  !                     various information, including scalar information (which would be
  !                     written in the header) as well a set of arrays written in columns
  !                     which the user would like to be able to plot.
  !
  ! CommentCharacters - String of characters, any of which, when found at the beginning of a line,
  !                     will cause the line to be interpreted as a header. The default is '#!cC*'
  !
  ! NewSection - If a new section is found and this argument is present, data will not be read,
  !              and the routine will return NewSection = .TRUE.
  !              This helps when reading an unknown number of data lines from a section. 
  SUBROUTINE ReadData(FileName,                                      &
       & Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, Int9, Int10, &
       & Real1, Real2, Real3, Real4, Real5, Real6, Real7, Real8,      &
       & Real9, Real10,                                               &
       & Double1, Double2, Double3, Double4, Double5, Double6,        &
       & Double7, Double8, Double9, Double10,                         &
       & Complex1, Complex2, Complex3, Complex4, Complex5, Complex6,  &
       & Complex7, Complex8, Complex9, Complex10,                     &
       & DComplex1, DComplex2, DComplex3, DComplex4, DComplex5,       &
       & DComplex6, DComplex7, DComplex8, DComplex9, DComplex10,      &
       & String1, String2, String3, String4, String5, String6,        &
       & String7, String8, String9, String10,                         &
       & SectionNumber, Headers, ReadDataFromHeader, CommentCharacters, &
       & NewSection, FileType, ExpectNewSection, ErrorMessage)
    CHARACTER*(*), INTENT(IN) :: FileName
    CHARACTER*(*), INTENT(IN), OPTIONAL :: ErrorMessage(:)
    INTEGER, INTENT(OUT), OPTIONAL :: Int1, Int2, Int3, Int4, Int5,    &
         & Int6, Int7, Int8, Int9, Int10
    REAL, INTENT(OUT), OPTIONAL :: Real1, Real2, Real3, Real4, Real5,  &
         & Real6, Real7, Real8, Real9, Real10
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: Double1, Double2,       &
         & Double3, Double4, Double5, Double6, Double7, Double8,      &
         & Double9, Double10
    COMPLEX, INTENT(OUT), OPTIONAL :: Complex1, Complex2, Complex3,    &
         & Complex4, Complex5,                                        &
         & Complex6, Complex7, Complex8, Complex9, Complex10
    COMPLEX*16, INTENT(OUT), OPTIONAL :: DComplex1, DComplex2,         &
         & DComplex3, DComplex4, DComplex5, DComplex6, DComplex7,     &
         & DComplex8, DComplex9, DComplex10
    CHARACTER*(*), INTENT(OUT), OPTIONAL :: String1, String2, String3, &
         & String4, String5, String6, String7, String8, String9,      &
         & String10

    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: CommentCharacters, FileType
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader, ExpectNewSection
    LOGICAL,INTENT(OUT),OPTIONAL :: NewSection
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber

    ! Local variables
    ! NumArgs - Number of data arguments passed
    ! iArg    - The current argument.
    ! iUnit   - unit number associated with FileName
    INTEGER NumArgs, iArg, iUnit
    ! IsPresent  - If IsPresent(i) is true then the ith argument has already been defined.
    ! HeaderData - If true, interpret lines beginning with #HD# as data.
    ! txt, pad  - These are here in anticipation of enabling various formats,
    !              including text, paded ascii to start, and possibly further formats 
    !              such as xml.
    ! NewSection - If a new section is found, return set to true and return without reading.
    LOGICAL IsPresent(10), HeaderData, txt, pad, ExpNewSect
    ! DataArray - Array of strings to hold the data. Data line is split into words contained 
    !             in DataArray, then the arguments are read from DataArray(iArg)
    CHARACTER(MaxStrLen) DataArray(10)
    ! DataTypeLine - The data type line lists the types of data that have been specified
    !                in the file on a line beginning with #DT#. The arguments specified by
    !                the user will be checked against the types defined in the file.
    CHARACTER(120) DataTypeLine
    ! DataTypes - Array of types. 'Int', 'Real', 'Double' ...
    CHARACTER(20) DataTypes(10)
    ! ! FileAction - String variable that holds the action associated with the IOFile
    !                specified by FileName.
    CHARACTER(80) Args(20), FileAction
    ! String of comment characters.
    CHARACTER(10) CmtChars, FlType
    ! String holding the current line.
    CHARACTER(1000) Line
    ! NSections  - number of current section.
    ! NDataTypes - number of data types specified on the data type line in the file.
    ! NFields - number of columns in the current line.
    ! iField - current column.
    INTEGER i1, NSections, NDataTypes, NFields, iField

    ! Initialization
    HeaderData   = .FALSE.
    IsPresent(:) = .FALSE.
    NumArgs      = 0
    txt          = .TRUE.
    pad         = .FALSE.
    DataTypeLine = ''
    DataTypes(:) = 'Unknown'
    NDataTypes = 0
    iField = 0
    FlType = 'TXT'
    ExpNewSect = .FALSE.

    IF(PRESENT(ExpectNewSection)) ExpNewSect = ExpectNewSection
    ! Set the type of file. Default is txt.
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          txt = .TRUE.
          pad = .FALSE.
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          pad = .TRUE.
          txt = .FALSE.
       ELSE
          CALL Error('Error: Illegal file type passed to WriteData')
       END IF
    END IF

    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = '#!cC*'
    END IF

    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    ! Open the file.
    IF(PRESENT(ErrorMessage)) THEN
       CALL OpenFl(FileName, FileStatus = 'OLD', FileAction = 'READ', ErrorMessage = ErrorMessage)
    ELSE
       CALL OpenFl(FileName, FileStatus = 'OLD', FileAction = 'READ')
    END IF

    ! Get the unit number, the section number, and the fileaction (Read or Write).
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit,NSections = NSections, FileAction = FileAction)    
    ! If this is the first call, set NSections to 1
    IF(NSections.eq.0) NSections = 1

    ! Check that the FileAction is 'READ' not 'WRITE'.
    IF(FileAction(1:5).eq.'WRITE') CALL Error('Error in ReadData. File ' // FileName // &
         & ' is already open for writing.')

    ! If section number has been supplied, and is not equal to the current section number
    ! read to the specified section. 
    IF(PRESENT(SectionNumber)) THEN
       IF(SectionNumber.ne.NSections) CALL ReadToSection(FileName,SectionNumber)
       NSections = SectionNumber
       ! If we have reached the end of the file, return.
       IF(EndOfFile(FileName)) RETURN
    END IF

    ! Read the headers, keeping them in Headers array if supplied. This will also update the
    ! DataTypeLine and section number if they are found.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName,Line,CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName,Line,CommentCharacters = CmtChars)
    END IF

    ! If end of file or end of section, return.
    IF(EndOfFile(FileName).or.ReadingNewSection(FileName,NSections)) THEN
       IF(.NOT.ExpNewSect) THEN
          IF(PRESENT(NewSection)) NewSection = .TRUE.
          RETURN
       END IF
    END IF

    ! Get the DataTypeLine
    CALL GetIOFileInfo(FileName, DataTypeLine = DataTypeLine, NSections = NSections)
    ! If the DataTypeLine is not defined in the file, ignore.
    IF(TRIM(ADJUSTL(DataTypeLine)).ne.'Unknown') THEN
       NDataTypes = 10
       CALL bwords(DataTypeLine,NDataTypes,DataTypes)
    END IF

    ! Read the Line.
    Line = TRIM(ADJUSTL(Line))

    ! If Line starts with #HD#, this is header data. Ignore #HD# and 
    ! break the rest of the line.
    IF(Line(1:4).eq.'#HD#') Line = Line(5:LEN(Line))
    NFields = 20
    CALL bwords2(Line,NFields, Args)


    ! Find the 1st argument.
    IF(PRESENT(Int1)) CALL GetArgInt(Int1, 1)
    IF(PRESENT(Real1)) CALL GetArgReal(Real1, 1)
    IF(PRESENT(Double1)) CALL GetArgDouble(Double1, 1)
    IF(PRESENT(Complex1)) CALL GetArgComplex(Complex1, 1)
    IF(PRESENT(DComplex1)) CALL GetArgDComplex(DComplex1, 1)
    IF(PRESENT(String1)) CALL GetArgString(String1, 1)

    ! The GetArgType(Arg,ColumnNumber) subroutines checks the presence of Arg. If Arg
    ! is present, the data is read into DataArray(ColumnNumber). Also
    ! iField is incremented.
    ! In addition, various errors are checked.
    ! Find the 2nd argument.
    IF(PRESENT(Int2)) CALL GetArgInt(Int2, 2)
    IF(PRESENT(Real2)) CALL GetArgReal(Real2, 2)
    IF(PRESENT(Double2)) CALL GetArgDouble(Double2, 2)
    IF(PRESENT(Complex2)) CALL GetArgComplex(Complex2, 2)
    IF(PRESENT(DComplex2)) CALL GetArgDComplex(DComplex2, 2)
    IF(PRESENT(String2)) CALL GetArgString(String2, 2)

    ! Find the 3rd argument.
    IF(PRESENT(Int3)) CALL GetArgInt(Int3, 3)
    IF(PRESENT(Real3)) CALL GetArgReal(Real3, 3)
    IF(PRESENT(Double3)) CALL GetArgDouble(Double3, 3)
    IF(PRESENT(Complex3)) CALL GetArgComplex(Complex3, 3)
    IF(PRESENT(DComplex3)) CALL GetArgDComplex(DComplex3, 3)
    IF(PRESENT(String3)) CALL GetArgString(String3, 3)

    ! Find the 4th argument.
    IF(PRESENT(Int4)) CALL GetArgInt(Int4, 4)
    IF(PRESENT(Real4)) CALL GetArgReal(Real4, 4)
    IF(PRESENT(Double4)) CALL GetArgDouble(Double4, 4)
    IF(PRESENT(Complex4)) CALL GetArgComplex(Complex4, 4)
    IF(PRESENT(DComplex4)) CALL GetArgDComplex(DComplex4, 4)
    IF(PRESENT(String4)) CALL GetArgString(String4, 4)

    ! Find the 5th argument.
    IF(PRESENT(Int5)) CALL GetArgInt(Int5, 5)
    IF(PRESENT(Real5)) CALL GetArgReal(Real5, 5)
    IF(PRESENT(Double5)) CALL GetArgDouble(Double5, 5)
    IF(PRESENT(Complex5)) CALL GetArgComplex(Complex5, 5)
    IF(PRESENT(DComplex5)) CALL GetArgDComplex(DComplex5, 5)
    IF(PRESENT(String5)) CALL GetArgString(String5, 5)

    ! Find the 6th argument.
    IF(PRESENT(Int6)) CALL GetArgInt(Int6, 6)
    IF(PRESENT(Real6)) CALL GetArgReal(Real6, 6)
    IF(PRESENT(Double6)) CALL GetArgDouble(Double6, 6)
    IF(PRESENT(Complex6)) CALL GetArgComplex(Complex6, 6)
    IF(PRESENT(DComplex6)) CALL GetArgDComplex(DComplex6, 6)
    IF(PRESENT(String6)) CALL GetArgString(String6, 6)

    ! Find the 7th argument.
    IF(PRESENT(Int7)) CALL GetArgInt(Int7, 7)
    IF(PRESENT(Real7)) CALL GetArgReal(Real7, 7)
    IF(PRESENT(Double7)) CALL GetArgDouble(Double7, 7)
    IF(PRESENT(Complex7)) CALL GetArgComplex(Complex7, 7)
    IF(PRESENT(DComplex7)) CALL GetArgDComplex(DComplex7, 7)
    IF(PRESENT(String7)) CALL GetArgString(String7, 7)

    ! Find the 8th argument.
    IF(PRESENT(Int8)) CALL GetArgInt(Int8, 8)
    IF(PRESENT(Real8)) CALL GetArgReal(Real8, 8)
    IF(PRESENT(Double8)) CALL GetArgDouble(Double8, 8)
    IF(PRESENT(Complex8)) CALL GetArgComplex(Complex8, 8)
    IF(PRESENT(DComplex8)) CALL GetArgDComplex(DComplex8, 8)
    IF(PRESENT(String8)) CALL GetArgString(String8, 8)

    ! Find the 9th argument.
    IF(PRESENT(Int9)) CALL GetArgInt(Int9, 9)
    IF(PRESENT(Real9)) CALL GetArgReal(Real9, 9)
    IF(PRESENT(Double9)) CALL GetArgDouble(Double9, 9)
    IF(PRESENT(Complex9)) CALL GetArgComplex(Complex9, 9)
    IF(PRESENT(DComplex9)) CALL GetArgDComplex(DComplex9, 9)
    IF(PRESENT(String9)) CALL GetArgString(String9, 9)

    ! Find the 10th argument.
    IF(PRESENT(Int10)) CALL GetArgInt(Int10, 10)
    IF(PRESENT(Real10)) CALL GetArgReal(Real10, 10)
    IF(PRESENT(Double10)) CALL GetArgDouble(Double10, 10)
    IF(PRESENT(Complex10)) CALL GetArgComplex(Complex10, 10)
    IF(PRESENT(DComplex10)) CALL GetArgDComplex(DComplex10, 10)
    IF(PRESENT(String10)) CALL GetArgString(String10, 10)

  CONTAINS

    ! The GetArgType(Arg,ColumnNumber) subroutines checks the presence of Arg. If Arg
    ! is present, the data is read into DataArray(ColumnNumber). Also
    ! iField is incremented.
    ! In addition, various errors are checked.    
    SUBROUTINE GetArgInt(Arg, iArg)
      INTEGER,INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Int')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iField),*) Arg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField),Arg)
      END IF
    END SUBROUTINE GetArgInt

    SUBROUTINE GetArgReal( Arg, iArg)
      REAL,INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Real')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iArg),*) Arg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadrDefault, Arg)
      END IF
    END SUBROUTINE GetArgReal

    SUBROUTINE GetArgDouble( Arg, iArg)
      DOUBLE PRECISION,INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Double')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iArg),*) Arg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadxDefault, Arg)
      END IF
    END SUBROUTINE GetArgDouble

    SUBROUTINE GetArgComplex( Arg, iArg)
      COMPLEX,INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg
      REAL ReArg, ImArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'Complex')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iField),*) ReArg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadrDefault, ReArg)
      END IF

      ! Check that the number of Fields in the data file is >= to iField. This was checked in
      ! CheckForErrors, but we need to check it again since we are reading again.
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iField),*) ImArg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadrDefault, ImArg)
      END IF
      Arg = ReArg + (0.0, 1.0)*ImArg
    END SUBROUTINE GetArgComplex

    SUBROUTINE GetArgDComplex( Arg, iArg)
      COMPLEX*16,INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg

      DOUBLE PRECISION ReArg, ImArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'DComplex')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iField),*) ReArg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadxDefault, ReArg)
      END IF

      ! Check that the number of Fields in the data file is >= to iField. This was checked in
      ! CheckForErrors, but we need to check it again since we are reading again.
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)
      iField = iField + 1
      IF(txt) THEN
         READ(Args(iField),*) ImArg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField), npadxDefault, ImArg)
      END IF
      Arg = ReArg + (0.d0, 1.d0)*ImArg
    END SUBROUTINE GetArgDComplex

    SUBROUTINE GetArgString( Arg, iArg)
      CHARACTER*(*),INTENT(OUT) :: Arg
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,'String')

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
      READ(Args(iField),*) Arg
      ELSEIF(pad) THEN
         CALL ReadPAD(Args(iField),Arg)
      END IF
    END SUBROUTINE GetArgString

    ! CheckForErrors checks for the following errors and stops if it finds one:
    ! 
    ! 1. 'Error in ReadData multiple arguments in place' - This error occurs 
    !    when the user has specified that two arguments should be written into
    !    one column, i.e. WriteData(FileName, Int1 = I, Real1 = R, ...)
    !
    ! 2. 'Too many arguments passed to ReadData. Max is 10.' - user has specified
    !    more than 10 aguments. This one is redundant, and error 1. should always 
    !    occur if this happens.
    !
    ! 3. 'Too many arguments passed to ReadData when compared to
    !     number of fields in file FileName' - user has specified
    !     more (or different) arguments than are available in file.
    !
    ! 4. 'Error: Number of arguments passed to ReadData does not
    !     match that defined in the file FileName' - user has
    !     specified more arguments than are defined in the data type
    !     line in the file.  
    !
    ! 5.  'Error: Arguments passed to ReadData are of different type
    !      than those defined in the file FileName.' - user is trying
    !      to read different types than those defined by the data
    !      type line in the file.
    !
    ! In addition, CheckForError increments the number of arguments.
    SUBROUTINE CheckForErrors(iArg,DTL)
      INTEGER,INTENT(IN) :: iArg
      CHARACTER*(*) DTL
      CHARACTER(2) ArgNum
      INTEGER iAllocateError
      ArgNum=''
      WRITE(ArgNum,'(I2)') iArg

      IF(IsPresent(iArg)) THEN
         CALL Error('Error in ReadData multiple arguments in place ' // TRIM(ArgNum) // '.')
      END IF

      ! Set IsPresent(iArg)
      IsPresent(iArg) = .TRUE.

      ! Check that number of arguments is not too large. Redundant, but oh well.
      IF(NumArgs.ge.10) CALL Error('Too many arguments passed to ReadData. Max is 10.')
      NumArgs = NumArgs + 1

      ! Check that the number of Fields in the data file is >= to iField
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)

      ! If data types have been defined in this file, check for data type errors.
      IF(NDataTypes.gt.0) THEN
         ! Check that there are not more arguments than data types.
         IF(iArg.gt.NDataTypes) THEN
            CALL Error('Error: Number of arguments passed to ReadData does not match ' // &
                 & 'that defined in the file ' // FileName)
         END IF

         ! Check that this argument is compatible with the datatype described in DataTypeLine.
         IF( TRIM(ADJUSTL(DataTypes(iArg))).ne.TRIM(ADJUSTL(DTL)) ) THEN
            CALL Error('Error: Arguments passed to ReadData are of different type than those' // &
                 & ' defined in the file ' // FileName // '.')
         END IF
      END IF
    END SUBROUTINE CheckForErrors

  END SUBROUTINE ReadData


  ! SUBROUTINE ReadArrayData reads up to 10 arrays out of the file
  ! FileName. This routine is almost identical to ReadData except
  ! that the data arguments are arbitrary length arrays. See
  ! description of ReadData for more information. The following
  ! information gives differences to ReadData.
  ! 
  ! ReadArrayData will read multiple lines into the given arrays
  ! until it reaches the end of the section, or the end of the file.
  ! Each column of data in the file is interpreted as a single array.
  ! All arrays passed to the routine must be the same length. If the
  ! arrays are not long enough to hold the data specified in the
  ! section to be read, an error will occur.
  !
  ! It has one additional optional argument when compared to ReadData.
  ! NumElements - ReadArrayData will return the number of elements it
  ! read from the specified file.
  SUBROUTINE ReadArrayData(FileName,                                      &
       & Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, Int9, Int10, &
       & Real1, Real2, Real3, Real4, Real5, Real6, Real7, Real8,      &
       & Real9, Real10,                                               &
       & Double1, Double2, Double3, Double4, Double5, Double6,        &
       & Double7, Double8, Double9, Double10,                         &
       & Complex1, Complex2, Complex3, Complex4, Complex5, Complex6,  &
       & Complex7, Complex8, Complex9, Complex10,                     &
       & DComplex1, DComplex2, DComplex3, DComplex4, DComplex5,       &
       & DComplex6, DComplex7, DComplex8, DComplex9, DComplex10,      &
       & String1, String2, String3, String4, String5, String6,        &
       & String7, String8, String9, String10,                         &
       & SectionNumber, Headers, ReadDataFromHeader, CommentCharacters, &
       & NewSection, NumElements, FileType, ExactLength)
    CHARACTER*(*), INTENT(IN) :: FileName
    INTEGER, INTENT(OUT), OPTIONAL :: Int1(:), Int2(:), Int3(:), Int4(:), Int5(:),    &
         & Int6(:), Int7(:), Int8(:), Int9(:), Int10(:)
    REAL, INTENT(OUT), OPTIONAL :: Real1(:), Real2(:), Real3(:), Real4(:), Real5(:),  &
         & Real6(:), Real7(:), Real8(:), Real9(:), Real10(:)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: Double1(:), Double2(:),       &
         & Double3(:), Double4(:), Double5(:), Double6(:), Double7(:), Double8(:),      &
         & Double9(:), Double10(:)
    COMPLEX, INTENT(OUT), OPTIONAL :: Complex1(:), Complex2(:), Complex3(:),    &
         & Complex4(:), Complex5(:),                                        &
         & Complex6(:), Complex7(:), Complex8(:), Complex9(:), Complex10(:)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: DComplex1(:), DComplex2(:),         &
         & DComplex3(:), DComplex4(:), DComplex5(:), DComplex6(:), DComplex7(:),     &
         & DComplex8(:), DComplex9(:), DComplex10(:)
    CHARACTER*(*), INTENT(OUT), OPTIONAL :: String1(:), String2(:), String3(:), &
         & String4(:), String5(:), String6(:), String7(:), String8(:), String9(:),      &
         & String10(:)

    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: CommentCharacters, FileType
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader, ExactLength
    LOGICAL,INTENT(OUT),OPTIONAL :: NewSection
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    INTEGER,INTENT(OUT),OPTIONAL :: NumElements

    INTEGER NumArgs, iArg, iUnit, NumData, NColumnLabels
    LOGICAL IsPresent(10), HeaderData, txt, pad, ExactLen
    CHARACTER(MaxStrLen),ALLOCATABLE :: DataArray(:,:)
    CHARACTER(120) DataTypeLine
    CHARACTER(20) DataTypes(10)
    CHARACTER(20) IntFormat, RealFormat, DoubleFormat, ComplexFormat, &
         & DComplexFormat, StringFormat, FileAction
    CHARACTER(80) Args(20)
    CHARACTER(10) CLFormats(30), CmtChars
    CHARACTER(1000) Line
    CHARACTER(10) TmpStr, FlType

    INTEGER i1, NSections, NDataTypes, NFields, iField, MaxData, iAllocateError

    ! Initialization
    TmpStr = ' '
    HeaderData   = .FALSE.
    IsPresent(:) = .FALSE.
    NumArgs      = 0
    txt          = .TRUE.
    pad          = .FALSE.
    DataTypeLine = ''
    DataTypes(:) = 'Unknown'
    CLFormats(:) = ''
    NColumnLabels = 0
    NDataTypes = 0
    iField = 0
    MaxData = 0
    FlType = 'TXT'
    ExactLen = .TRUE.

    IF(PRESENT(ExactLength)) ExactLen = ExactLength

    ! Set the type of file. Default is txt.
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          txt = .TRUE.
          pad = .FALSE.
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          pad = .TRUE.
          txt = .FALSE.
       ELSE
          CALL Error('Error: Illegal file type passed to WriteData')
       END IF
    END IF

    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    ! Open the file.
    CALL OpenFl(FileName, FileStatus = 'OLD', FileAction = 'READ')

    ! Get the unit number, and the section number.
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit,NSections = NSections, FileAction = FileAction)
    IF(NSections.eq.0) NSections = 1

    ! Check that the FileAction is read.
    IF(FileAction(1:5).eq.'WRITE') CALL Error('Error in ReadData. File ' // FileName // &
         & ' is already open for writing.')

    ! If section number has been supplied, read to that section. 
    IF(PRESENT(SectionNumber)) THEN
       IF(SectionNumber.ne.NSections) CALL ReadToSection(FileName,SectionNumber)
       NSections = SectionNumber
       IF(EndOfFile(FileName)) RETURN
    END IF

    ! Read the headers, keeping them in Headers if supplied. This will also update the
    ! DataTypeLine and section number if they are found.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName,Line, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName,Line,CommentCharacters = CmtChars)
    END IF
    IF(EndOfFile(FileName)) RETURN

    ! Get the DataTypeLine
    CALL GetIOFileInfo(FileName, DataTypeLine = DataTypeLine, NSections = NSections)
    IF(TRIM(ADJUSTL(DataTypeLine)).ne.'Unknown') THEN
       NDataTypes = 10
       CALL bwords(DataTypeLine,NDataTypes,DataTypes)
    END IF

    ! Get the maximum number of data lines.
    CALL GetMaxData

    ! Allocate space for string data array.
    ALLOCATE(DataArray(MaxData,20), STAT = iAllocateError)
    CALL CheckAllocation(iAllocateError, 'Error: Failed to allocate space for DataArray' //&
         & ' in subroutine WriteArrayData.')
    ! Read all the lines in this section, and split them into words.
    ! ReadHeaders reads the line.

    DO i1 = 1, MaxData 
       NFields = 20
       ! Clean up the line.
       Line = ADJUSTL(Line)
       ! If this is header data, remove the #HD# marker.
       IF(Line(1:4).eq.'#HD#') Line = Line(5:LEN(Line))
       ! Split line into words and copy to DataArray.
       CALL bwords2(Line,NFields,Args)
       DO iField = 1, NFields
          DataArray(i1,iField) = Args(iField)
       END DO

       ! If we have reached the end of the section or the end of the 
       ! file, return.
       IF(EndOfFile(FileName).or.ReadingNewSection(FileName,NSections)) THEN
          IF(PRESENT(NewSection)) NewSection = .TRUE.
          EXIT
       ELSE
          NumData = i1
       END IF
       IF(i1.eq.MaxData) EXIT
       ! Read past any comments and output next line.
       CALL ReadHeaders(FileName,Line,CommentCharacters = CmtChars)
    END DO
    ! Check that the next line is not data. If it is. The arrays that were passed
    ! are not long enough to hold all of the data in the file. Error.
    IF(.not.ExactLen) THEN
       IF(NextLineIsData(FileName)) THEN
          WRITE(TmpStr,'(I4)') NSections
          CALL Error('Error: Arrays passed to ReadArrayData are not long enough to hold', StopProgram = .FALSE.)
          CALL Error('the data defined in section ' // TRIM(ADJUSTL(TmpStr)) // ' of ' // TRIM(FileName), StopProgram = .FALSE.)
       END IF
    END IF

    iField = 0
    ! Find the 1st argument and read all elements into DataArray
    IF(PRESENT(Int1)) CALL GetArrayArgInt(Int1, 1)
    IF(PRESENT(Real1)) CALL GetArrayArgReal(Real1, 1)
    IF(PRESENT(Double1)) CALL GetArrayArgDouble(Double1, 1)
    IF(PRESENT(Complex1)) CALL GetArrayArgComplex(Complex1, 1)
    IF(PRESENT(DComplex1)) CALL GetArrayArgDComplex(DComplex1, 1)
    IF(PRESENT(String1)) CALL GetArrayArgString(String1, 1)

    ! Find the 2nd argument.
    IF(PRESENT(Int2)) CALL GetArrayArgInt(Int2, 2)
    IF(PRESENT(Real2)) CALL GetArrayArgReal(Real2, 2)
    IF(PRESENT(Double2)) CALL GetArrayArgDouble(Double2, 2)
    IF(PRESENT(Complex2)) CALL GetArrayArgComplex(Complex2, 2)
    IF(PRESENT(DComplex2)) CALL GetArrayArgDComplex(DComplex2, 2)
    IF(PRESENT(String2)) CALL GetArrayArgString(String2, 2)

    ! Find the 3rd argument.
    IF(PRESENT(Int3)) CALL GetArrayArgInt(Int3, 3)
    IF(PRESENT(Real3)) CALL GetArrayArgReal(Real3, 3)
    IF(PRESENT(Double3)) CALL GetArrayArgDouble(Double3, 3)
    IF(PRESENT(Complex3)) CALL GetArrayArgComplex(Complex3, 3)
    IF(PRESENT(DComplex3)) CALL GetArrayArgDComplex(DComplex3, 3)
    IF(PRESENT(String3)) CALL GetArrayArgString(String3, 3)

    ! Find the 4th argument.
    IF(PRESENT(Int4)) CALL GetArrayArgInt(Int4, 4)
    IF(PRESENT(Real4)) CALL GetArrayArgReal(Real4, 4)
    IF(PRESENT(Double4)) CALL GetArrayArgDouble(Double4, 4)
    IF(PRESENT(Complex4)) CALL GetArrayArgComplex(Complex4, 4)
    IF(PRESENT(DComplex4)) CALL GetArrayArgDComplex(DComplex4, 4)
    IF(PRESENT(String4)) CALL GetArrayArgString(String4, 4)

    ! Find the 5th argument.
    IF(PRESENT(Int5)) CALL GetArrayArgInt(Int5, 5)
    IF(PRESENT(Real5)) CALL GetArrayArgReal(Real5, 5)
    IF(PRESENT(Double5)) CALL GetArrayArgDouble(Double5, 5)
    IF(PRESENT(Complex5)) CALL GetArrayArgComplex(Complex5, 5)
    IF(PRESENT(DComplex5)) CALL GetArrayArgDComplex(DComplex5, 5)
    IF(PRESENT(String5)) CALL GetArrayArgString(String5, 5)

    ! Find the 6th argument.
    IF(PRESENT(Int6)) CALL GetArrayArgInt(Int6, 6)
    IF(PRESENT(Real6)) CALL GetArrayArgReal(Real6, 6)
    IF(PRESENT(Double6)) CALL GetArrayArgDouble(Double6, 6)
    IF(PRESENT(Complex6)) CALL GetArrayArgComplex(Complex6, 6)
    IF(PRESENT(DComplex6)) CALL GetArrayArgDComplex(DComplex6, 6)
    IF(PRESENT(String6)) CALL GetArrayArgString(String6, 6)

    ! Find the 7th argument.
    IF(PRESENT(Int7)) CALL GetArrayArgInt(Int7, 7)
    IF(PRESENT(Real7)) CALL GetArrayArgReal(Real7, 7)
    IF(PRESENT(Double7)) CALL GetArrayArgDouble(Double7, 7)
    IF(PRESENT(Complex7)) CALL GetArrayArgComplex(Complex7, 7)
    IF(PRESENT(DComplex7)) CALL GetArrayArgDComplex(DComplex7, 7)
    IF(PRESENT(String7)) CALL GetArrayArgString(String7, 7)

    ! Find the 8th argument.
    IF(PRESENT(Int8)) CALL GetArrayArgInt(Int8, 8)
    IF(PRESENT(Real8)) CALL GetArrayArgReal(Real8, 8)
    IF(PRESENT(Double8)) CALL GetArrayArgDouble(Double8, 8)
    IF(PRESENT(Complex8)) CALL GetArrayArgComplex(Complex8, 8)
    IF(PRESENT(DComplex8)) CALL GetArrayArgDComplex(DComplex8, 8)
    IF(PRESENT(String8)) CALL GetArrayArgString(String8, 8)

    ! Find the 9th argument.
    IF(PRESENT(Int9)) CALL GetArrayArgInt(Int9, 9)
    IF(PRESENT(Real9)) CALL GetArrayArgReal(Real9, 9)
    IF(PRESENT(Double9)) CALL GetArrayArgDouble(Double9, 9)
    IF(PRESENT(Complex9)) CALL GetArrayArgComplex(Complex9, 9)
    IF(PRESENT(DComplex9)) CALL GetArrayArgDComplex(DComplex9, 9)
    IF(PRESENT(String9)) CALL GetArrayArgString(String9, 9)

    ! Find the 10th argument.
    IF(PRESENT(Int10)) CALL GetArrayArgInt(Int10, 10)
    IF(PRESENT(Real10)) CALL GetArrayArgReal(Real10, 10)
    IF(PRESENT(Double10)) CALL GetArrayArgDouble(Double10, 10)
    IF(PRESENT(Complex10)) CALL GetArrayArgComplex(Complex10, 10)
    IF(PRESENT(DComplex10)) CALL GetArrayArgDComplex(DComplex10, 10)
    IF(PRESENT(String10)) CALL GetArrayArgString(String10, 10)

    IF(PRESENT(NumElements)) NumElements = NumData
    DEALLOCATE(DataArray)
  CONTAINS

    SUBROUTINE GetArrayArgInt(Arg, iArg)
      INTEGER,INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      INTEGER iData

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'Int')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in ReadArrayData: All arrays must be same size.')
      END IF
      ! Read the next column.
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField),*) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1,iField),Arg(i1))
         END DO
      END IF
    END SUBROUTINE GetArrayArgInt

    SUBROUTINE GetArrayArgReal( Arg, iArg)
      REAL,INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'Real')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in writearray: All arrays must be same size.')
      END IF

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField),*) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1, iField), npadrDefault, Arg(i1))
         END DO
      END IF
    END SUBROUTINE GetArrayArgReal

    SUBROUTINE GetArrayArgDouble( Arg, iArg)
      DOUBLE PRECISION,INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'Double')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in writearray: All arrays must be same size.')
      END IF

      ! Read the next argument
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField),*) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1, iField), npadrDefault, Arg(i1))
         END DO
      END IF
    END SUBROUTINE GetArrayArgDouble

    SUBROUTINE GetArrayArgComplex( Arg, iArg)
      COMPLEX,INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      REAL ReArg, ImArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'Complex')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in writearray: All arrays must be same size.')
      END IF

      ! Read the next argument
      iField = iField + 1
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField-1),*) ReArg
            READ(DataArray(i1,iField),*) ImArg
            Arg(i1) = ReArg + (0.0,1.0)*ImArg
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1, iField-1), npadrDefault, ReArg)
            CALL ReadPAD(DataArray(i1, iField), npadrDefault, ImArg)
            Arg(i1) = ReArg + (0.0,1.0)*ImArg
         END DO
      END IF
    END SUBROUTINE GetArrayArgComplex

    SUBROUTINE GetArrayArgDComplex( Arg, iArg)
      COMPLEX*16,INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg

      DOUBLE PRECISION ReArg, ImArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'DComplex')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in writearray: All arrays must be same size.')
      END IF

      ! Read the next column
      iField = iField + 1
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField-1),*) ReArg
            READ(DataArray(i1,iField),*) ImArg
            Arg(i1) = ReArg + (0.0,1.0)*ImArg
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1, iField-1), npadxDefault, ReArg)
            CALL ReadPAD(DataArray(i1, iField), npadxDefault, ImArg)
            Arg(i1) = ReArg + (0.d0,1.d0)*ImArg
         END DO
      END IF
    END SUBROUTINE GetArrayArgDComplex

    SUBROUTINE GetArrayArgString( Arg, iArg)
      CHARACTER*(*),INTENT(OUT) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckArrayForErrors(iArg,'String')

      IF(MaxData.ne.SIZE(Arg)) THEN
         CALL Error('Error in writearray: All arrays must be same size.')
      END IF

      ! Read the next column
      iField = iField + 1
      IF(txt) THEN
         DO i1 = 1, NumData
            READ(DataArray(i1,iField),*) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL ReadPAD(DataArray(i1, iField),Arg(i1))
         END DO
      END IF
    END SUBROUTINE GetArrayArgString
      
    SUBROUTINE CheckArrayForErrors(iArg,DTL)
      INTEGER,INTENT(IN) :: iArg
      CHARACTER*(*) DTL
      CHARACTER(2) ArgNum
      INTEGER iAllocateError
      ArgNum=''
      WRITE(ArgNum,'(I2)') iArg

      IF(IsPresent(iArg)) THEN
         CALL Error('Error in ReadData: multiple arguments in place ' // TRIM(ArgNum) // '.')
      END IF

      ! Set IsPresent(iArg)
      IsPresent(iArg) = .TRUE.

      ! Check that number of arguments is not too large. Redundant, but oh well.
      IF(NumArgs.ge.10) CALL Error('Too many arguments passed to ReadData. Max is 10.')
      NumArgs = NumArgs + 1
      ! Check that the number of Fields in the data file is >= to iField
      IF(iField.ge.NFields) CALL Error('Too many arguments passed to ReadData when compared ' //&
           & 'number of fields in file ' // FileName)

      ! If data types have been defined in this file, check for data type errors.
      IF(DataTypes(1).ne.'Unknown') THEN
         ! Check that there are not more arguments than data types.
         IF(iArg.gt.NDataTypes) THEN
            CALL Error('Error: Number of arguments passed to ReadData does not match ' // &
                 & 'that defined in the file ' // FileName)
         END IF

         ! Check that this argument is compatible with the datatype described in DataTypeLine.
         IF( TRIM(ADJUSTL(DataTypes(iArg))).ne.TRIM(ADJUSTL(DTL)) ) THEN
            CALL Error('Error: Arguments passed to ReadData are of different type than those' // &
                 & ' defined in the file ' // FileName // '.')
         END IF
      END IF
    END SUBROUTINE CheckArrayForErrors

    SUBROUTINE GetMaxData
      ! Find the 1st argument.
      IF(PRESENT(Int1)) MaxData = SIZE(Int1)
      IF(PRESENT(Real1)) MaxData = SIZE(Real1)
      IF(PRESENT(Double1)) MaxData = SIZE(Double1)
      IF(PRESENT(Complex1)) MaxData =SIZE(Complex1)
      IF(PRESENT(DComplex1)) MaxData =SIZE(DComplex1)
      IF(PRESENT(String1)) MaxData = SIZE(String1)
    END SUBROUTINE GetMaxData

  END SUBROUTINE ReadArrayData

  ! Subroutine WriteArrayData writes scalar data to the file filename.
  ! Up to 10 arrays of data can be written each specified by using
  ! the optional arguments Int1 - Int10, Real1 - Real10, etc.  
  ! WriteArrayData is very similar to WriteData, except that the data
  ! arguments are arbitrary length arrays instead of scalar valued
  ! variables. All data arrays passed to the routine must be the same
  ! length.
  SUBROUTINE WriteArrayData(FileName,                                      &
       & Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, Int9, Int10, &
       & Real1, Real2, Real3, Real4, Real5, Real6, Real7, Real8,      &
       & Real9, Real10,                                               &
       & Double1, Double2, Double3, Double4, Double5, Double6,        &
       & Double7, Double8, Double9, Double10,                         &
       & Complex1, Complex2, Complex3, Complex4, Complex5, Complex6,  &
       & Complex7, Complex8, Complex9, Complex10,                     &
       & DComplex1, DComplex2, DComplex3, DComplex4, DComplex5,       &
       & DComplex6, DComplex7, DComplex8, DComplex9, DComplex10,      &
       & String1, String2, String3, String4, String5, String6,        &
       & String7, String8, String9, String10,                         &
       & Headers, ColumnLabels, WriteDataInHeader, ForceNewSection,   &
       & FileType)
    CHARACTER*(*), INTENT(IN) :: FileName
    INTEGER, INTENT(IN), OPTIONAL :: Int1(:), Int2(:), Int3(:), Int4(:),  Int5(:),    &
         & Int6(:), Int7(:), Int8(:), Int9(:), Int10(:)
    REAL, INTENT(IN), OPTIONAL :: Real1(:), Real2(:), Real3(:), Real4(:), Real5(:),  &
         & Real6(:), Real7(:), Real8(:), Real9(:), Real10(:)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: Double1(:), Double2(:),       &
         & Double3(:), Double4(:), Double5(:), Double6(:), Double7(:), Double8(:),      &
         & Double9(:), Double10(:)
    COMPLEX, INTENT(IN), OPTIONAL :: Complex1(:), Complex2(:), Complex3(:),    &
         & Complex4(:), Complex5(:),                                        &
         & Complex6(:), Complex7(:), Complex8(:), Complex9(:), Complex10(:)
    COMPLEX*16, INTENT(IN), OPTIONAL :: DComplex1(:), DComplex2(:),         &
         & DComplex3(:), DComplex4(:), DComplex5(:), DComplex6(:), DComplex7(:),     &
         & DComplex8(:), DComplex9(:), DComplex10(:)
    CHARACTER*(*), INTENT(IN), OPTIONAL :: String1(:), String2(:), String3(:), &
         & String4(:), String5(:), String6(:), String7(:), String8(:), String9(:),      &
         & String10(:)

    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), ColumnLabels(:), FileType
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader, ForceNewSection

    INTEGER NumArgs, NumData, iUnit
    LOGICAL IsPresent(10), HeaderData, txt, pad, NewSection
    CHARACTER(80),ALLOCATABLE :: DataArray(:,:)
    CHARACTER(120) DataTypeLine
    CHARACTER(20) IntFormat, RealFormat, DoubleFormat, ComplexFormat, &
         & DComplexFormat, StringFormat, FileAction
    CHARACTER(80) ColumnLabelFormat
    CHARACTER(10) CLFormats(20)
    CHARACTER(4) TmpStr, FlType

    INTEGER i1, i2, NSections, NColumnLabels
    
    ! Initialization
    HeaderData   = .FALSE.
    IsPresent(:) = .FALSE.
    NumArgs      = 0
    txt          = .TRUE.
    pad         = .FALSE.
    DataTypeLine = ''
    CLFormats(:) = ''
    NColumnLabels = 0
    FlType = 'TXT'

    ! Set the type of file. Default is txt.
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          txt = .TRUE.
          pad = .FALSE.
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          pad = .TRUE.
          txt = .FALSE.
       ELSE
          CALL Error('Error: Illegal file type passed to WriteData')
       END IF
    END IF

    IF(PRESENT(ForceNewSection)) THEN
       NewSection = ForceNewSection
    ELSE
       NewSection = .FALSE.
    END IF

    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader

    ! Open the file.
    CALL OpenFl(FileName, FileStatus = 'REPLACE', FileAction = 'WRITE')

    ! Get the formats for each data type, and the FileAction.
    CALL GetIOFileInfo(FileName, IFormat = IntFormat, RFormat = RealFormat, &
         & DFormat = DoubleFormat, CFormat = ComplexFormat, &
         & DCFormat = DComplexFormat, SFormat = StringFormat, FileAction = FileAction)

    ! Check that this file was not already opened for READ.
    IF(TRIM(FileAction).eq.'READ') THEN 
       CALL Error('Error: file ' // FileName // &
            & ' was already opened for reading.', StopProgram = .FALSE.)
       CALL Error("Please use CALL CloseFl('" // FileName // "') before writing.")
    END IF

    ! Find the 1st argument.
    IF(PRESENT(Int1)) CALL SetArgInt(Int1, 1)
    IF(PRESENT(Real1)) CALL SetArgReal(Real1, 1)
    IF(PRESENT(Double1)) CALL SetArgDouble(Double1, 1)
    IF(PRESENT(Complex1)) CALL SetArgComplex(Complex1, 1)
    IF(PRESENT(DComplex1)) CALL SetArgDComplex(DComplex1, 1)
    IF(PRESENT(String1)) CALL SetArgString(String1, 1)

    ! Find the 2nd argument.
    IF(PRESENT(Int2)) CALL SetArgInt(Int2, 2)
    IF(PRESENT(Real2)) CALL SetArgReal(Real2, 2)
    IF(PRESENT(Double2)) CALL SetArgDouble(Double2, 2)
    IF(PRESENT(Complex2)) CALL SetArgComplex(Complex2, 2)
    IF(PRESENT(DComplex2)) CALL SetArgDComplex(DComplex2, 2)
    IF(PRESENT(String2)) CALL SetArgString(String2, 2)

    ! Find the 3rd argument.
    IF(PRESENT(Int3)) CALL SetArgInt(Int3, 3)
    IF(PRESENT(Real3)) CALL SetArgReal(Real3, 3)
    IF(PRESENT(Double3)) CALL SetArgDouble(Double3, 3)
    IF(PRESENT(Complex3)) CALL SetArgComplex(Complex3, 3)
    IF(PRESENT(DComplex3)) CALL SetArgDComplex(DComplex3, 3)
    IF(PRESENT(String3)) CALL SetArgString(String3, 3)

    ! Find the 4th argument.
    IF(PRESENT(Int4)) CALL SetArgInt(Int4, 4)
    IF(PRESENT(Real4)) CALL SetArgReal(Real4, 4)
    IF(PRESENT(Double4)) CALL SetArgDouble(Double4, 4)
    IF(PRESENT(Complex4)) CALL SetArgComplex(Complex4, 4)
    IF(PRESENT(DComplex4)) CALL SetArgDComplex(DComplex4, 4)
    IF(PRESENT(String4)) CALL SetArgString(String4, 4)

    ! Find the 5th argument.
    IF(PRESENT(Int5)) CALL SetArgInt(Int5, 5)
    IF(PRESENT(Real5)) CALL SetArgReal(Real5, 5)
    IF(PRESENT(Double5)) CALL SetArgDouble(Double5, 5)
    IF(PRESENT(Complex5)) CALL SetArgComplex(Complex5, 5)
    IF(PRESENT(DComplex5)) CALL SetArgDComplex(DComplex5, 5)
    IF(PRESENT(String5)) CALL SetArgString(String5, 5)

    ! Find the 6th argument.
    IF(PRESENT(Int6)) CALL SetArgInt(Int6, 6)
    IF(PRESENT(Real6)) CALL SetArgReal(Real6, 6)
    IF(PRESENT(Double6)) CALL SetArgDouble(Double6, 6)
    IF(PRESENT(Complex6)) CALL SetArgComplex(Complex6, 6)
    IF(PRESENT(DComplex6)) CALL SetArgDComplex(DComplex6, 6)
    IF(PRESENT(String6)) CALL SetArgString(String6, 6)

    ! Find the 7th argument.
    IF(PRESENT(Int7)) CALL SetArgInt(Int7, 7)
    IF(PRESENT(Real7)) CALL SetArgReal(Real7, 7)
    IF(PRESENT(Double7)) CALL SetArgDouble(Double7, 7)
    IF(PRESENT(Complex7)) CALL SetArgComplex(Complex7, 7)
    IF(PRESENT(DComplex7)) CALL SetArgDComplex(DComplex7, 7)
    IF(PRESENT(String7)) CALL SetArgString(String7, 7)

    ! Find the 8th argument.
    IF(PRESENT(Int8)) CALL SetArgInt(Int8, 8)
    IF(PRESENT(Real8)) CALL SetArgReal(Real8, 8)
    IF(PRESENT(Double8)) CALL SetArgDouble(Double8, 8)
    IF(PRESENT(Complex8)) CALL SetArgComplex(Complex8, 8)
    IF(PRESENT(DComplex8)) CALL SetArgDComplex(DComplex8, 8)
    IF(PRESENT(String8)) CALL SetArgString(String8, 8)

    ! Find the 9th argument.
    IF(PRESENT(Int9)) CALL SetArgInt(Int9, 9)
    IF(PRESENT(Real9)) CALL SetArgReal(Real9, 9)
    IF(PRESENT(Double9)) CALL SetArgDouble(Double9, 9)
    IF(PRESENT(Complex9)) CALL SetArgComplex(Complex9, 9)
    IF(PRESENT(DComplex9)) CALL SetArgDComplex(DComplex9, 9)
    IF(PRESENT(String9)) CALL SetArgString(String9, 9)

    ! Find the 10th argument.
    IF(PRESENT(Int10)) CALL SetArgInt(Int10, 10)
    IF(PRESENT(Real10)) CALL SetArgReal(Real10, 10)
    IF(PRESENT(Double10)) CALL SetArgDouble(Double10, 10)
    IF(PRESENT(Complex10)) CALL SetArgComplex(Complex10, 10)
    IF(PRESENT(DComplex10)) CALL SetArgDComplex(DComplex10, 10)
    IF(PRESENT(String10)) CALL SetArgString(String10, 10)

    ! Get the unit number
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit,NSections = NSections)

    ! If data type line is different that last, write a new section header and 
    ! data type line.
    IF(WritingNewSection(FileName,DataTypeLine).or.NewSection) THEN
       NSections = NSections + 1

       WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections
       WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
       WRITE(iUnit,'(A)') '#H#'

       WRITE(iUnit,'(A)') '#H# The following data types are written in this section.'
       WRITE(iUnit,'(2A)') '#DT# ', TRIM(DataTypeLine)
       CALL SetIOFileInfo(FileName,NSections = NSections, DataTypeLine = TRIM(DataTypeLine))

       ! If headers are present, write them
       IF(PRESENT(Headers)) THEN
          DO i1 = 1, SIZE(Headers)
             IF((LEN_TRIM(Headers(i1)).gt.0).and.(ICHAR(Headers(i1)(1:1)).ne.0)) &
                  & WRITE(iUnit,'(2A)') '#H# ', TRIM(Headers(i1))
          END DO
       END IF

       ! If column labels are present, get the columnlable format and write them.
       IF(PRESENT(ColumnLabels)) THEN
          ColumnLabelFormat = '(A5, '
          DO i1 = 1, NColumnLabels - 1
             ColumnLabelFormat = TRIM(ColumnLabelFormat) // TRIM(CLFormats(i1)) // ', '
          END DO
          ColumnLabelFormat = TRIM(ColumnLabelFormat) // TRIM(CLFormats(NColumnLabels)) // ')'
          WRITE(iUnit,ColumnLabelFormat) '#CL# ', (TRIM(ColumnLabels(i1)) // ' ', i1 = 1, NColumnLabels)
       END IF
       ! If only headers are present, don't update the section or write a data type line.
    ELSEIF(PRESENT(Headers).and.(NumArgs.eq.0)) THEN
       DO i1 = 1, SIZE(Headers)
          IF((LEN_TRIM(Headers(i1)).gt.0).and.(ICHAR(Headers(i1)(1:1)).ne.0)) &
               & WRITE(iUnit,'(2A)') '#H# ', TRIM(Headers(i1))
       END DO
    END IF

    ! Write data line.
    IF(HeaderData) THEN
       DO i1 = 1, NumData
          WRITE(iUnit, '(11A)') '#HD# ', (TRIM(DataArray(i1,i2)) // ' ', i2 = 1, NumArgs)
       END DO
    ELSE
       DO i1 = 1, NumData
          WRITE(iUnit, '(11A)') (TRIM(DataArray(i1,i2)) // ' ', i2 = 1, NumArgs)
       END DO
    END IF

    DEALLOCATE(DataArray)

  CONTAINS

    SUBROUTINE SetArgInt(Arg, iArg)
      INTEGER,INTENT(IN) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg

      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'Int')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            WRITE(DataArray(i1, iArg), IntFormat) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL WritePAD(DataArray(i1, iArg),Arg(i1))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 1
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) - 4
         ELSE
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) + 1
         END IF
         CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
      END IF
    END SUBROUTINE SetArgInt

    SUBROUTINE SetArgReal( Arg, iArg)
      REAL,INTENT(IN) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'Real')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            WRITE(DataArray(i1, iArg), RealFormat) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL WritePAD(DataArray(i1, iArg), npadrDefault, Arg(i1))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 1
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) - 4
         ELSE
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) + 1
         END IF
         CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
      END IF
    END SUBROUTINE SetArgReal

    SUBROUTINE SetArgDouble( Arg, iArg)
      DOUBLE PRECISION,INTENT(IN) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'Double')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            WRITE(DataArray(i1, iArg), DoubleFormat) Arg(i1)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL WritePAD(DataArray(i1, iArg), npadxDefault, Arg(i1))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 1
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) - 4
         ELSE
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) + 1
         END IF
         CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
      END IF
    END SUBROUTINE SetArgDouble

    SUBROUTINE SetArgComplex( Arg, iArg)
      COMPLEX,INTENT(IN) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'Complex')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            WRITE(DataArray(i1, iArg), ComplexFormat) REAL(Arg(i1)), IMAG(Arg(i1))
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL WritePAD(DataArray(i1, iArg), npadrDefault, Arg(i1))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 2
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) - 4
            CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) + 1
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         ELSE
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) + 1
            CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF
    END SUBROUTINE SetArgComplex

    SUBROUTINE SetArgDComplex( Arg, iArg)
      COMPLEX*16,INTENT(IN) :: Arg(:)
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'DComplex')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            WRITE(DataArray(i1, iArg), DComplexFormat) DBLE(Arg(i1)), DIMAG(Arg(i1))
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            CALL WritePAD(DataArray(i1, iArg), npadxDefault, Arg(i1))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 2
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) - 4
            CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) + 1
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         ELSE
            WRITE(TmpStr,'(I3)') INT(LEN_TRIM(DataArray(1, iArg))/2) + 1
            CLFormats(NColumnLabels - 1) = 'A' // TRIM(ADJUSTL(TmpStr))
            CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
         END IF
      END IF
    END SUBROUTINE SetArgDComplex

    SUBROUTINE SetArgString( Arg, iArg)
      CHARACTER*(*),INTENT(IN) :: Arg(:)
      CHARACTER(MaxStrLen) TmpStr
      INTEGER,INTENT(IN) :: iArg
      ! Check for errors, increment number of args,
      ! set numdata if this is the first argument, and write data to string
      ! array.
      CALL CheckForErrors(iArg,SIZE(Arg),'String')

      ! Write the data to the character array with the right format.
      ! For now fileformat is txt only.
      IF(txt) THEN
         DO i1 = 1, NumData
            ! Check that the string is not too long.
            TmpStr = TRIM(ADJUSTL(Arg(i1)))
            IF(LEN_TRIM(Arg(i1)).gt.MaxStrLen) CALL Error('Warning: String argument is too long in WriteArrayData. ' // &
                 & 'Will be truncated.' , StopProgram = .FALSE.)
            WRITE(DataArray(i1, iArg), StringFormat) TRIM(TmpStr)
         END DO
      ELSEIF(pad) THEN
         DO i1 = 1, NumData
            ! Check that the string is not too long.
            TmpStr = TRIM(ADJUSTL(Arg(i1)))
            IF(LEN_TRIM(Arg(i1)).gt.MaxStrLen) CALL Error('Warning: String argument is too long in WriteArrayData. ' // &
                 & 'Will be truncated.' , StopProgram = .FALSE.)

            CALL WritePAD(DataArray(i1,iArg),TRIM(TmpStr))
         END DO
      END IF

      ! Set the columnlabel format string
      IF(PRESENT(ColumnLabels)) THEN
         NColumnLabels = NColumnLabels + 1
         IF((iArg.eq.1).and.txt) THEN
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) - 4
         ELSE
            WRITE(TmpStr,'(I3)') LEN_TRIM(DataArray(1,iArg)) + 1
         END IF
         CLFormats(NColumnLabels) = 'A' // TRIM(ADJUSTL(TmpStr))
      END IF
    END SUBROUTINE SetArgString

    SUBROUTINE CheckForErrors(iArg,ArgSize,DTL)
      INTEGER,INTENT(IN) :: iArg, ArgSize
      CHARACTER*(*) DTL
      CHARACTER(2) ArgNum
      INTEGER iAllocateError
      ArgNum=''
      WRITE(ArgNum,'(I2)') iArg

      IF(IsPresent(iArg)) THEN
         CALL Error('Error in WriteData multiple arguments in place ' // TRIM(ArgNum) // '.')
      END IF

      ! Set IsPresent(iArg)
      IsPresent(iArg) = .TRUE.

      ! Check that number of arguments is not too large. Redundant, but oh well.
      IF(NumArgs.ge.iArg) CALL Error('Too many arguments passed to WriteData. Max is 10.')
      NumArgs = NumArgs + 1

      ! If this is the first data argument, set the size of the arrays and allocate space, 
      ! otherwise check that this array is the same size as the first array.        
      IF(iArg.eq.1) THEN
         NumData = ArgSize
         ALLOCATE(DataArray(NumData,10), STAT = iAllocateError)
         CALL CheckAllocation(iAllocateError, 'Error: Failed to allocate space for DataArray' //&
              & ' in subroutine WriteArrayData.')
      ELSE
         IF(NumData.ne.ArgSize) THEN
            CALL Error('Error in writearray: All arrays must be same size.')
         END IF
      END IF

      ! Set data type line
      DataTypeLine = TRIM(DataTypeLine) // ' ' // TRIM(DTL)
    END SUBROUTINE CheckForErrors

  END SUBROUTINE WriteArrayData


  ! SUBROUTINE WriteInt2D writes a 2D integer array to file.
  SUBROUTINE WriteInt2D(FileName,IntArray,Headers, FileType, WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    INTEGER,INTENT(IN) :: IntArray(:,:)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections, iFlType
    CHARACTER(100) IntFormat, DataStr
    CHARACTER(4) FlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, IFormat = IntFormat, &
         & NSections = NSections)

    N1 = SIZE(IntArray,1)
    N2 = SIZE(IntArray,2)

    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
    WRITE(iUnit,'(A)') '#H#'

    ! Write info about this array
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D integer array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write data to file.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          IF(iFlType.eq.itxt) THEN
             WRITE(iUnit,IntFormat, ADVANCE = 'NO',IOSTAT = iError) IntArray(i1,i2)
             IF(iError.ne.0) CALL Error('Error in WriteInt2D while trying to write to file: ' //&
                  & FileName)
          ELSEIF(iFlType.eq.ipad) THEN
             CALL WritePAD(DataStr,IntArray(i1,i2))
             WRITE(iUnit,'(A)', ADVANCE = 'NO',IOSTAT = iError) TRIM(DataStr) // ' '
          END IF
       END DO
       WRITE(iUnit,*, IOSTAT = iError)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteInt2D

  ! SUBROUTINE WriteReal2D writes a 2D integer array to file.
  SUBROUTINE WriteReal2D(FileName,RealArray,Headers,FileType,WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    REAL,INTENT(IN) :: RealArray(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections, iFlType
    CHARACTER(100) RealFormat, IntFormat, DataStr
    CHARACTER(4) FlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)

    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader
    
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, RFormat = RealFormat, &
         & NSections = NSections)

    N1 = SIZE(RealArray,1)
    N2 = SIZE(RealArray,2)

    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
    WRITE(iUnit,'(A)') '#H#'

    ! Write some other information
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D real array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write the data.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          IF(iFlType.eq.itxt) THEN
             WRITE(iUnit,RealFormat, ADVANCE = 'NO',IOSTAT = iError) RealArray(i1,i2)
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          ELSEIF(iFlType.eq.ipad) THEN
             CALL WritePAD(DataStr,npadrDefault,RealArray(i1,i2))
             WRITE(iUnit,'(A)', ADVANCE = 'NO',IOSTAT = iError) TRIM(DataStr) // ' '
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          END IF
       END DO
       WRITE(iUnit,*, IOSTAT = iError)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteReal2D

  ! SUBROUTINE WriteReal2D writes a 2D integer array to file.
  SUBROUTINE WriteDouble2D(FileName,DoubleArray,Headers,FileType,WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    DOUBLE PRECISION,INTENT(IN) :: DoubleArray(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections, iFlType
    CHARACTER(100) DoubleFormat, IntFormat, DataStr
    CHARACTER(4) FlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)

    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader

    iFlType = 1
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, DFormat = DoubleFormat, &
         & NSections = NSections)

    N1 = SIZE(DoubleArray,1)
    N2 = SIZE(DoubleArray,2)

    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
    WRITE(iUnit,'(A)') '#H#'

    ! Write some other information
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D double array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write the data.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          IF(iFlType.eq.itxt) THEN
             WRITE(iUnit,DoubleFormat,IOSTAT = iError, ADVANCE = 'NO') DoubleArray(i1,i2)
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          ELSEIF(iFlType.eq.ipad) THEN
             CALL WritePAD(DataStr,npadxDefault,DoubleArray(i1,i2))
             WRITE(iUnit,'(A)',IOSTAT = iError, ADVANCE = 'NO') TRIM(DataStr) // ' '
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          END IF
       END DO
       WRITE(iUnit,*)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteDouble2D

  ! SUBROUTINE WriteComplex2D writes a 2D integer array to file.
  SUBROUTINE WriteComplex2D(FileName,ComplexArray,Headers,FileType,WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    COMPLEX,INTENT(IN) :: ComplexArray(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections, iFlType
    CHARACTER(100) ComplexFormat, IntFormat, DataStr
    CHARACTER(4) FlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)

    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader
    
    iFlType = itxt
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, CFormat = ComplexFormat, &
         & NSections = NSections)

    N1 = SIZE(ComplexArray,1)
    N2 = SIZE(ComplexArray,2)
    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
    WRITE(iUnit,'(A)') '#H#'

    ! Write some other information
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D complex array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write the data.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          IF(iFlType.eq.itxt) THEN
             WRITE(iUnit,ComplexFormat, ADVANCE = 'NO',IOSTAT = iError) REAL(ComplexArray(i1,i2)),  &
                  & IMAG(ComplexArray(i1,i2))
             IF(iError.ne.0) CALL Error('Error in WriteComplex2D while trying to write to file: ' //&
                  & FileName)
          ELSEIF(iFlType.eq.ipad) THEN
             CALL WritePAD(DataStr,npadrDefault,ComplexArray(i1,i2))
             WRITE(iUnit,'(A)', ADVANCE = 'NO',IOSTAT = iError) TRIM(DataStr) // ' '
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          END IF
       END DO
       WRITE(iUnit,*, IOSTAT = iError)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteComplex2D

  ! SUBROUTINE WriteDComplex2D writes a 2D integer array to file.
  SUBROUTINE WriteDComplex2D(FileName,DComplexArray,Headers,FileType,WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    COMPLEX*16,INTENT(IN) :: DComplexArray(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections, iFlType
    CHARACTER(100) DComplexFormat, DataStr
    CHARACTER(4) FlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)

    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader
    
    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, DCFormat = DComplexFormat, &
         & NSections = NSections)

    N1 = SIZE(DComplexArray,1)
    N2 = SIZE(DComplexArray,2)

    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#DF# This section written in ' // FlType // '.'
    WRITE(iUnit,'(A)') '#H#'

    ! Write some other information
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D double complex array with sizes ', N1, N2
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D complex array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write the data.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          IF(iFlType.eq.itxt) THEN
             WRITE(iUnit,DComplexFormat, ADVANCE = 'NO',IOSTAT = iError) DBLE(DComplexArray(i1,i2)),  &
                  & DIMAG(DComplexArray(i1,i2))
             IF(iError.ne.0) CALL Error('Error in WriteDComplex2D while trying to write to file: ' // &
                  & FileName)
          ELSEIF(iFlType.eq.ipad) THEN
             CALL WritePAD(DataStr,npadrDefault,DComplexArray(i1,i2))
             WRITE(iUnit,'(A)', ADVANCE = 'NO',IOSTAT = iError) TRIM(DataStr) // ' '
             IF(iError.ne.0) CALL Error('Error in WriteReal2D while trying to write to file: ' //&
                  & FileName)
          END IF
       END DO
       WRITE(iUnit,*, IOSTAT = iError)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteDComplex2D

  ! SUBROUTINE WriteString2D writes a 2D integer array to file.
  SUBROUTINE WriteString2D(FileName,StringArray,Headers,FileType,WriteDataInHeader)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: Headers(:), FileType
    CHARACTER*(*),INTENT(IN) :: StringArray(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: WriteDataInHeader

    INTEGER iUnit, N1, N2, i1, i2, iError, NSections
    CHARACTER(100) StringFormat
    CHARACTER(4) FlType
    INTEGER iFlType
    LOGICAL HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)

    iFlType = 1
    HeaderData = .FALSE.
    IF(PRESENT(WriteDataInHeader)) HeaderData = WriteDataInHeader
    
    CALL OpenFl(FileName)
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, SFormat = StringFormat, &
         & NSections = NSections)

    N1 = SIZE(StringArray,1)
    N2 = SIZE(StringArray,2)

    ! Write the section number.
    WRITE(iUnit,'(A,I4)') '#SN#   Section: ', NSections + 1
    WRITE(iUnit,'(A)') '#H#'

    ! Write some other information
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D double complex array with sizes ', N1, N2
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D double complex array with sizes ', N1, N2
    WRITE(iUnit,'(A,I4,I4)') '#DT# 2D complex array with sizes ', N1, N2
    WRITE(iUnit,'(A)') '#H# File is organized as follows:  Array(1,i)     Array(1,i+1)    Array(1,i+2)  . . .'
    WRITE(iUnit,'(A)') '#H#                                Array(2,i)'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'
    WRITE(iUnit,'(A)') '#H#                                     .'

    ! Write headers if they exist
    IF(PRESENT(Headers)) THEN
       CALL WriteData(FileName, Headers = Headers)
    END IF

    ! Write the data.
    DO i1 = 1, N1
       ! If User has specified that data should be written in the header, do so
       IF(HeaderData.and.(iFlType.eq.itxt)) WRITE(iUnit, '(A)', ADVANCE = 'NO', IOSTAT = iError) '#HD# '
       DO i2 = 1, N2
          WRITE(iUnit,StringFormat, ADVANCE = 'NO',IOSTAT = iError) TRIM(ADJUSTL(StringArray(i1,i2)))
          IF(iError.ne.0) CALL Error('Error in WriteString2D while trying to write to file: ' //&
               & FileName)

          ! Write some space between elements
          WRITE(iUnit,'(A)', ADVANCE = 'NO') '  '
       END DO
       WRITE(iUnit,*, IOSTAT = iError)
    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE WriteString2D

  ! SUBROUTINE ReadInt2D reads a 2D integer array from file.
  SUBROUTINE ReadInt2D(FileName,IntArray, N1, N2, SectionNumber, Headers, ReadDataFromHeader, &
       & FileType, CommentCharacters)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input: 
    ! FileName          - name of file to read.
    ! FileType          - 'PAD' or 'TXT'
    ! CommentCharacters - string of characters, any of which can signify a header line. By default
    !                     'TXT' headers start with '#', '!', 'c', 'C', or '*'. When reading 'PAD' data
    !                     comment character is '#' only.
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output:
    ! Headers           - Read header lines into this variable.
    ! IntArray       - Holds data read from file.
    ! N1, N2            - Dimensions of the array as read from file.
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    INTEGER,INTENT(OUT) :: IntArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local Variables:
    ! iUnit        - unit number of the file.
    ! iError       - error code passed back from read.
    ! NSections    - Current section of file.
    ! iFlType      - index of file type. (1 = txt, 2 = pad)
    ! NE1, NE2     - used to hold the dimensions of the array as read from file.
    ! Length       - Length of the current word.
    ! IntFormat - Format to use for double precision data.
    ! CmtChars     - Holds string of possible comment characters.
    ! DataTypeLine - Holds the type of data being read from file, if defined.
    ! DataStr      - Holds the current word as read from file.
    ! FlType       - Holds the file type ('TXT' or 'PAD').
    ! EOF/EOR      - True if end of file/record has been reached
    ! IsHead       - True if a comment character has been found at beginning of line.
    INTEGER iUnit, iError, NSections, iFlType, NE1, NE2, Length
    CHARACTER(100) IntFormat, CmtChars, DataTypeLine
    CHARACTER(40) DataStr
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData
 
    ! Loop Variables: 
    INTEGER i1, i2

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    ! Initialization
    NE1 = 0
    NE2 = 0
    
    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          ! If filetype is 'PAD' comment characters must be '#' only.
          iFlType = ipad
          CmtChars='#'
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF
    ! End Initialization

    ! Open the file for reading.
    CALL OpenFl(FileName, FileAction = 'READ')

    ! Get some info about the file.
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, IFormat = IntFormat, &
         & NSections = NSections)

    ! Set the maximum size of the array to be read.
    N1 = SIZE(IntArray,1)
    N2 = SIZE(IntArray,2)

    ! Read headers and save if Headers has been passed.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
    END IF

    ! Backspace since ReadHeaders reads one line passed the last header.
    BACKSPACE(iUnit)

    ! Loop over dimensions of the array.
    DO i1 = 1, N1
       DO i2 = 1, N2

          ! Read the next word from the file and save in DataStr
          CALL ReadNextWord(iUnit,DataStr,Length,IsHeader,EOR,EOF,CmtChars)

          ! If end of file, set EOF and return.
          IF(EOF) THEN
             CALL SetIOFileInfo(FileName, EOF = EOF)
             RETURN
          END IF

          ! If this is a header line backspace and read headers.
          IF(IsHeader) THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
             
             ! Check if we are reading a new section. If so, return.
             IF(ReadingNewSection(FileName, NSections)) RETURN
          END IF

          ! If line is not blank, set the data array element.
          IF(Length.gt.0) THEN
             IF(iFlType.eq.itxt) THEN
                READ(DataStr(1:Length),FMT=*,IOSTAT = iError) IntArray(i1,i2)
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to read from file: ' //&
                     & FileName)
             ELSEIF(iFlType.eq.ipad) THEN
                IF(Length.gt.0) CALL ReadPAD(DataStr(1:Length),IntArray(i1,i2))
             END IF

             ! Set dimension 2.
             NE2 = i2
          END IF

          ! If end of record reached, exit this loop.
          IF(EOR) EXIT

       END DO

       ! Advance to the next line.
       IF(.not.EOR) READ(iUnit,*,END=20)              

       ! Set dimension 1
       N1 = i1

    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)

    ! End of file reached.
20  CONTINUE
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)
  END SUBROUTINE ReadInt2D


  ! SUBROUTINE ReadReal2D reads a 2D integer array from file.
  SUBROUTINE ReadReal2D(FileName,RealArray, N1, N2, SectionNumber, Headers, ReadDataFromHeader, &
       & FileType, CommentCharacters)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input: 
    ! FileName          - name of file to read.
    ! FileType          - 'PAD' or 'TXT'
    ! CommentCharacters - string of characters, any of which can signify a header line. By default
    !                     'TXT' headers start with '#', '!', 'c', 'C', or '*'. When reading 'PAD' data
    !                     comment character is '#' only.
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output:
    ! Headers           - Read header lines into this variable.
    ! RealArray       - Holds data read from file.
    ! N1, N2            - Dimensions of the array as read from file.
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    REAL,INTENT(OUT) :: RealArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local Variables:
    ! iUnit        - unit number of the file.
    ! iError       - error code passed back from read.
    ! NSections    - Current section of file.
    ! iFlType      - index of file type. (1 = txt, 2 = pad)
    ! NE1, NE2     - used to hold the dimensions of the array as read from file.
    ! Length       - Length of the current word.
    ! RealFormat - Format to use for double precision data.
    ! CmtChars     - Holds string of possible comment characters.
    ! DataTypeLine - Holds the type of data being read from file, if defined.
    ! DataStr      - Holds the current word as read from file.
    ! FlType       - Holds the file type ('TXT' or 'PAD').
    ! EOF/EOR      - True if end of file/record has been reached
    ! IsHead       - True if a comment character has been found at beginning of line.
    INTEGER iUnit, iError, NSections, iFlType, NE1, NE2, Length
    CHARACTER(100) RealFormat, CmtChars, DataTypeLine
    CHARACTER(40) DataStr
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData
 
    ! Loop Variables: 
    INTEGER i1, i2

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    ! Initialization
    NE1 = 0
    NE2 = 0

    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          ! If filetype is 'PAD' comment characters must be '#' only.
          iFlType = ipad
          CmtChars='#'
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF
    ! End Initialization

    ! Open the file for reading.
    CALL OpenFl(FileName, FileAction = 'READ')

    ! Get some info about the file.
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, RFormat = RealFormat, &
         & NSections = NSections)

    ! Set the maximum size of the array to be read.
    N1 = SIZE(RealArray,1)
    N2 = SIZE(RealArray,2)

    ! Read headers and save if Headers has been passed.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
    END IF

    ! Backspace since ReadHeaders reads one line passed the last header.
    BACKSPACE(iUnit)

    ! Loop over dimensions of the array.
    DO i1 = 1, N1
       DO i2 = 1, N2

          ! Read the next word from the file and save in DataStr
          CALL ReadNextWord(iUnit,DataStr,Length,IsHeader,EOR,EOF,CmtChars)

          ! If end of file, set EOF and return.
          IF(EOF) THEN
             CALL SetIOFileInfo(FileName, EOF = EOF)
             RETURN
          END IF

          ! If this is a header line backspace and read headers.
          IF(IsHeader) THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
             
             ! Check if we are reading a new section. If so, return.
             IF(ReadingNewSection(FileName, NSections)) RETURN
          END IF

          ! If line is not blank, set the data array element.
          IF(Length.gt.0) THEN
             IF(iFlType.eq.itxt) THEN
                READ(DataStr(1:Length),FMT=*,IOSTAT = iError) RealArray(i1,i2)
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to write to file: ' //&
                     & FileName)
             ELSEIF(iFlType.eq.ipad) THEN
                IF(Length.gt.0) CALL ReadPAD(DataStr(1:Length),npadrDefault,RealArray(i1,i2))
             END IF

             ! Set dimension 2.
             NE2 = i2
          END IF

          ! If end of record reached, exit this loop.
          IF(EOR) EXIT

       END DO

       ! Advance to the next line.
       IF(.not.EOR) READ(iUnit,*,END=20)              

       ! Set dimension 1
       N1 = i1

    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)

    ! End of file reached.
20  CONTINUE
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)
  END SUBROUTINE ReadReal2D


  ! SUBROUTINE ReadDouble2D reads a 2D integer array from file.
  SUBROUTINE ReadDouble2D(FileName,DoubleArray, N1, N2, SectionNumber,  Headers, ReadDataFromHeader,&
       & FileType, CommentCharacters)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input: 
    ! FileName          - name of file to read.
    ! FileType          - 'PAD' or 'TXT'
    ! CommentCharacters - string of characters, any of which can signify a header line. By default
    !                     'TXT' headers start with '#', '!', 'c', 'C', or '*'. When reading 'PAD' data
    !                     comment character is '#' only.
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output:
    ! Headers           - Read header lines into this variable.
    ! DoubleArray       - Holds data read from file.
    ! N1, N2            - Dimensions of the array as read from file.
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    DOUBLE PRECISION,INTENT(OUT) :: DoubleArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local Variables:
    ! iUnit        - unit number of the file.
    ! iError       - error code passed back from read.
    ! NSections    - Current section of file.
    ! iFlType      - index of file type. (1 = txt, 2 = pad)
    ! NE1, NE2     - used to hold the dimensions of the array as read from file.
    ! Length       - Length of the current word.
    ! DoubleFormat - Format to use for double precision data.
    ! CmtChars     - Holds string of possible comment characters.
    ! DataTypeLine - Holds the type of data being read from file, if defined.
    ! DataStr      - Holds the current word as read from file.
    ! FlType       - Holds the file type ('TXT' or 'PAD').
    ! EOF/EOR      - True if end of file/record has been reached
    ! IsHead       - True if a comment character has been found at beginning of line.
    INTEGER iUnit, iError, NSections, iFlType, NE1, NE2, Length
    CHARACTER(100) DoubleFormat, CmtChars, DataTypeLine
    CHARACTER(40) DataStr
    CHARACTER(200) line
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData
 
    ! Loop Variables: 
    INTEGER i1, i2

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    ! Initialization
    NE1 = 0
    NE2 = 0

    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          ! If filetype is 'PAD' comment characters must be '#' only.
          iFlType = ipad
          CmtChars='#'
       ELSE
          CALL Error('Illegal FileType passed to Read2D.')
       END IF
    END IF
    ! End Initialization

    ! Open the file for reading.
    CALL OpenFl(FileName, FileAction = 'READ')

    ! Get some info about the file.
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, DFormat = DoubleFormat, &
         & NSections = NSections)

    ! Set the maximum size of the array to be read.
    N1 = SIZE(DoubleArray,1)
    N2 = SIZE(DoubleArray,2)

    ! Read headers and save if Headers has been passed.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, line, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, line, CommentCharacters = CmtChars)
    END IF

    ! Backspace since ReadHeaders reads one line passed the last header.
    BACKSPACE(iUnit)

    ! Loop over dimensions of the array.
    DO i1 = 1, N1
       DO i2 = 1, N2

          ! Read the next word from the file and save in DataStr
          CALL ReadNextWord(iUnit,DataStr,Length,IsHeader,EOR,EOF,CmtChars)

          ! If end of file, set EOF and return.
          IF(EOF) THEN
!             CALL SetIOError('Array_Too_Small')
             CALL SetIOFileInfo(FileName, EOF = EOF)
             RETURN
          END IF
          
          ! If this is a header line backspace and read headers.
          IF(IsHeader) THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, line, CommentCharacters = CmtChars)
             
             ! Check if we are reading a new section. If so, return.
             IF(ReadingNewSection(FileName, NSections)) THEN
                RETURN
             ELSE
                ! Set Comment_Error to true. By default, this will produce a warning.
                CALL Error('Unexpected comment character found in file ' // FileName // '.')
                CALL Error(line)
             END IF
          END IF

          ! If line is not blank, set the data array element.
          IF(Length.gt.0) THEN
             IF(iFlType.eq.itxt) THEN
                READ(DataStr(1:Length),FMT=*,IOSTAT = iError) DoubleArray(i1,i2)
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to write to file: ' //&
                     & FileName)
             ELSEIF(iFlType.eq.ipad) THEN
                IF(Length.gt.0) CALL ReadPAD(DataStr(1:Length),npadrDefault,DoubleArray(i1,i2))
             END IF

             ! Set dimension 2.
             NE2 = i2
          END IF

          ! If end of record reached, exit this loop.
          IF(EOR) EXIT

       END DO

       ! Advance to the next line if needed.
       IF(.not.EOR) READ(iUnit,*,END=20)              

       ! Set dimension 1
       N1 = i1

    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)

    RETURN
    ! End of file reached.
20  CONTINUE
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)
  END SUBROUTINE ReadDouble2D


  ! SUBROUTINE ReadComplex2D reads a 2D integer array from file.
  SUBROUTINE ReadComplex2D(FileName,ComplexArray, N1, N2, SectionNumber, Headers, ReadDataFromHeader, &
       & FileType, CommentCharacters)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input: 
    ! FileName          - name of file to read.
    ! FileType          - 'PAD' or 'TXT'
    ! CommentCharacters - string of characters, any of which can signify a header line. By default
    !                     'TXT' headers start with '#', '!', 'c', 'C', or '*'. When reading 'PAD' data
    !                     comment character is '#' only.
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output:
    ! Headers           - Read header lines into this variable.
    ! ComplexArray       - Holds data read from file.
    ! N1, N2            - Dimensions of the array as read from file.
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    COMPLEX,INTENT(OUT) :: ComplexArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local Variables:
    ! iUnit        - unit number of the file.
    ! iError       - error code passed back from read.
    ! NSections    - Current section of file.
    ! iFlType      - index of file type. (1 = txt, 2 = pad)
    ! NE1, NE2     - used to hold the dimensions of the array as read from file.
    ! Length       - Length of the current word.
    ! ComplexFormat - Format to use for double precision data.
    ! CmtChars     - Holds string of possible comment characters.
    ! DataTypeLine - Holds the type of data being read from file, if defined.
    ! DataStr      - Holds the current word as read from file.
    ! FlType       - Holds the file type ('TXT' or 'PAD').
    ! EOF/EOR      - True if end of file/record has been reached
    ! IsHead       - True if a comment character has been found at beginning of line.
    REAL ReData, ImData
    INTEGER iUnit, iError, NSections, iFlType, NE1, NE2, ReLength, ImLength
    CHARACTER(100) ComplexFormat, CmtChars, DataTypeLine
    CHARACTER(40) ReDataStr, ImDataStr
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData
 
    ! Loop Variables: 
    INTEGER i1, i2

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    ! Initialization
    NE1 = 0
    NE2 = 0

    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          ! If filetype is 'PAD' comment characters must be '#' only.
          iFlType = ipad
          CmtChars='#'
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF
    ! End Initialization

    ! Open the file for reading.
    CALL OpenFl(FileName, FileAction = 'READ')

    ! Get some info about the file.
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, CFormat = ComplexFormat, &
         & NSections = NSections)

    ! Set the maximum size of the array to be read.
    N1 = SIZE(ComplexArray,1)
    N2 = SIZE(ComplexArray,2)

    ! Read headers and save if Headers has been passed.
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars)
    END IF

    ! Backspace since ReadHeaders reads one line passed the last header.
    BACKSPACE(iUnit)

    ! Loop over dimensions of the array.
    DO i1 = 1, N1
       DO i2 = 1, N2
          ! Read the next word from the file and save in ReDataStr
          CALL ReadNextWord(iUnit,ReDataStr,ReLength,IsHeader,EOR,EOF,CmtChars)
          IF(EOF.or.EOR) CALL Error('Error: Unexpected end of record while reading from ' // FileName // '.')
          ! Read the next word from the file and save in ImDataStr
          CALL ReadNextWord(iUnit,ImDataStr,ImLength,IsHeader,EOR,EOF,CmtChars)

          ! If end of file, set EOF and return.
          IF(EOF) THEN
             CALL SetIOFileInfo(FileName, EOF = EOF)
             RETURN
          END IF

          ! If this is a header line backspace and read headers.
          IF(IsHeader) THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars)
             
             ! Check if we are reading a new section. If so, return.
             IF(ReadingNewSection(FileName, NSections)) RETURN
          END IF

          ! If line is not blank, set the data array element.
          IF(iFlType.eq.itxt) THEN
             IF(ReLength.gt.0) THEN
                READ(ReDataStr(1:ReLength),FMT=*, iostat = iError) ReData
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to read to file: ' //&
                     & FileName)
             END IF
             IF(ImLength.gt.0) THEN
                READ(ImDataStr(1:ImLength),FMT=*, iostat = iError) ImData
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to read to file: ' //&
                     & FileName)
             ENDIF
          ELSEIF(iFlType.eq.ipad) THEN
             IF(ReLength.gt.0) CALL ReadPAD(ReDataStr(1:ReLength),npadrDefault,ReData)
             IF(ImLength.gt.0) CALL ReadPAD(ImDataStr(1:ImLength),npadrDefault,ImData)
          END IF

          ! Set the data array.
          ComplexArray(i1,i2) = ReData + (0.0,1.0)*ImData

          ! Set dimension 2.
          NE2 = i2
          
          ! If end of record reached, exit this loop.
          IF(EOR) EXIT

       END DO

       ! Advance to the next line.
       IF(.not.EOR) READ(iUnit,*,END=20)              

       ! Set dimension 1
       N1 = i1

    END DO

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)

    ! End of file reached.
20  CONTINUE
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)
  END SUBROUTINE ReadComplex2D


  ! SUBROUTINE ReadInt2D reads a 2D integer array from file.
  SUBROUTINE ReadDComplex2D(FileName,DComplexArray, N1, N2, SectionNumber, Headers, ReadDataFromHeader, &
       & FileType, CommentCharacters)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    COMPLEX*16,INTENT(OUT) :: DComplexArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader

    INTEGER iUnit, i1, i2, iError, NSections, iFlType, NWords, NE1, NE2, Length
    CHARACTER(100) DComplexFormat, CmtChars, DataTypeLine
    CHARACTER(100) ReDataStr, ImDataStr
    CHARACTER(20) Words(20)
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData
    DOUBLE PRECISION ReArg, ImArg

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    NE1 = 0
    NE2 = 0

    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Write2D.')
       END IF
    END IF

    CALL OpenFl(FileName, FileAction = 'READ')
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, DCFormat = DComplexFormat, &
         & NSections = NSections)

    N1 = SIZE(DComplexArray,1)
    N2 = SIZE(DComplexArray,2)

    ! Read headers and save if they exist
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars)
    END IF
    BACKSPACE(iUnit)

    ! Read data from file.
    DO i1 = 1, N1
       DO i2 = 1, N2
          CALL ReadNextWord(iUnit,ReDataStr,Length,IsHeader,EOR,EOF,CmtChars)
          CALL ReadNextWord(iUnit,ImDataStr,Length,IsHeader,EOR,EOF,CmtChars)
          IF(EOF) EXIT
          ! If this is a header line backspace and read headers.
          IF(ReDataStr(1:1).eq.'#') THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, ReDataStr, CommentCharacters = CmtChars)
             IF(ReadingNewSection(FileName, NSections)) RETURN
          END IF
          IF((LEN_TRIM(ReDataStr).gt.0).and.(LEN_TRIM(ImDataStr).gt.0)) THEN
             IF(iFlType.eq.itxt) THEN
                READ(ReDataStr,FMT=*,IOSTAT = iError) ReArg
                READ(ImDataStr,*,IOSTAT = iError) ImArg
                IF(iError.ne.0) CALL Error('Error in Read2D while trying to read from file: ' //&
                     & FileName)
             ELSEIF(iFlType.eq.ipad) THEN
                CALL ReadPAD(ReDataStr,npadrDefault,ReArg)
                CALL ReadPAD(ImDataStr,npadrDefault,ImArg)
             END IF
             DComplexArray(i2,i1) = ReArg + (0.0,1.0)*ImArg
             NE1 = NE1 + 1
          END IF
          IF(EOR) EXIT
       END DO
       IF(EOF) EXIT
       NE2 = NE2 + 1
    END DO
    EOF = .FALSE.
5   CONTINUE

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE ReadDComplex2D

  ! SUBROUTINE ReadInt2D reads a 2D integer array from file.
  SUBROUTINE ReadString2D(FileName, StringArray, N1, N2, SectionNumber, Headers, ReadDataFromHeader, &
       & FileType, CommentCharacters)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileType, CommentCharacters
    CHARACTER*(*),INTENT(OUT) :: StringArray(:,:)
    INTEGER,INTENT(OUT) :: N1, N2
    INTEGER,INTENT(IN),OPTIONAL :: SectionNumber
    LOGICAL,INTENT(IN),OPTIONAL :: ReadDataFromHeader

    INTEGER iUnit, i1, i2, iError, NSections, iFlType, NWords, NE1, NE2, Length
    CHARACTER(100) StringFormat, CmtChars, DataTypeLine
    CHARACTER(MaxStrLen) DataStr
    CHARACTER(20) Words(20)
    CHARACTER(4) FlType
    LOGICAL EOF, EOR, IsHeader, HeaderData

    ! Define file types.
    INTEGER itxt, ipad
    PARAMETER(itxt = 1, ipad = 2)
    
    NE1 = 0
    NE2 = 0

    HeaderData = .FALSE.
    IF(PRESENT(ReadDataFromHeader)) HeaderData = ReadDataFromHeader

    IF(PRESENT(SectionNumber)) THEN
       CALL ReadToSection(FileName, SectionNumber)
    END IF

    iFlType = iFileTypeDefault
    IF(PRESENT(CommentCharacters)) THEN 
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    IF(PRESENT(FileType)) THEN
       FlType = TRIM(ADJUSTL(FileType))
       CALL Upper(FlType)
       IF(TRIM(FlType).eq.'TXT') THEN
          iFlType = itxt
       ELSEIF(TRIM(FlType).eq.'PAD') THEN
          iFlType = ipad
       ELSE
          CALL Error('Illegal FileType passed to Read2D.')
       END IF
    END IF

    CALL OpenFl(FileName, FileAction = 'READ')
    CALL GetIOFileInfo(FileName, UnitNumber = iUnit, SFormat = StringFormat, &
         & NSections = NSections)

    N1 = SIZE(StringArray,1)
    N2 = SIZE(StringArray,2)

    ! Read headers and save if they exist
    IF(PRESENT(Headers)) THEN
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars, Headers = Headers)
    ELSE
       CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
    END IF
    BACKSPACE(iUnit)

    ! Read data to file.
    DO i1 = 1, N1
       DO i2 = 1, N2
          CALL ReadNextWord(iUnit,DataStr,Length,IsHeader,EOR,EOF,CmtChars)
          IF(EOF) EXIT
          ! If this is a header line backspace and read headers.
          IF(DataStr(1:1).eq.'#') THEN
             BACKSPACE(iUnit)
             CALL ReadHeaders(FileName, DataStr, CommentCharacters = CmtChars)
             IF(ReadingNewSection(FileName, NSections)) RETURN
          END IF
          IF(LEN_TRIM(DataStr).gt.0) THEN
             StringArray(i2,i1) = TRIM(ADJUSTL(DataStr))
             NE1 = NE1 + 1
          END IF
          IF(EOR) EXIT
       END DO
       IF(EOF) EXIT
       NE2 = NE2 + 1
    END DO
10  CONTINUE

    ! Increment the section number
    CALL SetIOFileInfo(FileName, NSections = NSections + 1)
  END SUBROUTINE ReadString2D

  ! Reads through the file FileName until it finds the section specified by
  ! SectionNumber. 
  SUBROUTINE ReadToSection(FileName,SectionNumber)
    CHARACTER*(*),INTENT(IN) :: FileName
    INTEGER,INTENT(IN) :: SectionNumber

    INTEGER iUnit, NSections, iError
    CHARACTER(10) ErrorNum
    CHARACTER(20) TmpStr
    CHARACTER(1000) Line

    ! Get the unit number.
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit, NSections = NSections)

    ! Rewind to start of file
    REWIND(iUnit)

    ! Read until section.
    NSections = 0
    DO
       READ(iUnit,'(A)',IOSTAT = iError, END = 5) Line

       ! If error stop.          
       IF(iError.ne.0) THEN
          WRITE(ErrorNum,'(I5)') iError
          CALL Error('Error ocurred in ReadToSection', StopProgram = .FALSE.)
          CALL Error('while reading from file ' // TRIM(ADJUSTL(FileName)) // &
               & '.',StopProgram = .FALSE.)
          CALL Error('Error: ' // ErrorNum)
       END IF

       TmpStr = ADJUSTL(Line)
       IF(TmpStr(1:4).eq.'#SN#') NSections = NSections + 1
       IF(NSections.eq.SectionNumber) EXIT
    END DO

    ! Set the section number.
    CALL SetIOFileInfo(FileName, NSections = NSections)

    RETURN

5   CONTINUE ! End of file reached. Set EOF and return.
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)

  END SUBROUTINE ReadToSection

  ! Reads lines until it finds one that is not a header. The next non
  !-header (non-blank) line will be passed back in Line.
  ! FileName - name of file to read from.
  ! Line - next data line
  ! CommentCharacters - string of characters, any of which will be
  ! interpreted as comment characters when found at the beginning of
  ! a line.
  ! Headers - optional: will store the headers found in Headers.
  SUBROUTINE ReadHeaders(FileName,Line,CommentCharacters,Headers)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(OUT) :: Line
    CHARACTER*(*),INTENT(IN),OPTIONAL :: CommentCharacters
    CHARACTER*(*),INTENT(OUT),OPTIONAL :: Headers(:)

    INTEGER iUnit, iError, NHead, MaxNHead, NSections
    CHARACTER(10) ErrorNum, CmtChars
    CHARACTER(200) TmpStr

    LOGICAL IsHead

    INTEGER i1

    ! Initialization
    IsHead = .TRUE.
    IF(PRESENT(Headers)) MaxNHead = SIZE(Headers)
    NHead = 0
    IF(PRESENT(CommentCharacters)) THEN
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    ! Get the unit number.
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit, NSections = NSections)
    
    ! Read until next non-header line, saving header lines if Headers variable
    ! is provided.
    DO
       READ(iUnit,'(A)',IOSTAT = iError, END = 5) Line

       ! If error stop.          
       IF(iError.ne.0) THEN
          WRITE(ErrorNum,'(I5)') iError
          CALL Error('Error ocurred in ReadHeaders', StopProgram = .FALSE.)
          CALL Error('while reading from file ' // TRIM(ADJUSTL(FileName)) // &
               & '.',StopProgram = .FALSE.)
          CALL Error('Error: ' // ErrorNum)
       END IF

       TmpStr = TRIM(ADJUSTL(Line))

       ! First run through a few cases that look like headers but are not.
       IF(TmpStr(1:4).eq.'#SN#') THEN 
          ! New section has been reached. Increment section #
          CALL SetIOFileInfo(FileName,NSections = NSections + 1)
          NSections = NSections + 1

       ELSEIF(TmpStr(1:4).eq.'#DT#') THEN
          ! This is a data type line. Set DataTypeLine and continue.
          TmpStr = TRIM(ADJUSTL(Line(5:LEN(Line))))
          CALL SetIOFileInfo(FileName, DataTypeLine = TmpStr)

       ELSEIF(TmpStr(1:4).eq.'#DF#') THEN
          ! This is a data format line. Eventually add an automatic reader.
          CONTINUE

       ELSEIF(TmpStr(1:4).eq.'#HD#') THEN
          ! This is header data: return.
          IF(NSections.eq.0) NSections = 1
          CALL SetIOFileInfo(FileName, NSections = NSections)
          RETURN

       ELSEIF(LEN_TRIM(TmpStr).ne.0) THEN ! This is not a blank line, check for header line.
          ! Header lines are anything that start with any of the CommentCharacters
          ! if provided, otherwise they start with #, *, c, C, or ! 
          IsHead = .FALSE.
          DO i1 = 1, LEN_TRIM(CmtChars)
             IF(TmpStr(1:1).eq.CmtChars(i1:i1)) THEN
                ! if this is a header, set headers if specified
                IF(PRESENT(Headers)) THEN
                   NHead = NHead + 1
                   IF(MaxNHead.le.NHead) Headers(NHead) = TRIM(ADJUSTL(Line))
                END IF
                IsHead = .TRUE.
                EXIT
             END IF
          END DO

          ! If this line was not a header, backspace and return.
          IF(.not.IsHead) THEN
             IF(NSections.eq.0) NSections = 1
             CALL SetIOFileInfo(FileName, NSections = NSections)
             RETURN
          END IF
       END IF
    END DO

5   CONTINUE ! End of file reached. Set EOF for this file.
    CALL SetIOFileInfo(FileName, EOF = .TRUE.)
    
  END SUBROUTINE ReadHeaders

  ! Checks if the end of the file has been reached.
  LOGICAL FUNCTION EndOfFile(FileName)
    CHARACTER*(*) FileName

    CALL GetIOFileInfo(FileName, EOF = EndOfFile)
  END FUNCTION EndOfFile

  ! Checks if the next line is data (not a header), and backspaces.
  LOGICAL FUNCTION NextLineIsData(FileName, CommentCharacters)
    CHARACTER*(*),INTENT(IN) :: FileName
    CHARACTER*(*),INTENT(IN),OPTIONAL :: CommentCharacters

    INTEGER iUnit, iError, i1
    CHARACTER(1000) Line
    CHARACTER(20) TmpStr, CmtChars
    CHARACTER(5) ErrorNum
    LOGICAL IsHead, EOF

    NextLineIsData = .FALSE.
    IF(PRESENT(CommentCharacters)) THEN
       CmtChars = CommentCharacters
    ELSE
       CmtChars = DefaultCommentCharacters
    END IF

    ! Get the unit number.
    CALL GetIOFileInfo(FileName,UnitNumber = iUnit, EOF = EOF)
    
    ! If EOF, return false.
    IF(EOF) THEN
       NextLineIsData = .FALSE.
       RETURN
    END IF

    ! Read the next line.
    READ(iUnit,'(A)',IOSTAT = iError, END = 5) Line

    ! If error stop.          
    IF(iError.ne.0) THEN
       WRITE(ErrorNum,'(I5)') iError
       CALL Error('Error ocurred in NextLineIsData', StopProgram = .FALSE.)
       CALL Error('while reading from file ' // TRIM(ADJUSTL(FileName)) // &
            & '.',StopProgram = .FALSE.)
       CALL Error('Error: ' // ErrorNum)
    END IF

    TmpStr = TRIM(ADJUSTL(Line))

    ! First run through a few cases that look like headers but are not.
    IF(TmpStr(1:4).eq.'#SN#') THEN 
       NextLineIsData = .FALSE.
    ELSEIF(TmpStr(1:4).eq.'#DT#') THEN
       NextLineIsData = .FALSE.
    ELSEIF(TmpStr(1:4).eq.'#HD#') THEN
       ! Header Data 
       NextLineIsData = .TRUE.
    ELSEIF(TmpStr(1:4).eq.'#DF#') THEN
       NextLineIsData = .FALSE.
    ELSEIF(LEN_TRIM(TmpStr).ne.0) THEN ! This is not a blank line, check for header line.
       ! Header lines are anything that start with any of the CommentCharacters
       ! if provided, otherwise they start with #, *, c, C, or ! 
       IsHead = .FALSE.
       DO i1 = 1, LEN_TRIM(CmtChars)
          IF(TmpStr(1:1).eq.CmtChars(i1:i1)) THEN
             IsHead = .TRUE.
             EXIT
          END IF
       END DO

       ! If this line was not a header, backspace and return.
       IF(.not.IsHead) THEN          
          NextLineIsData = .TRUE.
       END IF
    END IF
5   CONTINUE
    BACKSPACE(iUnit)
    
  END FUNCTION NextLineIsData

  SUBROUTINE ReadNextWord(iUnit,Word,Length,IsHeader,EOR,EOF,CommentCharacters)
    INTEGER,INTENT(IN) :: iUnit
    CHARACTER,INTENT(IN) :: CommentCharacters
    INTEGER,INTENT(OUT) :: Length
    CHARACTER*(*),INTENT(OUT) :: Word
    LOGICAL,INTENT(OUT) :: EOR, EOF, IsHeader
    CHARACTER c
    CHARACTER(200) line
    INTEGER i1, clen

    EOR = .FALSE.
    EOF = .FALSE.
    IsHeader = .FALSE.
    Length = 0
    Word = ' '
    ! Read until we hit the first character. If we hit the end of the line, keep reading.
    DO
       READ(iUnit,'(A1)',ADVANCE = 'NO', EOR = 5, END = 10) c
       
       ! If c is not blank. Set Word(1:1) and exit
       clen = LEN_TRIM(c)
       IF(clen.gt.0) THEN
          Word(1:1) = c
          EXIT
       ENDIF          
    END DO

    ! Now check if this is a comment line.
    DO i1 = 1, LEN(CommentCharacters)
       IF(c.eq.CommentCharacters(i1:i1)) THEN
          IsHeader = .TRUE.
          READ(iUnit,*)
          RETURN
       END IF
    END DO

    Length = 1
    ! Now read until we hit a blank character.
    DO
       READ(iUnit,'(A)',ADVANCE = 'NO', EOR = 5, END = 10) c
       IF(LEN_TRIM(c).gt.0) THEN
          Word = TRIM(Word) // c
          Length = Length + 1
       ELSE
          EXIT
       END IF
    END DO
    RETURN
    
    ! End of record has been reached.
5   CONTINUE
    EOR = .TRUE.
    RETURN

    ! End of file has been reached.
10  CONTINUE
    EOF = .TRUE.
  END SUBROUTINE ReadNextWord
  
  ! Strips #HD# from the beginning of header data lines.
  SUBROUTINE StripHDLine(iUnit, IsHeader, EOR, EOF)
    INTEGER,INTENT(IN) :: iUnit
    INTEGER iError
    CHARACTER(4) HDString
    LOGICAL IsHeader, EOR, EOF

    EOR = .FALSE.
    EOF = .FALSE.
    IsHeader = .FALSE.

    READ(iUnit,'(A4)',ADVANCE = 'NO', EOR = 5, END = 10, IOSTAT = iError) HDString 
    CALL ReadError(iUnit = iUnit, Routine = 'StripHDLine', iError = iError)

    ! If HDString = #HD#, good, return
    IF(HDString.eq.'#HD#') THEN
       RETURN
    ELSEIF(TRIM(ADJUSTL(HDString(1:1))).eq.'#') THEN
       ! This is a header, but not data. 
       IsHeader = .TRUE.       
       RETURN
    ELSE
       ! This is not a header line or header data. Error
       CALL ReadError(iUnit = iUnit, Routine = 'StripHDLine')
    END IF

5   CONTINUE
    ! End of record found. Return
    EOR = .TRUE.
    RETURN

10  CONTINUE
    ! End of file found.
    EOF = .TRUE.
  END SUBROUTINE StripHDLine

  ! Handles read errors. Prints error number, file name, and calling routine.
  SUBROUTINE ReadError(FileName, iUnit, Routine, iError)
    CHARACTER*(*),INTENT(IN),OPTIONAL :: FileName, Routine
    INTEGER,INTENT(IN),OPTIONAL :: iUnit, iError
    CHARACTER(300) Messg(3), FlName, ErrorNum

    IF(PRESENT(FileName)) THEN
       FlName = FileName
    ELSEIF(PRESENT(iUnit)) THEN
       INQUIRE(UNIT = iUnit, NAME = FlName)
    ELSE
       CALL Error('ReadError called without either FileName or UnitNumber present.')
    END IF

    Messg(1) = 'READ ERROR:'

    IF(PRESENT(iError)) THEN
       WRITE(ErrorNum,*) iError
       Messg(1) = TRIM(Messg(1)) // ' ' // TRIM(ADJUSTL(ErrorNum))
    END IF

    Messg(2) = 'Read error occured while reading from file ' // TRIM(FlName) // '.'

    IF(PRESENT(Routine)) Messg(3) = 'This error occured inside the ' // TRIM(Routine) // ' routine.'

    CALL Error(Messg)

  END SUBROUTINE ReadError

  INTEGER Function NumberOfLines(FileName)
    CHARACTER*(*) FileName
    INTEGER iUnit, i1
    CHARACTER c
    CHARACTER(200) Line
    ! Find the number of lines that contain data
    CALL ReadHeaders(FileName,Line)
    CALL GetIOFileInfo(FileName,UnitNumber=iUnit)    
    BACKSPACE(iUnit)
    NumberOfLines = 0
    DO
       READ(iUnit,'(a)',END=10) c
       NumberOfLines = NumberOfLines + 1
    END DO
10  CONTINUE
    REWIND(iUnit)
  END Function NumberOfLines

END MODULE IOMod
