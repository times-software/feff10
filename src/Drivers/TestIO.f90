!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: TestIO.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM TestIO
  USE IOMod

  INTEGER iUnit, I, IA(100), i1, NumElements
  REAL R, RA(100)
  DOUBLE PRECISION D, DA(100)
  COMPLEX C, CA(100)
  COMPLEX*16 DC, DCA(100)
  CHARACTER(20) S, SA(100), ColumnLabels(8)
  CHARACTER(100) Headers(100), FileName

  ! This program illustrates a few simple examples of how to 
  ! use IOMod

  ! Define some data, including arrays.
  I = 1
  IA(:) = 1
  R = 2.1
  RA(:) = 2.1
  D = 3.2
  DA(:) = 3.2
  C = (2.1,2.1)
  CA(:) = (2.1,2.1)
  DC = (3.2d0,3.2d0)
  DCA(:) = (3.2d0,3.2d0)
  S = "Cu"
  SA(:) = "Cu"

  ! Initialize an array that will hold headers, and define some column labels.
  Headers(:) = ' '
  ColumnLabels(1) = 'ie'
  ColumnLabels(2) = 'omega'
  ColumnLabels(3) = 'e'
  ColumnLabels(4) = 'Re[p]'
  ColumnLabels(5) = 'Im[p]'
  ColumnLabels(6) = 'Re[em]'
  ColumnLabels(7) = 'Im[em]'
  ColumnLabels(8) = 'atom'

  ! Try writing just some headers first, i.e. CALL WriteData(Filename, Headers(:))
  Headers(1) = 'Writing headers is done by using the optional Headers argument as follows:'
  Headers(2) = '-'
  Headers(3) = "-   CALL WriteData('TestIOData/WriteData.dat', Headers = HeaderArray)"
  Headers(4) = '-'
  Headers(5) = 'Note that null or blank header strings are not written to the file.'
  
  CALL WriteData('TestIO.dat', Headers = Headers)

  ! Write the following variables to file WriteData.dat in a row: I, R, D, C, DC, S
  ! Note that a section number line (beginning with #SN#) has been written, followed by a data
  ! type line (beginning with #DT#), then the data on the final line. 
  CALL WriteData('TestIOData/WriteData.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S)

  ! Write a few more headers with no data. Note that the section number is not updated until the next write.
  Headers(:) = ''
  Headers(1) = 'If you want to write headers alone, do'
  Headers(2) = "   CALL WriteData('TestIOData/data1.dat', Headers = Headers)"
  Headers(3) = 'Note that this does not increment the section #.'
  CALL WriteData('TestIOData/WriteData.dat', Headers = Headers)

  ! Now lets call WriteData with the same arguments as the first time. This will increment
  ! The Section number, and print a data type line which begins with #DT#.
  I = 2
  CALL WriteData('TestIOData/WriteData.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S)

  ! Now if we write again with the same type and number of arguments, i.e. I, R, D, C, DC, S
  ! the line is printed directly below the previous line with no new section number, etc.
  CALL WriteData('TestIOData/WriteData.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S)

  ! Now let's add some columnlabels and write with a different set of arguments. This will increment
  ! the section number, and print column labels directly above the data.
  I = 3
  CALL WriteData('TestIOData/WriteData.dat', DComplex1 = DC, Int2 = I, Real3 = R, Double4 = D, Complex5 = C, &
       & DComplex6 = DC, ColumnLabels = ColumnLabels)

  ! We can also write data in the header, which will be marked with #HD# using the WriteDataInHeader optional argument.
  ! This is usefull if you want to pass various data by file but only forsee wanting to plot one of the arrays of data.
  CALL WriteData('TestIOData/WriteData.dat', DComplex1 = DC, Int2 = I, Real3 = R, Double4 = D, Complex5 = C, &
       & DComplex6 = DC, ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE.)

  Headers(:) = ''
  Headers(1) = 'Here we will write an array of data to file by looping over WriteData.'
  Headers(2) = 'Note that headers and columnlabels will only be written out on the first call.'
  ! Writing an array of data can be performed in two ways.
  !
  ! 1) Using WriteData inside a loop. Note that we can write with the Headers and ColumnLabels optional arguments,
  !    and they will only be written on the first call, before any of the data for this section. A new section occurs when
  !    only when the data types or number of arguments passed to the subroutine change. 
  DO i1 = 1, 20
     CALL WriteData('TestIOData/WriteData.dat', Int1 = i1, Real2 = RA(i1), Headers = Headers, ColumnLabels = ColumnLabels)
  END DO

  Headers(1) = 'Here we will write an array of data to file by calling WriteArrayData.'
  Headers(2) = 'Note that this is much faster for large arrays.'
  ! 2) We can also use the subroutine WriteDataArray which is very similar to WriteData, but takes arrays instead of
  !    scalars as arguments. This is much faster for large arrays.
  CALL WriteArrayData('TestIOData/WriteData.dat', Int1 = IA(1:20), Real2 = RA(1:20), Headers = Headers, &
       & ColumnLabels = ColumnLabels, ForceNewSection = .TRUE.)

  ! Now we can try reading some data using ReadData
  ! First we have to close the file using CloseFl
  CALL CloseFl('TestIOData/WriteData.dat')

  ! Read the first section.
  I = 0; R = 0.0; D = 0.d0; C = 0.0; DC = 0.d0; S = ' '
  CALL ReadData('TestIOData/WriteData.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S)
  PRINT*, 'Section 1'
  PRINT*, I, R, D, C, DC, S

  ! Now read from section 3 using the optional SectionNumber argument.
  I = 0; R = 0.0; D = 0.d0; C = 0.0; DC = 0.d0; S = ' '
  CALL ReadData('TestIOData/WriteData.dat', DComplex1 = DC, Int2 = I, Real3 = R, Double4 = D, Complex5 = C, &
       & DComplex6 = DC, SectionNumber = 3)
  PRINT*, 'Section 3'
  PRINT*, I, R, D, C, DC, S
  
  ! Now we can read arrays by looping over ReadData.
  PRINT*, 'Section 4'
  DO i1 = 1, 10
     CALL ReadData('TestIOData/WriteData.dat', Int1 = IA(i1), Real2 = RA(i1), SectionNumber = 4)
     PRINT*, IA(i1), RA(i1)
  END DO

  IA(:) = 0
  RA(:) = 0.0
  ! Or by calling ReadArrayData. Note that in this case we want to read section 5, and we already
  ! read section 4, so we don't need to specify.
  CALL ReadArrayData('TestIOData/WriteData.dat', Int1 = IA, Real2 = RA, &
       & NumElements = NumElements)
  PRINT*, 'Section 5: With ReadArrayData'
  DO i1 = 1, 20
     PRINT*, IA(i1), RA(i1)
  END DO
  PRINT*, NumElements

  ! The following examples illustrate the most common ways you would use these routines.
  ! 
  ! Ex. 1: Writing an array to a file from the inside of a do loop.
  ! Declare some headers. Initialize the whole array to blank, then blank arrays
  ! will not be written. Alternatively you can pass Headers(1:n) instead of the 
  ! whole array.
  Headers(:) = ' '
  Headers(1) = 'Here we illustrate how to use WriteData to write an array'
  Headers(2) = 'to a file from inside a loop.'
  
  ! If we wanted we could write some data in the headers by calling WriteData before
  ! beginning the loop. Let's say we want to write E_Fermi in the header. We could 
  ! do this as follows:
  !
  ! CALL WriteData('IOExample1.dat', Real1 = EFermi, Headers = Headers, 
  ! ColumnLabels = ('E_Fermi'), WriteDataInHeaders = .TRUE.)
  !
  ! Now write the data arrays defined above to the file IOExamples.dat.
  ! Here we use six data variables. But you can use any number from one to ten. 
  ! The keywords are used to tell the routine what column you want the data to go in.
  ! Note that the headers, section #, data type line, and column labels are only written
  ! in the first call of the section. New sections occur when the data types or number of 
  ! data are changed. New sections can also be forced with the optional logical argument
  ! NewSection.
  DO i1 = 1, 20
     CALL WriteData('IOExample1.dat', Int1 = IA(i1), Real2 = RA(i1), Double3 = DA(i1), &
          & Complex4 = CA(i1), DComplex5 = DCA(i1), String6 = SA(i1), &
          & ColumnLabels = ColumnLabels, Headers = Headers)
  END DO

END PROGRAM TestIO
