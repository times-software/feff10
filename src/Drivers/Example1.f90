!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: Example1.f90,v $:
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

  ! The following example illustrate the most common ways you would use these routines.
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
     CALL WriteData('Example1a.dat', Int1 = IA(i1), Real2 = RA(i1), Double3 = DA(i1), &
          & Complex4 = CA(i1), DComplex5 = DCA(i1), &
          & ColumnLabels = ColumnLabels, Headers = Headers)
  END DO

  ! Close the file.
  CALL CloseFl('Example1a.dat')
  ! Reinitialize the arrays.
  IA(:) = 0; RA(:) = 0.0; DA(:) = 0.d0; CA(:) = (0.0,0.0); DCA(:) = (0.d0,0.d0)

  Headers(:) = ' '
  Headers(1) = 'Here we read the data out of Example1a.dat and wrote it into this file.'
  ! Now read the file back into the array within a do loop using ReadData. We will write
  ! it again to Example1b.dat to be sure that we read it in correctly.
  DO i1 = 1, 20
     CALL ReadData('Example1a.dat', Int1 = IA(i1), Real2 = RA(i1), Double3 = DA(i1), &
          & Complex4 = CA(i1), DComplex5 = DCA(i1))

     CALL WriteData('Example1b.dat', Int1 = IA(i1), Real2 = RA(i1), Double3 = DA(i1), &
          & Complex4 = CA(i1), DComplex5 = DCA(i1), &
          & ColumnLabels = ColumnLabels, Headers = Headers)
  END DO
END PROGRAM TestIO
