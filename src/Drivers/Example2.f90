!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: Example2.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Example
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
  Headers(1) = 'Here we illustrate how to use WriteArrayData to write an array'
  Headers(2) = 'to a file.'
  
  ! If we wanted we could write some data in the headers by calling WriteArrayData 
  ! before beginning the loop. Let's say we want to write E_Fermi in the header. 
  ! We could do this as follows:
  !
  ! CALL WriteData('IOExample1.dat', Real1 = EFermi, Headers = Headers, 
  ! ColumnLabels = ('E_Fermi'), WriteDataInHeaders = .TRUE.)
  !
  ! This would cause a new section to be defined for the array of data.
  !
  ! Now write the some of the data arrays defined above to the file Examplew.dat.
  ! Here we use three data variables. But you can use any number from one to ten. 
  ! The keywords are used to tell the routine what column you want the data to go in.
  ! Note that the headers, section #, data type line, and column labels are only written
  ! in the first call of the section. New sections occur when the data types or number of 
  ! data are changed. New sections can also be forced with the optional logical argument
  ! NewSection.
  CALL WriteArrayData('Example2a.dat', Real1 = RA, Complex2 = CA, &
       & DComplex3 = DCA, ColumnLabels = ColumnLabels, Headers = Headers)

  ! Now we could write some more headers if we want.
  Headers(:) = ' '
  Headers(1) = 'This file contains the arrays that define this example.'
  CALL WriteData('Example2a.dat', Headers = Headers)

  ! Close the file.
  CALL CloseFl('Example2a.dat')

  ! Now read the file and write the data to Example2b.dat
  ! Using the NumElements optional argument is usefull for getting
  ! the number of elements in the array.
  CALL ReadArrayData('Example2a.dat', Real1 = RA, Complex2 = CA, &
       & DComplex3 = DCA, NumElements = i1)

  PRINT*, 'The number of elements is ', i1
  Headers(:) = ' '
  Headers(1) = 'This is a copy of Example2a.dat'
  CALL WriteArrayData('Example2b.dat', Real1 = RA, Complex2 = CA, &
       & DComplex3 = DCA, ColumnLabels = ColumnLabels, Headers = Headers)
END PROGRAM Example
