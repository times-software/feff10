!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: Example4.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Example4
  USE IOMod
  INTEGER I
  REAL R
  DOUBLE PRECISION D
  COMPLEX C
  COMPLEX*16 DC
  CHARACTER(100) S, Headers(2)

  I  = 123456
  R  = 1.234567890
  D  = 1.23456789d0
  C  = (1.23456789,1.23456789)
  DC = (1.23456789d0,1.23456789d0)
  S  = 'S'
  Headers(:) = ' '

  PRINT*, I, R, D, C, DC, S

  Headers(1) = 'Using WriteData to write PAD to a file.'
  CALL WriteData('Example4.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S, Headers = Headers,FileType = 'PAD')
  CALL CloseFl('Example4.dat')
  I = 0
  R = 0.0
  D = 0.d0
  C = 0.0
  DC = 0.d0
  S = ' '
  ! Now we can use ReadData to read the PAD from the file.
  CALL ReadData('Example4.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S, Headers = Headers,FileType = 'PAD')
  PRINT*, I, R, D, C, DC, S
END PROGRAM Example4
  
