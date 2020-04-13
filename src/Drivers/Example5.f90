!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: Example5.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Example5
  USE IOMod
  INTEGER I(100000)
  REAL R(100000)
  DOUBLE PRECISION D(100000)
  COMPLEX C(100000)
  COMPLEX*16 DC(100000)
  CHARACTER(100) S(100000), Headers(2)

  I(:)  = 123456
  R(:)  = 1.234567890
  D(:)  = 1.23456789d0
  C(:)  = (1.23456789,1.23456789)
  DC(:) = (1.23456789d0,1.23456789d0)
  S(:)  = 'S'
  Headers(:) = ' '

  Headers(1) = 'Using WriteData to write PAD to a file.'
  PRINT*, 'PAD'
  CALL WriteArrayData('Example5a.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S, Headers = Headers,FileType = 'PAD')
  CALL WriteArrayData('Example5b.dat', Int1 = I, Real2 = R, Double3 = D, Complex4 = C, &
       & DComplex5 = DC, String6 = S, Headers = Headers,FileType = 'TXT')

  I = 0
  R = 0.0
  D = 0.d0
  C = 0.0
  DC = 0.d0
  S = ' '



END PROGRAM Example5
