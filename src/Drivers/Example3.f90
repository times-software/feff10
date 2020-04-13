!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: Example3.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Example3
  USE IOMod
  INTEGER IA(30,30), N1, N2
  REAL RA(30,30)
  DOUBLE PRECISION DA(30,30)
  COMPLEX CA(30,30)
  COMPLEX*16 DCA(2,300)
  CHARACTER(100) SA(30,30), Headers(2)

  IA(:,:)  = 1
  RA(:,:)  = 2.0
  DA(:,:)  = 3.d0
  CA(:,:)  = (4.0,5.0)
  DCA(:,:) = (6.d0,7.d0)
  SA(:,:)  = 'S'
  Headers(:) = ' '

  Headers(1) = 'Using Write2D to write an integer array.'
  CALL Write2D('Example3.dat', DCA, Headers = Headers, FileType = 'PAD')
  CALL CloseFl('Example3.dat')
  N2 = 0
  CALL Read2D('Example3.dat',DCA, N1, N2, FileType = 'PAD')
  PRINT*, N1, N2
  PRINT*, DCA
!  Headers(1) = 'Using Write2D to write a real array.'
!  CALL Write2D('Example3.dat', RA, Headers = Headers)

!  Headers(1) = 'Using Write2D to write a double array.'
!  CALL Write2D('Example3.dat', DA,Headers = Headers)

!  Headers(1) = 'Using Write2D to write a complex array.'
!  CALL Write2D('Example3.dat', CA,Headers = Headers)

!  Headers(1) = 'Using Write2D to write a double complex array.'
!  CALL Write2D('Example3.dat', DCA,Headers = Headers)

!  Headers(1) = 'Using Write2D to write a string array.'
!  CALL Write2D('Example3.dat', SA,Headers = Headers)
  
END PROGRAM Example3
  
