!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fegrid.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/04/07 19:35:57 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     ENERGY GRID ROUTINES
!=======================================================================            
! JK 4/2010
SUBROUTINE SetEGrid(emax, ne, em)
  REAL(8), INTENT(IN) :: emax
  INTEGER, INTENT(IN) :: ne
  COMPLEX*16, INTENT(OUT) :: em(ne)

  REAL(8) dx
  INTEGER ie

  ! This will make an exponential energy grid from E_Fermi to 
  ! E_Fermi + i*emax, i.e.
  ! E_i = exp[(N-1)dx] - 1; dx = ln[emax + 1]/(N-1)
  dx = LOG(emax+1.d0)/DBLE(ne-1)
  em(1) = 0.d0
  DO ie = 1, ne
     em(ie) = (EXP(DBLE(ne-ie)*dx) - 1)*(0.d0,1.d0)
  END DO
  RETURN
END SUBROUTINE SetEGrid
  
  
subroutine setegi(emin, emax, eimax, ermin, ner, nei, em, ne)

  use DimsMod, only: nex
  use constants
  implicit none
  !     makes a grid in complex e-plane like below
  !     en_grid is comlex energy in hartrees
  !
  !     ...................
  !     .                 .
  !     1                 n
  !     ---------------------------->Re (Im = 0)
  double precision emin, emax, eimax, ermin
  complex*16 em(nex)
  integer ner, nei, ne

  !     Local
  double precision s_real, csreal, csimag
  complex*16 de, dei, delta, em2(nex)
  integer i, neg

  !     never do calculations on real axis.
  if (ermin .le. 0.0) ermin = 0.05d0
  s_real = emax - emin

  !     differential energy step along real axis
  de  = (emax - emin) / (ner - 1)
  !     differential energy step along imaginary axis
  dei = coni*(eimax - ermin) / (nei - 1)

  neg = 1

  !=================================================================
  !     Prepare Grids
  !=================================================================      
  em2(neg) = emax + coni*ermin
  csreal  = 0.0d0
  csimag  = ermin
  delta   = dei

  do neg = 2, nex**2
     !        neg = neg + 1

     if (dreal(em2(neg-1)) .lt. emin) then
        !         going down along the imag axis
        delta = -dei
        if (dimag(em2(neg-1)) .le. ermin) then
           goto 60
        end if
     else
        if (abs(csimag) .ge. eimax) then
           !           change direction, steps along the real axis
           delta  = -de
           csimag = 0.0d0
        end if
     end if

     csimag = csimag + abs(dimag(delta))
     csreal = csreal + abs(dreal(delta))

     em2(neg) = em2(neg-1) + delta

  end do

60 continue

  if (dimag(em2(neg-1)) .le. 0.0) then
     neg = neg - 2
  else
     neg = neg - 1
  end if

  !     reverse the grid      
  do i = 1, neg
     em(i) = em2(neg+1-i)
  end do

  if (neg .gt. nex) then
     print *, 'Error: too many energy points:', neg
     stop
  end if

  ne = neg

  return

  !      end subroutine setegi
end subroutine setegi
!=======================================================================
!     END ENERGY GRID ROUTINES
!=======================================================================
