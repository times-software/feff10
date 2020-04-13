!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: frgrid.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     RADIAL GRID ROUTINES
!=======================================================================
subroutine setri(dx, x0, ilast, ri)
  use DimsMod, only: nrptx
  implicit none

  double precision ri(nrptx), dx, x0
  integer ilast, i

  do i = 1, ilast
     ri(i) = exp(-x0+(i-1)*dx)
  end do

  !      end subroutine setri
end subroutine setri

function getiat(x0, dx, rref)
  implicit none
  integer getiat
  double precision rref, x0, dx

  getiat = (log(rref) + x0) / dx + 1

  return
  !      end function getiat
end function getiat
!=======================================================================
!     END RADIAL GRID ROUTINES
!=======================================================================
