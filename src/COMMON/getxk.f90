!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getxk.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function getxk (e)
      implicit double precision (a-h, o-z)

!     Make xk from energy(in Hartrees) as
!          k =  sqrt(2*e)  for e > 0  (above the edge)
!          k = -sqrt(-2*e)  for e < 0  (below the edge)
      getxk = sqrt(abs(2*e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
