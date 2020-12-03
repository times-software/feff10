!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fdmocc.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function fdmocc (i,j)
!     product of the occupation numbers of the orbitals i and j

! Josh Kas - Changed array dimensions from 30 to 41 for high Z elements
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
 
      if (j.eq.i) then
         fdmocc=xnel(i)*(xnel(j)-1)
         a=2* abs(kap(i))
         fdmocc=fdmocc*a/(a-1.0)
      else
         fdmocc=xnel(i)*xnel(j)
      endif
      return
      end
