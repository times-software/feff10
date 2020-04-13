!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fdmocc.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function fdmocc (i,j)
!     product of the occupation numbers of the orbitals i and j

      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
 
      if (j.eq.i) then
         fdmocc=xnel(i)*(xnel(j)-1)
         a=2* abs(kap(i))
         fdmocc=fdmocc*a/(a-1.0)
      else
         fdmocc=xnel(i)*xnel(j)
      endif
      return
      end
