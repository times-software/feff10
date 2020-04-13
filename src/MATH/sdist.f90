!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sdist.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function sdist (r0, r1)
!     find distance squared between cartesian points r0 and r1
!     single precision
      dimension r0(3), r1(3)
      sdist = 0
      do 10  i = 1, 3
         sdist = sdist + (r0(i) - r1(i))**2
   10 continue
      sdist = sqrt(sdist)
      return
      end
