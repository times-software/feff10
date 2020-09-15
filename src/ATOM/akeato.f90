!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: akeato.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function akeato (i,j,k)
!     angular coefficient by the direct coulomb integral fk for orbitals
!     i and j
      implicit double precision (a-h,o-z)
      common/mulabk/afgk
      ! Dimensioning changed for high z elements. 
      dimension afgk(41,41,0:4) 
 
!     afgk angular coefficients by integrales fk and gk
!        coefficient of integral fk(i;j) is in  afgk(min,max)
!        and that of integral gk(i;j) is in  afgk(max,min)
!        max=max(i,j) min=min(i,j)
 
      if (i .le. j) then 
         akeato=afgk(i,j,k/2)
      else
         akeato=afgk(j,i,k/2)
      endif
      return
      end

      double precision function bkeato (i,j,k)
      implicit double precision (a-h,o-z)
      common/mulabk/afgk
      ! Dimensioning changed for high z elements. 
      dimension afgk(41,41,0:4)
!     angular coefficient at the exchange coulomb integral gk
 
      bkeato=0.0d00
      if (i .lt. j) then
         bkeato=afgk(j,i,k/2)
      elseif (i.gt.j) then
         bkeato=afgk(i,j,k/2)
      endif
      return
      end
