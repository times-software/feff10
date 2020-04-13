!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_euler.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine mdff_euler(a,b,g,E)

        implicit none
! IN/OUTPUT
        real*8,intent(in) :: a,b,g
        real*8,intent(out) :: E(3,3)
! LOCALS
      integer i,j
        real*8 det


      E(1,1)=dcos(a)*dcos(b)*dcos(g)-dsin(a)*dsin(g)
        E(2,1)=dsin(a)*dcos(b)*dcos(g)+dcos(a)*dsin(g)
        E(1,2)=-dcos(a)*dcos(b)*dsin(g)-dsin(a)*dcos(g)
        E(2,2)=-dsin(a)*dcos(b)*dsin(g)+dcos(a)*dcos(g)
        E(1,3)=dcos(a)*dsin(b)
        E(2,3)=dsin(a)*dsin(b)
        E(3,3)=dcos(b)
        E(3,1)=-dsin(b)*dcos(g)
        E(3,2)=dsin(b)*dsin(g)


       ! check :
        call mdff_determinant(E,det)
        if(dabs(det-1).gt.0.0001) then
         write(11,*) 'Warning : nasty Euler matrix.'
        do i=1,3
         write(11,*) (E(i,j),j=1,3)
        enddo
        write(11,*) 'det E = ',det
        endif


        return
        end






!BOP
! !ROUTINE: Determinant 
! !INTERFACE:
      SUBROUTINE mdff_Determinant(M, Det)
! !INPUT/OUTPUT PARAMETERS:
!     M  :   3\*3 real input matrix
!     Det :  determinant of M
! !DESCRIPTION:
!     calculate the determinant of a 3\*3 real matrix.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP      

      implicit none
          real*8,intent(in) :: M(3,3)
          real*8,intent(out) :: Det
            
      Det = M(1, 1) * M(2, 2) * M(3, 3)                                 &
     &      + M(1, 2) * M(2, 3) * M(3, 1)                               &
     &      + M(2, 1) * M(3, 2) * M(1, 3)                               &
     &      - M(3, 1) * M(2, 2) * M(1, 3)                               &
     &      - M(2, 1) * M(1, 2) * M(3, 3)                               &
     &      - M(1, 1) * M(3, 2) * M(2, 3)
      
      RETURN
      END

      
