!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: angularmesh.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/07/03 01:25:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !ROUTINE: Angularmesh
! !INTERFACE:
      SUBROUTINE AngularMesh (ThetaXCenter,ThetaYCenter)
! !USES:
        use eels_inp
          use qvectors,only : thxv,thyv
          use work
          use constants,only : pi
! !INPUT/OUTPUT PARAMETERS:
!     ThetaXCenter,ThetaYCenter : position of the detector in mrad
! !DESCRIPTION:
!     Creates the angular mesh : the detector is placed at 'Center'.
!     Now a circle of radius (collection + convergence semiangle) around
!     this center is sampled with NPos points.
!     Precise sampling depends on modus selected - Uniform, Logarithmic
!     or 1-dimensional.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!     Hacked for FEFF October 2005 (Kevin Jorissen)
!EOP



!     nqf points at the distance ThPart of the center, 3*nqf at the distance
!     3*ThPart, ... , (2*nqr-1)*nqf at the distance (2*nqr-1)*ThPart.
!     The coordinates of the center are (ThetaX, ThetaY).
!     The coordinates of the point i are (ThXV(i), ThYV(i)).
          implicit none
!   INPUT : Position of the center of the aperture:
      real*8,intent(in) :: ThetaXCenter,ThetaYCenter
!   LOCAL VARIABLES :
      integer IRay, ITour, IndexPos, NPresentTour
      real*8 InterAngle,dxx
 

      if(qmodus.eq.'L'.or.qmodus.eq.'1') then
          dxx= dlog((aconv+acoll)/th0)/dble(nqr-1)
      endif
      IF (npos.gt.1) THEN
         IndexPos = 0
         DO IRay = 1, nqr
            NPresentTour = nqf*(2*IRay-1)
            if(qmodus.eq.'1') NPresentTour=1
            InterAngle = dble(2) * PI / DBLE (NPresentTour)
            DO ITour = 1, NPresentTour
               IndexPos = IndexPos + 1
               if (qmodus.eq.'L'.or.qmodus.eq.'1') then
                 if(IRay.eq.1) then
                    ThXV(IndexPos)=ThetaXCenter+ Th0*DCOS(InterAngle*ITour)/dble(2)
                    ThYV(IndexPos)=ThetaYCenter+ Th0*DSIN(InterAngle*ITour)/dble(2)
                 else
                    ThXV(IndexPos)=ThetaXCenter+ Th0*DCOS(InterAngle*ITour)* dexp(dxx*dble(IRay-2))*(dble(1)+dexp(dxx))/dble(2)
                    ThYV(IndexPos)=ThetaYCenter+ Th0*DSIN(InterAngle*ITour)* dexp(dxx*dble(IRay-2))*(dble(1)+dexp(dxx))/dble(2)
                 endif

               else   
                 ThXV(IndexPos) = ThetaXCenter + DCOS(InterAngle*ITour) * (Thpart*(2*IRay-1))
                 ThYV(IndexPos) = ThetaYCenter + DSIN(InterAngle*ITour) * (Thpart*(2*IRay-1))
               endif
            ENDDO
         ENDDO      
      ELSE
         ThXV(1) = ThetaXCenter
         ThYV(1) = ThetaYCenter 
      ENDIF
      RETURN
      END


