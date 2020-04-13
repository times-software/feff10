!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: calculateweights.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/07/03 01:25:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOP
!     !ROUTINE: CalculateWeights
!     !INTERFACE:
      SUBROUTINE CalculateWeights
!     !USES:
      use work
      use qvectors
      use eels_inp,only : qmodus,nqr=>nqr,nqf=>nqf,aconv,acoll,th0,thetax,thetay
      use constants,only: pi
!     !DESCRIPTION:
!     Set up the mesh of (k,k') - pairs for which the double differential scattering cross-section
!     will be calculated.
!     The mesh is determined by collection and convergence semiangle, and by the type of mesh
!     chosen (Uniform, Logarithmic, 1-dimensional) and its specified size (nqr, nqf, NPos).
!     
!     A note : since we do not change the microscope and detector semiangles while we measure
!     a full edge, also the definition of the 'pixels' in our calculation, i.e. the *angular* mesh
!     used for sampling the Q vectors, should not change.
!     Since there is a relation connecting energy loss, the magnitude of impuls transfer, and the
!     scattering angle, keeping the latter fixed does imply that the vector Q and its length change
!     when energy changes.
!     Therefore, we run calculweight and newthxthy only once at the beginning of the calculation,
!     but we run qmesh for every energy in the edge.
!     !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     Hacked for FEFF October 2005 (Kevin Jorissen)
!     EOP



      implicit none
!     LOCAL variables
      integer IRay, ITour, IndexPos, NPresentTour
      real*8 ConvolValue, Lfactor, dxx
      real*8 Theta,sa,ca,p


      
      call make_Qvecs1(NPos)
      if(qmodus.eq.'L'.or.qmodus.eq.'1') then
         dxx= dlog((aconv+acoll)/th0)/dble(nqr-1)
      endif
      CALL AngularMesh (dble(0),dble(0))
      IF (npos.gt.1) THEN
         
         sa=acoll               ! for convenience, we use these abbreviations :
         ca=aconv
         IndexPos = 0
         DO IRay = 1, nqr
            if(qmodus.eq.'1') then
               NPresentTour=1
            else
               NPresentTour = nqf*(2*IRay-1)
            endif
            Theta = ThXV(IndexPos+NPresentTour) ! The last position in a ring has cos(2 pi) = 1; and here ThXCenter = 0
            
            IF     (Theta.LE.dabs(sa-ca)) THEN
               if (ca.gt.dble(0.000001).and.sa.gt.dble(0.000001)) then
                  ConvolValue = pi* min(ca,sa)**2  / (pi*ca**2)
               else
                  ConvolValue = dble(1)
               endif
            ELSEIF (Theta.GE.(sa+ca)) THEN
               ConvolValue = dble(0)
            ELSE                ! |a-b| < Theta < a+b
               p=(Theta**2+ca**2-sa**2)/(dble(2)*Theta)
! Modified by FDV
! Split line to compile with Solaris Studio
               convolvalue=pi/dble(2)*(ca**2+sa**2)- p*dsqrt(ca**2-p**2) - &
               (Theta-p)*dsqrt(sa**2-(Theta-p)**2)- sa**2*dasin((Theta-p)/sa) -ca**2*dasin(p/ca)
               convolvalue = convolvalue / (pi*ca**2) ! if ca=0, we wouldn't be here !
            ENDIF

            DO ITour = 1, NPresentTour
               IndexPos = IndexPos + 1
!     Weight is the surface around a point multiplied with the value at the
!     point of the convolution of acoll and Convergence distribution.

 !KJ calibration 12/2010 WeightV(IndexPos) = (1.D6*ThPart**2)/ DBLE(NPresentTour) &

               WeightV(IndexPos) = (ThPart**2)/ DBLE(NPresentTour) * PI * dble(4) * (2 * IRay - 1) * ConvolValue
!     Now correct the weights for nonuniform Q-meshes :
               if(qmodus.eq.'L'.or.qmodus.eq.'1') then
                  if (Iray.eq.1) then
                     Lfactor=(dble(nqr)*Th0/(sa+ca))**2*nqf/NPresentTour
                  else
                     Lfactor=(dble(nqr)*Th0*dexp(dxx*(iray-2)) /(sa+ca))**2*(dexp(2*dxx)-dble(1))* nqf/NPresentTour
                  endif
               else
                  Lfactor=dble(1)
               endif
               WeightV(IndexPos)=WeightV(IndexPos)*Lfactor

            ENDDO
         ENDDO
      ELSE
!     Only one spectrum
         WeightV(1) = dble(1)
!         write(11,'(4(f10.5,2x))') ThXV(1), ThYV(1),WeightV(1),dble(1)
      ENDIF

!     for the rest of the calculation:
      call AngularMesh(thetax,thetay) ! makes thx and thy for given alfa, beta
      RETURN
      END
      

