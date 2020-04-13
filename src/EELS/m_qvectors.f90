!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_qvectors.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module qvectors
!     This module contains the impuls transfer vectors and the final wave vectors.

!     Weights for integration over Q-vectors
          real*8, allocatable ::  WeightV(:)
!     Carthesian coordinates of final wave vectors in the 'detector plane' (careful : this is in a 
!     plane made by convoluting the beam cross section and the detector aperture)
      real*8, allocatable ::  ThXV(:), ThYV(:)
!     The Q-mesh for a particular energy
          real*8, allocatable ::  QV(:,:,:)
!     The length of each Q-vector, both relativistically and 'classically'
      real*8, allocatable ::  QLenV(:,:),QLenVClas(:,:)

          CONTAINS
           subroutine make_Qvecs1(nposmax)
!     allocate the final wave vectors
             integer nposmax
             allocate(WeightV(NPOSMAX), ThXV(NPOSMAX), ThYV(NPOSMAX))
                 weightv(:)=dble(0);thxv(:)=dble(0);thyv(:)=dble(0)
           end subroutine make_Qvecs1
           subroutine make_Qvecs2(nposmax,ndif)
!     allocate the Q-mesh
             integer nposmax,ndif
             allocate(QLenV(ndif,nposmax),qv(3,ndif,nposmax),           &
     &               QLenVClas(ndif,nposmax))	
                 qlenv(:,:)=dble(0);qv(:,:,:)=dble(0)
                 qlenvclas=dble(0)
           end subroutine make_Qvecs2
       subroutine destroy_qvecs1
!     deallocate final wave vectors and weights
             deallocate(weightv,thxv,thyv)
           end subroutine destroy_qvecs1
       subroutine destroy_qvecs2
!     deallocate Q-mesh
             deallocate(qlenv,qv,qlenvclas)
           end subroutine destroy_qvecs2
          end module qvectors
