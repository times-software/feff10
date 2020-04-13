!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: writeangulardependence1.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOP
!     !ROUTINE: WriteAngularDependence1
!     !INTERFACE:
      subroutine WriteAngularDependence1(qvs)
!     !USES:
      use eels_inp, only : npos, k0len, Energy
      use qvectors,only : WeightV
      use spectra_hyperfine,only : sdlm
      use constants, only : hbarc,mec2
      use energygrid,only : jemin,ene
      use program_control, only : relatQ,headers
!     !INPUT/OUTPUT PARAMETERS:
!     qvs   : set of Q-vectors in the laboratory frame, in spherical coordinates

!     !DESCRIPTION:
!     The angular differential partial cross sections are taken and the integration weights
!     for detector aperture and beam convergence are removed, so that a pure spectrum as a 
!     function of scattering vector is obtained.  These are translated into a scattering angle
!     (we have removed the collection*convergence convolution!).
!     Partial spectra as a function of scattering angle are written to file 60.

!     Only the first (j=l+1/2) edge is used.

!     !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     EOP

      implicit none

      real*8 Q(npos),pi(npos),sigma(npos),tot(npos),totdip(npos),       &
     &     sigmadip(npos),QVs(3,npos),theta(npos),mono(npos),quad(npos),&
     &     octo(npos)
      real*8 beta,k0local
      integer k

      do k=1,npos
         pi(k)=sdlm(3,k,1) /weightv(k)
         sigmadip(k)=(sdlm(2,k,1)+sdlm(4,k,1))/weightv(k)
         sigma(k)=(sdlm(2,k,1)+sdlm(4,k,1)+sdlm(1,k,1))/weightv(k)
         quad(k)=(sdlm(5,k,1)+sdlm(6,k,1)+sdlm(7,k,1)+sdlm(8,k,1)       &
     &        +sdlm(9,k,1))/weightv(k)
         octo(k)=sdlm(10,k,1)/weightv(k)
         mono(k)=sigma(k)-sigmadip(k)
         totdip(k)=pi(k)+sigmadip(k)
         tot(k)=pi(k)+sigma(k)+quad(k)
         Q(k)=QVs(1,k)
!     if(relatQ) then
!     beta=dsqrt((2+Energy/MeC2)/(2+Energy/MeC2+MeC2/Energy))
!     k0local=k0len-beta*(energy-ene(jemin))/hbarc
!     Note the sign : this is because of the stupid convention where the z-axis is antiparallel to the incident beam.
!     else
         k0local=k0len
!     endif
         theta(k)=Q(k)*dsin(QVs(2,k))*dble(-1000)/                      &
     &        dsqrt((k0len**2+Q(k)**2-2.0*k0len*Q(k)*dcos(QVs(2,k))))
!     This is a small angle approximation for the scattering angle theta(k) = angle < k0, k' >

!     From Q=k-k' (all vectors), we get cos (theta) = (Q^2 - k^2 -k'^2)/(-2kk')   (from now on, all scalars)
!     Now approximate cos theta = 1 - theta^2 / 2
!     And express k' as the dsqrt in the equation.
!     This gives theta^2 = ( (k-k')^2 - Q^2 )/ (-kk')    ( all scalars)
!     The first part equals k sin (ThetaQ).
!     This yields the equation used above.        

!     If there are relativistic corrections, Q(z) has been shortened by a vector parallel to k :
!     Then Q = k-k'-cr
!     Since cr || k , we can still use the above equation if we replace the length of k by the length of k minus the relativistic correction cr.
!     This is the purpose of the variable k0local.



      enddo


      if (headers) write(60,'(a)') '## theta  pi(1_0)  sigma(11 1-1 00) &
     &     total(**)  sigmadipole(11 1-1) totaldipole(1*) monopole(00) quadrupole(2*)'
      do k=1,npos
         write(60,10) theta(k),pi(k),sigma(k),tot(k),sigmadip(k),totdip(k),&
     &        mono(k),quad(k),octo(k)
      enddo
 10   format(9F17.11)

      return
      end

