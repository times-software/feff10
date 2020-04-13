!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: qmesh.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/07/03 01:25:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOP
!     !ROUTINE: QMesh
!     !INTERFACE:
      SUBROUTINE QMesh (Energy2,npos)
!     !USES:
      use program_control
      use qvectors
      use constants
      use eels_inp, energy=>ebeam
!     !INPUT/OUTPUT PARAMETERS:
!     Energy2 : Energy of the beam electron after scattering (in eV).
!     !DESCRIPTION:
!     The impuls transfer vector is calculated for one particular k0
!     and for a set of k' at energy Energy2 and corresponding to different
!     scattering directions, as expressed by ThXV and ThYV.
!     The Q-mesh will be used for evaluation and integration of the cross-section
!     over all initial and final states of the beam electron permitted by
!     collection and convergence semiangles.
!     
!     A relativistic correction is added to Q (=k0-k'-correction).

!     A little output is written to the master output file 6.
!     !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     EOP
      implicit none
      
!     INPUT (energy in eV):
      real*8,intent(in) ::  Energy2
      integer,intent(in) :: npos
!     LOCAL variables:
      real*8  Qrot(3), QLen, Q(3), KPrLen,beta,QLenClas
      integer IAtom, IPos
      real*8, external :: WaveLength
      real*8 K0Len              !KJ

      real*8 E(3,3)             !KJ Euler-rotation-matrix
      real*8 alfa1,alfa2,alfa3  ! 3 Euler angles KJ
      real*8 theta,phi,test
!	  real*8 q1(3),q2(3),q3(3),q4(3),q5(3) !testing
      integer natom,neqat,neqatdeb,i,j !KJ
!     KJ  set some variables :
      natom=1
      neqat=1
      neqatdeb=1


!     prepare matrix for rotation from laboratory frame to crystal frame :
      if(dabs(xivec(1)).lt.0.0001) then
	     if(xivec(2).gt.0.0001) then  !not that this matters too much, but it seems to make more sense to me
            alfa1=pi/dble(2)
		 else
		    alfa1=0
		 endif
      else
         alfa1=datan(xivec(2)/xivec(1))
      endif
      if(dabs(xivec(3)).lt.0.0001) then
         alfa2=pi/dble(2)
      else
         alfa2=datan(dsqrt(xivec(1)**2+xivec(2)**2)/xivec(3))
      endif
      alfa3=dble(0)             ! We don't care about an orientation of our mesh in the plane
                                ! perpendicular to the beam.  We only care that the beam direction is correct.
                                ! This may have to change later.
      call euler(alfa1,alfa2,alfa3,E)

      if(writeqmesh) then
         write(11,*) 'Euler angles :'
         write(11,'(3(f9.4,2x))') alfa1,alfa2,alfa3
         write(11,*) 'Rotation matrix for q-mesh :'
         do i=1,3
            write(11,'(3(f9.4,2x))') (E(i,j),j=1,3)
         enddo
         write(11,*)
         
      endif


      K0Len  = dble(2)*PI/WaveLength(Energy)
      KPrLen = dble(2)*PI/WaveLength(Energy2)

!     the vectors K0 and KPr are:
!     K0 = K0Len*(0, 0, -1)
!     KPr = KPrLen * (sin (ThX), sin(ThY), -SQRT(1 - sin^2(ThX) - sin^2(ThY)) )
!     Q = K0 - KPr 

      do IPos=1,NPos
!     The Q vectors are given in the basis of the observer:

         theta=(dsqrt(ThXV(IPos)**2+ThYV(IPos)**2))
!         theta=dasin(dsqrt(ThXV(IPos)**2+ThYV(IPos)**2))
         if(dabs(ThXV(IPos)).lt.0.000001) then
            if(ThYV(IPos).gt.0) then
               phi=pi/dble(2)
            else
               phi=-pi/dble(2)
            endif
         else
            phi=dabs(datan(ThYV(IPos)/ThXV(IPos)))
            if(ThYV(IPos).lt.0.and.ThXV(IPos).lt.0) then
        	   phi=pi+phi
            elseif(ThXV(IPos).lt.0) then
               phi=pi-phi
            elseif(ThYV(IPos).lt.0) then
               phi=-phi
            endif
         endif
!      open(59,file='angles.txt',form='formatted',position='append')
!	write(59,'(4e12.5)') thxv(ipos),thyv(ipos),theta,phi
!	close(59)

         Q(1) = - KPrLen * dsin(theta) * dcos(phi)
         Q(2) = - KPrLen * dsin(theta) * dsin(phi)
         Q(3) = KPrLen * dcos(theta) - K0Len
!		q1=q !test
!        qrot=q

!         Q(1) = - KPrLen * ThXV(IPos)
!         Q(2) = - KPrLen * ThYV(IPos)
!         Q(3) =   KPrLen * DSQRT(1-ThXV(IPos)**2-ThYV(IPos)**2)-K0Len
!		 q2=q !test
!      test=dsqrt((qrot(1)-q(1))**2+(qrot(2)-q(2))**2+(qrot(3)-q(3))**2)
!       if (test.gt.0.01 .and. .false.) then
!        open(59,file='qq.txt',position='append')
!        write(59,'(i5,10(e15.7,x))') Ipos,ThXV(IPos),THYv(Ipos),theta,phi
!        write(59,'(i5,10(e15.7,x))') IPos,qrot
!        write(59,'(i5,10(e15.7,x))') IPos,q
!        close(59)
!        endif
!      q=qrot

         QLenClas = DSQRT(Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3))
         beta=dsqrt((2+Energy/MeC2)/(2+Energy/MeC2+MeC2/Energy))
         if(RelatQ)                                                     &
!     1        Q(3)=Q(3)+(energy-energy2)/hbarc * beta
     &        Q(3)=Q(3)*(dble(1)- beta*beta)
!           write(17,*) beta*beta*q(3),(energy-energy2)/hbarc * beta
!     Note the sign : this is because of the stupid convention where the z-axis is antiparallel to the incident beam.

         QLen = DSQRT(Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3))
 !       q3=q !test

         DO IAtom = NEqAtDeb, NEqAt
!     For every equivalent atom of the unit cell, we put the Q vectors
!     in the local basis of the atom: Q -> Qrot 
!     IEqAtom numbers the atoms in the equivalency class (as defined in case.struct).

!     !KJ               CALL ProductMatVect(GeneralM(1,1,IAtom),Q,Qrot)
            CALL ProductMatVect(E,Q,Qrot)
!     !               Qrot=Q  !!KJ
!q4=qrot !test
            QV(1:3, IAtom, IPos) = Qrot(1:3)
            QLenV(IAtom, IPos)   = QLen
            QLenVClas(IAtom,IPos) = QLenClas
!        open(59,file='qq.txt',position='append')
!        write(59,'(5(3(e12.5,x),"-> "))') q1,q2,q3,q4
!        close(59)

         ENDDO
      enddo


      if(writeqmesh.and.verbosity.ge.2) then
!     print this only for the very first energy point (if present), otherwise file gets too large.
         write(11,'(/,a)') 'Output of subroutine QMesh:'
         write(11,'(a,i4,x,i4)') '# Qx      Qy        Qz  for atom '    &
     &        ,natom,neqatdeb
         do ipos=1,npos
            write(11,'(3(f10.5,x))') qv(1:3,neqatdeb,ipos)
         enddo
         writeqmesh=.false.
      endif

      RETURN
      END
