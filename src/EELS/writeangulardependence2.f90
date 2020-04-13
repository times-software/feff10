!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: writeangulardependence2.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/12/16 18:30:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOP
!     !ROUTINE: WriteAngularDependence2
!     !IntERFACE:
      subroutine WriteAngularDependence2
!     !USES:
      use eels_inp
      use work
      use qvectors
      use program_control
      use constants,only : pi,hbarc_eV
      use spectrum,only : s,ne
!     !DESCRIPTION:
!     The angular differential partial cross sections are integrated up to a variable
!     collection angle.  The resulting partial intensities are written, as a function
!     of collection semiangle, to file 59.

!     !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     Adapted for FEFF October 2005 (Kevin Jorissen)
!     EOP

      implicit none

!     LOCAL VARIABLES
      integer ncol,j,icol,imagic
      real*8 betastep,beta,dx,qfac,origin
      real*8, allocatable :: WeightVold(:),ints(:,:),sp2(:), collection(:),qve(:,:,:)
      real*8 aconvold,acollold
      integer nqrold,nqfold,nposold

!     Find energy index corresponding to energy emagic
      imagic=0
      origin=dble(-5)
      do j=1,ne
         if(s(j,2).gt.0.000001.and.origin.lt.0.0) origin=s(j,1)
         if (emagic.gt.(s(j,1)-origin).and.origin.ge.0.0) then
            imagic=j
            exit
         endif
      enddo
      if (j.gt.ne) then    
         write(11,'(/,a)') 'Emagic out of range, reset to last mesh point'
         imagic=ne
      else
         write(11,'(/,a,f12.5)') 'Magic angle at ',s(imagic,1)
      endif
      
      


      write(11,'(/,a)') 'Output from subroutine PlotAngularDependence:'
      write(11,'(/,a,f12.4,a)') 'Plotting at energy ',s(imagic,1),' eV.' 

!     We first save old data :
      nqrold=nqr
      nqfold=nqf
      aconvold=aconv
      acollold=acoll
      nposold=npos
      allocate(WeightVold(npos),qve(3,1,npos))
      WeightVold(:)=WeightV(:)


!     There are only beta/(alfa+beta) * nqr collection angles !!!
      if (qmodus.eq.'L'.or.qmodus.eq.'l'.or.qmodus.eq.'1') then
         dx=dlog((acoll+aconv)/th0)/dble(nqr-1)
         ncol=1+int(log(acoll/th0)/dx)
      else
         betastep=(aconv+acoll)/dble(nqr)
         ncol=int(dble(acoll/(acoll+aconv)) *dble(nqr))
      endif


      write(11,'(a,i6,a)') 'Using ',ncol,' collection semi-angles.'
      allocate(ints(ncol,3),sp2(ncol),collection(ncol))
      if(qmodus.eq.'U') then
         write(11,*) 'Using betastep = ',betastep
      else
         write(11,*) 'Using dx = ',dx
      endif
      ints(:,:)=dble(0)
      sp2(:)=dble(0)
      collection(:)=dble(0)
      if (headers) write(59,'(a)') '#    beta        sp2        pi        sigmadip        total'




      call qmesh(ebeam-s(imagic,1),nposold)
      qve=qv


      do icol=1,ncol
!     Loop over collection angles :      
         if (qmodus.eq.'L'.or.qmodus.eq.'l'.or.qmodus.eq.'1') then
            if(icol.eq.1) then
               acoll=th0
            else
               acoll=th0*dexp(dble(icol-1)*dx)
            endif
         else
            acoll=betastep*icol
         endif
         collection(icol)=acoll
         if(qmodus.eq.'1') then
            npos=icol
         else
            npos=icol*icol*nqf
         endif
         nqr=icol
         ThPart = (aconv + acoll) / DBLE(2*nqr)

!     For every collection angle, new integration weights have to be calculated :        
         call destroy_qvecs1
         call CalculateWeights
         
         if(npos.eq.1) then
            if(aconv.gt.0.00001) then
               weightv(1)=pi*((aconv+ acoll)*min(aconv,acoll)/aconv)**2  !12/10 took out 1000^2
            else
               weightv(1)=pi*((aconv+ acoll))**2 !id
            endif
            write(11,*) weightv(1)
         endif
!     Normally, we don't want weights for just one position - since, then, we feel the user does not want integration, but the DOUBLE dscs.
!     However, here we do need integration.  So we make the weight ourselves.

!     We integrate all contributions to the spectrum separately :
         do j=1,npos

            if(relatQ) then
               qfac=(QLenVClas(1,j)**2-(s(imagic,1)/hbarc_eV)**2)**2
            else
               qfac=QLenVClas(1,j)**4
            endif
            ints(icol,1)=ints(icol,1)+WeightV(j)/qfac   *   qve(3,1,j)*qve(3,1,j)*s(imagic,10) ! pi
            ints(icol,2)=ints(icol,2)+WeightV(j)/qfac   * ( qve(1,1,j)*qve(1,1,j)*s(imagic,2) +qve(2,1,j)*qve(2,1,j)*s(imagic,6) ) ! sigma dipole
         enddo                  !j
         ints(icol,3)=ints(icol,1)+ints(icol,2) ! total
         if(dabs(ints(icol,3)).gt.0)    sp2(icol)=ints(icol,1)/ints(icol,3) ! sp2 = pi/total

         write(59,'(5(F14.9,x),I7)') collection(icol),sp2(icol), (ints(icol,j),j=1,3),npos

      enddo                     ! icol


!     we restore old data :
      aconv=aconvold
      acoll=acollold
      nqr=nqrold
      nqf=nqfold
      npos=nposold
      ThPart = (aconv + acoll) / DBLE(2*nqr) !KJ
      call destroy_qvecs1
      call destroy_qvecs2       ! It is unlikely that they would be needed at the exact same energy, and dangerous to leave them lying around.
      call CalculateWeights
!     clean up :
      deallocate(ints,sp2,collection,qve,weightvold)


      return
      end
