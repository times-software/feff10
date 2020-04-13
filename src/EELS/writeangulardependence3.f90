!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: writeangulardependence3.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/12/16 18:30:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOP
!     !ROUTINE: WriteAngularDependence3
!     !IntERFACE:
      subroutine WriteAngularDependence3
!     !USES:
      use eels_inp
      use work
      use qvectors
      use program_control
      use constants,only : pi,hbarc_eV
      use spectrum,only : s,ne
!     !DESCRIPTION:
!
!     !REVISION HISTORY:
!     Created November 2004 (Kevin Jorissen)
!     Adapted for FEFF October 2005 (Kevin Jorissen)
!     EOP

      implicit none

!     LOCAL VARIABLES
      integer i,iq,info2_3,info1_4,nqq
      real*8 a0,info1_1,info1_2,info1_3,info2_1,info2_2,qmin,qmax,qfac
      real*8 prefac,gamma
      real*8,allocatable :: qq(:),xq(:,:)
      


      write(11,'(/,a)') 'Output from subroutine PlotAngularDependence3:'
      if(aver.ne.1) then
         call wlog ('Plotangulardependence3 requires averaging!!')
         return
      endif

! For now, start from info-variables
      info1_1=dble(0.44)
      info1_2=dble(0.1950)
      info1_3=dble(100)  !I have no idea what this is
      info1_4=20
      info2_1=dble(100) !for energy mesh
      info2_2=dble(10) !id
      info2_3=12 !id
      info2_3=ne  !my hack for now
      a0=dble(0.529177)


      qmin=info1_1*(dexp(info1_2)-dble(1))*a0
      nqq=info1_4
      qmax=info1_1*(dexp(dble(nqq)*info1_2)-dble(1))*a0

! Prepare to later start from given qmin,qmax,nq :

      info1_2=dlog((dble(1)+qmax)/(dble(1)+qmin))/dble(nqq-1)
      info1_1=qmin/(a0*(dexp(info1_2)-dble(1)))
      
      allocate(qq(nqq))
      do i=1,nqq
         qq(i)=info1_1*(dexp(dble(i)*info1_2)-dble(1))*a0
      enddo

      allocate(xq(nqq,ne))
      xq=dcmplx(0)
 
      
      call wlog ('Entering big loop over energy.')            
!     now, we have a loop over energy loss :
      do i=1,ne
         gamma=(dble(1) + ebeam/dble(511004))
         prefac=s(i,1)*ebeam*(dble(1)+gamma)/(dble(2)*gamma**2)  /(dble(4)*pi*dble(13.6)**2) * dble(1000)
         do iq=1,nqq
!     in dipole approximation, now a q-dependent factor (q'_i q'_j) / (q^2-(E/hbar c)^2 has to be added
!     Spectrum(iq,E) =     * SUM (i=1,3 ; j=1,3)  Spectrum(iq;i,j) * q-dependent_factor(iq;i,j2)

            if(relatQ) then
               qfac=(qq(iq)**2-(s(i,1)/hbarc_eV)**2)**2
            else
               qfac=qq(iq)**4
            endif
            xq(iq,i)=xq(iq,i)+qq(iq)**2/qfac *s(i,2)*prefac
         enddo

      enddo


      open(101,file='gos1.txt')
      do iq=1,nqq
         write(101,*) qq(iq),xq(iq,ne/2+1)
      enddo
      close(101)
      
      open(102,file='gos2.txt')
      write(102,'(a4)') 'OXYG'
      write(102,'(a6,2f7.4,f6.1,i3)') ' 1S1/2',info1_1,info1_2,info1_3, info1_4
      write(102,'(2f8.2,i3)') info2_1,info2_2,info2_3
      do i=1,ne
         write(102,'(5e16.8)') xq(:,i)
      enddo
      close(102)



      return
      end
