!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: structurefactor.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2012/02/17 07:39:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine structurefactor(p,k,gk)

        use struct,only:nats,alat
        use trafo
      use boundaries,only:nkkrmax,msize,mls,maxl
        use constants,only: pi

      implicit none
! in/output
        real*8,intent(in) :: k(3)  !  in a.u. and carthesian coordinates
        complex*16,intent(in) :: p
        complex,intent(out) :: gk(msize,msize)
! locals
        integer a,na,ia,b,nb,ib,j1,j2,joff1,joff2
      complex*16 maux(nkkrmax,nkkrmax),tauk(nkkrmax,nkkrmax),w1(mls,mls)
        complex*16 w2(mls,mls),pinv
        complex*16,parameter :: coni=(0.0,1.0)
        real*8 kfac,klocal(3)
		logical,parameter :: debug=.false.
		logical,parameter :: feff_basis=.true.

     !   write(*,*) 'entering structurefactor'
	!	write(*,*) 'p',p
	!	write(*,*) 'k',k
	!	call writematrix(gk,msize,'gk_1.txt')

        kfac= alat(1) /(dble(2)*pi)

      maux=dcmplx(0,0)
        tauk=dcmplx(0,0)
        gk=cmplx(0,0)
        pinv=dcmplx(1,0)/p

        klocal=k*kfac  ! Input is k-vector in carthesian coordinates and a.u. ; we rescale here by a/2pi

      if(debug) call wlog('calling strset')
      call strset(klocal(1),klocal(2),klocal(3),maux,tauk,p)


	if(debug) then
       open(12,file='hahaha.txt')
	   write(12,*) 'k= ',k
	   write(12,*) 'p= ',p
	   do a=1,msize
	      write(12,'(100e12.4)') tauk(a,1:msize)
	   enddo
	   close(12)
!	   stop
    endif


	   if(.not.feff_basis) then
          gk(:,:)=tauk(1:msize,1:msize)
          return
       endif


      do j1=1,nats
      do j2=1,nats
         joff1=(j1-1)*mls
         joff2=(j2-1)*mls      
         w1=dcmplx(0,0)
!        CAREFUL : In the Rehr-Albers convention, structure factors and t-matrices are dimensionless.
!                  Therefore, we remove the SPRKKR-factor of p.
           w2=tauk(1+joff1:mls+joff1,1+joff2:mls+joff2)*pinv
!        Convert from SPRKKR basis to FEFF basis
!        For nsp=1 : from real to spherical harmonics
         call changerep(w2,'RLM>CLM',w1,mls,mls,rc,crel,rrel,'          ',0)
!        To get the normalization of the spherical harmonics right, add factor i^l :
           ia=0;ib=0
           do a=0,maxl;do na=-a,a;ia=ia+1;ib=0
           do b=0,maxl;do nb=-b,b;ib=ib+1
              w1(ia,ib)=w1(ia,ib)*coni**(a-b)
           enddo;enddo
           enddo;enddo
        gk(1+joff1:mls+joff1,1+joff2:mls+joff2)=cmplx(w1)
      enddo
      enddo

!call writematrix(gk,msize,'gk.txt')
!stop
      return
      end
