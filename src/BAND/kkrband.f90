!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: kkrintegral.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/03/16 22:55:02 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine kkrband(G,Tmat,p,bk,freeprop,eryd)

  !  * create Gamma(k) * T (matrices of dimension (lx+1)**1 * natlat
  !  * invert (1-Gamma(k)T) 

  use struct,only: nats,alat,absorber,ppot
  use boundaries,only: maxl,msize,mls

  implicit none

  ! INPUT
  complex,intent(in) ::  Tmat(msize,msize)   ! scattering matrix
  complex*16,intent(in) :: p,eryd !energy at which we work - one point, not a mesh
  real*8, intent(in) ::  bk(3)
  logical,intent(in) :: freeprop
  ! OUTPUT
  complex,intent(out) :: G(msize,msize)  ! The fms-matrix G_fms - also used as a work array in this routine
  ! LOCAL
  complex Mmat(msize,msize),TT(msize,msize),SS(msize,msize)
  integer info,i,j,ik,ipiv(msize)
  character*13 trans
  integer n,m
  real*8,parameter :: pi = 3.1415926535897932384626433d0
  complex,parameter :: ie=(0.0,1.0),c0=(0.0,0.0),c1=(1.0,0.0)  ! careful : these are defined as complex*16 in other routines!
  logical,parameter :: debug=.false.


  !***************** PREPARATION **************************************
  !   Process input and initialize variables.
  G=cmplx(0,0)
  trans = 'NotTransposed'
  n=mls
  m=msize

        !        prepare Gamma(k)
        !note : G contains (-1) * structure factor
        call structurefactor(p,bk(:),G)
		if(freeprop)  return

        !       matrix inversion of Tmat  into Mmat
        call cgetrf(msize,msize,Tmat,msize,ipiv,info)
        if (info.ne.0)  call par_stop('Error in cgetrf when processing T')
        call unity(Mmat,msize)
        call cgetrs(trans,msize,msize,Tmat,msize,ipiv,Mmat,msize,info)
        if (info.lt.0)  call par_stop('*** Error in cgetrs when computing T^-1')

        !        now add Mmat        
        TT = G - Mmat
        G=TT
         if(debug) write(44,440) eryd,p,Mmat(1,1),G(1,1)
440      format(100(f14.7,2x))

		
!        !        matrix inversion of TT=((T^-1)-Gamma(k))  into SS=((T^-1)-Gamma(k))^-1
!        call cgetrf(msize,msize,TT,msize,ipiv,info)
!        if (info.ne.0)  call par_stop('Error in cgetrf when processing 1-Gamma(k)T')
!        call unity(SS,msize)
!        call cgetrs(trans,msize,msize,TT,msize,ipiv,SS,msize,info)
!        if (info.lt.0)  call par_stop('*** Error in cgetrs when computing (1-Gamma(k)T)^-1')
!
!        G=SS


  return
end subroutine kkrband 





!**********************************************************************************************************
!      subroutine unity(a,n)
!! Initialize a to the unit matrix of dimension n.
!        implicit none
!        integer n,i
!        complex a(n,n)
!        complex,parameter :: nul=(0.0,0.0),een=(1.0,0.0)
!
!      a=nul
!        do i=1,n
!           a(i,i)=een
!        enddo
!
!        return
!        end   ! unity


