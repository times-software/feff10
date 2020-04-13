!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_workstrfacs2.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2012/03/27 22:46:31 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       MORE WORK ARRAYS FOR THE KKR STRUCTURE FACTORS
!*****************************************************************************

!KJ  Note: I'm generally a proponent of zeroeing everything explicitly.
! However, some of the arrays here are huge, and it takes time.
! So, I'm checking carefully whether it's really necessary given the way the code is currently (2-2012)
! and eliminating it where possible (i.e. where entire array is explicitly initialized before being used).

      module workstrfacs2
! replaces COMMON /STRACD/ , /STRCOM/ , /STRETA/ , /STRG/, /STRTSC/ , /STRHP/

      REAL*8 ALPHA0,ETA,ETA0
        real*8,allocatable :: BGX(:),BGY(:),BGZ(:),BRX(:),BRY(:),BRZ(:),&
     &       CQMLTAB(:,:),GGJLRS(:,:,:,:),HP(:),QJLTAB(:,:),            &
     &       QQPX(:),QQPY(:),QQPZ(:),QX(:),QY(:),QZ(:),T(:,:)
      COMPLEX*16 d300
        complex*16,allocatable :: D1TERM3(:),EXPGNQ(:,:),QQMLRS(:,:,:)
      INTEGER,allocatable :: G1(:),G2(:),G3(:),INDR(:,:),R1(:),R2(:),R3(:)
      integer g123max,r123max

      CONTAINS
          subroutine init_workstrfacs2
           use boundaries,only: llarr,nllmmmax,nqmax
           implicit none
           real*8, parameter :: rnul=0.d0
           complex*16,parameter :: cnul=(0.d0,0.d0)
		   ! Very small or "exactly the right size" allocations:
           allocate(BGX(3),BGY(3),BGZ(3),BRX(3),BRY(3),BRZ(3), QX(NQMAX),       &
             QY(NQMAX),QZ(NQMAX),QJLTAB(0:LLARR,0:LLARR),D1TERM3(0:LLARR),         &
             CQMLTAB(-LLARR:LLARR,0:LLARR),T(0:LLARR,0:LLARR), HP(NLLMMMAX)  )

           bgx=rnul;bgy=rnul;bgz=rnul;brx=rnul;bry=rnul;brz=rnul
           d300=cnul;d1term3=cnul;eta0=rnul
           !g1=0;g2=0;g3=0;r1=0;r2=0;r3=0;indr=0
           cqmltab=rnul;hp=rnul;qjltab=rnul;qx=rnul;qy=rnul;qz=rnul;t=rnul
          end subroutine init_workstrfacs2

          subroutine init_workstrfacs2_b(n)
		  integer,intent(in) :: n
           real*8, parameter :: rnul=0.d0
	       ! Allocations containing nqqpmax:
		   if(n.lt.1) stop 'error calling init_workstrfacs2_b too soon'
		   allocate(QQPX(n),QQPY(n),QQPZ(n))
!		   qqpx=rnul;qqpy=rnul;qqpz=rnul;
		  end subroutine init_workstrfacs2_b
		  
		  subroutine exit_workstrfacs2_b
		  deallocate(qqpx,qqpy,qqpz)
		  end subroutine exit_workstrfacs2_b
		  
		  
		  subroutine init_workstrfacs2_c(nr0,ndr,n) !for nrdlmax0,nrdlmax,nqqpmax
           use boundaries,only: nllmmmax,llarr,j22max
           implicit none
		   integer,intent(in) :: nr0,ndr,n
           real*8, parameter :: rnul=0.d0
           complex*16,parameter :: cnul=(0.d0,0.d0)
		   ! Allocations containing nrdlmax/nrdlmax0:  (very large)
			 if(n.lt.1.or.ndr.lt.1.or.nr0.lt.1) stop 'error calling init_workstrfacs2_c too soon'
		   allocate( GGJLRS(0:J22MAX,0:LLARR,ndr,N),QQMLRS(NLLMMMAX,ndr,N),INDR(nr0,N),R1(nr0),R2(nr0),R3(nr0))
!           qqmlrs=rnul;ggjlrs=rnul
           r1=0;r2=0;r3=0;indr=0  !these arrays are not fully initialized because smax(:) depends on iqqp
		   ! they are only used for the initialized fields - still, reason for caution.
          end subroutine init_workstrfacs2_c

		  subroutine init_workstrfacs2_d(ng,n)  ! for ngrlmax,nqqpmax
           implicit none
		   integer,intent(in):: ng,n
           complex*16,parameter :: cnul=(0.d0,0.d0)
		   ! Allocations containing ngrlmax:  (very large)
			 if(n.lt.1.or.ng.lt.1) stop 'error calling init_workstrfacs2_d too soon'
		   allocate( EXPGNQ(ng,N),G1(ng),G2(ng),G3(ng))
!           expgnq=cnul
!           g1=0;g2=0;g3=0
          end subroutine init_workstrfacs2_d
		  
          subroutine exit_workstrfacs2
           use boundaries
           implicit none
           deallocate(BGX,BGY,BGZ,BRX,BRY,BRZ,CQMLTAB,GGJLRS,HP,QJLTAB, &
     &       QQPX,QQPY,QQPZ,QX,QY,QZ,T,D1TERM3,EXPGNQ,QQMLRS,G1,G2,G3,  &
     &       INDR,R1,R2,R3)
          end subroutine exit_workstrfacs2


        end module workstrfacs2
