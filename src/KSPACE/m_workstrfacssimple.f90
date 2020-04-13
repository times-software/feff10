!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_workstrfacssimple.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module workstrfacssimple
      
      REAL ETA
      real,allocatable :: BGX(:),BGY(:),BGZ(:),BRX(:),BRY(:),BRZ(:),    &
     &       QQPX(:),QQPY(:),QQPZ(:)
      COMPLEX d300
        complex,allocatable :: D1TERM3(:),EXPGNQ(:,:),QQMLRS(:,:,:)
      INTEGER,allocatable :: G1(:),G2(:),G3(:),R1(:),R2(:),R3(:)
      integer g123max,r123max
        complex edu
        real GMAXSQ
      INTEGER,allocatable :: SMAX(:)
      integer llmax,nmax,nqqp,nrtab,nl,nlm

      CONTAINS
          subroutine init_workstrfacssimple
           use boundaries
           use workstrfacs,only: smaxd=>smax,llmaxd=>llmax,             &
     & gmaxsqd=>gmaxsq,nmaxd=>nmax,nqqpd=>nqqp,                         &
     & nrtabd=>nrtab,nld=>nl,nlmd=>nlm	   
           use workstrfacs2,only: bgxd=>bgx,bgyd=>bgy,bgzd=>bgz,        &
     & brxd=>brx,bryd=>bry,brzd=>brz,qqpxd=>qqpx,                       &
     & qqpyd=>qqpy,qqpzd=>qqpz,g1d=>g1,g2d=>g2,g3d=>g3,                 &
     & r1d=>r1,r2d=>r2,r3d=>r3,g123maxd=>g123max,r123maxd=>r123max,     &
     & etad=>eta
           implicit none
           allocate(BGX(3),BGY(3),BGZ(3),BRX(3),BRY(3),BRZ(3),          &
     &       QQPX(NQQPMAX),QQPY(NQQPMAX),QQPZ(NQQPMAX),                 &
     &       D1TERM3(0:LLARR),EXPGNQ(NGRLMAX,NQQPMAX),                  &
     &       QQMLRS(NLLMMMAX,NRDLMAX,NQQPMAX),                          &
     &       G1(NGRLMAX),G2(NGRLMAX),G3(NGRLMAX),                       &
     &       R1(NRDLMAX0),R2(NRDLMAX0),R3(NRDLMAX0))
            allocate(SMAX(NQQPMAX))
! Copy only arrays that are fixed after initialization :	    	
            bgx=real(bgxd);bgy=real(bgyd);bgz=real(bgzd)
            brx=real(brxd);bry=real(bryd);brz=real(brzd)  
            qqpx=real(qqpxd);qqpy=real(qqpyd);qqpz=real(qqpzd)
            g1=g1d;g2=g2d;g3=g3d;r1=r1d;r2=r2d;r3=r3d
            g123max=g123maxd;r123max=r123maxd
            gmaxsq=real(gmaxsqd)
            eta=real(etad)
            smax=smaxd
            llmax=llmaxd;nmax=nmaxd;nqqp=nqqpd;nrtab=nrtabd
            nl=nld;nlm=nlmd

           end subroutine init_workstrfacssimple

          subroutine copy_workstrfacssimple
           use boundaries
           use workstrfacs,only: edud=>edu
           use workstrfacs2,only: d300d=>d300,d1term3d=>d1term3,        &
     & expgnqd=>expgnq,qqmlrsd=>qqmlrs
           implicit none
! Copy arrays that change for every energy	
            d300=cmplx(d300d)
            d1term3=cmplx(d1term3d)
            expgnq=cmplx(expgnqd)
            qqmlrs=cmplx(qqmlrsd)
            edu=cmplx(edud)
           end subroutine copy_workstrfacssimple
      
      end module workstrfacssimple
