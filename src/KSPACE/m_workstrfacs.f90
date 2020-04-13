!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_workstrfacs.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       WORK ARRAYS FOR THE KKR STRUCTURE FACTORS
!*****************************************************************************
      module workstrfacs
! replaces COMMON /STRQQP/ , /STRLIM/ , /STRGNT/ , /CSRREL/ , /STRENE/
! plus variables nl,nlm,nkmq,ind0q
      COMPLEX*16,allocatable :: CIPWL(:),SRREL(:,:,:),IILERS(:,:,:)
        complex*16 edu
      REAL*8,allocatable :: gnt(:)
        real*8 GMAX,GMAXSQ,RMAX
      INTEGER,allocatable :: IGNT(:),IJQ(:,:),IRREL(:,:,:),             &
     &        NGNT(:),NIJQ(:),NRREL(:,:),SMAX(:),NKMQ(:),IND0Q(:)
      integer llmax,mmllmax,nmax,nqqp,nrtab,nl,nlm

      CONTAINS
          subroutine init_workstrfacs(nq)
            use boundaries
            implicit none
            integer nq
            allocate(CIPWL((2*NLMAX)**2),SRREL(2,2,NKMMAX),             &
     &      GNT(LGNT123),IGNT(LGNT123),IRREL(2,2,NKMMAX),NGNT(LGNT12),               &
     &        NRREL(2,NKMMAX),nkmq(nq),ind0q(nq)          &
     &        ) !,IILERS(0:LLARR,NRDLMAX,NQQPMAX)  )
          end subroutine init_workstrfacs
		  
          subroutine exit_workstrfacs
            implicit none
            deallocate(CIPWL,SRREL,GNT,IGNT,IJQ,IRREL,NGNT,NIJQ,        &
     &        NRREL,SMAX,nkmq,ind0q,IILERS)
          end subroutine exit_workstrfacs
		  
		  subroutine init_workstrfacs_b
		     use boundaries
			 if(nqqp.lt.1 .or. nrdl.lt.1) stop 'error calling init_workstrfacs_b too soon'
             allocate(IILERS(0:LLARR,NRDL,NQQP))
		  end subroutine init_workstrfacs_b

		  subroutine init_workstrfacs_c
		     use boundaries
			 if(nqqp.lt.1) stop 'error calling init_workstrfacs_c too soon'
             allocate(IJQ(NQQP,NQQP),NIJQ(NQQP),SMAX(NQQP))
		  end subroutine init_workstrfacs_c
          
        end module workstrfacs
