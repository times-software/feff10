!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_kgenwork.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/01/30 06:01:58 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module kgenwork
!    Number of k-points asked for
      integer nka
!    Number of k-points in full mesh
      integer nkf
!    Number of k-points in irreducible mesh
      integer nki
!    Number of k-points in work arrays
      integer nkw
!    The mesh of k-points (full,irreducible,work)
      real*8,allocatable :: bkf(:,:),bki(:,:),bkw(:,:)
!    The mesh of k-points in sublattice units
      integer,allocatable :: bkwi(:,:)
!    A spare copy (different units ...)
      real*8,allocatable :: bki2(:,:)
!    The weights of k-points (full,irreducible,work)
      real*8,allocatable :: wf(:),wi(:),ww(:)
!    Correspondence of full to irr mesh
      integer,allocatable :: linkf(:),linkw(:)
!    Same correspondence - symmetry used
      integer,allocatable :: lsymf(:),lsymw(:)

        CONTAINS
          subroutine init_fullmesh(n)
            integer n
            allocate(bkf(3,n),wf(n),linkf(n),lsymf(n))
          end subroutine init_fullmesh

          subroutine init_irrmesh(n)
            integer n
            allocate(bki(3,n),wi(n))
            allocate(bki2(3,n))
          end subroutine init_irrmesh
        
          subroutine init_workmesh(n)
            integer n
            allocate(bkw(3,n),bkwi(3,n),ww(n),linkw(n),lsymw(n))
          end subroutine init_workmesh

          subroutine destroy_workmesh
            deallocate(bkw,bkwi,ww,linkw,lsymw)
          end subroutine destroy_workmesh

          subroutine destroy_meshes
            if(allocated(bkw)) deallocate(bkw)
            if(allocated(bkwi)) deallocate(bkwi)
            if(allocated(ww)) deallocate(ww)
            if(allocated(linkw)) deallocate(linkw)
            if(allocated(lsymw)) deallocate(lsymw)
            if(allocated(bki)) deallocate(bki)
            if(allocated(bki2)) deallocate(bki2)
            if(allocated(wi)) deallocate(wi)
            if(allocated(bkf)) deallocate(bkf)
            if(allocated(wf)) deallocate(wf)
            if(allocated(linkf)) deallocate(linkf)
            if(allocated(lsymf)) deallocate(lsymf)
          end subroutine destroy_meshes	


      end module kgenwork
