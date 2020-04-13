!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_controls.f90,v $:
! $Revision: 1.7 $
! $Author: bmattern $
! $Date: 2012/02/09 18:04:57 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       CONTROL THE WAY FEFF WORKS
!*****************************************************************************
      module controls
!    Switch to 1 for real space, 0 for reciprocal space
      integer ispace
!    Read sprkkr-structure file
      integer sprkkrstruct
!    Read sprkkr-potential file
      integer sprkkrpot
!    Read sprkkr-klist file
      integer sprkkrklist
!    Switch between real and complex spherical harmonics for the KKR structure factors
!    F for real, T for complex
      logical,parameter :: cplxylm=.false.
!    Arrays are allocated via kprep
      logical allocated
!    Use spin/relativistic matrices or not (LM basis) :
      integer irel
!    Set verbosity of SPRKKR subroutines
      integer iprint
!    Set up k-mesh in ffmod3
      logical makekmeshnow
!    Use a core hole or not
      logical corehole
!    Strength of the core hole - 1 is normal, 0 is nohole
      real*8 cholestrength       ! multiply core hole t-matrix by this number.  Currently strongly suggested to stay away from it!
!    Use single precision in strbbdd
      logical,parameter :: singleprec=.false.
!    Use full potential (t-matrix) or muffin tin potential (phases)
      logical fullpot      
      logical gglu_save_slice !BAM 2/2012

        CONTAINS
        subroutine init_controls
        ispace=1  ! real space
        sprkkrstruct=0
        sprkkrpot=0
        sprkkrklist=0
        makekmeshnow=.false. !use k-mesh from file
        allocated=.false.  ! not yet been in kprep
        irel=1  ! work in LM-basis
        iprint=0
        corehole=.false. !no core hole
        cholestrength=dble(1)
        open(96,file='ini.inp',status='old',err=2341)
        read(96,*,err=2341,end=2341) sprkkrstruct,sprkkrpot,sprkkrklist
        close(96)
2341    continue
        fullpot=(sprkkrpot.eq.1)
        gglu_save_slice = .false. ! BAM - save slice of g with n=0 in gglu
        return	          	
        end subroutine init_controls

        end module controls
