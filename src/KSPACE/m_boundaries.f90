!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_boundaries.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/04/03 22:39:49 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
!       DIMENSIONS FOR THE KKR STRUCTURE FACTOR ROUTINES
!*****************************************************************************
      module boundaries
      integer nlmax,nkmax,nkmmax,nkmpmax,nmuemax,nkkrmax
      integer nqmax
      integer NLLMMMAX,NQQPMAX,LGNT12,LLARR
	  integer NRDL,NG,NR ! The actualized values for NRDLMAX,NGRLMAX,NRDLMAX0 !KJ 2-2012

! summation limit >> setting up the structure constants
      integer,parameter :: J13MIN    = 113     
      integer,parameter :: J22MAX    = 100 !200 !120 !30  !Increase to 200 for energies up to E=40 Ry (series will peak at j=100)      
      !integer,parameter :: NMARR     = 15 !12   
	  integer  :: NMARR = 10  !Will be set to appropriate value in straa   
! store all non-zero gaunt coefficients:
! set dynamically upon initialization
!                 LGNT123 >=   100  for  l_max = 2
!                 LGNT123 >=   400  for  l_max = 3
!                 LGNT123 >=  1200  for  l_max = 4
      integer :: LGNT123 =  1200    

! The following 3 are no longer as important as before.
! Dynamical values are now used for allocation of arrays throughout the program.
! The 3 below are only used inside strvecgen for local work arrays.
! As such, however, they still function as maximum values.
	  
! G-vectors in reciprocal lattice
      integer,parameter :: ngrlmax =  5000
! R-vectors in direct lattice
      integer,parameter :: nrdlmax =  5000
! R-vectors in direct lattice -- initially
      integer,parameter :: nrdlmax0=  10000

!      integer :: ngrlmax =  5000
!      integer :: nrdlmax =  5000
!      integer :: nrdlmax0=  10000

!  FINALLY, these are MINE  (KJ) :
      integer maxl  ! maximal angular momentum index in the crystal
      integer msize ! size of most matrices in reciprocal space
      integer mls   ! size of matrices that do not depend on position


        CONTAINS
          subroutine init_boundaries(lin,nq)
            use controls,only : cplxylm
            implicit none
            integer,intent(in) :: lin,nq
! atom positions
            nqmax = nq
! gaunt coefficients:
		    if(lin.le.2) then
			   lgnt123=100
			elseif(lin.eq.3) then
			   lgnt123=400
			elseif(lin.ge.4) then
			   lgnt123=1200  !Will probably crash if lin>4 - but then so will much of the FEFF code ...
			endif
! angular momentum expansion of wave functions
            nlmax=lin+1
            nkmax=2*nlmax-1
            nkmmax=2*nlmax**2
            nkmpmax=2*nlmax**2+2*nlmax
            nmuemax=2*nlmax
            nkkrmax = nqmax * 2 * nlmax**2 
            if(cplxylm) then
               lgnt12=nlmax**4
            else   ! matrix is symmetric
             lgnt12=(nlmax**2*(nlmax**2+1)/2)
            endif
! inequivalent combination of lattice sites in structure constants matrix
            nqqpmax = nqmax*(nqmax-1)+1
            llarr=2*lin 
            nllmmmax=(2*lin+1)**2 

        end subroutine init_boundaries
        end module boundaries
