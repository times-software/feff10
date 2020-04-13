!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_mtdp.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Mtdp

! Muffin-Tin density and potential

  use m_Kinds

  type :: Mtdp_Data_Type
    integer(kind=i08)                            :: nR
    integer(kind=i08)                            :: nAt
    integer(kind=i08), dimension(:), allocatable :: At_AN
    real(kind=r08), dimension(:,:), allocatable  :: At_XYZ
    real(kind=r08), dimension(:), allocatable    :: At_R
    integer(kind=i08), dimension(:), allocatable :: At_iR
    real(kind=r08), dimension(:,:), allocatable  :: At_Den
    real(kind=r08), dimension(:,:), allocatable  :: At_Pot
    integer(kind=i08)                            :: nESph
    real(kind=r08), dimension(:,:), allocatable  :: ESph_XYZ
    real(kind=r08), dimension(:), allocatable    :: ESph_R
    integer(kind=i08), dimension(:), allocatable :: ESph_iR
    real(kind=r08), dimension(:,:), allocatable  :: ESph_Den
    real(kind=r08), dimension(:,:), allocatable  :: ESph_Pot
    real(kind=r08)                               :: V_Int
    real(kind=r08)                               :: V_HOMO
    real(kind=r08)                               :: V_LUMO
  end type

  private
  public :: Mtdp_Data_Type
  public :: Write_Mtdp
  public :: Read_Mtdp

contains

  subroutine Write_Mtdp(Mtdp_Data)

    implicit none

    type(Mtdp_Data_Type), intent(in) :: Mtdp_Data

! Number of points in the radial grid
    write(6,fmt='(i12)') Mtdp_Data%nR

! Number of atoms
    write(6,fmt='(i12)') Mtdp_Data%nAt

! Atomic numbers
    write(6,fmt='(i12)') Mtdp_Data%At_AN

! Atomic coordinates
    write(6,fmt='(1p,e20.10)') Mtdp_Data%At_XYZ

! Muffin-Tin Radii
    write(6,fmt='(1p,e20.10)') Mtdp_Data%At_R

! Index of the radii on the radial grid
    write(6,fmt='(i12)') Mtdp_Data%At_iR

! Electron density inside each Muffin-Tin
    write(6,fmt='(1p,e20.10)') Mtdp_Data%At_Den

! Potential (H+XC) inside each Muffin-Tin
    write(6,fmt='(1p,e20.10)') Mtdp_Data%At_Pot

! Number of empty spheres
    write(6,fmt='(i12)') Mtdp_Data%nESph

! Coordinates of the empty sphere's centers
    write(6,fmt='(1p,e20.10)') Mtdp_Data%ESph_XYZ

! Empty sphere's radii
    write(6,fmt='(1p,e20.10)') Mtdp_Data%ESph_R

! Index of the empty sphere's radii on the radial grid
    write(6,fmt='(i12)') Mtdp_Data%ESph_iR

! Electron density inside each empty sphere
    write(6,fmt='(1p,e20.10)') Mtdp_Data%ESph_Den

! Potential (H+XC) inside each empty sphere
    write(6,fmt='(1p,e20.10)') Mtdp_Data%ESph_Pot

! Interstitial potential
    write(6,fmt='(1p,e20.10)') Mtdp_Data%V_Int

! HOMO energy
    write(6,fmt='(1p,e20.10)') Mtdp_Data%V_HOMO

! LUMO energy
    write(6,fmt='(1p,e20.10)') Mtdp_Data%V_LUMO

  end subroutine Write_Mtdp

  subroutine Read_Mtdp(iunit,Mtdp_Data)

    implicit none

    integer, intent(in) :: iunit
    type(Mtdp_Data_Type), intent(out) :: Mtdp_Data

! Number of points in the radial grid
    read(iunit,*) Mtdp_Data%nR

! Number of atoms
    read(iunit,*) Mtdp_Data%nAt

! Atomic numbers
    allocate(Mtdp_Data%At_AN(Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_AN

! Atomic coordinates
    allocate(Mtdp_Data%At_XYZ(3,Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_XYZ

! Muffin-Tin Radii
    allocate(Mtdp_Data%At_R(Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_R

! Index of the radii on the radial grid
    allocate(Mtdp_Data%At_iR(Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_iR

! Electron density inside each Muffin-Tin
    allocate(Mtdp_Data%At_Den(Mtdp_Data%nR,Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_Den

! Potential (H+XC) inside each Muffin-Tin
    allocate(Mtdp_Data%At_Pot(Mtdp_Data%nR,Mtdp_Data%nAt))
    read(iunit,*) Mtdp_Data%At_Pot

! Number of empty spheres
    read(iunit,*) Mtdp_Data%nESph

! Coordinates of the empty sphere's centers
    allocate(Mtdp_Data%ESph_XYZ(3,Mtdp_Data%nESph))
    read(iunit,*) Mtdp_Data%ESph_XYZ

! Empty sphere's radii
    allocate(Mtdp_Data%ESph_R(Mtdp_Data%nESph))
    read(iunit,*) Mtdp_Data%ESph_R

! Index of the empty sphere's radii on the radial grid
    allocate(Mtdp_Data%ESph_iR(Mtdp_Data%nESph))
    read(iunit,*) Mtdp_Data%ESph_iR

! Electron density inside each empty sphere
    allocate(Mtdp_Data%ESph_Den(Mtdp_Data%nR,Mtdp_Data%nESph))
    read(iunit,*) Mtdp_Data%ESph_Den

! Potential (H+XC) inside each empty sphere
    allocate(Mtdp_Data%ESph_Pot(Mtdp_Data%nR,Mtdp_Data%nESph))
    read(iunit,*) Mtdp_Data%ESph_Pot

! Interstitial potential
    read(iunit,*) Mtdp_Data%V_Int

! HOMO energy
    read(iunit,*) Mtdp_Data%V_HOMO

! LUMO energy
    read(iunit,*) Mtdp_Data%V_LUMO

  end subroutine Read_Mtdp

end module Mtdp
