!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_ptable.f90,v $:
! $Revision: 1.3 $
! $Author: fer $
! $Date: 2010/09/20 18:24:13 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_PTable

! Define a periodic table and procedure to access its information

  use m_Kinds
  use m_Strings

! type :: PTable_Info
! end type PTable_Info

! type(PTable_Info), parameter :: PTable

  integer(kind=i08), parameter :: PTable_nElem = 103
  real(kind=r08), parameter ::    PTable_AM(PTable_nElem) = (/ &
     1.00798d0,   4.0026d0,    6.941d0,     9.012182d0, 10.811d0, &
    12.011d0,    14.00672d0,  15.9994d0,   18.9984d0,   20.18d0, &
    22.989767d0, 24.3051d0,   26.981538d0, 28.0855d0,   30.973762d0, &
    32.064d0,    35.4527d0,   39.948d0,    39.0983d0,   40.078d0, &
    44.95591d0,  47.88d0,     50.9415d0,   51.9961d0,   54.93805d0, &
    55.844d0,    58.9332d0,   58.69d0,     63.546d0,    65.4d0, &
    69.723d0,    72.59d0,     74.92159d0,  78.96d0,     79.904d0, &
    83.8d0,      85.4678d0,   87.62d0,     88.90585d0,  91.224d0, &
    92.90638d0,  95.93d0,     98.0d0,     101.07d0,    102.9055d0, &
   106.42d0,    107.8681d0,  112.412d0,   114.82d0,    118.71d0, &
   121.76d0,    127.6d0,     126.90447d0, 131.29d0,    132.90544d0, &
   137.327d0,   138.9054d0,  140.115d0,   140.90765d0, 144.24d0, &
   145.0d0,     150.36d0,    151.965d0,   157.25d0,    158.92534d0, &
   162.5d0,     164.93032d0, 167.26d0,    168.93421d0, 173.04d0, &
   174.967d0,   178.49d0,    180.9479d0,  183.85d0,    186.207d0, &
   190.2d0,     192.22d0,    195.08d0,    196.96654d0, 200.6d0, &
   204.3833d0,  207.2d0,     208.98037d0, 209.0d0,     210.0d0, &
   222.0d0,     223.0d0,     226.0d0,     227.0d0,     232.0381d0, &
   231.04d0,    238.0289d0,  237.0d0,     244.0d0,     243.0d0, &
   247.0d0,     247.0d0,     251.0d0,     252.0d0,     257.0d0, &
   258.0d0,     259.0d0,     262.0d0   /)
 
  integer, parameter :: PTable_Symb_Len = 2

! Periodic table in our standard format:
!   Symbols are left-justified
!   First char is uppercase
!   Second char is lowercase
  character(len=PTable_Symb_Len), parameter :: &
    PTable_Symb(PTable_nElem) = (/ &
    'H ',    'He',    'Li',    'Be',    'B ', &
    'C ',    'N ',    'O ',    'F ',    'Ne', &
    'Na',    'Mg',    'Al',    'Si',    'P ', &
    'S ',    'Cl',    'Ar',    'K ',    'Ca', &
    'Sc',    'Ti',    'V ',    'Cr',    'Mn', &
    'Fe',    'Co',    'Ni',    'Cu',    'Zn', &
    'Ga',    'Ge',    'As',    'Se',    'Br', &
    'Kr',    'Rb',    'Sr',    'Y ',    'Zr', &
    'Nb',    'Mo',    'Tc',    'Ru',    'Rh', &
    'Pd',    'Ag',    'Cd',    'In',    'Sn', &
    'Sb',    'Te',    'I ',    'Xe',    'Cs', &
    'Ba',    'La',    'Ce',    'Pr',    'Nd', &
    'Pm',    'Sm',    'Eu',    'Gd',    'Tb', &
    'Dy',    'Ho',    'Er',    'Tm',    'Yb', &
    'Lu',    'Hf',    'Ta',    'W ',    'Re', &
    'Os',    'Ir',    'Pt',    'Au',    'Hg', &
    'Tl',    'Pb',    'Bi',    'Po',    'At', &
    'Rn',    'Fr',    'Ra',    'Ac',    'Th', &
    'Pa',    ' U',    'Np',    'Pu',    'Am', &
    'Cm',    'Bk',    'Cf',    'Es',    'Fm', &
    'Md',    'No',    'Lr'  /)

  private
  public :: PTable_nElem , PTable_Symb_Len
  public :: PTable_AM, PTable_Symb, PTable_AN
  public :: Find_AN

contains

! Find the atomic number stored in a label or the atomic number
! associated with the symbol stored in the label.
  function Find_AN(Label) result(an)

    implicit none

    character(len=PTable_Symb_Len+1) :: Label
    integer :: an

    integer :: iElem
    character(len=PTable_Symb_Len)   :: Symb

! Check if we have all numbers
    if ( Is_AllNum(Label) ) then
      read(Label,*) an
      return
    end if

! Check if we have all characters
    if ( Is_AllAlpha(Label) ) then

! Adjust the label to get it in our standard symbol format (see above)
      Symb = adjustl(Label)
      Symb(1:1) = Uppercase(Symb(1:1))
      Symb(2:2) = Lowercase(Symb(2:2))

! Search for the symbol in the PTable
      do iElem=1,PTable_nElem
        if ( Symb == PTable_Symb(iElem) ) then
          an = iElem
          return
        end if
      end do

! If we get here, then we couldn't find the element
! NOTE: Printing from a function doesn't produce a <cr>, need to fix
      write(6,fmt='(a)') 'Element not found in Find_AN'
      stop
    end if

! If neither, then we have a problem
! NOTE: Printing from a function doesn't produce a <cr>, need to fix
    write(6,fmt='(a)') ' Error in Find_AN'
    stop

  end function Find_AN

! Find the atomic number associated with an atomic symbol
  function PTable_AN(Symb_Label) result(an)

    implicit none

    character(len=*) :: Symb_Label
    integer          :: an

    integer :: iElem
    character(len=PTable_Symb_Len)   :: Symb

! The first two non-space characters of a symbol label should correspond
! to the symbol
    Symb = adjustl(Symb_Label)

! Check if we have all alphabetic characters
    if ( Is_AllAlpha(Symb) ) then

! Adjust the symbol to get it in our standard symbol format (see above)
      Symb(1:1) = Uppercase(Symb(1:1))
      Symb(2:2) = Lowercase(Symb(2:2))

! Search for the symbol in the PTable
      do iElem=1,PTable_nElem
        if ( Symb == PTable_Symb(iElem) ) then
          an = iElem
          return
        end if
      end do

! If we get here, then we couldn't find the element
! NOTE: Printing from a function doesn't produce a <cr>, need to fix
      write(6,fmt='(a)') 'Error in PTable_AN: Element not found'
      stop

    else

! Throw error if we have non-alphabetic characters
! NOTE: Printing from a function doesn't produce a <cr>, need to fix
      write(6,fmt='(a)') 'Error in PTable_AN: Symbol label has wrong format'
      stop

    end if

  end function PTable_AN

end module m_PTable


