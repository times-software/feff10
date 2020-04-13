program pot_generator_test

! Written by FDV, Aug 2014

! Simple wrapper to test the smart potential type generator module

  use m_Kinds
  use m_Pot_Generator
  use m_PTable

  implicit none

! Test the potential generator from an XYZ format structure

  call Test_Gen_Pot_From_XYZ

contains

  subroutine Test_Gen_Pot_From_XYZ

    type(Pot_Gen_Rules_Format) :: Pot_Gen_Rule
    integer(kind=i08) :: Abs_Index
    type(XYZ_Format) :: Struc
    type(Pot_List_Type) :: Pot_List
    
! Auxiliary stuff
    integer(kind=i08) :: iAt
    character(len=2),  dimension(:),   allocatable :: Symb
    
! Since we are testing, we keep the input very simple:
! Read the name of the rule used for potential type generation
    read(5,*) Pot_Gen_Rule%Name
    
! Debug
!   write(6,*) Pot_Gen_Rule%Name
    
! Read the absorber index
    read(5,*)  Abs_Index
    
! Debug
!   write(6,*) Abs_Index
    
! From here to the end it is a standard xyz file
! Read the number of atoms
    read(5,*) Struc%nAt
    
! Debug
!   write(6,*) Struc%nAt
    
! Allocate the arrays to store the structural info
    allocate(Symb(Struc%nAt),Struc%AN(Struc%nAt),Struc%XYZ(Struc%nAt,3))
    
! Read the XYZ title line
    read(5,fmt='(a80)') Struc%Title
    
! Debug
!   write(6,*) Struc%Title
    
! Read the XYZ structure: Atomic symbol, x, y, z
    do iAt=1,Struc%nAt
      read(5,*) Symb(iAt), &
                Struc%XYZ(iAt,1), Struc%XYZ(iAt,2), Struc%XYZ(iAt,3)
    end do
    
! Debug
!   do iAt=1,Struc%nAt
!     write(6,*) Symb(iAt), &
!                Struc%XYZ(iAt,1), Struc%XYZ(iAt,2), Struc%XYZ(iAt,3)
!   end do
    
! Convert the symbols vector into atomic numbers
    
! Check each atom to find the atomic number, if the symbol is incorrect the
! function should throw an error and stop.
    do iAt=1,Struc%nAt
      Struc%AN(iAt) = PTable_AN(Symb(iAt))
    end do
    
! Generate the potential types
    call Gen_Pot_From_XYZ(Pot_Gen_Rule,Abs_Index,Struc,Pot_List,.true.)
    
! Simple print out of the generated potentials section
    call Print_Pot_List(Pot_List)
    
! Simple print out of the generated atoms section
    call Print_Struc(Struc)
    
  end subroutine Test_Gen_Pot_From_XYZ
    
  subroutine Print_Pot_List(Pot_List)
    
! Simple subroutine to write the generated potentials list.
! NOTE: This is not exactly a feff.inp POTENTIALS section
    
    implicit none
    
    type(Pot_List_Type),  intent(in)    :: Pot_List
    
    integer(kind=i08) :: iPot
    
    write(6,fmt='(a)') '# Generated potentials list:'
    write(6,fmt='(a)') '# Index  AN  Symb  Label (Num)  Label (Str)'
    do iPot=0,Pot_List%nPot
      write(6,fmt='(3x,i2,3x,i3,3x,a,4x,i3,10x,a)') &
        iPot, &
        Pot_List%AN(iPot), &
        PTable_Symb(Pot_List%AN(iPot)), &
        Pot_List%NumLab(iPot), &
        Pot_List%StrLab(iPot)
    end do

  end subroutine Print_Pot_List

  subroutine Print_Struc(Struc)

! Simple subroutine to write the generated atoms list.
! NOTE: This is not exactly a feff.inp ATOMS section

    implicit none

    type(XYZ_Format),                    intent(in)  :: Struc

    integer(kind=i08) :: iAt

    write(6,fmt='(a)') '# Generated atoms list:'
    write(6,fmt='(a)') &
      '#       X           Y          Z      Pot  Symb Label (Num)  Label (Str)'
    do iAt=1,Struc%nAt
      write(6,fmt='(3f12.6,3x,i2,3x,a,3x,i3,10x,a)') &
        Struc%XYZ(iAt,:), &
        Struc%PotInd(iAt), &
        PTable_Symb(Struc%AN(iAt)), &
        Struc%PotNumLab(iAt), &
        Struc%PotStrLab(iAt)
    end do

  end subroutine Print_Struc

end program pot_generator_test

