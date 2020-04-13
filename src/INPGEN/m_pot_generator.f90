module m_Pot_Generator

! Written by FDV, Aug 2014

! Get the numerical type kinds
  use m_Kinds

! Access some useful string functions
  use m_Strings

! Access the periodic table data and functions
  use m_PTable

  implicit none

! Maximum number of potentials that we consider "safe" for the SCF
! NOTE: We should get some experience with this to get a better number.
!       Things have worked OK for me up to 40-45, so I'm guessing 50 is a good
!       starting number.
  integer, parameter :: mx_Safe_nPot = 50

! Rule structure type.
! NOTE: This should make it easy to expand the rule syntax in the future
!       For now we only store the rule name
  integer, parameter :: mx_Rule_Name_Len = 8
  integer, parameter :: Pot_Label_Len = 8
  type :: Pot_Gen_Rules_Format
    character(len=mx_Rule_Name_Len) :: Name 
  end type

! The XYZ format, extended to also hold the feff potential type information:
!    The potential index, a numerical label and a string label.
! NOTE: This type is a bit redundant, at least when taken together with the
!       potential list format defined below. But I prefer to keep all the
!       information about the system in one structure, so you don't have to
!       go looking for the corresponding potential list. If this is to much
!       we can always simplify it later.
  type :: XYZ_Format
    integer(kind=i08)                                       :: nAt
    character(len=80)                                       :: Title
    integer(kind=i08), dimension(:),   allocatable          :: AN
    real(kind=r08),    dimension(:,:), allocatable          :: XYZ
    integer(kind=i08), dimension(:),   allocatable          :: PotInd
    integer(kind=i08), dimension(:),   allocatable          :: PotNumLab
    character(len=Pot_Label_Len), dimension(:), allocatable :: PotStrLab
  end type

! The feff potential list, generalized to include both numerical and string
! labels.
  type :: Pot_List_Type
    integer(kind=i08)                                         :: nPot
    integer(kind=i08), dimension(:),   allocatable            :: AN
    integer(kind=i08), dimension(:),   allocatable            :: NumLab
    character(len=Pot_Label_Len), dimension(:),   allocatable :: StrLab
  end type

  private
  public :: mx_Rule_Name_Len
  public :: Pot_Gen_Rules_Format, XYZ_Format
  public :: Pot_List_Type
  public :: Gen_Pot_From_XYZ

contains

  subroutine Gen_Pot_From_XYZ(Pot_Gen_Rule,Abs_Index,Struc, &
                              Pot_List,Verbose_Flag)

    implicit none

    type(Pot_Gen_Rules_Format),          intent(in)     :: Pot_Gen_Rule
    integer(kind=i08),                   intent(in)     :: Abs_Index
    type(XYZ_Format),                    intent(inout)  :: Struc
    type(Pot_List_Type),                 intent(inout)  :: Pot_List
    logical,                   optional, intent(in)     :: Verbose_Flag

    logical           :: Verbose
    integer(kind=i08) :: iAt

! Setup the level of verbosity of the output: By default we run quietly, but
! the user might get extra information as to how we got the potentials
    if ( .not. present(Verbose_Flag) ) then
      Verbose = .false.
    else
      Verbose = Verbose_Flag
    end if

    if (Verbose) then
      write(6,fmt='(a)') '# Generating potentials from XYZ input'
    end if

! Sanity checks on the input
! Check that we have the atomic numbers and coordinates
    if ( .not. allocated(Struc%AN) ) then
      write(6,fmt='(a)') &
        'Error in Gen_Pot_From_XYZ call:', &
        'Missing atomic numbers in structure.'
        stop
    end if

    if ( .not. allocated(Struc%XYZ) ) then
      write(6,fmt='(a)') &
        'Error in Gen_Pot_From_XYZ call:', &
        'Missing cartesian coordinates in structure.'
        stop
    end if

! Check that all sizes are consistent
    if ( Struc%nAt /= size(Struc%AN,1) ) then
      write(6,fmt='(a)') &
        'Error in Gen_Pot_From_XYZ call:', &
        'Inconsistent sizes in structure (AN).'
        stop
    end if

    if ( Struc%nAt /= size(Struc%XYZ,1) .or. &
                 3 /= size(Struc%XYZ,2)      ) then
      write(6,fmt='(a)') &
        'Error in Gen_Pot_From_XYZ call:', &
        'Inconsistent sizes in structure (XYZ).'
        stop
    end if

! Check that the absorber index is in the right range
    if ( Abs_index < 1 .or. &
         Abs_Index > Struc%nAt ) then
      write(6,fmt='(a)') &
        'Error in Gen_Pot_From_XYZ call:', &
        'Absorber index is out of range.'
        stop
    end if

    if (Verbose) then
      write(6,fmt='(a,i3)') '# Will assign absorber to center ', Abs_Index
    end if

! Now we proceed according to the rule name, the default is to use the
! "atomnum" rule.
! NOTE: If you add a rule, please document it briefly right before the call, and
! more extensively inside the handling procedure.
    select case (Lowercase(adjustl(Pot_Gen_Rule%Name)))

! "atomnum" rule: Give all atoms with the same atomic number the same potential
!                 label (with the exception of the absorber)
      case ( "atomnum" )

        if (Verbose) then
          write(6,fmt='(a,a)') '# Using rule: ', &
                               trim(Lowercase(adjustl(Pot_Gen_Rule%Name)))
        end if

        call Rule_atomnum(Abs_Index,Struc,Pot_List,Verbose)

! The default is to use the "atomnum"
      case default

        if (Verbose) then
          write(6,fmt='(a,a,a)') '# Warning: rule ', &
                               trim(Lowercase(adjustl(Pot_Gen_Rule%Name))), &
                               ' not defined.'
          write(6,fmt='(a)') '# Using default atomnum instead.'
        end if

        call Rule_atomnum(Abs_Index,Struc,Pot_List,Verbose)

    end select

    if (Verbose .and. (Pot_List%nPot > mx_Safe_nPot) ) then
      write(6,fmt='(a)') &
        'WARNING:', &
        '  Gen_Pot_From_XYZ generated a large number of unique potentials.', &
        '  This could cause issues in the SCF module.'
    end if

! Before we return, we make sure that the output structures have all the data
! by generating any missing information
! NOTE: This code is tentative until we decide which of the labels is the
!       dominant one used to determine all the rest. For now I will use the
!       potential index/type to create the numerical and string labels. Later
!       on we probably want to do it the other way around.
!       Also, the code below should be made into subroutines.
    if (Verbose) then
      write(6,fmt='(a)') '# Enforcing labels:', &
        '  Making sure that the potentials types and labels are consistent'
    end if

    if ( .not. allocated(Struc%PotNumLab) ) then
      allocate(Struc%PotNumLab(Struc%nAt))
    end if
    if ( .not. allocated(Struc%PotStrLab) ) then
      allocate(Struc%PotStrLab(Struc%nAt))
    end if
    do iAt=1,Struc%nAt
      Struc%PotNumLab(iAt) = Pot_List%NumLab(Struc%PotInd(iAt))
      Struc%PotStrLab(iAt) = Pot_List%StrLab(Struc%PotInd(iAt))
    end do

  end subroutine Gen_Pot_From_XYZ

  subroutine Rule_atomnum(Abs_Index,Struc,Pot_List,Verbose)

    implicit none

    integer(kind=i08),                   intent(in)     :: Abs_Index
    type(XYZ_Format),                    intent(inout)  :: Struc
    type(Pot_List_Type),                 intent(out)    :: Pot_List
    logical,                             intent(in)     :: Verbose

    integer(kind=i08), dimension(:),   allocatable :: Unique_atomnum
    integer(kind=i08)                              :: iUnique, iAt, iPot
    integer(kind=i08), dimension(:),   allocatable :: Indices

! First we allocate the PotTyp vector where the main output will go
    allocate(Struc%PotInd(Struc%nAt))

! Initialize to -1 for sanity checks later
    Struc%PotInd = -1

! Set the absorber to 0
    Struc%PotInd(Abs_Index) = 0

! Get the unique atomic numbers to generate potential type list
    Unique_atomnum = Vec_Unique_Elem(Struc%AN)

! Debug
!   print *, Unique_atomnum

    if (Verbose) then
      write(6,fmt='(a,i3,a)') '# atomnum: Found ', &
                              size(Unique_atomnum), ' unique atomic numbers'
    end if

! Check if there is only one atom of the absorber type to dimension the
! potential list correctly
    Indices = pack( (/ (iAt, iAt=1,Struc%nAt) /), &
                    (Struc%AN == Struc%AN(Abs_Index)) )

    if ( size(Indices) > 1 ) then
      Pot_List%nPot = size(Unique_atomnum)
      if (Verbose) then
        write(6,fmt='(a)') &
          '# atomnum: Multiple instances of absorber type.', &
          '           Consider averaging over absorbers.'
                                 
      end if
    else
      Pot_List%nPot = size(Unique_atomnum) - 1
    end if

    allocate(Pot_List%AN(0:Pot_List%nPot), &
             Pot_List%NumLab(0:Pot_List%nPot), &
             Pot_List%StrLab(0:Pot_List%nPot))

! Initialize the absorber potential in the list
! Set the atomic number corresponding to the absorber type
    Pot_List%AN(0)        = Struc%AN(Abs_Index)

! Set the numeric potential type of the absorber
! NOTE: For now we assign the same as the index
    Pot_List%NumLab(0)   = 0

! Set the label for absorber type
! NOTE: For now we set to the atomic symbol + the "_Abs" string
    Pot_List%StrLab(0) = PTable_Symb(Struc%AN(Abs_Index)) // '_Abs'

! Debug
!   print *, Struc%AN
!   print *, pack( (/ (iAt, iAt=1,Struc%nAt) /), (Struc%AN == 17) )

! Loop over the unique atomic numbers and assign the potentials
    iPot = 0
    do iUnique=1,size(Unique_atomnum)
! Get the indices of the atoms with this unique atomic number
      Indices = pack( (/ (iAt, iAt=1,Struc%nAt) /), &
                      (Struc%AN == Unique_atomnum(iUnique) .and. &
                      (Struc%PotInd /= 0) ) )
      if ( size(Indices) > 0 ) then
        iPot = iPot + 1

! Assign this potential to all the atoms with this atomic number
        Struc%PotInd(Indices) = iPot

! Set the atomic number corresponding to this potential type
        Pot_List%AN(iPot)        = Struc%AN(Indices(1))

! Set the numeric potential type
! NOTE: For now we assign the same as the index
        Pot_List%NumLab(iPot)   = iPot

! Set the label for this potential type
! NOTE: For now we set to the atomic symbol
        Pot_List%StrLab(iPot) = PTable_Symb(Struc%AN(Indices(1)))
      end if
    end do

! Debug
!   print *, Struc%PotInd
!   print *, Pot_List%NumLab
!   print *, Pot_List%AN
!   print *, Pot_List%StrLab

! Sanity check: At this point all values in Struc%PotTyp should have been
!               set to something different than -1, but we check anyway just
!               in case. We might remove this check later on.
    if ( any(Struc%PotInd < 0) ) then
      write(6,fmt='(a)') &
        'Error in Rule_atomnum:', &
        'Not all potentials were properly set.'
        stop
    end if

  end subroutine Rule_atomnum

  function Vec_Unique_Elem(a) result(b)

! Take a vector of integers a and return a vector with the unique elements
! a:  1d array of size n_a
! b:  1d array of size to be determined

! NOTE: This is probably not the most efficient implementation. Should only
!       be used for small vectors.

    implicit none

    integer(kind=i08), dimension(:), intent(in)  :: a
    integer(kind=i08), dimension(:), allocatable :: b

    integer(kind=i08)                     :: nUnique
    logical,           dimension(size(a)) :: Unique_Mask
    integer(kind=i08), dimension(size(a)) :: a_sorted

    integer(kind=i08) :: iEl

! Sort the input vector
    a_sorted = Vec_Shell_Sort(a)

    Unique_Mask = .false.

! Do one pass to find the position of the unique elements
    Unique_Mask(1) = .true.
    do iEl=1,size(a)-1
      if ( a_sorted(iEl) /= a_sorted(iEl+1) ) then
        Unique_Mask(iEl+1) = .true.
      end if
    end do

    b = pack(a_sorted,Unique_Mask)

  end function Vec_Unique_Elem

! Simple Shellsort implementation for vectors of integers
  function Vec_Shell_Sort(Data_Unsorted) result(Data)

    implicit none

    integer(kind=i08), dimension(:), intent(in)           :: Data_Unsorted
    integer(kind=i08), dimension(size(Data_Unsorted))     :: Data

    integer(kind=i08), dimension(size(Data_Unsorted))     :: Swap
    integer(kind=i08) :: iData, jData, nData, Inc
    integer(kind=i08) :: iTemp
    integer(kind=i08)    :: Temp

! Number of data to sort
    nData = size(Data_Unsorted)

! Initialize the array to be sorted
    Data = Data_Unsorted

! Initialize the index swap array
    Swap = (/ (iData,iData=1,nData) /)

    Inc = nint(nData/2.0_r08)
    do
      if ( Inc <= 0 ) exit
      do iData=Inc,nData
        Temp  = Data(iData)
        iTemp = Swap(iData)
        jData = iData
        do
          if ( jData <= Inc ) exit
          if ( Data(jData-Inc) <= Temp ) exit
          Data(jData) = Data(jData-Inc)
          Swap(jData) = Swap(jData-Inc)
          jData = jData - Inc
        end do
        Data(jData) = Temp
        Swap(jData) = iTemp
      end do
      Inc = nint(Inc/2.2_r08)
    end do

  end function Vec_Shell_Sort

end module m_Pot_Generator

