module m_Strings

! Define some useful string manipulation routines

  private
  public :: Uppercase, Lowercase, Is_AllNum, Is_AllAlpha
  public :: pad

contains

  function Uppercase(String1) result(String2)

    implicit none

    character(len=*)            :: String1
    character(len=len(String1)) :: String2

    integer :: iChar

    do iChar=1,len(String1)
      String2(iChar:iChar) = String1(iChar:iChar)
      if ( 'a' <= String1(iChar:iChar) .and. String1(iChar:iChar) <= 'z' ) then
        String2(iChar:iChar) = achar(iachar(String1(iChar:iChar)) - 32)
      end if
    end do

  end function Uppercase

  function Lowercase(String1) result(String2)

    implicit none

    character(len=*)            :: String1
    character(len=len(String1)) :: String2

    integer :: iChar

    do iChar=1,len(String1)
      String2(iChar:iChar) = String1(iChar:iChar)
      if ( 'A' <= String1(iChar:iChar) .and. String1(iChar:iChar) <= 'Z' ) then
        String2(iChar:iChar) = achar(iachar(String1(iChar:iChar)) + 32)
      end if
    end do

  end function Lowercase

  function Is_AllNum(String) result(Answer)
  
    implicit none

    character(len=*) :: String
    logical :: Answer
    
    integer :: iChar

    Answer = .true.
    do iChar=1,len(String)
      if ( .not. ( string(iChar:iChar) == ' ' .or. &
                   ( '0' <= String(iChar:iChar) .and. &
                            String(iChar:iChar) <= '9' ) ) ) then
        Answer = .false.
      end if
    end do

  end function Is_AllNum

  function Is_AllAlpha(String) result(Answer)
  
    implicit none

    character(len=*) :: String
    logical :: Answer

    integer :: iChar

    Answer = .true.
    do iChar=1,len(String)
      if ( .not. ( String(iChar:iChar) == ' ' .or. &
                   ( 'a' <= String(iChar:iChar) .and. &
                            String(iChar:iChar) <= 'z' ) .or. &
                   ( 'A' <= String(iChar:iChar) .and. &
                            String(iChar:iChar) <= 'Z' ) &
                   ) ) then
        Answer = .false.
      end if
    end do

  end function Is_AllAlpha

  function pad(nPad,Num) result(String)

    implicit none

    integer, intent(in) :: nPad, Num
    character(len=nPad) :: String

    integer :: NumLen, iChar
! Define a temp string to hold the working padded number. Since in principle
! it has to hold a very large integer, we make it the maximum length of an
! integer + the padding
    character(len=32+nPad) :: String_Tmp

    NumLen = int(log10(real(Num))) + 1
    if ( NumLen > nPad ) then
      print *, 'Error in pad: Number is bigger than padding size'
      stop
    end if
    write(String_Tmp,*) Num
    String_Tmp = adjustl(String_Tmp)
    do iChar=1,nPad-NumLen
      String_Tmp = "0" // adjustl(String_Tmp)
    end do
    String = adjustl(String_Tmp)

  end function pad

end module m_Strings

