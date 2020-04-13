!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_cmdline.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_CmdLine

! Command line operations

! Characters that determine an option
  character(len=*), parameter :: Opt_Flag = '--'

!! Program usage
!  character(len=*), parameter :: Usage = &
!     'Usage: pden [options] Label' // achar(10) // &
!     ' where Label is the label of the density files' // achar(10) // &
!     '' // achar(10) // &
!     'Options:' // achar(10) // &
!     ''
  integer, parameter :: mx_FName_Len = 50
  integer, parameter :: mx_Cmd_Len = 20
  integer, parameter :: mx_Opt_Len = 10
  integer, parameter :: mx_Val_Len = 100
  integer, parameter :: mx_Arg_Len = 100
  integer, parameter            :: mx_Line_Len = 300
  character(len=*), parameter   :: Fmt_Line = '(A300)'
  integer, parameter :: mx_Field_Len = 50
  integer, parameter :: Sys_Cmd_Buffer_Len = 100
  integer, parameter :: Label_Len = 40
  integer, parameter :: mx_Fields_per_Line = 40

  type :: Cmd_Line
    character(len=mx_Cmd_Len)                            :: Cmd
    character(len=mx_Opt_Len), dimension(:), allocatable :: Opt
    character(len=mx_Val_Len), dimension(:), allocatable :: Val
    character(len=mx_Val_Len), dimension(:), allocatable :: Val_NoOpt
    character(len=mx_Val_Len)                            :: Last
  end type

  private
  public :: mx_FName_Len
  public :: Cmd_Line
  public :: Parse_Cmd_Line
  public :: Print_Cmd_Line

contains

  subroutine Parse_Cmd_Line(Line)

! Parse the command line
! Line: parsed command line

    implicit none

    type(Cmd_Line), intent(out) :: Line

    integer :: iArg, nArg, iOpt, nOpt, nNoOpt
    character(len=mx_Arg_Len), dimension(:), allocatable :: Arg
    logical, dimension(:), allocatable :: Mask_Opt, Mask_Val, Mask_NoOpt_Val
    integer, dimension(:), allocatable :: Ind_Opt, Ind_Opt_Val, Ind_NoOpt_Val
    integer, dimension(:), allocatable :: Ind_Tmp
    integer :: i

   !interface
   !   integer function iargc()
   !   end function iargc
   !end interface

! Number of command line arguments
! This function is not f90 standard
    !nArg = iargc()
    nArg = COMMAND_ARGUMENT_COUNT()

! Allocate and read the arguments
    allocate(Arg(0:nArg))
    do iArg=0,nArg
      !call getarg(iArg,Arg(iArg))
      call GET_COMMAND_ARGUMENT(iArg,Arg(iArg))
    end do

! Start parsing the command line
! Command name
    Line%Cmd = Arg(0)

! If the last argument is not an option, make it special
    Line%Last = ''
    if ( (Arg(nArg)(1:len(Opt_Flag)) /= Opt_Flag) .and. (nArg > 0) ) then
      Line%Last = Arg(nArg)
      nArg = nArg - 1
    end if

! Find number of options
    allocate(Mask_Opt(nArg))
    Mask_Opt = (Arg(1:nArg)(1:len(Opt_Flag)) == Opt_Flag)
    nOpt = count(Mask_Opt)

! Find position of the options
    allocate(Ind_Opt(nOpt))
    Ind_Opt = pack((/(i,i=1,nArg)/),Mask_Opt)

! Find position of the option values
    allocate(Ind_Tmp(nOpt))
    Ind_Tmp = Ind_Opt + 1
    allocate(Mask_Val(nOpt))
    Mask_Val = ( (Ind_Tmp>=1)    .and. &
                 (Ind_Tmp<=nArg) .and. &
                 .not. Mask_Opt(Ind_Tmp) )
    allocate(Ind_Opt_Val(count(Mask_Val)))
    Ind_Opt_Val = pack(Ind_Tmp,Mask_Val)

! Find the values not associated with an option
    allocate(Mask_NoOpt_Val(nArg))
    Mask_NoOpt_Val = .true.
    Mask_NoOpt_Val(Ind_Opt) = .false.
    Mask_NoOpt_Val(Ind_Opt_Val) = .false.
    nNoOpt = count(Mask_NoOpt_Val)
    allocate(Ind_NoOpt_Val(nNoOpt))
    Ind_NoOpt_Val = pack((/(i,i=1,nArg)/),Mask_NoOpt_Val)

! Debug
!   write(6,FMT='(I5)') Ind_Opt
!   write(6,FMT='(A)') ''
!   write(6,FMT='(I5)') Ind_Opt_Val
!   write(6,FMT='(A)') ''
!   write(6,FMT='(I5)') Ind_NoOpt_Val
    
! Store the options and values
    allocate(Line%Opt(nOpt))
    Line%Opt = ''
    do iOpt=1,size(Ind_Opt)
      Line%Opt(iOpt) = Arg(Ind_Opt(iOpt))
    end do
    allocate(Line%Val(nOpt))
    Line%Val = ''
    do iOpt=1,size(Ind_Opt_Val)
      Line%Val(iOpt) = Arg(Ind_Opt_Val(iOpt))
    end do
    allocate(Line%Val_NoOpt(nNoOpt))
    Line%Val_NoOpt = ''
    do iOpt=1,size(Ind_NoOpt_Val)
      Line%Val_NoOpt(iOpt) = Arg(Ind_NoOpt_Val(iOpt))
    end do

  end subroutine Parse_Cmd_Line

  subroutine Print_Cmd_Line(Line)

! Line: parsed command line

    implicit none

    type(Cmd_Line), intent(in) :: Line
    integer :: iOpt

! Command name
    write(6,FMT='(A,A)')  'Command:       ', trim(Line%Cmd)

! Options with their values
    do iOpt=1,size(Line%Opt)
      write(6,FMT='(A,A,3X,A)') &
            'Option, Value: ', trim(Line%Opt(iOpt)), trim(Line%Val(iOpt))
    end do
! Values without an option
    do iOpt=1,size(Line%Val_NoOpt)
      write(6,FMT='(A,A)') &
            'Value (No Opt): ', trim(Line%Val_NoOpt(iOpt))
    end do

! Last value
    write(6,FMT='(A,A)')  'Last Arg:      ', trim(Line%Last)

  end subroutine Print_Cmd_Line

end module m_CmdLine
