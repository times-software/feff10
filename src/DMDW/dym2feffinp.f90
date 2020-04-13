!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dym2feffinp.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program dym2feffinp

! Program to create a feff.inp file that matches perfectly with the
! structure stored in a .dym file

  use m_Kinds
  use m_Const_and_Conv
  use m_CmdLine
  use m_DMDW

  implicit none

! Program usage
  character(len=*), parameter :: Usage = &
     'Usage: dym2feffinp [Options] dymfile' // achar(10) // &
     ' where dymfile is the name of the file' // &
     ' containing the dynamical matrix and' &
     // achar(10) // &
     '' // achar(10) // &
     'Options:' // achar(10) // &
     '  --c  iAbs   Center the feff.inp file on atom iAbs' // achar(10) // &
!    '              If this option is not present, the program chooses' &
!    // achar(10) // &
     '  --f  fname  Write feff input to file fname' // achar(10) // &
     '  --d  dname  Write adjusted dym file to file dname' // achar(10) // &
     '  --j         Write only ATOMS and POTENTIALS (used with JFEFF)' // achar(10) // &
     '  --s  EXAFS/XANES  Set input to EXAFS or XANES' // achar(10) // &
     ''

  type :: Input_Param
    character(len=mx_FName_Len)    :: dymName
    integer(kind=i08)              :: Center
    character(len=mx_FName_Len)    :: feffName
    character(len=mx_FName_Len)    :: ndymName
    character(len=5)               :: spectrum
    logical                        :: writeHeader
  end type

  type(Cmd_Line)     :: Cmd
  type(Input_Param)  :: Inp_Par

  type(dym_Info)     :: dym_In

  integer(kind=i08), dimension(:), allocatable   :: Swap

! Parse the command line
  call Parse_Cmd_Line(Cmd)
 
! Debug
! call Print_Cmd_Line(Cmd)

! Set the input parameters for this run
  call Get_Input_Param(Cmd,Inp_Par)

! Debug
! call Print_Input_Param(Inp_Par)
! stop
  
! Read the dym file
  call Read_dym_Info(Inp_Par%dymName,dym_In)

! Debug
! call Print_dym_Info(dym_In)

! Check that we have a valid center number
  if (Inp_Par%Center > dym_In%nAt ) then
    write(6,FMT='(A)') Usage
    stop
  end if

! Allocate arrays
  allocate(Swap(dym_In%nAt))

! Write a feff input file based on the dym info and a modified version of
! the dym file that matches the dym file, and return the index swapping list
  call Write_Feffinp(dym_In,Inp_Par%Center,Inp_Par%feffName,Inp_Par%spectrum,Inp_Par%writeHeader, Swap)

  call Write_dym(dym_In,Inp_Par%ndymName,Swap)

  contains  

    subroutine Get_Input_Param(Cmd,Param)

! Convert the command line into the input parameter for the run
! Cmd: Parsed command line
! Param: Parameter for the run

      implicit none

      type(Cmd_Line), intent(in) :: Cmd
      type(Input_Param), intent(out) :: Param

      integer :: iOpt, nOpt, nNoOpt
      character(len=256) :: Buffer

! Must have a last argument, which is the dym filename
      if ( trim(Cmd%Last) == '' ) then
        write(6,FMT='(A)') Usage
        stop
      else
        Param%dymName = trim(Cmd%Last)
      end if

! Check how many options we have
      nOpt = size(Cmd%Opt)

! Check how values without options we have
      nNoOpt = size(Cmd%Val_NoOpt)

! If we have too many options or values without option, stop
!     if ( nOpt > 1 .or. nNoOpt > 0 ) then
!       write(6,FMT='(A)') Usage
!       stop
!     end if

! Set the defaults
! Set center to 1 for now. Need to do something better later
      Param%Center   = 1
      Param%feffName = 'feff.inp'
      Param%ndymName = 'feff.dym'
      Param%spectrum = 'EXAFS'
      Param%writeHeader = .true.

! Check which options we have
      do iOpt=1,nOpt
        select case (trim(Cmd%Opt(iOpt)))

! Option: --c. If present, it must have a value
          case ( '--c' )
! Check that the option has a value attached and read
            if ( trim(Cmd%Val(iOpt)) /= '' ) then
              Buffer = trim(Cmd%Val(iOpt))
              read(Buffer,*) Param%Center
            else
              write(6,FMT='(A)') Usage
              stop
            end if

! Option: --f. If present, it must have a value
          case ( '--f' )
! Check that the option has a value attached and read
            if ( trim(Cmd%Val(iOpt)) /= '' ) then
              Buffer = trim(Cmd%Val(iOpt))
              read(Buffer,*) Param%feffName
            else
              write(6,FMT='(A)') Usage
              stop
            end if

! Option: --j. If present, it must not have a value
          case ( '--j' )
            Param%writeHeader=.false.

! Option: --s. If present, it must have a value
          case ( '--s' )
! Check that the option has a value attached and read
            if ( trim(Cmd%Val(iOpt)) /= '' ) then
              Buffer = trim(Cmd%Val(iOpt))
              read(Buffer,*) Param%spectrum
            else
              write(6,FMT='(A)') Usage
              stop
            end if
            if (Param%spectrum.ne.'XANES') Param%spectrum='EXAFS'

! Option: --d. If present, it must have a value
          case ( '--d' )
! Check that the option has a value attached and read
            if ( trim(Cmd%Val(iOpt)) /= '' ) then
              Buffer = trim(Cmd%Val(iOpt))
              read(Buffer,*) Param%ndymName
            else
              write(6,FMT='(A)') Usage
              stop
            end if

! If we have an unknown option, print usage
          case default
            write(6,FMT='(A,A)') 'Unknown option: ', trim(Cmd%Opt(iOpt))
            stop
        end select
      end do

    end subroutine Get_Input_Param

    subroutine Print_Input_Param(Param)

! Print the input parameters for this run
! Param: Parameters for the run

      implicit none

      type(Input_Param), intent(in) :: Param

! Write out the parameters for the run
      write(6,FMT='(a,i4)')  ' Center:            ', Param%Center
      write(6,FMT='(a,a)')   ' Original dym file: ', trim(Param%dymName)
      write(6,FMT='(a,a)')   ' Feff input file:   ', trim(Param%feffName)
      write(6,FMT='(a,a)')   ' Adjusted dym file: ', trim(Param%ndymName)

    end subroutine Print_Input_Param

  end program dym2feffinp

