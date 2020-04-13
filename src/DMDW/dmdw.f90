!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dmdw.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program dmdw

! Test wrapper for the dmdw library routines

  use m_Kinds
  use m_Const_and_Conv
  use m_DMDW
#ifdef FEFF
 use ff2x_inp,only: ff2x_read,idwopt
 use errorfile
 use par
#endif

  implicit none

  type(Lanczos_Info) :: Lanc_In
  type(Paths_Info)   :: Paths_In
  type(dym_Info)     :: dym_In
! type(dym_Info)     :: dym_In_WEnv

  real(kind=r08), dimension(:), allocatable :: sigc
  real(kind=r08), dimension(:), allocatable :: vfe,mef

  type(Error_Info)   :: Err

#ifdef FEFF
 call par_begin
 if(worker) goto 400
 call OpenErrorfileAtLaunch('dmdw')
 call ff2x_read
 if(idwopt.ne.5) goto 402
 call wlog('Calculating Debye-Waller factors in the Dynamical Matrix approach ...')
!KJ 11-2011 Added following section to mimic CONTROL card:
#endif 

! Open IO units
  call DMDW_Open_I(Err)

  if ( Err%Flag ) then
    write(6,fmt='(a)') trim(Err%Message)
    stop
  end if
  
! call wlog('Using dynamical matrix DW factors.') !KJ

  call DMDW_Open_O
  call DMDW_Open_E

! Read the input file
  call Read_Lanczos_Info(Lanc_In)

! Read the paths info
! Modified by Fer to swith computation of the VFE from the old way to the
! new one.
! if ( (Lanc_In%RunTyp .ne. 1) .and. (Lanc_In%RunTyp .ne. 2) ) then
  if ( Lanc_In%RunTyp .ne. 2 ) then
    call Read_Paths_Info(Paths_In)
  else
! SE/SF option doesn't use paths --> fill with dummy path
    Paths_In%nDesc = 1

! Allocate the paths info arrays
    allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
             Paths_In%Desc(Paths_In%nDesc,0:2-1), &
             Paths_In%Desc_mxR(Paths_In%nDesc))

! Use single dummy path: "2 1 0     30.0"
    Paths_In%Desc_Len(1) = 2
    Paths_In%Desc(1,0) = 1
    Paths_In%Desc(1,1) = 0
    Paths_In%Desc_mxR(1) = 30.0*ang2au
  end if

! Read the dym file data
  call Read_dym_Info(Lanc_In%dym_file,dym_In)

! Debug
! call Print_dym_Info(dym_In)

! Debug
! Here we process the dym input, for now we will test some simple things
! like  adding a shell of cells by replication of the input cell, if we have
! the information.
! NOTE: This seems to be unnecesarry with the way VASP prints out its
! dynamical matrices. We leave the code here for now, for use in the future.
! call Add_Cell_Env(dym_In,dym_In_WEnv)
! dym_In = dym_In_WEnv

! Create the full dynamical matrix
  call Make_DM(dym_In)
! call Make_DM(dym_In_WEnv)

! Calculate the transformation to internal coordinates (eliminates
! rotations and translations)
! NOTE: The subroutine is still missing parts of the code and shouldn't be used.
!       If used, the TrfD is deallocated before return and any attempt to use it
!       will result in a segfault.
  call Make_TrfD(dym_In)

! Debug
! call Print_dym_Info(dym_In)
! call Print_dym_Info(dym_In_WEnv)
! stop

! Print out a header
  call Print_Header(Lanc_In,Paths_In,dym_In)

! Proceed depending on the type of run selected
  select case ( Lanc_In%RunTyp )

! EXAFS Debye-Waller factors (s2)
    case ( 0 )

      call RunTyp_S2(Lanc_In,dym_In,Paths_In)

! Crystallographic Debye-Waller factors (u2)
    case ( 3 )

      call RunTyp_U2(Lanc_In,dym_In,Paths_In)

! Vibrational Free Energy
    case ( 1 )

      call RunTyp_VFE(Lanc_In,dym_In,Paths_In)

! Self-Energy and related quantities
    case ( 2 )

      call RunTyp_SE(Lanc_In,dym_In,Paths_In)

! IR spectra based on Gaussian's dipole derivatives
    case ( 4 )

      call RunTyp_IR(Lanc_In,dym_In,Paths_In)

! Compute partial and total phonon densities of states
    case ( 5 )

! Debug
!     print *, Paths_In%nDesc
!     print *, allocated(Paths_In%Desc_Len), Paths_In%Desc_Len
!     print *, allocated(Paths_In%Desc), Paths_In%Desc
!     print *, allocated(Paths_In%Desc_mxR), Paths_In%Desc_mxR
      call RunTyp_PDOS(Lanc_In,dym_In,Paths_In)

! If we get here it means we don't have this property yet
    case default

      write(IO_Err,fmt='(a)') &
        'Unrecognized calculation type. Allowed values are:'
      write(IO_Err,fmt='(a)') &
        ' 0: EXAFS Debye-Waller factors (s^2)'
      write(IO_Err,fmt='(a)') &
        ' 1: Projected vibrational free energies'
      write(IO_Err,fmt='(a)') &
        ' 2: Phonon Self-Energy and related quantities'
      write(IO_Err,fmt='(a)') &
        ' 3: Crystallographic Debye-Waller factors (u^2)'
!     write(IO_Err,fmt='(a)') &
!       ' 4: IR spectra (no yet implemented)'
      write(IO_Err,fmt='(a)') &
        ' 5: Partial and total phonon density of states (beta)'

  end select

!#############################################################################

! Close the IO units
  call DMDW_Close_I
  call DMDW_Close_O
  call DMDW_Close_E
#ifdef FEFF
401 continue   ! exit point for master if no need to do DMDW
call wlog('Done with module: Debye-Waller factors (DMDW).'//char(13)//char(10))
402 continue   ! DMDW skipped
call WipeErrorfileAtFinish
400 continue !exit point for worker
call par_end
#endif
  stop

  end program dmdw

