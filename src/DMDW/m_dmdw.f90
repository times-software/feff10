! $Revision: 1.8 $
! $Author: fer $
! $Date: 2010/10/26 18:07:13 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_dmdw.f90,v $:
! $Revision: 1.8 $
! $Author: srw $
! $Date: 2011/9/29 $
! $Comments: Self energy added
!
! $Date: 2011/11/16
! $Comments: Spectral function added
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_DMDW

  use m_Kinds
  use m_Const_and_Conv
  use m_PTable
  use m_Math
  use m_Strings
  
  integer, parameter :: IO_In=80, IO_Out=81, IO_dym=89, IO_Feff=90
  integer, parameter :: IO_PDOS=82
  integer, parameter :: IO_Err=6, IO_Dbg=6
  
  type :: Lanczos_Info
      integer(kind=i08)           :: IOFlag
      integer(kind=i08)           :: nPoles
      integer(kind=i08)           :: nT
      integer(kind=i08)           :: E_k_opt
      real(kind=r08)              :: E_k
      real(kind=r08), allocatable :: T(:)
      integer(kind=i08)           :: RunTyp
      integer(kind=i08)           :: disp_opt !srw
      character(len=100)          :: dym_file
      character(len=100)          :: pds_file !sms: can use full path to files
      character(len=100)          :: a2f_file
      logical                     :: PDOS_Poles = .false.
      logical                     :: PDOS_Rect  = .false.
      logical                     :: PDOS_Gauss = .false.
      logical                     :: PDOS_DropL = .false.
      logical                     :: PDOS_Part  = .false.
      real(kind=r08)              :: PDOS_Broad
      real(kind=r08)              :: PDOS_Res
  end type Lanczos_Info
  
  type :: Paths_Info
      integer(kind=i08) :: nDesc
      integer(kind=i08), dimension(:), allocatable   :: Desc_Len
      integer(kind=i08), dimension(:,:), allocatable :: Desc
      real(kind=r08), dimension(:), allocatable      :: Desc_mxR
  end type Paths_Info
  
  type :: dym_Info
      integer(kind=i08)              :: Type
      integer(kind=i08)              :: nAt
      integer(kind=i08), allocatable :: an(:)
      real(kind=r08), allocatable    :: am(:)
      real(kind=r08), allocatable    :: xyz(:,:)
      real(kind=r08), allocatable    :: red(:,:)
      real(kind=r08), allocatable    :: cell(:,:)
      real(kind=r08), allocatable    :: R_CM(:)
      real(kind=r08), allocatable    :: ToI(:,:)
      real(kind=r08), allocatable    :: MoI(:)
      real(kind=r08), allocatable    :: PAoR(:,:)
      real(kind=r08), allocatable    :: dm_block(:,:,:,:)
      real(kind=r08), allocatable    :: dm(:,:)
      real(kind=r08), allocatable    :: TrfD(:,:)
      real(kind=r08), allocatable    :: dm_Int(:,:)
      real(kind=r08), allocatable    :: Dip_Der(:,:,:)
      integer(kind=i08), allocatable :: BasAtInd(:)
      real(kind=r08), allocatable    :: u_xyz(:,:,:)!srw
      integer(kind=i08)              :: unAt,natom !srw
      integer(kind=i08), allocatable :: udeg(:),utype(:)!srw 
      integer(kind=i08), allocatable :: centeratomindex(:,:)
  end type dym_Info
  
  type :: DW_Out_Info
      integer(kind=i08)              :: Path_nAt
      integer(kind=i08), allocatable :: Path(:)
      integer(kind=i08)              :: Path_Pert
      real(kind=r08)                 :: Path_Len
      real(kind=r08), allocatable    :: s2(:)
      real(kind=r08), allocatable    :: u2(:)
      real(kind=r08), allocatable    :: vfe(:)
      real(kind=r08), allocatable    :: Poles_Frq(:)
      real(kind=r08), allocatable    :: Poles_Wgt(:)
      real(kind=r08)                 :: SPole_Frq
      real(kind=r08)                 :: RedMass
  end type DW_Out_Info
  
  type :: Paths
      integer(kind=i08)                              :: N
      integer(kind=i08), dimension(:,:), allocatable :: Ind
      real(kind=r08),    dimension(:), allocatable   :: Len
  end type Paths
  
  integer, parameter :: Error_Message_Len = 200
  type :: Error_Info
      logical                          :: Flag
      character(len=Error_Message_Len) :: Message
  end type Error_Info
  
  character          :: u2Pert_Label(0:2) = (/ 'x', 'y', 'z' /)

  integer, parameter :: FName_Len = 300

  private
  public :: Lanczos_Info, Paths_Info, dym_Info, Paths
  public :: IO_In, IO_Out, IO_Err
  public :: DMDW_Open_I, DMDW_Open_O, DMDW_Open_E
  public :: DMDW_Close_I, DMDW_Close_O, DMDW_Close_E
  public :: Calc_DW
  public :: Read_Lanczos_Info
  public :: Read_Paths_Info
  public :: Read_dym_Info, Print_dym_Info
  public :: Print_Header
  public :: Paths_Init, Paths_DeInit
  public :: Make_DM, Make_TrfD
  public :: Write_Feffinp, Write_dym
  public :: Error_Info
  public :: Add_Cell_Env
  public :: RunTyp_S2, RunTyp_U2, RunTyp_VFE
  public :: RunTyp_SE, RunTyp_IR, RunTyp_PDOS
  public :: DW_Out_Info

contains

  subroutine DMDW_Open_I(Err)
      implicit none
      type(Error_Info)   :: Err
      logical :: File_Present

      inquire(file="dmdw.inp",exist=File_Present)
      if ( File_Present ) then
          open(unit=IO_In,file='dmdw.inp',status='unknown')
          Err%Flag = .false.
          Err%Message = ""
      else
          Err%Flag = .true.
          Err%Message = "dmdw.inp: File not found"
      end if
  end subroutine DMDW_Open_I

  subroutine DMDW_Open_O
      implicit none
      open(unit=IO_Out,file='dmdw.out',status='unknown')
  end subroutine DMDW_Open_O

  subroutine DMDW_Open_E
      implicit none
      if ( .not. (IO_Err .eq. 6) ) then
        open(unit=IO_Err)
      end if
  end subroutine DMDW_Open_E

  subroutine DMDW_Open_PDOS(Label)

      implicit none

      character(len=*),   intent(in) :: Label

      character(len=FName_Len) :: PDOS_FName

! Here we make the right filename to print out the PDOS
      PDOS_FName = 'dmdw_pdos.' // trim(Label) // '.dat'

      open(unit=IO_PDOS,file=trim(PDOS_FName),status='unknown')

  end subroutine DMDW_Open_PDOS

  subroutine DMDW_Close_I
      implicit none
      close(unit=IO_In)
  end subroutine DMDW_Close_I

  subroutine DMDW_Close_O
      implicit none
      close(unit=IO_Out)
  end subroutine DMDW_Close_O

  subroutine DMDW_Close_E
      implicit none
! NOTE: FDV
! This if ensures that the standard IO is not "closed" avoiding the generation
! of fort.6 files under certain circumstances
      if ( .not. (IO_Err .eq. 6) ) then
        close(unit=IO_Err)
      end if
  end subroutine DMDW_Close_E

  subroutine DMDW_Close_PDOS
      implicit none
      close(unit=IO_PDOS)
  end subroutine DMDW_Close_PDOS

  subroutine Calc_DW(Lanc_In, dym_In, lpath_in, nleg, itt, iu2Pert, DW_Out)
 
    implicit none

    integer(kind=i08), parameter :: numxp=100000

    type(Lanczos_Info), intent(in)  :: Lanc_In
    type(dym_Info),     intent(in)  :: dym_In
    integer(kind=i08),  intent(in)  :: lpath_in(:)
    integer(kind=i08),  intent(in)  :: nleg
    integer(kind=i08),  intent(in)  :: itt
    integer(kind=i08),  intent(in)  :: iu2Pert
    type(DW_Out_Info),  intent(out) :: DW_Out

    real(kind=r08), dimension(Lanc_In%nT)     :: sig2
    real(kind=r08), dimension(Lanc_In%nT)     :: t_dF
!   real(kind=r08), dimension(Lanc_In%nT)     :: mef

    integer(kind=i08) :: inull
    real(kind=r08), dimension(:), allocatable :: beta, cotarg, cth, sig
    real(kind=r08), dimension(Lanc_In%nT) :: s_dF
    real(kind=r08) ::    dr(3) 
    real(kind=r08) ::    mu_inv 
    real(kind=r08) ::    xp, pnx, ratio

    integer(kind=i08) :: m, l
    integer(kind=i08) :: ni, mnull

    integer(kind=i08) :: ip,iq,jq,iAt,iDeg

    real(kind=r08)    :: u2Pert(0:2)
    integer(kind=i08) :: iAt_Tgt, iVec
    real(kind=r08)    :: SsMassj
    real(kind=r08)    :: Dummy
    logical           :: Project = .true.
!   logical           :: Project = .false.

    integer(kind=i08) :: nn,mm,ii,jj
    integer(kind=i08), dimension(:),  allocatable :: lpath, lpath_p, lpath_m
    real(kind=r08)    :: mu
    real(kind=r08), dimension(:,:,:), allocatable :: rc
    real(kind=r08), dimension(:),     allocatable :: qj0
    real(kind=r08), dimension(:),     allocatable :: w_pole,wil
    real(kind=r08), dimension(:),     allocatable :: wil_tmp,w_pole_tmp
    real(kind=r08) :: SPole_EinsteinFreq, W_Scal

!srw    
    complex,dimension(Lanc_In%nT)                     :: SE_a2f
    real,dimension(0:Lanc_In%nPoles-1,2)              :: a2f
    real,dimension(2,400)                             :: a2,a2fall
    real(kind=r08)                                    :: cf,norm,check,Gam,betashift,tot,width,xx,yy,Zk,nk,wh,wp
    complex                                           :: omega,deli,coni,cheightL,cheightR,integral,atot
    real(kind=r08)                                    ::  lowE, highE,E0,dE,  Ek, w0, t0,dt,t
    integer                                           :: ilowE,ihighE,NE,iE, iEk,     Nt,it  
    real(kind=r08),dimension(-50000:50000)            :: w,reSE,imSE,betak       
    complex(kind=r08),dimension(-50000:50000)         :: cIntegrand,Ckt,Akt,Akw,AkwQP,AktS,AkwS

    character(len=25), parameter           :: FMT_A="(A)",FMT_AFAF="(A,F10.3,A,F10.3)",FMT_AIAI="(A,I6,A,I6)"
    character(len=25), parameter           :: FMT_AFAFAI="(A,F8.3,A,F8.3,A,I6)"

    integer(kind=i08) :: iStat

    logical,parameter                      :: DEBUG = .FALSE.
    ! unit numbers assigned to output files
    integer,parameter                      :: u_PSinfo=40,u_Einfo=41,u_ReSE=42,u_ImSE=43,u_Akw=44
    ! unit numbers assigned to debug files
    integer,parameter                      :: u_beta=50,u_Ckt=51,u_Akt=52,u_Aqp=53,u_Asat=54
    !constants 
    real                                   :: pi,kb
      parameter(pi=3.141592653589793e0, kb=8.6173342e-5)

    allocate(sig(Lanc_In%nT), &
        beta(Lanc_In%nT),     &
        cotarg(Lanc_In%nT),   &
        cth(Lanc_In%nT))
         
    allocate(lpath(0:nleg-1),                           &
             lpath_m(0:nleg-1),                         &
             lpath_p(0:nleg-1),                         &
             rc(0:dym_In%nAt-1,0:dym_In%nAt-1,0:2),     &
             qj0(0:3*dym_In%nAt-1),                     &  
             w_pole(0:Lanc_In%nPoles),                  &
             wil(0:Lanc_In%nPoles),                     &
             w_pole_tmp(0:Lanc_In%nPoles*dym_In%natom), &
             wil_tmp(0:Lanc_In%nPoles*dym_In%natom)       )  !*dym_In%natom or **dym_In%unAt ?     

! Debug
!   print *, 'Entering Calc_DW'

! Adjust the lpath_in indices to the ones used locally (from 0)
    lpath(0:nleg-1) = lpath_in-1

! Debug
!   print *, 'lpath: ', lpath(0:nleg-1)

    beta = 1.0_r08/(K_B*Lanc_In%T)
    cotarg = 0.5_r08*hbar*beta
!##############################################################################
! Initialize director cosines
  rc = 0.0_r08

! Calculate distances and director cosines
  do m = 0, dym_In%nAt-1
    do l = m+1, dym_In%nAt-1
      dr = dym_In%xyz(l+1,:) - dym_In%xyz(m+1,:)
      rc(l,m,:) = dr/sqrt(sum(dr**2))
      rc(m,l,:) = -rc(l,m,:)
    end do
  end do

! Create the shifted path vectors
  lpath_p = cshift(lpath,+1)
  lpath_m = cshift(lpath,-1)

! Modified by FDV
! If we have a single atom path, then we need a different reduced mass.
! (It should be noted that the actual value of the reduced mass on the
! determination of the q0 is irrelevant in this case.)

! Debug (FDV)
! print *, 'Got here 1'
  if ( nleg .eq. 1 ) then

! Debug
    mu_inv = 1.0_r08/dym_In%am(lpath(0)+1)

  else

! Calculate the reduced mass for the path
    mu_inv = 0.0_r08

! Debug FDV
! print *, 'nleg ', nleg
    do l=0,nleg-1   
!     print *, 'l         ', l
!     print *, 'lpath(l)  ', lpath(l)
!     print *, 'lpath_m(l)', lpath_m(l)
!     print *, 'lpath_p(l)', lpath_p(l)
      ! sum of (r_i_i-1 + r_i_i+1)**2
      mu_inv = mu_inv + &
               sum((rc(lpath(l),lpath_m(l),:) + &
               rc(lpath(l),lpath_p(l),:))**2)/(4.0_r08*dym_In%am(lpath(l)+1))
    end do 

  end if

! Debug (FDV)
! print *, 'Got here 2'
  if ( mu_inv <= 0.0_r08 ) then
    write(IO_Err,fmt='(a,e16.10)') ' Problem in inverse reduced mass, mu_inv = ', mu_inv
    close (unit=94)
    stop
  end if
  mu = 1.0_r08/mu_inv

! Debug (FDV)
! print *, 'Got here 3'
! construct path vector
  ! initialize the path vector; zero qj
  qj0 = 0.0_r08

  ! srw:get vector |0> = qj
  If( Lanc_In%RunTyp.Eq.0 ) then!you've elected to calculate sig2 
    do l=0,nleg-1
      qj0( (/ lpath(l),                &
              lpath(l) +   dym_In%nAt, &
              lpath(l) + 2*dym_In%nAt    /) ) = &
            0.5_r08*dsqrt(mu/dym_In%am(lpath(l)+1)) * &
              (rc(lpath(l),lpath_m(l),:) + rc(lpath(l),lpath_p(l),:))
    end do

! TEST CODE: Project out the
! Projecting out the translation and rotation components
!   open(1)
!   read(1,*) Project
!   write(6,fmt='(a,l)') ' Projecting?: ', Project
    if ( Project ) then
    do iVec=1,6
      qj0 = qj0 - sum(qj0*Dym_In%TrfD(:,iVec))*Dym_In%TrfD(:,iVec)
    end do
    end if
!   close(1)

! Normalize the seed
    qj0 = qj0/sqrt(sum(qj0**2))

! Debug
! call Print_dym_Info(dym_In)

!Get PHDOS
    Call Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                 w_pole,wil,mnull,SPole_EinsteinFreq)
     
! Added by FDV
! Adding more general support for u2 and pDOS than the stuff used below for the
! phonon spectral function stuff.
! NOTE: Might have to make this more compatible in the future
  else if ( Lanc_In%RunTyp == 1 .or. &
            Lanc_In%RunTyp == 3 .or. &
            Lanc_In%RunTyp == 5        ) then

! Debug
! Do some temporary tests to make sure that we have what wee need.
!   write(6,fmt='(3f12.8)') transpose(dym_In%red)
!   write(6,fmt='(3f12.8)') transpose(dym_In%xyz)
!   write(6,fmt='(f16.8)') dym_In%dm(0,:)
!   write(6,fmt='(3f12.8)') dym_In%R_CM
!   write(6,fmt='(3f16.8)') dym_In%ToI
!   print *, nleg
!   print *, size(lpath)
!   print *, lpath(0)
!   print *, allocated(dym_In%TrfD)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,1)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,2)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,3)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,4)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,5)
!   write(6,fmt='(f12.8)') dym_In%TrfD(:,6)
!   print *, sum(dym_In%TrfD(:,3)*dym_In%TrfD(:,4))
!   print *, sum(dym_In%TrfD(:,3)*dym_In%TrfD(:,5))
!   print *, sum(dym_In%TrfD(:,3)*dym_In%TrfD(:,6))
!   stop

! If we got to this point, we have just one atom in the path, so we define
! some variables to make things easier:
    iAt_Tgt = lpath(0)+1

! Debug
!   print *, dym_In%am(iAt_Tgt)
!   stop

! NOTE: Here we use the masses in AMU, rather than au. This is because it only
!       introduces a scaling factor that is irrelevant.
    SsMassj = 0.0_r08
    do iAt=1,dym_In%nAt
      if ( iAt /= iAt_Tgt ) then
        SsMassj = SsMassj + sqrt(dym_In%am(iAt))
      end if
    end do

! Set the tgt atom displacement to 1, as before
    u2Pert = (/ 0.0_r08, 0.0_r08, 0.0_r08 /)
    u2Pert(iu2Pert) = 1.0_r08
    qj0( (/ iAt_Tgt-1,                &
            iAt_Tgt-1 +   dym_In%nAt, &
            iAt_Tgt-1 + 2*dym_In%nAt    /) ) =  u2Pert

! To preserve the CM, also displace the rest of the atoms
!   u2Pert(iu2Pert) = -sqrt(dym_In%am(iAt_Tgt))/SsMassj
!   do iAt=1,dym_In%nAt
!     if ( iAt /= iAt_Tgt ) then
!       qj0( (/ iAt-1,                &
!               iAt-1 +   dym_In%nAt, &
!               iAt-1 + 2*dym_In%nAt    /) ) =  u2Pert
!     end if
!   end do

! TEST CODE: Project out the
! Projecting out the translation and rotation components
!   open(1)
!   read(1,*) Project
!   write(6,fmt='(a,l)') ' Projecting?: ', Project
    if ( Project ) then
    do iVec=1,6
      qj0 = qj0 - sum(qj0*Dym_In%TrfD(:,iVec))*Dym_In%TrfD(:,iVec)
    end do
    end if
!   close(1)

! Debug
!   do l=0,dym_In%nAt*3-1
!     print *, qj0(l)
!   end do
!   stop

! Debug
! Test if we are preserving the CM, just in case I screwed up in my derivation
! of the displacements
!   Dummy = 0.0_r08
!   do iAt=1,dym_In%nAt
!     Dummy = Dummy + qj0(iAt-1)*sqrt(dym_In%am(iAt))
!   end do
!   print *, Dummy
!   stop

! This is the old code that doesn't preserve the CM
! Define the direction of the perturbation
!   u2Pert = (/ 0.0_r08, 0.0_r08, 0.0_r08 /)
!   u2Pert(iu2Pert) = 1.0_r08
!   do l=0,nleg-1
!     qj0( (/ lpath(l),                &
!             lpath(l) +   dym_In%nAt, &
!             lpath(l) + 2*dym_In%nAt    /) ) =  u2Pert
!   end do

! This is probably a broken attempt at trying to preserve the CM. Please ignore
!   do iAt=1,dym_In%nAt
!     u2Pert = (/ 0.0_r08, 0.0_r08, 0.0_r08 /)
!     if ( iAt-1 == lpath(0) ) then
!       u2Pert(iu2Pert) = 1.0_r08/dym_In%am(iAt)
!     else
!       u2Pert(iu2Pert) = -1.0_r08/dym_In%am(iAt)/(dym_In%nAt-1)
!     end if
!     qj0( (/ iAt-1,                &
!             iAt-1 +   dym_In%nAt, &
!             iAt-1 + 2*dym_In%nAt    /) ) =  u2Pert
!   end do

! Normalize the seed
    qj0 = qj0/sqrt(sum(qj0**2))

! Debug
! This is overkill, since I konw it is ok, but: Test if we are preserving
! the CM, after the renormalization
!   Dummy = 0.0_r08
!   do iAt=1,dym_In%nAt
!     Dummy = Dummy + qj0(iAt-1)*sqrt(dym_In%am(iAt))
!   end do
!   print *, Dummy
!   stop

! Debug
!   do l=0,dym_In%nAt*3-1
!     print *, qj0(l)
!   end do
!   stop

!Get PHDOS
    Call Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                 w_pole,wil,mnull,SPole_EinsteinFreq)
     
! Added by FDV
! Using a different seed, based on the dipole derivatives, to compute the
! infrared spectrum.
  else If( Lanc_In%RunTyp.Eq.4 ) then

! Inside here, the coordinates of the q0 are ordered like:
! Atom:    0 1 2 3 ... 0 1 2 3 ... 0 1 2 3 ...
! Coord:   0 0 0 0 ... 1 1 1 1 ... 2 2 2 2 ...
! or:      x x x x ... y y y y ... z z z z ...
! So, we have to be very careful with the assignement of the dipole
! derivatives, which are ordered as:
! Atom:    0 0 0 1 1 1 2 2 2 ...
! Coord:   0 1 2 0 1 2 0 1 2 ...
! or:      x y z x y z x y z ...

      do iAt=1,dym_In%nAt
        do ip=1,3
          write(IO_Dbg,fmt='(2i5,1p,3e18.6)') iAt, ip, &
                (dym_In%Dip_Der(jq,iAt,ip),jq=1,3)
        end do
      end do

    do iAt=1,dym_In%nAt
      do ip=1,3
!       qj0(iAt-1+(ip-1)*dym_In%nAt) = &
!         dym_In%Dip_Der(1,iAt,ip)**2
!       qj0(iAt-1+(ip-1)*dym_In%nAt) = &
!         dsqrt(dym_In%am(iAt))*dym_In%Dip_Der(1,iAt,ip)**2
        qj0(iAt-1+(ip-1)*dym_In%nAt) = &
          dsqrt(dym_In%am(iAt))*dym_In%Dip_Der(2,iAt,ip)**2
!       qj0(iAt-1+(ip-1)*dym_In%nAt) = &
!         dsqrt(dym_In%am(iAt))*dym_In%Dip_Der(3,iAt,ip)**2
      end do
    end do

! Debug
    print *, sum(dym_In%Dip_Der(1,:,:)**2)
    print *, sum(dym_In%Dip_Der(2,:,:)**2)
    print *, sum(dym_In%Dip_Der(3,:,:)**2)

! Normalize the seed
    qj0 = qj0/sqrt(sum(qj0**2))

! Debug
!   do ip=0,3*dym_In%nAt-1
!     write(IO_Dbg,fmt='(i5,1p,e18.6)') ip, qj0(ip)
!   end do
!   stop

!Get PHDOS
    Call Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                 w_pole,wil,mnull,SPole_EinsteinFreq)
     

! Modified by FDV to separate the self-energy calculation from the VFE one.
! Else If( (Lanc_In%RunTyp.Eq.1).OR.(Lanc_In%RunTyp.Eq.2) ) then!you'll be getting vibrational contribution to free energy
  Else If( Lanc_In%RunTyp.Eq.2 ) then!you'll be getting vibrational contribution to free energy
  Else If( (Lanc_In%RunTyp.Eq.1).OR.(Lanc_In%RunTyp.Eq.2) ) then!you'll be getting vibrational contribution to free energy
    qj0(:) = 0.0D0                                               !or phonon contribution to self energy 
    norm=0 
!note:dym_In%xyz(ip,r) ip=atom label in cluster,r=x,y,z position => {1,2,3}
    open(unit=100,file="dmdw_lanc.info ",status="unknown")

    !Make the requested displacement
    IF(Lanc_In%disp_opt.EQ.0)THEN !displace along each atom for full PHDOS 
      !This is not implemented properly
      write(100,*) "natom=",dym_In%natom
      write(100,*) "disp_opt=",Lanc_In%disp_opt
      wil_tmp(:) = 0.0
      w_pole_tmp(:) = 0.0
      Do iAt=1,dym_In%unAt!for nAt>1 want contribution from each type of atom in cell
        Do iDeg=1,dym_In%udeg(iAt)
          ii = dym_In%centeratomindex(iAt,iDeg)
          jj = (iAt-1)*Lanc_In%nPoles
          DO nn=1,3
            ! |0> = qj0 :Define Lanczos seed vector and normalize
            qj0(:) = 0.0
            wil(:) = 0.0
            w_pole(:) = 0.0
            qj0(3*(ii-1)+(nn-1)) = 1.0D0 !distort nn coordinate of central atom of this type 
            qj0 = qj0/dsqrt(sum(qj0**2))
            Call Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                         w_pole,wil,mnull,SPole_EinsteinFreq)
!           Do iq=0,Lanc_In%nPoles-1! nPoles/ atom in cell  
            Do iq=0,mnull! nPoles/ atom in cell  
              ! wil/w_pole indexed from 0.  Also, weighted according to atomic mass (am/n = atomic mass/number)
              wil_tmp(jj+iq)=wil_tmp(jj+iq)+wil(iq)*dym_In%am(ii)
              w_pole_tmp(jj+iq)=w_pole_tmp(jj+iq)+w_pole(iq)/3.0
            End Do
            write(100,*)
          END DO
!         DO iq=0,Lanc_In%nPoles-1
          DO iq=0,mnull
            norm = norm + wil_tmp(jj+iq)
          END DO
        End Do
      End Do
      !norm = dsqrt(sum(wil_tmp**2))
      write(100,*) "Pre-re-normalized wil"!comes out of Lanczos routine normalized, so needs to be re-normalized when stitched together
      write(100,*) "wil norm",norm,"SHOULD EQUAL natom,nice place for a warning probably" 
      wil_tmp=wil_tmp/norm
      wil(:) = 0.0
      w_pole(:) = 0.0
      DO iAt=1,dym_In%unAt
        jj = (iAt-1)*Lanc_In%nPoles
!       DO iq=0,Lanc_In%nPoles-1
        DO iq=0,mnull
          wil(   iq) = wil(   iq) + wil_tmp(   jj+iq)
          w_pole(iq) = w_pole(iq) + w_pole_tmp(jj+iq)
        END DO
      END DO 
      wil=wil/dym_In%unAt
      w_pole=w_pole/dym_In%unAt

    Else !displace along x/y/z only
        
      write(100,*) "natom=",dym_In%natom
      write(100,*) "disp_opt=",Lanc_In%disp_opt
      Do iAt=1,dym_In%unAt!for nAt>1 want contribution from each type of atom in cell
        Do iDeg=1,dym_In%udeg(iAt)
          ii = dym_In%centeratomindex(iAt,iDeg)
          jj = (iAt-1)*Lanc_In%nPoles
          qj0(:) = 0.0
          norm = 0.0
          ! |0> = qj0 :Define Lanczos seed vector and normalize
          qj0(3*(ii-1)+(Lanc_In%disp_opt-1)) = 1.0D0 !distort x/y/z coordinate of central atom of this type 
          qj0 = qj0/dsqrt(sum(qj0**2))
          Call Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                       w_pole,wil,mnull,SPole_EinsteinFreq)
!         Do iq=0,Lanc_In%nPoles-1! nPoles/ atom in cell  
          Do iq=0,mnull
            ! wil/w_pole indexed from 0.  Also, weighted according to atomic mass (am/n = atomic mass/number)
            wil_tmp(jj+iq)=wil(iq)*dym_In%am(ii)
            w_pole_tmp(jj+iq)=w_pole(iq) 
            norm=norm+wil_tmp(jj+iq)
            write(100,*)  dym_In%an(ii),ii,w_pole_tmp(jj+iq),wil_tmp(jj+iq)
          End Do
          write(100,*)
        End Do
      End Do

      write(100,*) "Pre-re-normalized wil"!comes out of Lanczos routine normalized, so needs to be re-normalized when stitched together
      write(100,*) "wil norm",norm,"SHOULD EQUAL natom,nice place for a warning probably" 
      wil_tmp=wil_tmp/norm
      wil(:) = 0.0
      w_pole(:) = 0.0
      DO iAt=1,dym_In%unAt
        jj = (iAt-1)*Lanc_In%nPoles
!       DO iq=0,Lanc_In%nPoles-1
        DO iq=0,mnull
          wil(   iq) = wil(   iq) + wil_tmp(   jj+iq)
          w_pole(iq) = w_pole(iq) + w_pole_tmp(jj+iq)
        END DO
      END DO 
      wil=wil/dym_In%unAt
      w_pole=w_pole/dym_In%unAt

      !!!write(100,*) "Normalized PHDOS FOR ALL ATOMS"
      !!!DO iAt=1,dym_In%unAt
      !!!  DO iDeg=1,dym_In%udeg(iAt)
      !!!    write(100,*) iAt,w_pole_tmp((iAt-1)*Lanc_In%nPoles+iq-1),wil_tmp((iAt-1)*Lanc_In%nPoles+iq-1)
      !!!  End DO
      !!!End DO
      
    End If ! Lanc_In%disp_opt
  End If ! Land_In%RunTyp 

  !Lanczos done.  Give the people what they asked for!
  IF ( (Lanc_In%RunTyp .Eq. 0) .or. (Lanc_In%RunTyp .Eq. 3) ) Then
    sig = 0.0_r08
! Scale the low frequency pole to address w=0 issue
!   if ( Lanc_In%RunTyp .eq. 3 ) then
!     wil(0) = wil(1)/w_pole(1)*w_pole(0)
!   end if
    do nn = 0, mnull
      cth = 1.0_r08 / dtanh(cotarg * w_pole(nn))
! Added by FDV:
!   I removed the zeroing out of imag. freq poles in the Lanczos (to get more
!   flexibility in the code). This means that now all sums based on the poles
!   need to check if the freq. is imag. (i.e. negative, cause we encode them
!   that way) and ignore the pole
      if ( w_pole(nn) > 0.0_r08 ) then
        sig = sig + wil(nn) / w_pole(nn) * cth
      end if
! Debug
!     print *, 'sig: ', sig
    end do
    s_dF = 0.0_r08   
    t_dF = 0.0_r08
!   mef=0.0e0 

  else IF ( Lanc_In%RunTyp .Eq. 1 ) Then
    s_dF = 0.0_r08   
!   do nn = 0, Lanc_In%nPoles-1
! Added by FDV:
! This is temporary code: Doing a first pass to get the scaling parameter
! to account for the loss of weight to the imag. freqs.
    W_Scal = 1.0_r08
!   W_Scal = 0.0_r08
!   do nn = 0, mnull
!     if ( w_pole(nn) > 0.0_r08 ) then
!       W_Scal = W_Scal + wil(nn)
!     end if
!   end do
!   W_Scal = 1.0_r08/W_Scal
    do nn = 0, mnull
      cth = 1.0_r08 / dtanh(cotarg * w_pole(nn))
! Added by FDV:
!   I removed the zeroing out of imag. freq poles in the Lanczos (to get more
!   flexibility in the code). This means that now all sums based on the poles
!   need to check if the freq. is imag. (i.e. negative, cause we encode them
!   that way) and ignore the pole
      if ( w_pole(nn) > 0.0_r08 ) then
! Modified by FDV:
! The change below ensured taht we don't have overflows on the dsinh operation
        where ( cotarg*w_pole(nn) <= 50.0_r08 )
!         s_dF = s_dF + wil(nn)*dlog(2.0_r08*dsinh(cotarg*w_pole(nn)))
          s_dF = s_dF + W_Scal*wil(nn)*dlog(2.0_r08*dsinh(cotarg*w_pole(nn)))
        else where
!         s_dF = s_dF + wil(nn)*cotarg*w_pole(nn)
          s_dF = s_dF + W_Scal*wil(nn)*cotarg*w_pole(nn)
        end where
!       s_dF = s_dF + wil(nn)*dlog(2.0_r08*dsinh(cotarg*w_pole(nn))) !srw       
!       print *, nn, w_pole(nn), cotarg(1), dsinh(cotarg(1)*w_pole(nn))
!       print *, nn, 2.0_r08*dsinh(cotarg(1)*w_pole(nn)), dlog(2.0_r08*dsinh(cotarg(1)*w_pole(nn)))
!       print *, nn, 2.0_r08*dsinh(cotarg(1)*w_pole(nn)), dlog(2.0_r08)+dlog(dsinh(cotarg(1)*w_pole(nn)))
      end if
    end do
!       stop
! Changed by FDV
! We remove the atom degeneracy for now. Will deal with this later
!   s_dF=s_dF*dym_In%udeg(itt)
! Changed by FDV
! We remove the factor of 3 because we now calculate all the directional
! components individually and add them up.
!   t_dF = 3.0D0*8.314713470D0*Lanc_In%T*s_dF
    t_dF = 8.314713470D0*Lanc_In%T*s_dF
!   mef=0.0e0 

  Else If(Lanc_In%RunTyp.Eq.2 )Then !This is where self-energy/spectral function calculation is performed
  !Get A2F
    !Separate out the electron-phonon couplings from A2F and PHDOS
    !Note: names of these files are currently hard coded in this routine
    call phonon_coupling(a2fall,a2,norm,Lanc_In)

    !A number of tests for a2f/weight PHDOS with coupling matrix elements
    open(unit=60,file="dmdw_a2f.info",status="unknown")
    Write(60,*) "# DMDW Option",Lanc_In%RunTyp
    Write(60,*) "# Displacement Option",Lanc_In%disp_opt
    Write(60,*) "# Lanczos Order",Lanc_In%nPoles
    Write(60,*) "# "
    Write(60,*) "# Lanczos Pole in Thz/weight PHDOS"
!   Do inull=0,Lanc_In%nPoles-1
    Do inull=0,mnull
      write(60,*) w_pole(inull)/6.28,wil(inull) !Thz
    End Do
    write(60,*) "# norm",norm
    write(60,*)

    ii=0
    !Get as close as possible to ii+1'th pole and match with appropriate coupling
    DO inull=1,400
      check=w_pole(ii)/6.28-a2(1,inull)*27.211396/4.135667516*1000!variable check is in Thz

      IF(check.LT.0.0)THEN
!        write(60,*)  inull
!        write(60,*)  "w_pole in Thz, wil "
!        write(60,*)  w_pole(ii)/6.28,wil(ii) !Thz
!        write(60,*)  "w in Thz, a2"
!        write(60,*)  a2(1,inull)*27.211396/4.135667516*1000,a2(2,inull)!1000for Thz.  
!        write(60,*) "w in Thz, a2f"
!        write(60,*)  a2(1,inull)*27.211396/4.135667516*1000,a2(2,inull)*wil(ii)

        !get a2f 
        a2f(ii,1) =  w_pole(ii)*6.58211928*0.0001!last factor is 10-16 from hbar and 10-12 from angular Thz  
! Added by FDV:
!   I removed the zeroing out of imag. freq poles in the Lanczos (to get more
!   flexibility in the code). This means that now all sums based on the poles
!   need to check if the freq. is imag. (i.e. negative, cause we encode them
!   that way) and ignore the pole
!   The code below should result in the same behavior as before.
        if ( w_pole(ii) > 0.0_r08 ) then
          a2f(ii,2) = a2(2,inull)*wil(ii)*norm!this norm factor accounts for the fact that Lanczos weights of PHDOS were normalized
        else
          a2f(ii,2) = 0.0_r08
        end if
!        a2f(ii,2) = wil(ii)*norm

!        write(60,*) "w in eV, a2f in array"
!        write(60,*) a2f(ii,1),a2f(ii,2)
!        write(60,*) " "
        
        ii=ii+1
      End IF
      
      IF(ii.eq.Lanc_In%nPoles)THEN 
!        write(60,*) ii, "IM GONE!"
        goto 13
      end if
    End Do ! inull=1,400

    13  write(60,*)
!        write(60,*) "what leaves:eV,a2f"

    write(60,*) "Pole/weight a2f in eV/Arb"
    w0 = 0.0
    tot = 0.0
    xx = 0.0
!   DO inull=0,Lanc_In%nPoles-1
    DO inull=0,mnull
      write(60,*) a2f(inull,1),a2f(inull,2)
      w0 = w0 + a2f(inull,1)*a2f(inull,2)
      tot = tot + a2f(inull,2)
      xx = xx + 2.d0*a2f(inull,2)/a2f(inull,1)
    END DO
    w0 = w0/tot
    write(60,*) "lambda =", xx
    write(60,*) "w0 =", w0


! Open all info files
open(unit=u_PSinfo,file="dmdw_spectral.info",status="unknown")        !A general diagnostic for spectral function
open(unit=u_Einfo, file="dmdw_Egrid.info",status="unknown")      !Energy grid    

! Open All output files for Self Energy 
open(unit=u_ReSE,file="dmdw_reSE_a2F.dat",status="unknown") !Real part of Self-Energy
open(unit=u_ImSE,file="dmdw_imSE_a2F.dat",status="unknown") !Imaginary part of Self-Energy

! Open all outputs for Spectral function
open(unit=u_Akw,file="dmdw_Akw.dat",status="unknown")

! Debugging files
IF (DEBUG) THEN
  open(unit=u_beta,file="dmdw_beta.dat",status="unknown")     !Beta functions Beta^p(w),Beta^h(w)
  open(unit=u_Ckt,file="dmdw_Ckt.dat",status="unknown")         !!! REMOVE EVENTUALLY?
  open(unit=u_Akt,file="dmdw_Akt.dat",status="unknown")         !!! REMOVE EVENTUALLY?
  open(unit=u_Aqp,file="dmdw_Aqp.dat",status="unknown")         !!! REMOVE EVENTUALLY?
  open(unit=u_Asat,file="dmdw_Asat.dat",status="unknown")       !!! REMOVE EVENTUALLY?
END IF

        ! ******************************** !
        ! Calculate SE(omega) and A(omega) !
        ! ******************************** !

    coni = (0.0d0,1.0d0)
    deli = 1.0d-30*coni
    
    ! time resolution
    Nt = 40000 + 1
    t0 = 2000.0
    dt = 2.0*t0/DBLE(Nt-1)

    !electron energy at which spectral function will be calculated. Recall epk=0 => fermi level 
    Ek = Lanc_In%E_k
    IF (Lanc_In%E_k_opt.EQ.1) THEN ! E_k was given in meV
      Ek = Ek * 1.0e-3 / w0 ! convert to eV and divide by w0 --> now Ek in units of w0
    END IF
    call SelfEn(Lanc_In,dym_In,cmplx(0.0e0,0.0e0),SE_a2f,a2f)  ! CAREFUL: SelfEn expects omega in eV!!
    betashift = abs(aimag(SE_a2f(1))/w0) ! in units of w0
    Gam = max(betashift,0.005)           ! in units of w0
    call SelfEn(Lanc_In,dym_In,cmplx(Ek*w0,0.0e0),SE_a2f,a2f)  ! CAREFUL: SelfEn expects omega in eV!!
    write(u_PSinfo,*) "Gamma_k = ",Gam
    write(u_PSinfo,*) "epk = E_k - ReSE(E_k) =", Ek-real(SE_a2f(1))/w0

    !Define energy range over which self-energy and spectral function are computed
    NE = 10000 + 1
    E0 = 10.0d0
    lowE = Ek - E0   ! in w0
    highE= Ek + E0   !  "  "  
    dE = (highE-lowE)/DBLE(NE-1) 
    !Generate energy grid (in units of w0)
    ilowE = nint(lowE/dE)
    ihighE = nint(highE/dE)
    DO iE=ilowE,ihighE
      w(iE) = DBLE(iE)*dE
    END DO

    ! get index of and E_k
    iEk = 0
    DO iE = ilowE,ihighE-1
      IF ((w(iE).LE.Ek).AND.(Ek.LE.w(iE+1))) THEN
        IF (Ek-w(iE) < w(iE+1)-Ek) THEN
          iEk = iE
        ELSE
          iEk = iE+1
        END IF
      END IF
    END DO

    ! print out energy grid info
    write(u_Einfo,FMT_A)    "#  Energies printed in meV"
    write(u_Einfo,FMT_AFAF) "#  lowE ",lowE*w0*1.0e3," highE ",highE*w0*1.0e3
    write(u_Einfo,FMT_AFAF) "#  dE  = ",dE*w0*1.0e3," w0 = ",w0*1.0e3
    write(u_Einfo,FMT_AFAF) "#  Ek  = ",Ek*w0*1.0e3," --> E = ",w(iEk)*w0*1.0e3

    !Initialize some things 
    imSE = 0.d0
    reSE = 0.d0
    atot = 0.d0
    write(u_ReSE,FMT_A) "#  Real part of the Self-energy"
    write(u_IMSE,FMT_A) "#  Imaginary part of the Self-energy"
    Do iE=ilowE,ihighE
      !Calculate Self-Energy
      omega = cmplx(w(iE)*w0,0.0d0)               
      call SelfEn(Lanc_In,dym_In,omega,SE_a2f,a2f)
      reSE(iE) = real(SE_a2f(1))/w0
      imSE(iE) = -sign(1.0,real(omega))*abs(aimag(SE_a2f(1)))/w0
      ! write out SE vs w (eV/eV)
      write(u_ReSE,*) real(omega),reSE(iE)*w0
      write(u_ImSE,*) real(omega),imSE(iE)*w0
      
      ! calculate betap and betah
      omega = cmplx((Ek+w(iE))*w0,0.0d0)
      call SelfEn(Lanc_In,dym_In,omega,SE_a2f,a2f)
      betak(iE) = (abs(aimag(SE_a2f(1))/w0) - betashift)/pi
      IF (DEBUG) THEN
        write(u_beta,*) w(iE),betak(iE)
      END IF
      
      ! For QP peak weight
      IF (iE > ilowE) THEN
        width = w(iE) - w(iE-1)
        cheightL = betak(iE-1)/(w(iE-1)-deli)**2
        cheightR = betak( iE )/(w( iE )-deli)**2
        IF (iE.EQ.0) THEN
          width = w(iE-1) - w(iE-2)
          cheightL = betak(iE-2)/(w(iE-2)-deli)**2
          cheightR = betak(iE-1)/(w(iE-1)-deli)**2
        ELSE IF (iE.EQ.1) THEN
          width = w(iE+1) - w(iE)
          cheightL = betak( iE )/(w( iE )-deli)**2
          cheightR = betak(iE+1)/(w(iE+1)-deli)**2
        END IF
        atot = atot + 0.5*(cheightL+cheightR)*width
      END IF
 
    End DO
    Zk = exp(-atot)
    write(u_PSinfo,*) ""
    atot = -1.d0*(reSE(iEk+1)+coni*imSE(iEk+1)-reSE(iEk-1)-coni*imSE(iEk-1))/(w(iEk+1)-w(iEk-1))
    Zk = exp(-atot)
    write(u_PSinfo,*) "atot    = ", atot
    write(u_PSinfo,*) "Zk      = ", Zk
    write(u_PSinfo,*) ""
   
    ! Calculate cumulant and analytically continue for t<0
    cIntegrand = 0.d0
    DO it = 1,(Nt-1)/2
      t = REAL(it)*dt
      DO jj = 1,NE
        iE = ilowE+jj-1
        xx = cos(w(iE)*t)
        yy = sin(w(iE)*t)
        ! integrand for a_k(t) integral
        cIntegrand(jj) = betak(iE)*(xx-coni*yy-1.d0)/(w(iE)-deli*sign( 1.0_r08,t))**2
      END DO
      ! trap rule as vector operation
      Ckt( it) = dE*(sum(cIntegrand)-0.5*(cIntegrand(1)+cIntegrand(NE)))
      Ckt(-it) = CONJG(Ckt(it))
    END DO
    Ckt(0) = 0.d0
    Akt = exp(Ckt)
    AktS = Akt - Zk

    IF (DEBUG) THEN
      ! print Ck(t) and Ak(t)
      write(u_Ckt,*) "# t, mag, ph, re, im"
      write(u_Akt,*) "# t, mag, ph, re, im"
      DO it = -(Nt-1)/2,(Nt-1)/2
        t = DBLE(it)*dt                              
        xx =  REAL(Ckt(it))
        yy = AIMAG(Ckt(it))
        write(u_Ckt,'(20f20.10)') t,sqrt(xx*xx+yy*yy),atan2(yy,xx),xx,yy
        xx =  REAL(Akt(it))
        yy = AIMAG(Akt(it))
        write(u_Akt,'(20f20.10)') t,sqrt(xx*xx+yy*yy),atan2(yy,xx),xx,yy
      END DO
    END IF

    ! FT to get Ak(w) and calculate normalization factor
    integral = 0.0
    AkwQP = 0.0d0
    AkwS = 0.0d0
    DO iE = ilowE,ihighE
      DO it = -(Nt-1)/2,(Nt-1)/2-1
        t = DBLE(it)*dt
        width = dt
        ! common integrand factor
        cheightL = EXP(coni*(w(iE)-Ek)*( t  )-Gam*abs( t  ))/(2.d0*pi)
        cheightR = EXP(coni*(w(iE)-Ek)*(t+dt)-Gam*abs(t+dt))/(2.d0*pi)
        Akw(iE) = Akw(iE) + 0.5d0*(Akt(it)*cheightL + Akt(it+1)*cheightR)*width
        ! QP peak and just sats
        AkwQP(iE) = AkwQP(iE) + 0.5d0*Zk*(cheightL+cheightR)*width
        AkwS(iE) = AkwS(iE) + 0.5d0*(AktS(it)*cheightL + AktS(it+1)*cheightR)*width
      END DO
    END DO
    integral = dE*(sum(Akw) - 0.5*(Akw(ilowE)+Akw(ihighE)))
    norm = sqrt(integral**2)
    Akw = Akw/norm
    AkwQP = AkwQP/norm

    ! print Ak>(w) and Ak<(w) and sum
    write(u_Akw,*) "# norm = ", norm
    write(u_Akw,*) "# w [meV], mag, ph, re, im"
    DO iE = ilowE,ihighE
      xx =  REAL(Akw(iE))/(w0*1.0e3)
      yy = AIMAG(Akw(iE))/(w0*1.0e3)
      write(u_Akw,'(20f20.10)') w(iE)*w0*1.0e3,sqrt(xx*xx+yy*yy),atan2(yy,xx),xx,yy
      IF (DEBUG) THEN
        xx =  REAL(AkwQP(iE))/(w0*1.0e3)
        yy = AIMAG(AkwQP(iE))/(w0*1.0e3)
        write(u_Aqp,'(20f20.10)')  w(iE)*w0*1.0e3,sqrt(xx*xx+yy*yy),atan2(yy,xx),xx,yy
        xx =  REAL(AkwS(iE))/(w0*1.0e3)
        yy = AIMAG(AkwS(iE))/(w0*1.0e3)
        write(u_Asat,'(20f20.10)')  w(iE)*w0*1.0e3,sqrt(xx*xx+yy*yy),atan2(yy,xx),xx,yy
      END IF
    END DO

    ! Close all outputs for Self Energy and Spectral function
    Close(unit=u_PSinfo)
    Close(unit=u_Einfo)
    Close(unit=u_ReSE)
    Close(unit=u_ImSE)
    Close(unit=u_Akw)

    IF (DEBUG) THEN
      Close(unit=u_beta)
      Close(unit=u_Ckt)
      Close(unit=u_Akt)
      Close(unit=u_Aqp)
      Close(unit=u_Asat)
    END IF
  End If ! Lanc_In%RunTyp (==2)

! ********************************************* !
! End Self-Energy and Spectral function section !
! ********************************************* !
     
!! Debug
!   if ( q0_opt ) then
!     sig2 = hbarc * caps / (2.0_r08 * dym_In%am(1) * amu) * sig
!     if ( temp_opt ) then
!       t_dF = 3.0D0*8.314713470D0*Lanc_In%T*s_dF
!       write(IO_Err,FMT='(A,F20.10)') 'Free Energy:', t_dF
!     end if
!   else
! Debug
!               print *, ' hbarc', hbarc
!               print *, ' caps', caps
!               print *, ' amu', amu
!               print *, ' sig', sig
      sig2 = hbarc * caps / (2.0D0 * mu * amu) * sig
!   end if
! Debug
!   write(IO_Err,FMT='(A,2F12.8)') ' mu, am(1):', mu, am(1)

! Debug
! Added by Fer to create a table of vibrational free energies
! as a function of temperature

! Commenting out for now 

!    if ( temp_opt ) then
!      write(IO_Err,FMT='(A)') "# Vibrational free energy (J/mol-c)"
!      write(IO_Err,FMT='(A)') "#  T          F(J/mol-c)"
!!     do Tempi=1.0D0,temp,1.0D0
!      do Tempi=1.0D0,700.0D0,1.0D0
!        s_dF = 0.0D0
!       tk = 8.617385d-5 * Tempi
!        beta = 1.0_r08/(K_B*Tempi)
!        cotarg = 0.5d0 * hbar * beta
!        do nn = 0, mnull
!          ome = dsqrt(xnull(nn))
!          s_dF = s_dF + wil(nn)*dlog(2.0D0*dsinh(cotarg*ome))
!        end do
!        t_dF = 3.0D0*8.314713470D0*Tempi*s_dF
!        write(IO_Err,FMT='(F6.2,F20.10)') Tempi, t_dF
!      end do
!    end if

! NOTE: This code is somewhat temporary. Trying to clean up the interface to
! Calc_DW. This code is an intermediate step keeping the variables used
! internally the same, and copying the results to the actual output variables.
! We make sure that things are allocated, or (re)allocate approprietly.

! Force deallocate all the output arrays, skip errors using "stat"
  deallocate(DW_Out%s2, DW_Out%u2, DW_Out%vfe, &
             DW_Out%Poles_Frq, DW_Out%Poles_Wgt, stat=iStat)

! We include information about the path for which these results were computed.
! This makes it easier later on when we want to print stuff.
  DW_Out%Path_nAt  = nleg
  allocate(DW_Out%Path(DW_Out%Path_nAt))
  DW_Out%Path      = lpath_in
  DW_Out%Path_Pert = iu2Pert
! Calculate the path length
  DW_Out%Path_Len  = 0.0_r08
  do iAt=1,DW_Out%Path_nAt-1
    DW_Out%Path_Len  = DW_Out%Path_Len + &
      sqrt(sum((dym_In%xyz(lpath_in(iAt),:)-dym_In%xyz(lpath_in(iAt+1),:))**2))
  end do

! The poles freqs and weights are always returned
  allocate(DW_Out%Poles_Frq(mnull+1), &
           DW_Out%Poles_Wgt(mnull+1) )
  DW_Out%Poles_Frq = w_pole(0:mnull)/(2.0D0*pi)
  DW_Out%Poles_Wgt = wil(0:mnull)
  DW_Out%SPole_Frq = SPole_EinsteinFreq
  DW_Out%RedMass   = mu

! Depending on the type of calculation we allocate only certain things
  select case ( Lanc_In%RunTyp )

! EXAFS Debye-Waller factors (s2)
    case ( 0 )

      allocate(DW_Out%s2(Lanc_In%nT))
      DW_Out%s2 = sig2

! Crystallographic Debye-Waller factors (u2)
    case ( 3 )

      allocate(DW_Out%u2(Lanc_In%nT))
      DW_Out%u2 = sig2

! Vibrational Free Energy
    case ( 1 )

      allocate(DW_Out%vfe(Lanc_In%nT))
      DW_Out%vfe = t_dF

! Self-Energy and related quantities
    case ( 2 )

! NOTE: No output here. This should be changed since the SE part does
! return stuff

! IR spectra based on Gaussian's dipole derivatives
    case ( 4 )

! NOTE: No output for this yet.

! Compute partial and total phonon densities of states
    case ( 5 )

! NOTE: No output for this yet.

! If we get here it means we don't have this property yet
    case default

      write(IO_Err,fmt='(a)') &
        'Internal error in CalcDW: Unrecognized calculation type'
      stop

  end select

! Deallocate the arrays
    deallocate(sig,beta,cotarg,cth)
    deallocate(lpath_m, lpath_p, rc, w_pole, wil, qj0)

  end subroutine Calc_DW

!##############################################################################


  function Poly_Y(Poly_Type,ni,xp,anj,bnj) result(Y)

  implicit none

  character                     :: Poly_Type
  integer(kind=i08)             :: ni
  real(kind=r08)                :: xp
  real(kind=r08), dimension(0:) :: anj, bnj
  real(kind=r08)                :: Y

  integer(kind=i08) :: n
  real(kind=r08)    :: Y_m1, Y_m2

    select case (Poly_Type )
      case ( 'S' )
        Y_m1 = 1.0_r08
        Y    = xp - anj(0)
      case ( 'R' )
        Y_m1 = 0.0_r08
        Y    = 1.0_r08
      case ( 'P' )
        Y_m1 = 0.0_r08
        Y    = 1.0_r08
      case default
        write(IO_Err,fmt='(a)') ' Error in Poly_Y: '
        write(IO_Err,fmt='(a,a)') ' Poly_Type = ', Poly_Type
        stop
    end select

    do n=2,ni
      Y_m2 = Y_m1
      Y_m1 = Y
      Y = (xp - anj(n-1))*Y_m1 - bnj(n-1)**2*Y_m2
    end do

  end function Poly_Y

  function PolyD_Y(Poly_Type,ni,xp,anj,bnj) result(Yp)

  implicit none

  character                     :: Poly_Type
  integer(kind=i08)             :: ni
  real(kind=r08)                :: xp
  real(kind=r08), dimension(0:) :: anj, bnj
  real(kind=r08)                :: Yp

  integer(kind=i08) :: n
  real(kind=r08)    :: Y_m1, Y_m2, Y_m3, Yp_m1, Yp_m2

    select case (Poly_Type )
      case ( 'S' )
        Y_m2 = 0.0_r08
        Y_m1 = 1.0_r08 
        Yp_m1 = 0.0_r08
        Yp    = 1.0_r08
      case default
        write(IO_Err,fmt='(a)') ' Error in Poly_Y: '
        write(IO_Err,fmt='(a,a)') ' Poly_Type = ', Poly_Type
        stop
    end select

    do n=2,ni
      Y_m3 = Y_m2
      Y_m2 = Y_m1
      Y_m1 = (xp - anj(n-2))*Y_m2 - bnj(n-2)**2*Y_m3
      Yp_m2 = Yp_m1
      Yp_m1 = Yp
      Yp = Y_m1 + (xp - anj(n-1))*Yp_m1 - bnj(n-1)**2*Yp_m2
    end do

  end function PolyD_Y

!#######################################################################
  logical function Are_Close(d,x,y)

  implicit none

  real(kind=r08) :: d, x, y

  if ( dabs(x-y) > d ) then
    Are_Close = .false.
  else
    Are_Close = .true.
  end if 

  end function Are_Close
!#######################################################################

  subroutine Make_DM(dym_In)

    implicit none

    type(dym_Info), intent(inout) :: dym_In
    
    integer(kind=i08) ::    iAt, jAt
    integer(kind=i08) ::    ip, jq
 
    real(kind=r08)    :: Avg_Val, Max_Val, Avg_Asym
    real(kind=r08)    :: Asym_T1, Asym_T2

! Allocate the full dynamical matrix
    if ( .not. allocated(dym_In%dm) ) then
      allocate(dym_In%dm(0:3*dym_In%nAt-1,0:3*dym_In%nAt-1))
    end if

! Convert the block dynamical matrix into a full one
    do iAt=1,dym_In%nAt
      do jAt=1,dym_In%nAt

        do ip=1,3
          do jq=1,3
            dym_In%dm((iAt-1)+dym_In%nAt*(ip-1),(jAt-1)+dym_In%nAt*(jq-1)) = &
              auf2npm*npm2amups2* &
              dym_In%dm_block(iAt,jAt,ip,jq) / &
              sqrt(dym_In%am(iAt)*dym_In%am(jAt))
          end do
        end do

      end do
    end do

! Debug
! stop

! Check that the dynamical matrix is symmetric within a certain threshold
    Avg_Val = sum(abs(dym_In%dm))/(3*dym_In%nAt)**2
    Max_Val = maxval(abs(dym_In%dm))
    Avg_Asym = sum(abs(dym_In%dm - transpose(dym_In%dm)))/(3*dym_In%nAt)**2
    Asym_T1 = (Avg_Asym/Avg_Val)*100
    Asym_T2 = (Avg_Asym/Max_Val)*100
    if ( Asym_T1 > 50 .or. Asym_T2 > 5 ) then
      write(IO_Err,fmt='(a)') 'Warning: The dynamical matrix is not symmetric'
      write(IO_Err,fmt='(a)') &
        '         All results should be carefully examined'
    end if

  end subroutine Make_DM

  subroutine Read_dym_Info(dym_file,dym_In)

! Read the dym data file

    implicit none

    character(len=*)               :: dym_file
    type(dym_Info), intent(out)    :: dym_In

    integer(kind=i08)              :: iAt, jAt, ip, jq, iPair
    integer(kind=i08)              :: ideg,itype,xx

    open(unit=IO_dym,file=dym_file,status='old')

! Read the dym type flag
    read(IO_dym,*) dym_In%Type

! Debug
!   write(6,*) dym_In%Type

! Read the number of atoms
    read(IO_dym,*) dym_In%nAt

! Allocate the arrays right before they are filled with their data.
! NOTE: Before we used to allocate everything here, but now I have decided
!       that unless an array has its corresponding data, it won't be
!       allocated. This way we can check its allocation status to see if it
!       makes sense to try and use the data.

! Read the atomic number for each atom
    allocate(dym_In%an(dym_In%nAt))
    read(IO_dym,*) (dym_In%an(iAt),iAt=1,dym_In%nAt)

! Read the atomic mass for each atom
    allocate(dym_In%am(dym_In%nAt))
    read(IO_dym,*) (dym_In%am(iAt),iAt=1,dym_In%nAt)

! Fix the atomic numbers and masses in case either is 0
    call Fix_AN_AM(dym_In)

! Depending on the type of dym file, we either read cartesian coordinates
! or reduced coordinates

    if      ( dym_In%Type == 1 .or. &
              dym_In%Type == 2 .or. &
              dym_In%Type == 3      ) then

! Read the atomic coordinates
      allocate(dym_In%xyz(dym_In%nAt,3))
      do iAt=1,dym_In%nAt
        read(IO_dym,*) dym_In%xyz(iAt,:)
      end do

    else if ( dym_In%Type == 4 ) then

! Read the reduced coordinates
      allocate(dym_In%red(dym_In%nAt,3))
      do iAt=1,dym_In%nAt
        read(IO_dym,*) dym_In%red(iAt,:)
      end do

    end if

! NOTE: Removed this, they are already in au
! Convert to atomic units
!   dym_In%xyz = dym_In%xyz/ang2au

! Read the block dynamical matrix
    allocate(dym_In%dm_block(dym_In%nAt,dym_In%nAt,3,3))
    do iPair=1,dym_In%nAt**2
      read(IO_dym,*) iAt, jAt
      do ip=1,3
        read(IO_dym,*) (dym_In%dm_block(iAt,jAt,ip,jq),jq=1,3)        
      end do
    end do

! Debug (FDV)
! Adding a quick on the fly conversion from eV/Ang^2 to au of force constant
! This is done only to check the Gilles Hug results
!   dym_In%dm_block = dym_In%dm_block*0.0102908487617_r08

! Add check to read extra info only when the dym file type is correct
    if ( dym_In%Type .eq. 2 ) then

      open(unit=97,file="dmdw_dymtest.info",status="unknown")

!srw:read additions to .dym
      Read(IO_dym,*) dym_In%unAt,dym_In%natom
!Now that we know how many unique types there are (unAt) allocate things...
      allocate(dym_In%udeg(dym_In%unAt),dym_In%utype(dym_In%unAt), &
               dym_In%u_xyz(dym_In%unAt,dym_In%nAt,3), &
               dym_In%centeratomindex(dym_In%unAt,dym_In%natom))
!... and read the rest
      dym_In%u_xyz(:,:,:)=-1.0d0
      Do itype=1,dym_In%unAt
        Read(IO_dym,*) dym_In%utype(itype),dym_In%udeg(itype)
        Do ideg=1,dym_In%udeg(itype)
          Read(IO_dym,*) dym_In%centeratomindex(itype,ideg), &
                         xx,dym_In%u_xyz(itype,ideg,:)  !is this right?     
          write(97,*) dym_In%centeratomindex(itype,ideg), xx,dym_In%u_xyz(itype,ideg,:)
        End Do
      End Do

      write(97,*) ""
      write(97,*) dym_In%unAt,dym_In%natom 
      Do itype=1,dym_In%unAt
        DO ideg=1,dym_In%udeg(itype)
          write(97,*) itype,ideg,dym_In%centeratomindex(itype,ideg)
        End Do
      End DO

      close(97)

    end if

! Added by FDV to read the dipole derivatives for IR spectrum
    if ( dym_In%Type .eq. 3 ) then

! Empty separator line
      read(IO_dym,*)

! Read the dipole derivatives and store in the right order
! This order should allow us to use the dip. derivs directly as Lanczos seed.

      allocate(dym_In%Dip_Der(3,dym_In%nAt,3))
      do iAt=1,dym_In%nAt
        do ip=1,3
          read(IO_dym,*) (dym_In%Dip_Der(jq,iAt,ip),jq=1,3)        
        end do
      end do

    end if

! Added by FDV to read the cell vectors, if we have the right dym type
    if ( dym_In%Type == 4 ) then

! Empty separator line
      read(IO_dym,*)

! Read the three cell vectors, in atomic units
      allocate(dym_In%cell(3,3))
      do ip=1,3
        read(IO_dym,*) (dym_In%cell(ip,jq),jq=1,3)        
      end do

    end if

! Now that we have read the information, we do some post-processing, if needed

! For dym type 4, we read reduced coordinates and cell vectors, so we generate
! the cartesian coordinates internally
    if ( dym_In%Type == 4 ) then

      allocate(dym_In%xyz(dym_In%nAt,3))
      do iAt=1,dym_In%nAt
        dym_In%xyz(iAt,:) = matmul(dym_In%red(iAt,:),dym_In%cell)
      end do

    end if

! Now we allocate and fill the array with the indices of the basis atoms.
! These are the atoms that must be used to obtain a complete picture of the
! system. For instance, they would be the unique atoms in a unit cell,
! embedded in a cluster.
! NOTE: In the future we will add a weight or repetition array to do the
!       averaging even simpler.
! NOTE: Since we dont have anything intelligent to do with this yet, we simply
!       allocate and fill the array with a list of all indices for the atoms
!       in the system.
    allocate(dym_In%BasAtInd(dym_In%nAt))
    dym_In%BasAtInd = (/ (iAt,iAt=1,dym_In%nAt) /)

  end subroutine Read_dym_Info

  subroutine Read_Lanczos_Info(Lanc_In)

    implicit none

    type(Lanczos_Info), intent(out) :: Lanc_In
    real(kind=r08)    :: T_Min, T_Max, Temp
    integer(kind=i04) :: iT
    integer, parameter :: Buffer_Length = 256
    character(len=*), parameter :: Buffer_Format = '(A256)'
    character(len=Buffer_Length) :: Buffer
    integer :: iRStat, PDOS_Fmt = 0

! Read the level of output we want
! More general than before, to get information on printout style of pDOS
    read(IO_In,*) Lanc_In%IOFlag
    IF(Lanc_In%IOFlag.EQ.-999) STOP

! Read from standard input
    read(IO_In,*) Lanc_In%nPoles

! NOTE: This code is probably temporary
! Used to be able to always have a dmdw.inp file even when the file is empty
    if ( Lanc_In%nPoles == -999 ) then
      stop
    end if

! Read the temperature grid parameters
! Modifying this to improve input file format
    read(IO_In,fmt=Buffer_Format) Buffer
    read(Buffer,*) Lanc_In%nT
    if ( Lanc_In%nT .eq. 1 ) then
      read(Buffer,*) Lanc_In%nT, T_Min
      T_Max = T_Min
    else
      read(Buffer,*,end=10) Lanc_In%nT, T_Min, T_Max
      if ( T_Max < T_Min ) then
        Temp  = T_Max
        T_Max = T_Min
        T_Min = Temp
      end if
    end if

! Create the temperature grid
    allocate(Lanc_In%T(Lanc_In%nT))
    if (Lanc_In%nT == 1 ) then
      Lanc_In%T = T_Min
    else
      Lanc_In%T = (/ &
                    (T_Min + &
                    (T_Max-T_Min)/(Lanc_In%nT-1)*iT,iT=0,Lanc_In%nT-1) &
                  /)
    end if

! Read the type of perturbation
    read(IO_In,fmt=Buffer_Format) Buffer
    read(Buffer,*) Lanc_In%RunTyp
    select case ( Lanc_In%RunTyp )
      case ( 0 )
! Do nothing (for now)
      case ( 1 )
! Do nothing (for now)
      case ( 2 )
! Do nothing (for now)
      case ( 3 )
! Do nothing (for now)
      case ( 4 )
! Do nothing (for now)
      case ( 5 )
! Read options for pDOS printing
        read(Buffer,*,iostat=iRStat) Lanc_In%RunTyp, PDOS_Fmt
        if ( iRStat < 0 ) then
          PDOS_Fmt = 0
        end if
        select case ( PDOS_Fmt )
          case ( 0 )
! This is the default and has fixed settings
            Lanc_In%PDOS_Poles = .true.
            read(Buffer,*,iostat=iRStat) &
              Lanc_In%RunTyp, PDOS_Fmt, Lanc_In%PDOS_Part
            if ( iRStat < 0 ) then
              Lanc_In%PDOS_Part  = .false.
            end if
          case ( 1 )
            Lanc_In%PDOS_Rect  = .true.
            read(Buffer,*,iostat=iRStat) &
              Lanc_In%RunTyp, PDOS_Fmt, Lanc_In%PDOS_Part, Lanc_In%PDOS_DropL
            if ( iRStat < 0 ) then
              Lanc_In%PDOS_DropL = .false.
              Lanc_In%PDOS_Part  = .false.
            end if
          case ( 2 )
            Lanc_In%PDOS_Gauss = .true.
            read(Buffer,*,iostat=iRStat) &
              Lanc_In%RunTyp, PDOS_Fmt, Lanc_In%PDOS_Part, &
              Lanc_In%PDOS_Broad, Lanc_In%PDOS_Res
            if ( iRStat < 0 ) then
              Lanc_In%PDOS_Broad = 0.500_r08
              Lanc_In%PDOS_Res   = 0.001_r08
            end if
          case ( 10 )
            Lanc_In%PDOS_Poles = .true.
            Lanc_In%PDOS_Rect  = .true.
            Lanc_In%PDOS_Gauss = .true.
            read(Buffer,*,iostat=iRStat) &
              Lanc_In%RunTyp, PDOS_Fmt, Lanc_In%PDOS_Part, &
              Lanc_In%PDOS_DropL, Lanc_In%PDOS_Broad, Lanc_In%PDOS_Res
            if ( iRStat < 0 ) then
              Lanc_In%PDOS_DropL = .false.
              Lanc_In%PDOS_Part  = .false.
              Lanc_In%PDOS_Broad = 0.500_r08
              Lanc_In%PDOS_Res   = 0.001_r08
            end if
          case default
            write(IO_Err,fmt='(a)') 'Error: Unrecognized PDOS format'
            stop
        end select
      case default
        write(IO_Err,fmt='(a)') 'Error: Unrecognized run type'
        stop
    end select

! Debug
!   print *, 'PDOS_Poles', Lanc_In%PDOS_Poles
!   print *, 'PDOS_Rect ', Lanc_In%PDOS_Rect
!   print *, 'PDOS_Gauss', Lanc_In%PDOS_Gauss
!   print *, 'PDOS_DropL', Lanc_In%PDOS_DropL
!   print *, 'PDOS_Part ', Lanc_In%PDOS_Part
!   print *, 'PDOS_Broad', Lanc_In%PDOS_Broad
!   print *, 'PDOS_Res  ', Lanc_In%PDOS_Res  

!   read(IO_In,*) Lanc_In%RunTyp

! Modified by Fer to swith computation of the VFE from the old way to the
! new one.
!   if ( ( Lanc_In%RunTyp .eq. 1 ) .or. &
!        ( Lanc_In%RunTyp .eq. 2 )        ) then
    if ( Lanc_In%RunTyp .eq. 2 ) then

      read(IO_In,*) Lanc_In%disp_opt !srw

! Read E_k
      read(IO_In,*) Lanc_In%E_k_opt,Lanc_In%E_k!sms

    end if

! Read the dym filename
    read(IO_In,*) Lanc_In%dym_file

! Modified by Fer to swith computation of the VFE from the old way to the
! new one.
! Read the ABINIT PDS,A2F filenames
!   if ( ( Lanc_In%RunTyp .eq. 1 ) .or. &
!        ( Lanc_In%RunTyp .eq. 2 )        ) then
    if ( Lanc_In%RunTyp .eq. 2 ) then
      read(IO_In,*) Lanc_In%pds_file
      read(IO_In,*) Lanc_In%a2f_file
    end if

! If everything work
    return

! Error handling code

10  write(IO_Err,fmt='(a)') 'Error while reading temperature grid'
    write(IO_Err,fmt='(a)') 'Did you include both the min and max temperatures?'
    stop

20  write(IO_Err,fmt='(a)') 'Error while reading IO information'
    write(IO_Err,fmt='(a)') 'Please check the manual for detailed options'
    stop

  end subroutine Read_Lanczos_Info

  subroutine Read_Paths_Info(Paths_In)

    implicit none

    type(Paths_Info), intent(out) :: Paths_In

    integer, parameter :: Buffer_Length = 80
    character(len=*), parameter :: Buffer_Format = '(A80)'
    character(len=Buffer_Length), dimension(:), allocatable :: Buffer
    integer(kind=i04) :: iDesc, iAtom
    integer(kind=i04) :: jleg, mxjleg

! Read the number of path descriptors
    read(IO_In,*) Paths_In%nDesc

! Debug
!   write(IO_Out,*) Paths_In%nDesc

! If we don't have any descriptors, we return without allocating anything
! and we let each type of calculation decide what to do for paths
    if ( Paths_In%nDesc < 1 ) then
      return
    end if

! Allocate and read to the buffer
    allocate(Buffer(Paths_In%nDesc))
    do iDesc=1,Paths_In%nDesc
      read(IO_In,fmt=Buffer_Format) Buffer(iDesc)
! Debug
!     write(IO_Out,*) iDesc, Buffer(iDesc)
    end do

! Find the highest degree of multiple-scattering
    mxjleg = -1
    do iDesc=1,Paths_In%nDesc
      read(Buffer(iDesc),*) jleg
      if ( jleg > mxjleg ) then
        mxjleg = jleg
      end if
    end do

! Debug
!   write(IO_Out,*) mxjleg

! Allocate the paths info arrays
    allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
             Paths_In%Desc(Paths_In%nDesc,0:mxjleg-1), &
             Paths_In%Desc_mxR(Paths_In%nDesc))

! Read the complete path information
    do iDesc=1,Paths_In%nDesc
      read(Buffer(iDesc),*) Paths_In%Desc_Len(iDesc), &
                            (Paths_In%Desc(iDesc,iAtom), &
                               iAtom=0,Paths_In%Desc_Len(iDesc)-1), &
                            Paths_In%Desc_mxR(iDesc)
    end do

! Debug
!   write(IO_Out,*) Paths_In%nDesc
!   write(IO_Out,*) Paths_In%Desc_Len
!   write(IO_Out,*) Paths_In%Desc
!   write(IO_Out,*) Paths_In%Desc_mxR

! Convert the maximum path length to atomic units
    Paths_In%Desc_mxR = Paths_In%Desc_mxR*ang2au

  end subroutine Read_Paths_Info

  subroutine Print_Header(Lanc_In,Paths_In,dym_In)

    implicit none
  
    type(Lanczos_Info), intent(in) :: Lanc_In
    type(Paths_Info),   intent(in) :: Paths_In
    type(dym_Info),     intent(in) :: dym_In

! Print an output header

    write(IO_Out,fmt='(a,i4)')   '# Lanczos recursion order: ', Lanc_In%nPoles
    if ( Lanc_In%nT == 1 ) then
      write(IO_Out,fmt='(a,f7.2)') '# Temperature: ', Lanc_In%T
    else
      write(IO_Out,fmt='(a,f7.2)') '# Temperature: (See list Below)'
    end if
    write(IO_Out,fmt='(a,a)') &
      '# Dynamical matrix file: ', trim(Lanc_In%dym_file)

If( Lanc_In%RunTyp.Eq.0 ) then
    write(IO_Out,fmt='(a)') &
      ''
Else IF( Lanc_In%RunTyp.Eq.1 ) then
! Changed by FDV
! Eliminating this for now to modify the output for the VFE
!   write(IO_Out,fmt='(a)') &
!     '#    Z        Degeneracy    Vibrational Free Energy'
Else IF ( Lanc_In%RunTyp.Eq.2 ) then
    write(IO_Out,fmt='(a)') &
      'Mass Enchancement Factor  '
End If 

  end subroutine Print_Header

  subroutine Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)

    implicit none

    type(dym_Info), intent(in)     :: dym_In
    type(Paths_Info), intent(in)   :: Paths_In
    integer(kind=i08), intent(in)  :: iDesc
    type(Paths), intent(out)       :: Desc_Paths

    integer(kind=i08), dimension(:), allocatable :: Path_Begin, Path_End, &
                                                    Path_N
    integer(kind=i08)              :: nPaths_Unpruned
    integer(kind=i08), dimension(:,:), allocatable :: Paths_Unpruned
    real(kind=r08), dimension(:), allocatable :: Paths_Len_Unpruned
    integer(kind=i08)              :: iPath, jPath, iAtom, jAtom
    logical, dimension(:), allocatable :: Paths_Keep
    logical :: Rep_At
    real(kind=r08)  :: Length

! Until we find a more intelligent way of doing this, we use two different
! pieces of code to treat the single (for u2) and multi-atom (for s2) paths
    if ( Paths_In%Desc_Len(iDesc) == 1 ) then

! If the descriptor is 0, then we create a list of all atoms in the cluster
      if ( Paths_In%Desc(iDesc,0) == 0 ) then
        
        Desc_Paths%N = dym_In%nAt
        allocate(Desc_Paths%Ind(Desc_Paths%N,0:0), &
                 Desc_Paths%Len(Desc_Paths%N))
        Desc_Paths%Ind(:,0) = (/ (iAtom,iAtom=1,dym_In%nAt) /)
        Desc_Paths%Len(Desc_Paths%N)   = 0.0_r08
      else
! If the descriptor is not 0, then we just use that atom number
        Desc_Paths%N = 1
        allocate(Desc_Paths%Ind(Desc_Paths%N,0:0), &
                 Desc_Paths%Len(Desc_Paths%N))
        Desc_Paths%Ind(Desc_Paths%N,0) = Paths_In%Desc(iDesc,0)
        Desc_Paths%Len(Desc_Paths%N)   = 0.0_r08
      end if

    else

! Allocate auxiliary variables
      allocate(Path_Begin(0:Paths_In%Desc_Len(iDesc)-1), &
               Path_End(0:Paths_In%Desc_Len(iDesc)-1),   &
               Path_N(0:Paths_In%Desc_Len(iDesc)-1))

      Path_Begin = merge(1, &
                         Paths_In%Desc(iDesc,0:Paths_In%Desc_Len(iDesc)-1), &
                         Paths_In%Desc(iDesc,0:Paths_In%Desc_Len(iDesc)-1) == 0)
      Path_End   = merge(dym_In%nAt, &
                         Paths_In%Desc(iDesc,0:Paths_In%Desc_Len(iDesc)-1), &
                         Paths_In%Desc(iDesc,0:Paths_In%Desc_Len(iDesc)-1) == 0)
      Path_N     = Path_End - Path_Begin + 1_i08

! Debug
!     print *, Path_Begin
!     print *, Path_End
!     print *, Path_N

! Calculate the number of unpruned paths for this descriptor
      nPaths_Unpruned = product(Path_N)

! Debug
!     print *, nPaths_Unpruned

      allocate(Paths_Unpruned(nPaths_Unpruned,0:Paths_In%Desc_Len(iDesc)-1), &
               Paths_Len_Unpruned(nPaths_Unpruned), &
               Paths_Keep(nPaths_Unpruned))

! Generate all possible paths for this descriptor and
! create a mask for prunning
      Paths_Keep = .true.
      do iPath=0,nPaths_Unpruned-1
        Paths_Unpruned(iPath+1,:) = Path_Begin + Loop_Ind(iPath,Path_N)

! Check that there aren't any repeated atoms
        Rep_At = .false.
! NOTE: This piece of code removes paths with any repeated atom. Some people
!       are interested in paths with non-consecutive repeated atoms, so we
!       have modified the code to include those.
!       We also make sure that the first and last atoms are not the same.
!ou1:   do iAtom=0,Paths_In%Desc_Len(iDesc)-2
!         do jAtom=iAtom+1,Paths_In%Desc_Len(iDesc)-1
!           if ( Paths_Unpruned(iPath+1,iAtom) == &
!                Paths_Unpruned(iPath+1,jAtom) ) then
!             Rep_At = .true.
!             exit ou1
!           end if
!         end do
!       end do ou1
        do iAtom=0,Paths_In%Desc_Len(iDesc)-2
          if ( Paths_Unpruned(iPath+1,iAtom) == &
               Paths_Unpruned(iPath+1,iAtom+1) ) then
            Rep_At = .true.
            exit
          end if
        end do
        if ( Paths_Unpruned(iPath+1,0) == &
             Paths_Unpruned(iPath+1,Paths_In%Desc_Len(iDesc)-1) ) then
          Rep_At = .true.
        end if
        if ( Rep_At ) then
          Paths_Keep(iPath+1) = .false.
        else

! Check that the path is not too long
          Length = 0.0_r08
          do iAtom=0,Paths_In%Desc_Len(iDesc)-2
            Length = Length + sqrt(sum( &
                       (dym_In%xyz(Paths_Unpruned(iPath+1,iAtom),:) - &
                        dym_In%xyz(Paths_Unpruned(iPath+1,iAtom+1),:))**2 &
                     ))
          end do
! NOTE: This piece of code removes paths where the accumulated length
!       between the first and last atom is below a threshold. Some users
!       want the path length to be the effective path length, that is, the
!       accumulated length between the start atom until the start atom, divided
!       by 2. We have made this modification in the code.
!         Paths_Len_Unpruned(iPath+1) = Length
          Length = Length + &
            sqrt(sum( &
              (dym_In%xyz(Paths_Unpruned(iPath+1,0),:) - &
               dym_In%xyz(Paths_Unpruned(iPath+1,Paths_In%Desc_Len(iDesc)-1),:)&
              )**2 &
               ))
          Paths_Len_Unpruned(iPath+1) = Length/2.0_r08
          if ( Paths_Len_Unpruned(iPath+1) > Paths_In%Desc_mxR(iDesc) ) then
            Paths_Keep(iPath+1) = .false.
          end if

        end if

      end do

! Debug
!     print *, Paths_Keep
!     stop

! Debug
! Horrible hack to make single atom paths available
! These will be used to compute the u2's
!     Paths_Keep = .true.

! Calculate the number of pruned paths
!     nPaths = count(Paths_Keep)
      Desc_Paths%N = count(Paths_Keep)

! Allocate the pruned paths
      allocate(Desc_Paths%Ind(Desc_Paths%N,0:Paths_In%Desc_Len(iDesc)-1), &
               Desc_Paths%Len(Desc_Paths%N))

! Select the pruned paths using the mask
      jPath = 1
      do iPath=1,nPaths_Unpruned
        if ( Paths_Keep(iPath) ) then
          Desc_Paths%Ind(jPath,:) = Paths_Unpruned(iPath,:)
          Desc_Paths%Len(jPath) = Paths_Len_Unpruned(iPath)
          jPath = jPath + 1
        end if
      end do

      deallocate(Path_Begin,Path_End,Path_N, &
                 Paths_Unpruned,Paths_Len_Unpruned,Paths_Keep)

    end if

  end subroutine Paths_Init

  subroutine Paths_DeInit(Desc_Paths)

    implicit none

    type(Paths), intent(inout)     :: Desc_Paths

    deallocate(Desc_Paths%Ind, &
               Desc_Paths%Len)

  end subroutine Paths_DeInit

  function Loop_Ind(i,n) result(Ind)

    implicit none

    integer(kind=i08), intent(in)               :: i
    integer(kind=i08), dimension(:), intent(in) :: n
    integer(kind=i08), dimension(size(n))       :: Ind

    integer(kind=i08), dimension(size(n))       :: cn
    integer(kind=i08)                           :: ic, ii, p

! Define the nesting order
    p = size(n)

! Compute cumulative products
    cn(p) = 1
    do ii=p-1,1,-1
      cn(ii) = n(ii+1)*cn(ii+1)
    end do

! Copy i to work on it
    ic = i

! Compute indices
    Ind = 0_i08
    do ii=1,p
      Ind(ii) = int(ic/cn(ii))
      ic = mod(ic,cn(ii));
    end do

  end function Loop_Ind

  subroutine Print_dym_Info(dym_In)

! Print the dym data

    implicit none

    type(dym_Info), intent(in)    :: dym_In

    integer(kind=i08) :: iAt, jAt, ip, jq

! Debug
    real(kind=r08) :: Pair_EVec(3,3), Pair_EVal(3)

! Print the dym type flag
    write(IO_Dbg,fmt='(a,i4)') ' dym Type: ', dym_In%Type
    write(IO_Dbg,fmt='(a)') ''

! Print the number of atoms
    write(IO_Dbg,fmt='(a,i4)') ' nAtoms:   ', dym_In%nAt
    write(IO_Dbg,fmt='(a)') ''

! Print the index, atomic number, atomic mass and cartesian coordinates
! If present we also print the cell vectors and the reduced coordinates
! for each reference atom
    if ( allocated(dym_In%an) .and. &
         allocated(dym_In%am) .and. &
         allocated(dym_In%xyz)        ) then

! If we have the reduced coordinates we print in a slightly different way
      if ( .not. allocated(dym_In%red) .and. &
           .not. allocated(dym_In%cell)      ) then

        write(IO_Dbg,fmt='(a)') '      At    At'
        write(IO_Dbg,fmt='(a)') &
          ' Idx Num   Mass                  XYZ'
        do iAt=1,dym_In%nAt
          write(IO_Dbg,fmt='(i4,1x,i3,1x,f8.4,1x,3f10.6)') &
                iAt, dym_In%an(iAt), dym_In%am(iAt), dym_In%xyz(iAt,:)
        end do
        write(IO_Dbg,fmt='(a)') ''

      else

        write(IO_Dbg,fmt='(a)') ' Cell vectors: '
        do ip=1,3
          write(IO_Dbg,fmt='(3f10.6)') dym_In%cell(ip,:)
        end do
        write(IO_Dbg,fmt='(a)') ''

        write(IO_Dbg,fmt='(a)') '      At    At'
        write(IO_Dbg,fmt='(a)') &
          ' Idx Num   Mass                  XYZ                            Red'
        do iAt=1,dym_In%nAt
          write(IO_Dbg,fmt='(i4,1x,i3,1x,f8.4,1x,3f10.6,1x,3f10.6)') &
                iAt, dym_In%an(iAt), &
                dym_In%am(iAt), dym_In%xyz(iAt,:), dym_In%red(iAt,:)
        end do
        write(IO_Dbg,fmt='(a)') ''

      end if

    end if

! Print the block dynamical matrix
    if ( allocated(dym_In%dm_block) ) then

      print *, 'size 1: ', size(dym_In%dm_block,1)
      print *, 'size 2: ', size(dym_In%dm_block,2)
      print *, 'size 3: ', size(dym_In%dm_block,3)
      print *, 'size 4: ', size(dym_In%dm_block,4)
      write(IO_Dbg,fmt='(a)') ' Block dynamical matrix:'
      do iAt=1,dym_In%nAt
        do jAt=1,dym_In%nAt
          write(IO_Dbg,fmt='(2i5)') iAt, jAt
          do ip=1,3
            write(IO_Dbg,fmt='(1p,3e14.6)') &
                  (dym_In%dm_block(iAt,jAt,ip,jq),jq=1,3)
          end do
        end do
      end do
      write(IO_Dbg,fmt='(a)') ''

! Debug
! More detailed printout of the block dym, for testing purposes
!     write(IO_Dbg,fmt='(a)') ' Block dynamical matrix:'
!     do iAt=1,dym_In%nAt
!       do jAt=1,dym_In%nAt
!         write(IO_Dbg,fmt='(2i5)') iAt, jAt
!         do ip=1,3
!           write(IO_Dbg,fmt='(1p,3e14.6)') &
!                 (dym_In%dm_block(iAt,jAt,ip,jq),jq=1,3)
!         end do
!         write(IO_Dbg,fmt='(a)') ''
!         write(IO_Dbg,fmt='(6f12.6)') dym_In%red(iAt,:), dym_In%red(jAt,:)
!         write(IO_Dbg,fmt='(3f12.6,5x,f12.6)') &
!           PBCD(dym_In%red(jAt,:)-dym_In%red(iAt,:)), &
!           sqrt(sum(PBCD(dym_In%red(jAt,:)-dym_In%red(iAt,:))**2))
!         call DSYEVH3(dym_In%dm_block(iAt,jAt,:,:),Pair_EVec,Pair_EVal)
!         write(IO_Dbg,fmt='(a)') ''
!         write(IO_Dbg,fmt='(3f14.6)') Pair_EVal
!         do ip=1,3
!           write(IO_Dbg,fmt='(3f14.6)') &
!                 (Pair_EVec(ip,jq),jq=1,3)
!         end do
!       end do
!     end do
!     write(IO_Dbg,fmt='(a)') ''

    end if

! Print the center of mass of the system
    if ( allocated(dym_In%R_CM) ) then
      write(IO_Dbg,fmt='(a)') ' Center of mass:'
      write(IO_Dbg,fmt='(3e14.6)') dym_In%R_CM
      write(IO_Dbg,fmt='(a)') ''
    end if

! Print the tensor of inertia
    if ( allocated(dym_In%ToI) ) then
      write(IO_Dbg,fmt='(a)') ' Tensor of Inertia:'
      do ip=1,3
        write(IO_Dbg,fmt='(3e14.6)') dym_In%ToI(ip,:)
      end do
      write(IO_Dbg,fmt='(a)') ''
    end if

! Print the moments of inertia
    if ( allocated(dym_In%MoI) ) then
      write(IO_Dbg,fmt='(a)') ' Moments of Inertia:'
      write(IO_Dbg,fmt='(3e14.6)') dym_In%MoI
      write(IO_Dbg,fmt='(a)') ''
    end if

! Print the principal axes or rotation
    if ( allocated(dym_In%PAoR) ) then
      write(IO_Dbg,fmt='(a)') ' Principal Axes of Rotation:'
      do ip=1,3
        write(IO_Dbg,fmt='(3e14.6)') dym_In%PAoR(ip,:)
      end do
      write(IO_Dbg,fmt='(a)') ''
    end if

! Print the whole dynamical matrix
!   if ( allocated(dym_In%dm) ) then

!     write(IO_Dbg,fmt='(a)') &
!           ' Dynamical matrix:'
!     call Prn_Mat(3*dym_In%nAt,3*dym_In%nAt,dym_In%dm)

!   end if

! Print the transformation matrix to remove rotations and translations
    if ( allocated(dym_In%TrfD) ) then

      write(IO_Dbg,fmt='(a)') &
            ' Transformation matrix:'
      call Prn_Mat(3*dym_In%nAt,3*dym_In%nAt,dym_In%TrfD)

    end if

! Print the dipole derivatives (if present)
    if ( allocated(dym_In%Dip_Der) ) then

      write(IO_Dbg,fmt='(a)') ''
      write(IO_Dbg,fmt='(a)') &
            ' Dipole derivatives:'
      write(IO_Dbg,fmt='(a)') &
            '   iAt   p       dmux/dR_ip        dmuy/dR_ip        dmuz/dR_ip'
      do iAt=1,dym_In%nAt
        do ip=1,3
          write(IO_Dbg,fmt='(2i5,1p,3e18.6)') iAt, ip, &
                (dym_In%Dip_Der(jq,iAt,ip),jq=1,3)
        end do
      end do

    end if

! Print out some debugging info
!     write(IO_Dbg,*) 'centeratomindex: ', allocated(dym_In%centeratomindex)
!     write(IO_Dbg,*) 'u_xyz:           ', allocated(dym_In%u_xyz)
!     write(IO_Dbg,*) 'udeg:            ', allocated(dym_In%udeg)
!     write(IO_Dbg,*) 'utype:           ', allocated(dym_In%utype)
!     write(IO_Dbg,*) 'dm:              ', allocated(dym_In%dm)
!     write(IO_Dbg,*) 'TrfD:            ', allocated(dym_In%TrfD)
!     write(IO_Dbg,*) 'dm_Int:          ', allocated(dym_In%dm_Int)

!     integer(kind=i08), allocatable :: centeratomindex(:,:)
!     real(kind=r08), allocatable    :: u_xyz(:,:,:)!srw
!     integer(kind=i08), allocatable :: udeg(:),utype(:)!srw 
!     real(kind=r08), allocatable    :: dm(:,:)
!     real(kind=r08), allocatable    :: TrfD(:,:)
!     real(kind=r08), allocatable    :: dm_Int(:,:)

!     integer(kind=i08)              :: unAt,natom !srw
!     real(kind=r08)                 :: R_CM(3)
!     real(kind=r08)                 :: ToI(3,3)
!     real(kind=r08)                 :: MoI(3)
!     real(kind=r08)                 :: PAoR(3,3)

  end subroutine Print_dym_Info

  subroutine Fix_AN_AM(dym_In)

! Fix the atomic numbers and atomic masses when either one is 0

    implicit none

    type(dym_Info), intent(inout)    :: dym_In

    integer(kind=i08) :: iAt, iElem

! Check if we have any atom where both atomic numers and masses are undefined
    if ( any((dym_In%an <= 0) .and. (dym_In%am < 0.2_r08)) ) then
      write(IO_Err,fmt='(a)') &
            'Error: An atom has undefined atomic number and mass'
      stop
    endif

! Check each atom and determine the undefined quantity
    do iAt=1,dym_In%nAt
      if ( dym_In%an(iAt) <= 0 ) then 
        do iElem=1,PTable_nElem
          if ( ( dym_In%am(iAt) < (PTable_AM(iElem) + 0.2_r08) ) .and. &
               ( dym_In%am(iAt) > (PTable_AM(iElem) - 0.2_r08) ) ) then
            dym_In%an(iAt) = iElem
          end if
        end do
      end if
      if ( dym_In%am(iAt) < 0.2_r08 ) then 
        dym_In%am(iAt) = PTable_AM(dym_In%an(iAt))
      end if
    end do

  end subroutine Fix_AN_AM

  subroutine Write_Feffinp(dym_In,Center,feffName,spectrum,writeHeader,Swap)

! Write a feff input file based on the dym info and a modified version of
! the dym file that matches the dym file

    implicit none

    type(dym_Info), intent(in)    :: dym_In
    integer(kind=i08), intent(in) :: Center
    character(len=*), intent(in)  :: feffName
	character*5, intent(in) :: spectrum
	logical,intent(in) :: writeHeader
    integer(kind=i08), dimension(dym_In%nAt)     :: Swap

    real(kind=r08), dimension(:,:), allocatable  :: xyz
    real(kind=r08), dimension(:), allocatable    :: Dist, Dist_Sorted
    integer(kind=i08)                            :: iAt, nTypes, iType, i
    integer(kind=i08), dimension(:), allocatable :: At_Type, Types
    character(len=2), dimension(:), allocatable  :: At_Symb

! Allocate local variables
    allocate(xyz(dym_In%nAt,3), &
             Dist(dym_In%nAt), &
             At_Type(dym_In%nAt), &
             At_Symb(dym_In%nAt))

! Copy the coordinate to a local variable
    xyz = dym_In%xyz

! Shift the coordinates to the cluster center
    xyz = xyz - spread(xyz(Center,:),1,dym_In%nAt)

! Create a list of distances to the cluster center
    Dist = sqrt(sum(xyz**2,2))

! Allocate arrays
    allocate(Dist_Sorted(dym_In%nAt))

! Sort the distances and get an index swap list
    call Shell_Sort(Dist,Swap,Dist_Sorted)
	
! Set the atomic symbols
    At_Symb = PTable_Symb(dym_In%AN)

! Set the atomic types
    allocate(Types(PTable_nElem))
    Types = 0
    nTypes = 0
    do iAt=1,dym_In%nAt
      if ( Types(dym_In%AN(iAt)) == 0 .and. iAt /= Center ) then
        nTypes =  nTypes + 1
        Types(dym_In%AN(iAt)) = nTypes
      end if
    end do
    At_Type = Types(dym_In%AN)
    At_Type(Center) = 0

! Open the feff input file
    open(unit=IO_Feff,file=feffName,status='unknown')

! Print the header part of the feff.inp file
    if (writeHeader) &
       call Write_Feffinp_Header(IO_Feff,spectrum)

! Print the potentials part of the feff.inp file
    write(IO_Feff,fmt='(a)') 'POTENTIALS'
    write(IO_Feff,fmt='(2i5,3x,a2)') 0, dym_In%AN(Center), At_Symb(Center)
    do iType=1,nTypes
      do iAt=1,dym_In%nAt
        if ( At_Type(iAt) == iType ) then
          write(IO_Feff,fmt='(2i5,3x,a2)') &
                iType, dym_In%AN(iAt), At_Symb(iAt)
          exit
        end if
      end do
    end do
    write(IO_Feff,fmt='(a)') ''

! Print the structure part of the feff.inp file
    write(IO_Feff,fmt='(a)') 'ATOMS'
    do iAt=1,dym_In%nAt
      write(IO_Feff,fmt='(3f11.5,i5,3x,a2,f9.5,i5)') &
            xyz(Swap(iAt),:)/ang2au, At_Type(Swap(iAt)), &
            At_Symb(Swap(iAt)), Dist_Sorted(iAt)/ang2au, iAt-1
    end do
    write(IO_Feff,fmt='(a)') 'END'

! Close the feff input file
    close(IO_Feff)

! Deallocate local variables
    deallocate(xyz,Dist)

  end subroutine Write_Feffinp

  subroutine Shell_Sort(Data_Unsorted,Swap,Data)

! Simple Shellsort implementation
! There are better choices for the increments, could be improved in the future.

    implicit none

    real(kind=r08), dimension(:), intent(in)            :: Data_Unsorted
    integer(kind=i08), dimension(size(Data_Unsorted)), intent(out) :: Swap
    real(kind=r08), dimension(size(Data_Unsorted)), intent(out)    :: Data

    integer(kind=i08) :: iData, jData, nData, Inc
    integer(kind=i08) :: iTemp
    real(kind=r08)    :: Temp
	real(kind=r08), parameter    :: tolerance =1.d-04  !KJ July 2014 to avoid meaningless swapping around precision limit

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
          if ( Data(jData-Inc) <= (Temp + tolerance) ) exit  ! <--- KJ adding tolerance
          Data(jData) = Data(jData-Inc)
          Swap(jData) = Swap(jData-Inc)
          jData = jData - Inc
        end do
        Data(jData) = Temp
        Swap(jData) = iTemp
      end do
      Inc = nint(Inc/2.2_r08)
    end do

  end subroutine Shell_Sort

  subroutine Write_Feffinp_Header(Unit,spectrum)

! Write the header of a feff.inp file making educated guesses for some
! of the values needed. The output should be very similar to that of ATOMS.

    implicit none

    integer(kind=i08), intent(in) :: Unit
	character*5, intent(in) :: spectrum
!   type(dym_Info), intent(in)    :: dym_In
!   character(len=*), intent(in)  :: feffName

! File generation information
    write(Unit,fmt='(a)') &
          ' * This feff9 input file was generated by dym2feffinp'
    write(Unit,fmt='(a)')  ''

! Information for title, this will go into the feff output
    write(Unit,fmt='(a,a)') &
          ' TITLE dymfile name:  ', 'Need to fix'
    write(Unit,fmt='(a,i4)') &
          ' TITLE absorbing atom:', 0
    write(Unit,fmt='(a)')  ''

! Edge and s0^2 information, default to K edges and 1.0 for now
    write(Unit,fmt='(a,a)') &
          ' EDGE      ', 'L3'
    write(Unit,fmt='(a,f6.4)') &
          ' S02       ', 1.0_r08
    write(Unit,fmt='(a)')  ''

! Control and print cards. Default to run everything and extra output for
! the pot module
    write(Unit,fmt='(a)') &
          ' *              pot    xsph     fms   paths  genfmt  ff2chi'
    write(Unit,fmt='(a,6i8)') &
          ' CONTROL   ', 1, 1, 1, 1, 1, 1
    write(Unit,fmt='(a,6i8)') &
          ' PRINT     ', 1, 0, 0, 0, 0, 0
    write(Unit,fmt='(a)')  ''

! Exchange controls. Defaults to 0 (Hedin-Lundqvist)
    write(Unit,fmt='(a)') &
          ' *          ixc  [ Vr  Vi ]'
    write(Unit,fmt='(a,i4)') &
          ' EXCHANGE  ', 0
    write(Unit,fmt='(a)')  ''

! SCF radius. Defaults to 4.0
    write(Unit,fmt='(a)') &
          ' *            r_scf  [ l_scf   n_scf   ca ]'
    write(Unit,fmt='(a,f8.3)') &
          ' SCF       ', 4.0_r08
    write(Unit,fmt='(a)')  ''

! XANES upper limit. Defaults to 4.0
    write(Unit,fmt='(a)') &
          ' *             kmax   [ delta_k  delta_e ]'
	if(spectrum.eq.'XANES') then
    write(Unit,fmt='(a,f8.3)') &
          ' XANES     ', 4.0_r08
	else
    write(Unit,fmt='(a,f8.3)') &
          ' * XANES     ', 4.0_r08
    endif	
    write(Unit,fmt='(a)')  ''

! FMS radius and l_lms. Defaults to 6.0
    write(Unit,fmt='(a)') &
          ' *            r_fms     l_fms'
	if(spectrum.eq.'XANES') then
    write(Unit,fmt='(a,f8.3)') &
          ' FMS       ', 6.0_r08
    else
    write(Unit,fmt='(a,f8.3)') &
          ' * FMS       ', 6.0_r08
	endif
    write(Unit,fmt='(a)')  ''

! LDOS grid. Defaults to a grid from -30.0 eV to 20 eV with a
! broadening of 0.1 eV
    write(Unit,fmt='(a)') &
          ' *             emin    emax   eimag'
    write(Unit,fmt='(a,3f8.3)') &
          ' * LDOS      ', -30.0_r08, 20.0_r08,  0.1_r08
    write(Unit,fmt='(a)')  ''

! EXAFS information
	if(spectrum.eq.'XANES') then
    write(Unit,fmt='(a,f8.3)') &
          ' RPATH     ',  0.1_r08
    write(Unit,fmt='(a,f8.3)') &
          '*EXAFS     ', 20.0_r08
    else
    write(Unit,fmt='(a,f8.3)') &
          ' RPATH     ',  6.0_r08
    write(Unit,fmt='(a,f8.3)') &
          ' EXAFS     ', 20.0_r08
    write(Unit,fmt='(a,i8)') &
          ' NLEG      ', 3
	endif
    write(Unit,fmt='(a)')  ''
	
!  DEBYE information
    write(Unit,fmt='(a)') &
          ' *        Temp  Debye_Temp  DW_Opt  dymfile  DMDW_Order  DMDW_Type  DMDW_Route'
    write(Unit,fmt='(a,f8.3)') &
          ' DEBYE    300.0   315.0  5 feff.dym  16  0  1'	
    write(Unit,fmt='(a)')  ''

  end subroutine Write_Feffinp_Header

  subroutine Write_dym(dym_In,ndymName,Swap)

! Write a dym file where the atoms are swapped

    implicit none

    type(dym_Info), intent(in)                             :: dym_In
    character(len=*), intent(in)                           :: ndymName
    integer(kind=i08), dimension(dym_In%nAt),   intent(in) :: Swap

    integer(kind=i08) :: iAt, jAt, ip, jq

! Open the feff input file
    open(unit=IO_dym,file=ndymName,status='unknown')

! Print the dym type flag
    write(IO_dym,fmt='(i5)') dym_In%Type

! Print the number of atoms
    write(IO_dym,fmt='(i5)') dym_In%nAt

! Print the atomic numbers
    write(IO_dym,fmt='(i5)') dym_In%an(Swap)

! Print the atomic masses
    write(IO_dym,fmt='(f12.6)') dym_In%am(Swap)

! Print the coordinates in cartesian or reduced format, depending on dym type
    if ( dym_In%Type == 4 ) then

! Print the reduced coordinates
      do iAt=1,dym_In%nAt
        write(IO_dym,fmt='(3f14.8)') &
          dym_In%red(Swap(iAt),:) - dym_In%red(Swap(1),:)
      end do

    else

! Print the cartesian coordinates
      do iAt=1,dym_In%nAt
        write(IO_dym,fmt='(3f14.8)') &
          dym_In%xyz(Swap(iAt),:) - dym_In%xyz(Swap(1),:)
      end do

    end if

! Print the block dynamical matrix
    do iAt=1,dym_In%nAt
      do jAt=1,dym_In%nAt
        write(IO_dym,fmt='(2i5)') iAt, jAt
        do ip=1,3
          write(IO_dym,fmt='(1p,3e14.6)') &
                (dym_In%dm_block(Swap(iAt),Swap(jAt),ip,jq),jq=1,3)
        end do
      end do
    end do

! Added by FDV to write the cell vectors, if we have the right dym type
    if ( dym_In%Type == 4 ) then

! Empty separator line
      write(IO_dym,fmt='(a)') ''

! Read the three cell vectors, in atomic units
      do ip=1,3
        write(IO_dym,fmt='(3f14.8)') (dym_In%cell(ip,jq),jq=1,3)
      end do

    end if

! Close the dym file
    close(IO_Dym)

  end subroutine Write_dym

!#######################################################################

!NOTE: This subroutine is incomplete yet and shouldn't be used.

  subroutine Make_TrfD(dym_In)

    implicit none

    type(dym_Info), intent(inout) :: dym_In

    integer(kind=i08) :: iAt, jAt, iCoord, nCoord, iVec
    real(kind=r08), allocatable :: Mass(:)
    real(kind=r08), allocatable :: dm_Int_Full(:,:)

! Calculate the center of mass
    if ( .not. allocated(dym_In%R_CM) ) then
      allocate(dym_In%R_CM(3))
    end if
    call Calc_R_CM(dym_In%xyz,dym_In%am,dym_In%R_CM)

! Debug
!   write(IO_Dbg,fmt='(f14.8)') dym_In%R_CM/ang2au

! Shift origin of cartesian coordinates in dym_In to R_CM
    do iAt=1,dym_In%nAt
      dym_In%xyz(iAt,:) = dym_In%xyz(iAt,:) - dym_In%R_CM
    end do

! Debug
!   write(IO_Dbg,fmt='(3f14.8)') transpose(dym_In%xyz)/ang2au

! Calculate the Tensor of Inertia
    if ( .not. allocated(dym_In%ToI) ) then
      allocate(dym_In%ToI(3,3))
    end if
    call Calc_ToI(dym_In%xyz,dym_In%am,dym_In%ToI)

! Debug
!   write(IO_Dbg,fmt='(3f18.8)') dym_In%ToI*amu2au
!   write(IO_Dbg,fmt='(3f18.8)') dym_In%ToI
!   stop

! Calculate the eigenvectors (Principal Axes of Rotation) and
! eigenvalues (Moments of Inertia) of the system
    if ( .not. allocated(dym_In%PAoR) ) then
      allocate(dym_In%PAoR(3,3))
    end if
    if ( .not. allocated(dym_In%MoI) ) then
      allocate(dym_In%MoI(3))
    end if
! Changed by FDV:
! Changed the routine used here to get better accuracy in these results. It
! should have little effect on the calculation, this is moslty to get exact
! numerical agreement with the octave tools.
    call DSYEVQ3(dym_In%ToI,dym_In%PAoR,dym_In%MoI)
!   call DSYEVH3(dym_In%ToI,dym_In%PAoR,dym_In%MoI)

! Debug
!   write(IO_Dbg,fmt='(3f20.8)') dym_In%MoI*amu2au
!   write(IO_Dbg,fmt='(3f20.8)') dym_In%MoI
!   write(IO_Dbg,fmt='(3f20.8)') transpose(dym_In%PAoR)
!   stop

! Total number of XYZ coordinates
    nCoord = 3*dym_In%nAt

! Build the transformation matrix
    if ( .not. allocated(dym_In%TrfD) ) then
      allocate(dym_In%TrfD(nCoord,nCoord))
    end if

! Initialize all the zero elements
    dym_In%TrfD = 0.0_r08

! Translation vectors
! x coords
    dym_In%TrfD( (/ (iCoord+0*dym_In%nAt,iCoord=1,dym_In%nAt) /) ,1) = &
      sqrt(dym_In%am)
! y coords
    dym_In%TrfD( (/ (iCoord+1*dym_In%nAt,iCoord=1,dym_In%nAt) /) ,2) = &
      sqrt(dym_In%am)
! z coords
    dym_In%TrfD( (/ (iCoord+2*dym_In%nAt,iCoord=1,dym_In%nAt) /) ,3) = &
      sqrt(dym_In%am)
! NOTE: This was my original code. I believe it is wrong because of the way
!       coordinates are ordered in the dm and q vecrors. They should be all x
!       first, then all y and then all z, as done above.
! x coords
!   dym_In%TrfD( (/ (iCoord, iCoord=1,nCoord,3) /) ,1) = sqrt(dym_In%am)
! y coords
!   dym_In%TrfD( (/ (iCoord, iCoord=2,nCoord,3) /) ,2) = sqrt(dym_In%am)
! z coords
!   dym_In%TrfD( (/ (iCoord, iCoord=3,nCoord,3) /) ,3) = sqrt(dym_In%am)

! Debug
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') sqrt(amu2au)*dym_In%TrfD(iCoord,1:3)
!   end do
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,1:3)
!   end do

! Rotation vectors
    do iVec=1,3
      do iAt=1,dym_In%nAt
        dym_In%TrfD( (/ iAt, iAt+1*dym_In%nAt, iAt+2*dym_In%nAt /) ,3+iVec) = &
             Cross(dym_In%PAoR(:,iVec),dym_In%xyz(iAt,:))*sqrt(dym_In%am(iAt));
      end do
    end do

! NOTE: This was my original code. I believe it is wrong because of the way
!       coordinates are ordered in the dm and q vecrors. They should be all x
!       first, then all y and then all z, as done above.
! Rotation vectors
!   do iVec=1,3
!     do iAt=1,dym_In%nAt
!       dym_In%TrfD(3*(iAt-1)+1:3*(iAt-1)+3,3+iVec) = &
!            Cross(dym_In%PAoR(:,iVec),dym_In%xyz(iAt,:))*sqrt(dym_In%am(iAt));
!     end do
!   end do

! Debug
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') sqrt(amu2au)*dym_In%TrfD(iCoord,4:6)
!   end do
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,4:6)
!   end do
!   stop

! Add the rest of the transformation vectors: Use as many bond vectors as
! needed to complete the internal coordinates basis set.

!   iVec = 7
!   iAt  = 1
!   jAt  = 2
!   do

!     dym_In%TrfD(3*(iAt-1)+1:3*(iAt-1)+3,iVec) = &
!          ( dym_In%xyz(jAt,:)*sqrt(dym_In%am(jAt)) - &
!            dym_In%xyz(iAt,:)*sqrt(dym_In%am(iAt)) )

!     dym_In%TrfD(3*(jAt-1)+1:3*(jAt-1)+3,iVec) = &
!         -( dym_In%xyz(jAt,:)*sqrt(dym_In%am(jAt)) - &
!            dym_In%xyz(iAt,:)*sqrt(dym_In%am(iAt)) )

!     if (jAt == dym_In%nAt) then
!       iAt = iAt + 1
!       jAt = iAt + 1
!     else
!       jAt = jAt + 1
!     end if
!     iVec = iVec + 1

!     if ( iVec > 3*dym_In%nAt ) exit

!   end do

! Debug
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') sqrt(amu2au)*dym_In%TrfD(iCoord,7:9)
!   end do

! Normalize non-orthogonal vectors
    do iCoord=1,6 !nCoord JJK - commented out to avoid divide by zero.
         dym_In%TrfD(:,iCoord) = dym_In%TrfD(:,iCoord) / &
            sqrt(sum(dym_In%TrfD(:,iCoord)*dym_In%TrfD(:,iCoord)))
    end do

! Debug
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,1:3)
!   end do
!   write(IO_Dbg,fmt='(a)') ''
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,4:6)
!   end do
!   write(IO_Dbg,fmt='(a)') ''
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,7:9)
!   end do

!   Gram-Schmidt orthogonalization of the coordinate basis set in TrfD

!   call GramSchmidt(dym_In%TrfD,dym_In%TrfD)

! Debug
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,1:3)
!   end do
!   write(IO_Dbg,fmt='(a)') ''
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,4:6)
!   end do
!   write(IO_Dbg,fmt='(a)') ''
!   do iCoord=1,nCoord
!     write(IO_Dbg,fmt='(3f20.12)') dym_In%TrfD(iCoord,7:9)
!   end do

! Mass weights
    allocate(Mass(nCoord))
! x coords
    Mass( (/ (iCoord, iCoord=1,nCoord,3) /) ) = dym_In%am
! y coords
    Mass( (/ (iCoord, iCoord=2,nCoord,3) /) ) = dym_In%am
! z coords
    Mass( (/ (iCoord, iCoord=3,nCoord,3) /) ) = dym_In%am

! Allocate the transformed dynamical matrix
!   if ( .not. allocated(dym_In%dm_Int) ) then
!     allocate(dym_In%dm_Int(nCoord-6,nCoord-6))
!   end if

! Allocate temporary storage
!   allocate(dm_Int_Full(nCoord,nCoord))

! Debug
!   stop

!   do iCoord=1,nCoord
!     do iAt=1,dym_In%nAt
!       dym_In%TrfD(3*(iAt-1)+1:3*(iAt-1)+3,3+iVec) = &
!            Cross(dym_In%PAoR(:,iVec),dym_In%xyz(iAt,:))*sqrt(dym_In%am(iAt));
!     end do
!   end do
!   H_MW

! Debug: Check the E_Vals for a water (3 internal coords) example)

!   deallocate(Mass,dm_Int_Full)

! Debug
! Temporary code!!!!
! As a safety measure against using this routine yet, I deallocate the
! transformation matrix, so noone can use it.
!   deallocate(dym_In%TrfD)

! Debug
! I am attempting to fix problems with the use of this routine when running
! within feff: What happens is that the shift to CM coordinates messes up the
! coincidence between the xyz coords in dym_In and those in the feff input.
! Attempting an unshift to fix the problem
! Unshift origin of cartesian coordinates in dym_In to R_CM
    do iAt=1,dym_In%nAt
      dym_In%xyz(iAt,:) = dym_In%xyz(iAt,:) + dym_In%R_CM
    end do

  end subroutine Make_TrfD

  subroutine Calc_R_CM(xyz,M,R_CM)

    implicit none

    real(kind=r08), intent(in)  :: xyz(:,:)
    real(kind=r08), intent(in)  :: M(:)
    real(kind=r08), intent(out) :: R_CM(3)

    integer :: ip

! Check we have the right sizes
    if ( size(M) /= size(xyz,1) .or. size(xyz,2) /= 3 ) then
      write(IO_Err,fmt='(a)') 'Calc_R_CM: Arrays have inconsistent dimensions'
      stop
    end if

    do ip=1,3
      R_CM(ip) = sum(xyz(:,ip)*M)
    end do

    R_CM = R_CM/sum(M)

  end subroutine Calc_R_CM

  subroutine Calc_ToI(xyz,M,ToI)

    implicit none

    real(kind=r08), intent(in)  :: xyz(:,:)
    real(kind=r08), intent(in)  :: M(:)
    real(kind=r08), intent(out) :: ToI(3,3)

    integer(kind=i08) :: iAt

! Check we have the right sizes
    if ( size(M) /= size(xyz,1) .or. size(xyz,2) /= 3 ) then
      write(IO_Err,fmt='(a)') 'Calc_ToI: Arrays have inconsistent dimensions'
      stop
    end if

    ToI = 0.0_r08
    do iAt=1,size(xyz,1)
      ToI(1,1) = ToI(1,1) + M(iAt)*(xyz(iAt,2)**2+xyz(iAt,3)**2)
      ToI(2,2) = ToI(2,2) + M(iAt)*(xyz(iAt,1)**2+xyz(iAt,3)**2)
      ToI(3,3) = ToI(3,3) + M(iAt)*(xyz(iAt,1)**2+xyz(iAt,2)**2)
      ToI(2,1) = ToI(2,1) - M(iAt)*(xyz(iAt,2)*xyz(iAt,1))
      ToI(3,1) = ToI(3,1) - M(iAt)*(xyz(iAt,3)*xyz(iAt,1))
      ToI(3,2) = ToI(3,2) - M(iAt)*(xyz(iAt,3)*xyz(iAt,2))
    end do

    ToI(1,2) = ToI(2,1);
    ToI(1,3) = ToI(3,1);
    ToI(2,3) = ToI(3,2);

  end subroutine Calc_ToI

      Subroutine Lanczos(Lanc_In,dym_In,lpath,lpath_p,lpath_m,nleg,qj0, &
                         w_pole,wil,mnull,SPole_EinsteinFreq)
      implicit none

      type(Lanczos_Info),  intent(in)  :: Lanc_In
      type(dym_Info),      intent(in)  :: dym_In
      integer(kind=i08),   intent(in)  :: lpath(0:)
      integer(kind=i08),   intent(in)  :: lpath_p(0:nleg-1)
      integer(kind=i08),   intent(in)  :: lpath_m(0:nleg-1)
      integer(kind=i08),   intent(in)  :: nleg
      integer(kind=i08),   intent(out) :: mnull
      real(kind=r08),      intent(out) :: w_pole(0:Lanc_In%nPoles)
      real(kind=r08),      intent(out) :: wil(0:Lanc_In%nPoles)
      real(kind=r08),      intent(out) :: SPole_EinsteinFreq

      integer(kind=i08), parameter :: numxp=100000

      integer(kind=i08) :: nCoor1
      integer(kind=i08) :: inull

      real(kind=r08)    ::    dr(3)
      real(kind=r08)    ::    mu_inv
      real(kind=r08)    ::    adm, om2m, xpeps, xpm, pnxm
      real(kind=r08)    ::    xp, pnx, ratio


      integer(kind=i08) :: m, l
      integer(kind=i08) :: ni

      integer(kind=i08) :: ip,iq,       iAt,jAt,jq

      integer(kind=i08)                                            :: nn
      real(kind=r08)                                               :: mu
      real(kind=r08), dimension(0:dym_In%nAt-1,0:dym_In%nAt-1,0:2) :: rc
      real(kind=r08), dimension(0:3*dym_In%nAt-1),intent(IN)       :: qj0
      real(kind=r08), dimension(:),     allocatable                :: qj, qp, qm
      real(kind=r08), dimension(0:3*dym_In%nAt-1,0:3*dym_In%nAt-1) :: DM_Temp
      real(kind=r08), dimension(:),     allocatable                :: xnull
      real(kind=r08), dimension(:),     allocatable                :: anj, bnj


     allocate(xnull(0:Lanc_In%nPoles),                    &            
              anj(0:Lanc_In%nPoles),                      &
              bnj(0:Lanc_In%nPoles+1),                    &                           
              qj(0:3*dym_In%nAt-1),                       &
              qp(0:3*dym_In%nAt-1),                       &
              qm(0:3*dym_In%nAt-1) )

!intalize DM_Temp
      DM_Temp=dym_In%dm
! Number of coordinates minus 1
    nCoor1 = 3*dym_In%nAt-1

!!

!!
!      Write(100,*) "lpath",lpath
!      Write(100,*) "nleg", nleg 
!      Write(100,*) " "
!      Write(100,*) "qj0", qj0(:)
!      Write(100,*) " "
!      Write(100,*) "lpath_p",lpath_p
!      Write(100,*) "lpath_m",lpath_m
!      write(100,*)
!      write(100,*) "DM_Temp",DM_Temp(:,:)
!      write(100,*)
!      write(100,*) ".dym",dym_In%dm

! start correlated DM
! get sigma for each path
    qj = qj0

! start Lanczos iteration with n=0
! a_0 = <0|dm|0> ;  <0| = qj_trans ;  |0> = qj
    adm = 0.0_r08
    do m=0,nCoor1
      adm = adm + qj(m)*sum(dym_In%dm(:,m)*qj)
    end do
    anj(0) = adm
    bnj(0) = 0.0_r08
!write(100,*) adm,nCoor1

! Debug
!   write(IO_Err,FMT='(A,I3,2F20.10)') 'anj, bnj: ', 0, anj(0), bnj(0)

! Save the single pole frequency for printing later
    SPole_EinsteinFreq = sqrt(anj(0))/(2.0D0*pi)
! Some extra output
!   if ( Lanc_In%IOFlag .gt. 0 ) then
!     write(IO_Out,FMT='(A,$)') ' pDOS Einstein Freq (THz): '
!     write(IO_Out,FMT='(F8.3)') sqrt(anj(0))/(2.0D0*pi)
!   end if

! get |1> = qp = (dm - a0) * qj
    do m=0,nCoor1
      DM_Temp(m,m) = dym_In%dm(m,m) - adm
      qp(m) = sum(DM_Temp(:,m)*qj)
    end do

! normalize |1> by b_1
    bnj(1) = dsqrt(sum(qp**2))
    qp = qp/bnj(1)

! |1> = qp

    do ni=1,Lanc_In%nPoles
    
! iteration step :  |n-1> = |n>  ,   |n> = |n+1>
      qm = qj
      qj = qp

! get a_n ; n = 1, ...
      adm = 0.0_r08
      do m=0,nCoor1
        adm = adm + qj(m)*sum(dym_In%dm(:,m)*qj)
      end do
      anj(ni) = adm

! Debug
!     write(IO_Err,FMT='(A,I3,2F20.10)') 'anj, bnj: ', ni, anj(ni), bnj(ni)

! get |n+1> = qp = (dm - a_n) * |n> - b_n * |n-1>
      do m=0,nCoor1
        DM_Temp(m,m) = dym_In%dm(m,m) - adm
        qp(m) = sum(DM_Temp(:,m)*qj) - bnj(ni) * qm(m)
      end do

! normalize |n+1> by b_n+1
      bnj(ni+1) = dsqrt(sum(qp**2))
      qp = qp/bnj(ni+1)

    end do

! Debug
! Test polynomial derivative
!   xp = anj(0)
!   write(IO_Err,*) xp
!   write(IO_Err,*) PolyD_Y('S',5,xp,anj,bnj)
!   write(IO_Err,*) ( Poly_Y('S',5,xp+0.5e-04_r08,anj,bnj) - &
!                     Poly_Y('S',5,xp-0.5e-04_r08,anj,bnj)     ) / 1.0e-04_r08
!   stop

!   om2m = 3.0_r08 * anj(0)
!   om2m =1000.0_r08 * anj(0)
! This number corresponds to a freq of ~5000 cm-1. No system that I know
! of has a fundamental that is higher.
    om2m = 810000.0_r08


! Debug
!   write(IO_Err,fmt='(a,f20.10)') 'om2m: ', om2m

! Debug
!   xpeps = om2m/(Lanc_In%nPoles*numxp)
    xpeps = 2.0_r08*om2m/(Lanc_In%nPoles*numxp)

! Debug
!   xp = 0.0_r08
!   xp = -om2m
!   xpeps = 2.0_r08*om2m/(Lanc_In%nPoles*numxp)
!   do ip=-Lanc_In%nPoles*numxp,Lanc_In%nPoles*numxp
!     xp = xp + xpeps
!     write(IO_Err,fmt='(2e20.10)') xp, Poly_Y('S',Lanc_In%nPoles,xp,anj,bnj)
!   end do
!   stop

! get Pade polynoms
    inull = 0
    xpm = 0.0_r08
    pnxm = 0.0_r08
! Debug
!   xp = 0.0_r08
    xp = -om2m
    do ip = 1, Lanc_In%nPoles * numxp
      xp = xp + xpeps
      pnx = Poly_Y('S',Lanc_In%nPoles,xp,anj,bnj)

! find zero of pnx
      if (ip > 1) then
        if (pnx /= 0.0D0) then
          if (pnx*pnxm < 0.0_r08) then
            ratio = dabs(pnxm) / (dabs(pnxm) + dabs(pnx))
            xnull(inull) = ratio * (xp - xpm) + xpm
! Debug
!           write(IO_Err,fmt='(a,i5,f20.10)') 'xnull: ', inull, xnull(inull)
            inull = inull + 1
          end if
        else
          xnull(inull) = xp
          inull = inull + 1
        end if
      end if
      xpm = xp
      pnxm = pnx
    end do
! x loop finished
! Check that we have the right number of poles
    if ( inull /= Lanc_In%nPoles ) then
      write(IO_Err,fmt='(a)') &
            ' Warning: The number of poles found is not equal', &
            '          to the number of iterations'
      write(IO_Err,fmt='(a,i5)') 'nPoles = ', inull
!     stop
    end if

    mnull = inull-1


! for each xnull
! get p_l-1(xnull), q_l-1(xnull) and p'_l
    do nn = 0, mnull

      if ( xnull(nn) < 0.0_r08 ) then

!       xnull(nn) =  -xnull(nn)
        w_pole(nn) = -sqrt(-xnull(nn))
!       wil(nn) = 0.0_r08
        wil(nn) = Poly_Y('R',Lanc_In%nPoles,xnull(nn),anj,bnj) / &
                  PolyD_Y('S',Lanc_In%nPoles,xnull(nn),anj,bnj)

! Print error/warning depending on the "gravity" of the problem
! Since these problems can be pretty serious, we send to both the error and
! output streams.
! NOTE: Please do not remove these warnings, they are here for a reason and
!       users should be aware of the problems with their calculation.
        if ( wil(nn) >= 0.05 ) then
          write(IO_Err,FMT='(A)') &
            'Error: Found imaginary frequency with large weight.', &
            '       Will ignore this pole for all properties.',    &
            '       Results are likely wrong.',                    &
            '       Please check guide for possible solutions.'
          write(IO_Err,FMT='(A,I5,2F16.8)') &
            '       Pole, Frq, Wgt: ', nn, w_pole(nn)/(2.0D0*pi), wil(nn)
          write(IO_Out,FMT='(A)') &
            'Error: Found imaginary frequency with large weight.', &
            '       Will ignore this pole for all properties.',    &
            '       Results are likely wrong.',                    &
            '       Please check guide for possible solutions.'
          write(IO_Out,FMT='(A,I5,2F16.8)') &
            '       Pole, Frq, Wgt: ', nn, w_pole(nn)/(2.0D0*pi), wil(nn)
        else if ( (0.01 <= wil(nn)) .and. (wil(nn) <= 0.05) ) then
          write(IO_Err,FMT='(A)') &
            'Warning: Found imaginary frequency with small weight.', &
            '       Will ignore this pole for all properties.',      &
            '       Results should be carefully checked.',           &
            '       Please check guide for details.'
          write(IO_Err,FMT='(A,I5,2F16.8)') &
            '       Pole, Frq, Wgt: ', nn, w_pole(nn)/(2.0D0*pi), wil(nn)
          write(IO_Out,FMT='(A)') &
            'Warning: Found imaginary frequency with small weight.', &
            '       Will ignore this pole for all properties.',      &
            '       Results should be carefully checked.',           &
            '       Please check guide for details.'
          write(IO_Out,FMT='(A,I5,2F16.8)') &
            '       Pole, Frq, Wgt: ', nn, w_pole(nn)/(2.0D0*pi), wil(nn)
        end if

      else

        w_pole(nn) = sqrt(xnull(nn))
        wil(nn) = Poly_Y('R',Lanc_In%nPoles,xnull(nn),anj,bnj) / &
                  PolyD_Y('S',Lanc_In%nPoles,xnull(nn),anj,bnj)

      end if

!     write(IO_Err,FMT='(A,I3,2F20.10)') ' Pol:', &
!           nn, Poly_Y('R',Lanc_In%nPoles,xnull(nn),anj,bnj), &
!               PolyD_Y('S',Lanc_In%nPoles,xnull(nn),anj,bnj)
!     write(IO_Err,FMT='(A,2F20.10)') ' w:', &
!           sqrt(xnull(nn))/(2.0D0*pi), wil(nn)

    end do

! Debug
!   do nn =0,mnull
!     write(6,fmt='(i5,2f20.10)') nn, xnull(nn), wil(nn)
!   end do

!srw now poles and weights have been determined
! Debug
! Scale the low frequency pole
!   if ( q0_opt .and. .not. temp_opt ) then
!     wil(0) = wil(1)/xnull(1)*xnull(0)
! Renormalize (in case we want to)
!     wils = 0.0D0
!     do nn=0,mnull
!       wils = wils + wil(nn,Lanc_In%nPoles)
!     end do
!     do nn=0,mnull
!       wil(nn) = wil(nn)/wils
!       write(IO_Err,FMT='(A,2F20.10)') ' w:', xnull(nn), wil(nn)
!     end do
!   end if

! Debug
! Calculate PDOS moments
!   wmm2 = 0.0D0
!   wmm1 = 0.0D0
!   wm0  = 0.0D0
!   wmp1 = 0.0D0
!   wmp2 = 0.0D0
!   do nn=0,mnull
!     wmm2 = wmm2 + sqrt(xnull(nn))**(-2)*wil(nn)
!     wmm1 = wmm1 + sqrt(xnull(nn))**(-1)*wil(nn)
!     wm0  = wm0  + sqrt(xnull(nn))**( 0)*wil(nn)
!     wmp1 = wmp1 + sqrt(xnull(nn))**(+1)*wil(nn)
!     wmp2 = wmp2 + sqrt(xnull(nn))**(+2)*wil(nn)
!   end do
!   write(IO_Err,FMT='(A,I3,5F12.6)') &
!         ' nPoles, moms:', Lanc_In%nPoles, wmm2, wmm1, wm0, wmp1, wmp2
!
!
! Debug
! Print out the low and high temp limits
!     if ( q0_opt ) then
!       write(IO_Err,FMT='(A,I3,5F12.6)') ' nPoles, sig(0):', &
!             Lanc_In%nPoles, hbarc*caps/(two*am(1)*amu)*wmm1
!       write(IO_Err,FMT='(A,I3,5F12.6)') ' nPoles, sig(700):', &
!             Lanc_In%nPoles, hbarc/hbar*caps*8.617385D-5* &
!                 700.0D0/(am(1)*amu)*wmm2
!     else
!       write(IO_Err,FMT='(A,I3,5F12.6)') ' nPoles, sig(0):',
!    &        Lanc_In%nPoles, hbarc*caps/(two*mu*amu)*wmm1
!       write(IO_Err,FMT='(A,I3,5F12.6)') ' nPoles, sig(700):',
!    &        Lanc_In%nPoles, hbarc/hbar*caps*8.617385D-5*
!    &            700.0D0/(mu*amu)*wmm2
!     end if
!#######################################################################

      return
      End Subroutine

    
!!*** 
!!	added by srw
!!      Last Update 8/25/2011
!!
!!	Name SelfEn
!!
!!	Description
!!	  This subroutine calculates the real and imaginary parts of the phonon contribution to the electron self energy within the Einstein model.
!!       The method is described in Eiguren and Ambrosch-Draxl  PRL 101, 036402 (2008) Complex Quasiparticle Band Structure...
!!       The einstein electron phonon self energy is calculated from Eq. 3 and integrated against the Eliashberg function (A2F) as in Eq. 1 and
!!       the digamma function needed for its evaluation is contained inside this routine.  
!!       The Eliashberb function is calculated in ABINIT+ANADDB and is converted into a pole/weight representation which is read in from the file
!!       A2F.inp.  At prestnt this name is hard coded just prior to the call for this routine.  
!!	 In addition to calculating the self energy, the mass enhancment factor is also calculated and written to the main DMDW output.  This is
!!	 done from the first derivative of the real part of the self energy and is done as a consistency check since values of lambda for several
!!       systems are known experimentally.
!!
!!	Inputs
!!       Lanc_In = Lannczos input variables
!!	 dym_In  = dynamical matrix input variables
!!       a2f     = the Eliashberg function
!!       z       = photo-electron energy (measured down from the fermi level
!!
!!     Output
!!      SE_a2f   = Self Energy (real and imaginary parts)
!!
!!     Notes
!!      All energies are in ev, but are adjusted to meV when printed.  A2F MUST BE READ IN Ha FOR THE TIME BEING 
!!                                                                  (becuase i dont know how to use ABINIT for beans)
!!***
      subroutine SelfEn(Lanc_In,dym_In,z,SE_a2f,a2fpole)           
      implicit none

      type(Lanczos_Info), intent(in)        :: Lanc_In
      type(dym_Info), intent(in)            :: dym_In

      integer                               :: inull,ii
     
      real,dimension(0:Lanc_In%nPoles-1,2),intent(IN)    :: a2fpole
      real                                  :: check
      
      complex, intent(in)                   :: z
      complex                               :: Er0
      complex, dimension(Lanc_In%nT)        :: SE_a2f
!      complex, dimension(Lanc_In%nT)        :: argp,argm,digp,digm,SEei !for my digamma
      complex                               :: argp,argm,digp,digm,SEei,SEei_tmp  !for fancy digamma
      real                                  :: w,T,a2f
      complex                               :: w_cpx,T_cpx,n_cpx
      complex                               :: I,half,kb,pi,two,twopi,zero
       parameter (  I=cmplx(0.0e0,1.0e0),       &
    &            half=cmplx(0.5e0,0.0e0),       &
    &              kb=cmplx(8.6173342e-5,0.0e0),    &
    &              pi=cmplx(3.141592653589792e0,0.0e0), &
    &             two=cmplx(2.0e0,0.0e0),       &
    &           twopi=cmplx(6.283185307179584e0,0.0e0), &
    &            zero=cmplx(0.0e0,0.0e0) )     
      

      !!!open(unit=88,file="SEdiagnostic.info",status="unknown")
      
! Main loop, integral in Eq. 1
      SE_a2f=cmplx(0.0e0,0.0e0)
      SEei=cmplx(0.0e0,0.0e0)
      !!!write(88,*) "Number of poles in A2f rep.",Lanc_In%nPoles, "equals Lanc. order"
      !!!write(88,*) "Energy z",z
      !!!write(88,*)  "Temp", Lanc_In%T      
      !!!write(88,*) " "
!     !!! do inull=0,Lanc_In%nPoles-1
!     !!!  write(88,*) a2fpole(inull,1),a2fpole(inull,2)  
!     !!! End Do
      !!!write(88,*) " "


       !!! sms: absorb kb into T?
      T = Lanc_In%T(1)
      T_cpx = cmplx(T,0.0e0)
!    For use with Lanczos, a2fpole
      Do inull=0,Lanc_In%nPoles-1 !integrate over phonon energies

       !fancy digamma - Fancy digamma agrees with mine to 3 decimals, but is much... much faster.
       w = a2fpole(inull,1)
       w_cpx = cmplx(w,0.0e0)
       a2f = a2fpole(inull,2)
       n_cpx = cmplx(1.0e0/(Exp(w/(real(kb)*T)) - 1.0e0),0.0e0)
       argp = half + I*(w_cpx-z)/(twopi*kb*T_cpx)
       argm = half - I*(w_cpx+z)/(twopi*kb*T_cpx)
       SEei   = CPSI(argp) - CPSI(argm) - I*twopi*(n_cpx+half) 
       SE_a2f = SE_a2f +(SEei*cmplx(a2f,0.0e0)) !dw in weight, in this version

!       Self Energy Diagnostics 
        !IF(abs(real(z)-0.0177).LE.0.001) THEN
         !write(88,*)  "inull",inull,"z",real(z)      
         !write(88,*) aimag(CPSI(argp))+aimag(CPSI(argm))!-real(twopi*half)
         !write(88,*)  "arg/dig p",argp,CPSI(argp)!digp
         !write(88,*)  "arg/dig m",argm,CPSI(argm)!digm
         !write(88,*)  "w/A2F",a2fpole(inull,1),a2fpole(inull,2) 
         !write(88,*)  "Bose term of SE", 1.0e0/(Exp(a2fpole(inull,1)/(8.617e-5*Lanc_In%T)) - 1.0e0)
         !write(88,*)  "SEei",SEei
         !write(88,*)  "SE",SE_a2f
         !write(88,*)  " "
        !END IF

!         open(unit=89,file="reEImode_24meV.out",status="unknown")
!         open(unit=90,file="imEImode_24meV.out",status="unknown")

!          IF(inull.EQ.5)THEN
!           write(89,*) real(z), real(SEei)/1000  !z, which was read in, is in eV 
!           write(90,*) real(z), aimag(SEei)/1000 !SEei are arbitrary.  The 1/1000 is just to make the grid reasonable for when SEh... 
                                                 !report in normalized units or something.  THIS IS ONLY FOR WHEN RUNNING SEei THROUGH SPECTRAL FUNCTION
!           SEei_tmp=SEei/1000
!         End IF

     End Do 
         !!!write(88,*) " "

      !Only use if you want to run SEei through spectral function
!      SE_a2f=SEei_tmp
      return
      end subroutine

      subroutine digamma(dig,z,Lanc_In)
!     DiGamma evaluated by Ab.Stg. 6.3.21 with simpsons rule...
!     Added srw
      Implicit none

      type(Lanczos_Info), intent(in)             :: Lanc_In
      complex,dimension(Lanc_In%nT),intent(in)   :: z
      complex,dimension(Lanc_In%nT),intent(out)  :: dig
      Integer ii,n
      real b,dt,rp,ip
!      real, parameter ::  e= 2.718281828904590!,y=-0.577215664901532
      complex, parameter ::  e= (2.718281828904590,0.0e0)

      n=1.0e5 !choose this to be even;# of sampling points in (0,B)
      dt=5.0e-4
      b=2*n*dt !2 since each n=2 sample points,see below
      
      dig =Exp(-1.0*dt)/dt - Exp(-1.0*z*dt)/(1-e**(-1.0*dt))   &
     &   + Exp(-1.0*b )/b  - Exp(-1.0*z*b) /(1-e**(-1.0*b))     !endpoints
      Do  ii=1,n
       dig=dig +                                               &
     &    4*(Exp(-1.0*(2*ii-1)*dt)/((2*ii-1)*dt) -             &!Odd term
     &       Exp(-1.0*z*(2*ii-1)*dt)/(1-e**(-1.0*(2*ii-1)*dt)))&
     &  + 2*(Exp(-1.0*(2*ii*dt))/(2*ii*dt) -                   &!Even term
     &       Exp(-1.0*z*(2*ii*dt))/(1-e**(-1.0*(2*ii*dt))))

!      dig =e**cmplx((-1.0*dt),0.0e0)/cmplx(dt,0.0e0) - e**(cmplx((-1.0*dt),0.0e0)*z)/(cmplx(1.0e0,0.0e0) - e**cmplx((-1.0*dt),0.0e0) )   &
!     &   + e**(cmplx(-1.0*b,0.0e0) )/cmplx(b,0.0e0 ) - e**(cmplx((-1.0*dt),0.0e0)*z)/(cmplx(1.0e0,0.0e0) - e**cmplx((-1.0*b),0.0e0)  )     !endpoints
!      Do  ii=1,n
!       dig=dig +                                               &
!     &    4*(e**(cmplx(-1.0*(2*ii-1)*dt,0.0e0))  /(cmplx((2*ii-1)*dt,0.0e0)) -             &!Odd term
!     &       e**(cmplx(-1.0*(2*ii-1)*dt,0.0e0)*z)/(cmplx(1.0e0,0.0e0)-e**(cmplx(-1.0*(2*ii-1)*dt,0.0e0))) )&
!     &  + 2*(e**(cmplx(-1.0*(2*ii*dt),0.0e0))/(cmplx(2*ii*dt,0.0e0)) -                     &!Even term
!     &       e**(cmplx(-1.0*(2*ii*dt),0.0e0)*z)/(cmplx(1.0e0,0.0e0)-e**(cmplx(-1.0*(2*ii*dt),0.0e0))))  
      End DO
      dig=dig*dt/3.0e0
      Return
      End Subroutine

! phonon_coupling reads in PHDOS and A2F from ABINIT and determines the electron phonon matrix elements.  
! Note: These elements are integrated over the Fermi surface
! The method is to divide A2F/F = A2 = electron phonon matrix elements.  This is done so that we may use 
! our in house Lanczos determination of F to provide a pole-weight model of A2F  to improve computation time

      Subroutine phonon_coupling(a2fall,a2,norm,Lanc_In)
      implicit none

      real              :: phdos,eli,w,wprev,w0,z
      real(kind=r08), intent(out) :: norm
      real,dimension(2,400),intent(out) :: a2,a2fall
      integer           ::j
  
      !These lines are for plot of a2F*SEei
      type(Lanczos_Info), intent(in)        :: Lanc_In
      complex                               :: argp,argm,digp,digm,SEei  !for fancy digamma
      complex                               :: I,half,kb,pi,two,twopi,zero
       parameter (  I=cmplx(0.0e0,1.0e0),       &
    &            half=cmplx(0.5e0,0.0e0),       &
    &              kb=cmplx(8.617e-5,0.0e0),    &
    &              pi=cmplx(3.1415926e0,0.0e0), &
    &             two=cmplx(2.0e0,0.0e0),       &
    &           twopi=cmplx(6.2831853e0,0.0e0), &
    &            zero=cmplx(0.0e0,0.0e0) )


      Open(unit=57,file=Lanc_In%pds_file,status="old")
      Open(unit=58,file=Lanc_In%a2f_file,status="old")
      Open(unit=59,file="dmdw_A2.dat",status="unknown")

     ! Open(unit=69,file="reSEei.out") 
     ! Open(unit=70,file="imSEei.out") 
     ! Open(unit=71,file="a2f.out") 
     ! Open(unit=72,file="n_w0.out") 
    !read (skip over) headers
     Do j=1,10
      read(57,*)
      read(58,*)
     End Do

     norm=0.0e0
     wprev=0.0e0   
     !read phdos a2f and divide out electron-phonon coupling 
     Do j=1,400
      read(57,*) w, phdos
      read(58,*) w, eli
      norm=norm+phdos*(w-wprev)*27.211396!w's are in Ha from ABINIT
      wprev=w

      a2fall(1,j)=w*27.211396132
      a2fall(2,j)=eli

      a2(1,j)=w
      a2(2,j)=eli/phdos
      write(59,*) a2(1,j),a2(2,j)

      !For a2F*SEei
       !z=0=e_f   
!       SEei=cmplx(0.0e0,0.0e0)

!       z=-1.0e0*w*27.211396132 !now in eV, this is just a way to plot z for the same range of energy as a2f
!       w0=1.5e-2               !this is the chosen phonon energy (w in equation) 
!       IF(w.GT.0.0e0)THEN
!argp = half + I*(cmplx(w0-z,0.0e0))/(twopi*kb*cmplx(Lanc_In%T(1),0.0e0))
!argm = half - I*(cmplx(w0-z,0.0e0))/(twopi*kb*cmplx(Lanc_In%T(1),0.0e0))
!SEei = CPSI(argp) - CPSI(argm) - I*twopi*(cmplx(1.0e0/(Exp(w0/(8.617e-5*Lanc_In%T(1))) - 1.0e0),0.0e0)+half)
!        write(69,*) z*1.0e3, real(SEei)
!        write(70,*) z*1.0e3, aimag(SEei)*-1.0e0
!        write(71,*) a2fall(1,j),a2fall(2,j) 
!        write(72,*) z*1.0e3,real(twopi*(cmplx(1.0e0/(Exp(w0/(8.617e-5*Lanc_In%T(1))) - 1.0e0),0.0e0)+half)) 
!       End IF
     End Do

      close(unit=57)
      close(unit=58)
      close(unit=59)

     ! close(unit=69) 
     ! close(unit=70) 
     ! close(unit=71) 
     ! close(unit=72) 
      return
      End Subroutine



FUNCTION CPSI (ZIN)
!
!! CPSI computes the Psi (or Digamma) function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7C
!***TYPE      COMPLEX (PSI-S, DPSI-D, CPSI-C)
!***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! PSI(X) calculates the psi (or digamma) function of X.  PSI(X)
! is the logarithmic derivative of the gamma function of X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CCOT, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  CPSI
  COMPLEX CPSI
  COMPLEX ZIN, Z, Z2INV, CORR!, CCOT
  DIMENSION BERN(13)
  LOGICAL FIRST
!  EXTERNAL CCOT
  SAVE BERN, PI, NTERM, BOUND, DXREL, RMIN, RBIG, FIRST
  DATA BERN( 1) /   .83333333333333333E-1 /
  DATA BERN( 2) /  -.83333333333333333E-2 /
  DATA BERN( 3) /   .39682539682539683E-2 /
  DATA BERN( 4) /  -.41666666666666667E-2 /
  DATA BERN( 5) /   .75757575757575758E-2 /
  DATA BERN( 6) /  -.21092796092796093E-1 /
  DATA BERN( 7) /   .83333333333333333E-1 /
  DATA BERN( 8) /  -.44325980392156863E0 /
  DATA BERN( 9) /   .30539543302701197E1 /
  DATA BERN(10) /  -.26456212121212121E2 /
  DATA BERN(11) /   .28146014492753623E3 /
  DATA BERN(12) /  -.34548853937728938E4 /
  DATA BERN(13) /   .54827583333333333E5 /
  DATA PI / 3.141592653589793E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  CPSI
  if (FIRST) THEN
     NTERM = -0.30*LOG(R1MACH(3))
! MAYBE BOUND = N*(0.1*EPS)**(-1/(2*N-1)) / (PI*EXP(1))
     BOUND = 0.1171*NTERM*(0.1*R1MACH(3))**(-1.0/(2*NTERM-1))
     DXREL = SQRT(R1MACH(4))
     RMIN = EXP (MAX (LOG(R1MACH(1)), -LOG(R1MACH(2))) + 0.011 )
     RBIG = 1.0/R1MACH(3)
  end if
  FIRST = .FALSE.
!
  Z = ZIN
  X = REAL(Z)
  Y = AIMAG(Z)
  if (Y < 0.0) Z = CONJG(Z)
!
  CORR = (0.0, 0.0)
  CABSZ = ABS(Z)
  if (X >= 0.0 .AND. CABSZ > BOUND) go to 50
  if (X < 0.0 .AND. ABS(Y) > BOUND) go to 50
!
  if (CABSZ < BOUND) go to 20
!
! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
! ABS(AIMAG(Y)) SMALL.
!
  CORR = -PI*CCOT(PI*Z)
  Z = 1.0 - Z
  go to 50
!
! USE THE RECURSION RELATION FOR ABS(Z) SMALL.
!
 20   if (CABSZ  <  RMIN) STOP 'CPSI CALLED WITH Z SO NEAR 0 THAT CPSI OVERFLOWS' !call XERMSG ('SLATEC', 'CPSI', &
!     'CPSI CALLED WITH Z SO NEAR 0 THAT CPSI OVERFLOWS', 2, 2)
!
  if (X >= (-0.5) .OR. ABS(Y) > DXREL) go to 30
  if (ABS((Z-AINT(X-0.5))/X)  <  DXREL) STOP 'SLATEC CPSI ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER'
 
  if (Y  ==  0.0 .AND. X  ==  AINT(X)) STOP 'SLATEC CPSI Z IS A NEGATIVE INTEGER'
!
 30   N = SQRT(BOUND**2-Y**2) - X + 1.0
  DO 40 I=1,N
    CORR = CORR - 1.0/Z
    Z = Z + 1.0
 40   CONTINUE
!
! NOW EVALUATE THE ASYMPTOTIC SERIES FOR SUITABLY LARGE Z.
!
 50   if (CABSZ > RBIG) CPSI = LOG(Z) + CORR
  if (CABSZ > RBIG) go to 70
!
  CPSI = (0.0, 0.0)
  Z2INV = 1.0/Z**2
  DO 60 I=1,NTERM
    NDX = NTERM + 1 - I
    CPSI = BERN(NDX) + Z2INV*CPSI
 60   CONTINUE
  CPSI = LOG(Z) - 0.5/Z - CPSI*Z2INV + CORR
!
 70   if (Y < 0.0) CPSI = CONJG(CPSI)
!
  return
end function



FUNCTION CCOT (Z)
!
!! CCOT computes the cotangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      COMPLEX (COT-S, DCOT-D, CCOT-C)
!***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CCOT(Z) calculates the complex trigonometric cotangent of Z.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CCOT
  COMPLEX CCOT
  COMPLEX Z
  SAVE SQEPS
  DATA SQEPS /0./
!***FIRST EXECUTABLE STATEMENT  CCOT
  if (SQEPS == 0.) SQEPS = SQRT (R1MACH(4))
!
  X2 = 2.0*REAL(Z)
  Y2 = 2.0*AIMAG(Z)
!
  SN2X = SIN (X2)
  !call XERCLR
!
  DEN = COSH(Y2) - COS(X2)
  if (DEN  ==  0.) STOP 'SLATEC CCOT COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI AND Y IS 0)'
!
  if (ABS(DEN) > MAX(ABS(X2),1.)*SQEPS) go to 10
  !call XERCLR
  STOP 'SLATEC CCOT ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR 0 OR PI'
!
 10   CCOT = CMPLX (SN2X/DEN, -SINH(Y2)/DEN)
!
  return
end function
FUNCTION CDCDOT (N, CB, CX, INCX, CY, INCY)
!
!! CDCDOT computes the inner product of two vectors with extended ...
!            precision accumulation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      COMPLEX (SDSDOT-S, CDCDOT-C)
!***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       CB  complex scalar to be added to inner product
!       CX  complex vector with N elements
!     INCX  storage spacing between elements of CX
!       CY  complex vector with N elements
!     INCY  storage spacing between elements of CY
!
!     --Output--
!   CDCDOT  complex dot product (CB if N  <=  0)
!
!     Returns complex result with dot product accumulated in D.P.
!     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY)
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CDCDOT
  COMPLEX CDCDOT
  INTEGER N, INCX, INCY, I, KX, KY
  COMPLEX CX(*), CY(*), CB
  DOUBLE PRECISION DSDOTR, DSDOTI, DT1, DT2, DT3, DT4
!***FIRST EXECUTABLE STATEMENT  CDCDOT
  DSDOTR = DBLE(REAL(CB))
  DSDOTI = DBLE(AIMAG(CB))
  if (N  <=  0) go to 10
  KX = 1
  KY = 1
  if ( INCX < 0) KX = 1+(1-N)*INCX
  if ( INCY < 0) KY = 1+(1-N)*INCY
  DO 5 I = 1,N
    DT1 = DBLE(REAL(CX(KX)))
    DT2 = DBLE(REAL(CY(KY)))
    DT3 = DBLE(AIMAG(CX(KX)))
    DT4 = DBLE(AIMAG(CY(KY)))
    DSDOTR = DSDOTR+(DT1*DT2)-(DT3*DT4)
    DSDOTI = DSDOTI+(DT1*DT4)+(DT3*DT2)
    KX = KX+INCX
    KY = KY+INCY
    5 CONTINUE
   10 CDCDOT = CMPLX(REAL(DSDOTR),REAL(DSDOTI))
  return
end function



FUNCTION R1MACH (I)
!
!! R1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  R1MACH
!
  real r1mach
  INTEGER SMALL(2)
  INTEGER LARGE(2)
  INTEGER RIGHT(2)
  INTEGER DIVER(2)
  INTEGER LOG10(2)
!
  REAL RMACH(5)
  INTEGER IMACH(5)
  SAVE RMACH
!
  EQUIVALENCE (RMACH(1),SMALL(1))
  EQUIVALENCE (RMACH(2),LARGE(1))
  EQUIVALENCE (RMACH(3),RIGHT(1))
  EQUIVALENCE (RMACH(4),DIVER(1))
  EQUIVALENCE (RMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7EFFFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1) / 16#00800000 /
!     DATA LARGE(1) / 16#7FFFFFFF /
!     DATA RIGHT(1) / 16#33800000 /
!     DATA DIVER(1) / 16#34000000 /
!     DATA LOG10(1) / 16#3E9A209B /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA RMACH(1) / Z400800000 /
!     DATA RMACH(2) / Z5FFFFFFFF /
!     DATA RMACH(3) / Z4E9800000 /
!     DATA RMACH(4) / Z4EA800000 /
!     DATA RMACH(5) / Z500E730E8 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
!
!     DATA RMACH(1) / O1771000000000000 /
!     DATA RMACH(2) / O0777777777777777 /
!     DATA RMACH(3) / O1311000000000000 /
!     DATA RMACH(4) / O1301000000000000 /
!     DATA RMACH(5) / O1157163034761675 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA RMACH(1) / Z"3001800000000000" /
!     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
!     DATA RMACH(3) / Z"3FD2800000000000" /
!     DATA RMACH(4) / Z"3FD3800000000000" /
!     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7FFFFFFF' /
!     DATA RMACH(3) / Z'34800000' /
!     DATA RMACH(4) / Z'35000000' /
!     DATA RMACH(5) / Z'3F9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 OR -pd8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CC0000000000000' /
!     DATA RMACH(4) / Z'3CD0000000000000' /
!     DATA RMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC RMACH(5)
!
!     DATA SMALL /    20K,       0 /
!     DATA LARGE / 77777K, 177777K /
!     DATA RIGHT / 35420K,       0 /
!     DATA DIVER / 36020K,       0 /
!     DATA LOG10 / 40423K,  42023K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA RMACH(1) / '00000080'X /
!     DATA RMACH(2) / 'FFFF7FFF'X /
!     DATA RMACH(3) / '00003480'X /
!     DATA RMACH(4) / '00003500'X /
!     DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /

      DATA IMACH(1) / Z'00800000' /
      DATA IMACH(2) / Z'7F7FFFFF' /
      DATA IMACH(3) / Z'33800000' /
      DATA IMACH(4) / Z'34000000' /
      DATA IMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*4 IS THE DEFAULT REAL)
!
!     DATA SMALL(1) / '00800000'X /
!     DATA LARGE(1) / '7F7FFFFF'X /
!     DATA RIGHT(1) / '33800000'X /
!     DATA DIVER(1) / '34000000'X /
!     DATA LOG10(1) / '3E9A209B'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
!     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA RMACH(1) / O402400000000 /
!     DATA RMACH(2) / O376777777777 /
!     DATA RMACH(3) / O714400000000 /
!     DATA RMACH(4) / O716400000000 /
!     DATA RMACH(5) / O776464202324 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1) / 00004000000B /
!     DATA LARGE(1) / 17677777777B /
!     DATA RIGHT(1) / 06340000000B /
!     DATA DIVER(1) / 06400000000B /
!     DATA LOG10(1) / 07646420233B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA SMALL(1) / 1.18E-38      /
!     DATA LARGE(1) / 3.40E+38      /
!     DATA RIGHT(1) / 0.595E-07     /
!     DATA DIVER(1) / 1.19E-07      /
!     DATA LOG10(1) / 0.30102999566 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1) /    8388608 /
!     DATA LARGE(1) / 2147483647 /
!     DATA RIGHT(1) /  880803840 /
!     DATA DIVER(1) /  889192448 /
!     DATA LOG10(1) / 1067065499 /
!
!     DATA RMACH(1) / O00040000000 /
!     DATA RMACH(2) / O17777777777 /
!     DATA RMACH(3) / O06440000000 /
!     DATA RMACH(4) / O06500000000 /
!     DATA RMACH(5) / O07746420233 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /   128,     0 /
!     DATA LARGE(1), LARGE(2) / 32767,    -1 /
!     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
!     DATA DIVER(1), DIVER(2) / 13568,     0 /
!     DATA LOG10(1), LOG10(2) / 16282,  8347 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
!     DATA DIVER(1), DIVER(2) / O032400, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020233 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CA0000000000000' /
!     DATA RMACH(4) / Z'3CB0000000000000' /
!     DATA RMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     DATA RMACH(1) / O000400000000 /
!     DATA RMACH(2) / O377777777777 /
!     DATA RMACH(3) / O146400000000 /
!     DATA RMACH(4) / O147400000000 /
!     DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA SMALL(1), SMALL(2) /     0,    256/
!     DATA LARGE(1), LARGE(2) /    -1,   -129/
!     DATA RIGHT(1), RIGHT(2) /     0,  26880/
!     DATA DIVER(1), DIVER(2) /     0,  27136/
!     DATA LOG10(1), LOG10(2) /  8347,  32538/
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    STOP 'SLATEC R1MACH I OUT OF BOUNDS'
  end if

! Changed by FDV
! R1MACH = RMACH(I)
  R1MACH = REAL(IMACH(I))

  return
end function

      subroutine Prn_Mat(Size,mxSize,A)

!-------------------- Description ------------------------------------
! Prints a general array in sets of nCol columns.
!---------------------------------------------------------------------

      implicit none

      character    Blank
      parameter(Blank     = ' ')

      integer*4 Size,mxSize
      real*8    A(mxSize,mxSize)
      integer*4 i,j,k
      integer*4 nCol
      integer*4 nl,nr

      nCol = 5
      nl = Size/nCol
      nr = mod(Size,nCol)

      do k=1,nl
        write(6,FMT='(6X,A,$)') Blank
        do j=(k-1)*nCol+1,k*nCol
          write(6,FMT='(I7,6X,A,$)') j, Blank
        end do
        write(6,FMT='(A)') Blank
        do i=1,Size
          write(6,FMT='(I5,$)') i
          do j=(k-1)*nCol+1,k*nCol
            write(6,FMT='(1X,E13.6,$)') A(i,j)
          end do
          write(6,FMT='(A)') Blank
        end do
        write(6,FMT='(A)') Blank
      end do

      if (nr .ne. 0) then
        write(6,FMT='(6X,A,$)') Blank
        do j=nl*nCol+1,Size
          write(6,FMT='(I7,6X,A,$)') j, Blank
        end do
        write(6,FMT='(A)') Blank
        do i=1,Size
          write(6,FMT='(I5,$)') i
          do j=nl*nCol+1,Size
            write(6,FMT='(1X,E13.6,$)') A(i,j)
          end do
          write(6,FMT='(A)') Blank
        end do
        write(6,FMT='(A)') Blank
      end if

      end subroutine Prn_Mat

!#######################################################################

  subroutine Add_Cell_Env(dym_In,dym_In_WEnv)

    implicit none

    type(dym_Info), intent(in)  :: dym_In
    type(dym_Info), intent(out) :: dym_In_WEnv

! Check that we have the cell info
    if ( .not. allocated(dym_In%cell) ) then
      write(IO_Err,fmt='(a)') &
        'Warning: Not enough information to add cell environment.', &
        '         Will proceed with dym information in input'
    end if

! To add the environment we use the cell replication routine, adding replicas
! of the input cell in each direction, along each of the axes. The final
! system should be 27 times larger than the original.

! NOTE: this is a tentative interface for this routine, not totally sold on
!       it yet, but it might work for now.
    call Replicate_dym(dym_In,(/1,1,1/),'s',dym_In_WEnv)

  end subroutine Add_Cell_Env

  subroutine Replicate_dym(dym_In,Rep,Symm,dym_Out)

    implicit none

    type(dym_Info),   intent(in)  :: dym_In
    integer,          intent(in)  :: Rep(3)
    character(len=1), intent(in)  :: Symm
    type(dym_Info),   intent(out) :: dym_Out

! This routine takes the information in dym_In and replicates it along
! the a, b anc c axes, according to the repetitions Rep_[]. The Symm
! variable controls if the repetitions are synmmetric around the origin.

    integer(kind=i08)              :: Ind_lb(3), Ind_ub(3)
    integer(kind=i08)              :: nInd, iInd, Ind0, Ind1, Ind2, Ind3
    integer(kind=i08)              :: iAt, jAt, iVec
    integer(kind=i08), allocatable :: Indices(:,:), RepInd(:,:)
    integer(kind=i08), allocatable :: Rep2Ref(:)
    real(kind=r08)                 :: Cell_Red_Sca(3)

! Check we have a dym file that can be repeated
    if ( dym_In%Type /= 4 ) then
      write(IO_Err,fmt='(a)') &
        ' Error in Replicate_dym: wrong type of dym data'
      stop
    end if
 
! Check that the repetitions are meaningful
    if ( any(Rep < 1) ) then
      write(IO_Err,fmt='(a)') &
        ' Error in Replicate_dym: Number of replicas must be larger than 1'
      stop
    end if

! Generate the upper and lower bounds for the indices of replication,
! depending on what symmetry was required.
! NOTE: for now we only go with the simplest possibilities.
    Ind_ub = Rep
    select case ( Symm )
      case ( 's', 'S' )
        Ind_lb = -Ind_ub
      case ( 'n', 'N' )
        Ind_lb = (/ 0, 0, 0 /)
      case default
        write(IO_Err,fmt='(a)') &
          ' Error in Replicate_dym: Unrecognized symmetry directive'
        stop
    end select

! Generate a list of the indices, sorted the "right" way (that is, with the
! (/ 0, 0, 0 /) triad at the top. This way we always end up with the atoms
! of the original structure at the top of the cell. This minimizes any
! confusions.
    nInd = product(Ind_ub-Ind_lb+1)
    allocate(Indices(nInd,3))
    Indices = 0
    iInd = 1
    do Ind1=Ind_lb(1),Ind_ub(1)
      do Ind2=Ind_lb(2),Ind_ub(2)
        do Ind3=Ind_lb(3),Ind_ub(3)
          Indices(iInd,:) = (/ Ind1, Ind2, Ind3 /)
          if ( all(Indices(iInd,:) == (/ 0, 0, 0 /)) ) then
            Ind0 = iInd
          end if
          iInd = iInd + 1
        end do
      end do
    end do

! Debug
!   write(6,fmt='(3i4)') transpose(Indices)
!   write(6,fmt='(a)') ''

! Shuffle things to put the (/ 0, 0, 0 /) at the beginning
    do iInd=Ind0-1,1,-1
      Indices(iInd+1,:) = Indices(iInd,:)
    end do
    Indices(1,:) = (/ 0, 0, 0 /)

! Debug
!   write(6,fmt='(3i4)') transpose(Indices)
!   stop

! Create a vector of reference indices, to simplify things
    allocate(Rep2Ref(dym_In%nAt*nInd))
    do iInd=1,nInd
      Rep2Ref((iInd-1)*dym_In%nAt+1:iInd*dym_In%nAt) = &
        (/ (iAt, iAt=1,dym_In%nAt) /)
    end do

! Debug
!   write(6,fmt='(i4)') Rep2Ref
!   stop

! Now that we have the indices we start by setting up all data we can of the
! replicated dym file

! Set the dym type flag
    dym_Out%Type = dym_In%Type

! Set the number of atoms
    dym_Out%nAt = nInd*dym_In%nAt

! Set the atomic number for each atom
    allocate(dym_Out%an(dym_Out%nAt))
    dym_Out%an = dym_In%an(Rep2Ref)

! Debug
!   write(6,fmt='(i5)') dym_Out%an
!   stop

! Set the atomic mass for each atom
    allocate(dym_Out%am(dym_Out%nAt))
    dym_Out%am = dym_In%am(Rep2Ref)

! Debug
!   write(6,fmt='(f12.5)') dym_Out%am
!   stop

! Total scaling in each directions, to use for rescaling the cell
  Cell_Red_Sca = real(Ind_ub - Ind_lb + 1)

! Set the reduced coordinates and also save an array with the cell indices
! for each replicated atom. This simplifies things later on.
    allocate(dym_Out%red(dym_Out%nAt,3))
    allocate(RepInd(dym_Out%nAt,3))
    do iInd=1,nInd
      do iAt=1,dym_In%nAt
        dym_Out%red((iInd-1)*dym_In%nAt+iAt,:) = &
          (dym_In%red(iAt,:) + real(Indices(iInd,:)))/Cell_Red_Sca
        RepInd((iInd-1)*dym_In%nAt+iAt,:) = Indices(iInd,:)
      end do
    end do

! Scale the cell vectors accordingly
    allocate(dym_Out%cell(3,3))
    do iVec=1,3
      dym_Out%cell(iVec,:) = Cell_Red_Sca(iVec)*dym_In%cell(iVec,:)
    end do

! Set the cartesian coordinates
    allocate(dym_Out%xyz(dym_Out%nAt,3))
    do iAt=1,dym_Out%nAt
      dym_Out%xyz(iAt,:) = matmul(dym_Out%red(iAt,:),dym_Out%cell)
    end do

! Debug
!   write(6,fmt='(3f12.5)') (dym_In%cell(iVec,iVec),iVec=1,3)
!   write(6,fmt='(3f12.5)') (/ 90.0, 90.0, 90.0 /)
!   write(6,fmt='(3f12.5)') transpose(dym_In%red)
!   write(6,fmt='(a)') ''
!   write(6,fmt='(3f12.5)') (dym_Out%cell(iVec,iVec),iVec=1,3)
!   write(6,fmt='(3f12.5)') (/ 90.0, 90.0, 90.0 /)
!   write(6,fmt='(3f12.5)') transpose(dym_Out%red)
!   stop

! Set the block dynamical matrix
    allocate(dym_Out%dm_block(dym_Out%nAt,dym_Out%nAt,3,3))
! NOTE: This is the naive way of generating the replicated block matrix
!   do iAt=1,dym_Out%nAt
!     do jAt=1,dym_Out%nAt
!       dym_Out%dm_block(iAt,jAt,:,:) = &
!         dym_In%dm_block(Rep2Ref(iAt),Rep2Ref(jAt),:,:)
!     end do
!   end do
! Here we try a slightly different way, in a more real-space vein
    dym_Out%dm_block = 0.0_r08
!   do iAt=1,dym_Out%nAt
!     do jAt=1,dym_Out%nAt
!       if ( all(abs(RepInd(jAt,:)-RepInd(iAt,:)) < 2) ) then
!         dym_Out%dm_block(iAt,jAt,:,:) = &
!           dym_In%dm_block(Rep2Ref(iAt),Rep2Ref(jAt),:,:)
!       end if
!     end do
!   end do

! Debug
! Replicate only in the "cell-diagonal", this should give exactly the same
! result as the unreplicated cell
    do iAt=1,dym_Out%nAt
      do jAt=1,dym_Out%nAt
        if ( all(RepInd(jAt,:) == RepInd(iAt,:)) ) then
          dym_Out%dm_block(iAt,jAt,:,:) = &
            dym_In%dm_block(Rep2Ref(iAt),Rep2Ref(jAt),:,:)
        end if
      end do
    end do

! Debug
!   stop

! Now we allocate and fill the array with the indices of the basis atoms.
! These are the atoms that must be used to obtain a complete picture of the
! system. For instance, they would be the unique atoms in a unit cell,
! embedded in a cluster.
! NOTE: In the future we will add a weight or repetition array to do the
!       averaging even simpler.
! NOTE: Since we don't have anything intelligent to do with this yet, we simply
!       allocate and fill the array with a list of all indices for the atoms
!       in the system.
! NOTE: In the case of the replicated cell here is that we keep the indices
!       of the origianl atoms in the cell.
    allocate(dym_Out%BasAtInd(dym_In%nAt))
    dym_Out%BasAtInd = (/ (iAt,iAt=1,dym_In%nAt) /)

  end subroutine Replicate_dym

  function PBC(v) Result(u)

! Take a reduced coordinate and adjust the components so they all lie in the
! [0,1] range (i.e. "standard" reduced coordinates).

    implicit none

    real(kind=r08) :: v(3), u(3)

    u = v - dint(v)
    u = u - merge(1.0_r08,0.0_r08,(u >= 1.0_r08)) + &
            merge(1.0_r08,0.0_r08,(u <  0.0_r08))

  end function PBC

  function PBCD(v) Result(u)

! Take a vector in reduced coordinates and adjust it so it corresponds to the
! minimal reduced distance between images.

    implicit none

    real(kind=r08) :: v(3), u(3)

    u = v - dint(v)
    u = u - merge(1.0_r08,0.0_r08,(u > +0.5_r08)) + &
            merge(1.0_r08,0.0_r08,(u < -0.5_r08))

  end function PBCD

  function PBC_Dist(a,v1,v2) Result(D)

    implicit none

    real(kind=r08) :: a(3), v1(3), v2(3), u(3)
    real(kind=r08) :: D

    u = abs(v1-v2)
    u = u - merge(1.0_r08,0.0_r08,(u >= 0.5_r08))
    D = sqrt(sum((a*u)**2))

  end function PBC_Dist

  subroutine RunTyp_S2(Lanc_In,dym_In,Paths_In)

    implicit none

    type(dym_Info),       intent(in) :: dym_In
    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(Paths_Info),     intent(in) :: Paths_In

    type(DW_Out_Info)                :: DW_Out

    type(Paths)        :: Desc_Paths
    integer(kind=i08)  :: iDesc
    integer(kind=i08)  :: iT
    integer(kind=i08) :: iPath, iAtom

! Debug
    type(Lanczos_Info) :: Lanc_In_Tmp

! Process all the DW factors generated by each descriptor
      do iDesc=1,Paths_In%nDesc

! Added by FDV
! If this is a single atom path, skip. Single atom paths are used for the
! calculation of u2.
        if ( Paths_In%Desc_Len(iDesc) .eq. 1 ) then
          write(IO_Out,fmt='(a)') ' Skipping single atom path descriptor'
          cycle
        end if

! Initialize the paths list for this descriptor
        call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)

! Debug
!       write(IO_Out,*) iDesc
!       write(IO_Out,*) Desc_Paths%N
!       write(IO_Out,*) size(Desc_Paths%Ind,1), size(Desc_Paths%Ind,2)
!       write(IO_Out,*) Desc_Paths%Ind(1,1), Desc_Paths%Ind(1,2)
!       write(IO_Out,*) Desc_Paths%Ind(2,1), Desc_Paths%Ind(2,2)
!       write(IO_Out,*) Desc_Paths%Len

! Process all paths for this descriptor
        do iPath=1,Desc_Paths%N

          write(IO_Out,fmt='(a)') &
                '--------------------------------------------------------------'

          call Calc_DW(Lanc_In, dym_In, &
                       Desc_Paths%Ind(iPath,:), Paths_In%Desc_Len(iDesc), &
                       0,0,DW_Out)

! Debug
! Testing the projected PDOS generation
!         Lanc_In_Tmp = Lanc_In
!         Lanc_In_Tmp%PDOS_Poles = .true.
!         Lanc_In_Tmp%PDOS_Part  = .true.
!         call PrnPlot_PDOS_Poles(Lanc_In_Tmp,DW_Out)

          call Print_DW_Out(Lanc_In,DW_Out)

        end do

! Deallocate stuff
        call Paths_DeInit(Desc_Paths)

        write(IO_Out,fmt='(a)') &
              '--------------------------------------------------------------'
      end do

  end subroutine RunTyp_S2

  subroutine RunTyp_U2(Lanc_In,dym_In,Paths_In)

    implicit none

    type(dym_Info),       intent(in) :: dym_In
    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(Paths_Info),     intent(in) :: Paths_In

    type(DW_Out_Info)                :: DW_Out

    type(Paths)        :: Desc_Paths
    integer(kind=i08)  :: iDesc
    integer(kind=i08)  :: iPath, iAtom
    integer(kind=i08)  :: iu2Pert
    integer(kind=i08)  :: iT

! Process all the u2's generated by each descriptor
! We need to call it three times, for x, y, and z.
    do iDesc=1,Paths_In%nDesc

! Added by FDV
! If this is a scattering path, skip. Scattering paths are used for the
! calculation of s2.
      if ( Paths_In%Desc_Len(iDesc) .gt. 1 ) then
        write(IO_Out,fmt='(a)') ' Skipping scattering path descriptor'
        cycle
      end if

! Initialize the paths list for this descriptor
      call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)

! Debug
!     write(IO_Out,*) iDesc
!     write(IO_Out,*) Desc_Paths%N
!     write(IO_Out,*) size(Desc_Paths%Ind,1), size(Desc_Paths%Ind,2)
!     write(IO_Out,*) Desc_Paths%Ind(:,:)
!     write(IO_Out,*) Desc_Paths%Len

! Process all paths for this descriptor
      do iPath=1,Desc_Paths%N

        write(IO_Out,fmt='(a)') &
              '=============================================================='

! Loop over all three perturbations (x, y, z, or 0, 1, 2)
        do iu2Pert=0,2

          call Calc_DW(Lanc_In, dym_In, &
                       Desc_Paths%Ind(iPath,:), Paths_In%Desc_Len(iDesc), &
                       0,iu2Pert,DW_Out)

          call Print_DW_Out(Lanc_In,DW_Out)

        end do

      end do

! Deallocate stuff
      call Paths_DeInit(Desc_Paths)

    end do
    write(IO_Out,fmt='(a)') &
          '=============================================================='

  end subroutine RunTyp_U2

  subroutine RunTyp_SE(Lanc_In,dym_In,Paths_In)

    implicit none

    type(dym_Info),       intent(in) :: dym_In
    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(Paths_Info),     intent(in) :: Paths_In

    type(DW_Out_Info)                :: DW_Out

    type(Paths)        :: Desc_Paths
    real(kind=r08), dimension(:), allocatable :: vfe, mef, sigc

    call Paths_Init(dym_In,Paths_In,1,Desc_Paths)

    call Calc_DW(Lanc_In, dym_In, &
                 Desc_Paths%Ind(1,:), Paths_In%Desc_Len(1), &
                 0,0,DW_Out)

!   write(IO_Out,*) mef

  end subroutine RunTyp_SE

  subroutine RunTyp_IR(Lanc_In,dym_In,Paths_In)

    implicit none

    type(dym_Info),       intent(in)  :: dym_In
    type(Lanczos_Info),   intent(in)  :: Lanc_In
    type(Paths_Info),     intent(in)  :: Paths_In

    type(DW_Out_Info)                 :: DW_Out

    type(Paths)        :: Desc_Paths
    real(kind=r08), dimension(:), allocatable :: vfe, mef, sigc

! Debug
! Initialize the Paths as for the self-energy (aka, dummy)
    call Paths_Init(dym_In,Paths_In,1,Desc_Paths)

    call Calc_DW(Lanc_In, dym_In, &
                 Desc_Paths%Ind(1,:), Paths_In%Desc_Len(1), &
                 0,0,DW_Out)

    call Print_DW_Out(Lanc_In,DW_Out)

  end subroutine RunTyp_IR

  subroutine RunTyp_VFE(Lanc_In,dym_In,Paths_In_)

    implicit none

    type(dym_Info),       intent(in) :: dym_In
    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(Paths_Info),     intent(in) :: Paths_In_

    type(DW_Out_Info)                :: DW_Out

    type(DW_Out_Info),   allocatable :: All_DW_Outs(:)
    type(Paths_Info)   :: Paths_In
    integer(kind=i08)  :: iPath, iAtom
    integer(kind=i08)  :: iStat
    integer(kind=i08)  :: nTot_VFE
    type(Paths)        :: Desc_Paths
    integer(kind=i08)  :: iDesc
    integer(kind=i08)  :: iu2Pert
    integer(kind=i08)  :: iDW_Out, iT
    real(kind=r08)     :: Tot_VFE(Lanc_In%nT)

! Debug
!   print *, Paths_In_%nDesc
!   print *, allocated(Paths_In_%Desc_Len), Paths_In_%Desc_Len
!   print *, allocated(Paths_In_%Desc), Paths_In_%Desc
!   print *, allocated(Paths_In_%Desc_mxR), Paths_In_%Desc_mxR
!   stop

! Initialize the Paths for the calculation of the VFE:
    if ( Paths_In_%nDesc == 0 ) then

! Deallocate everything, just in case
      deallocate(Paths_In%Desc_Len, &
                 Paths_In%Desc,     &
                 Paths_In%Desc_mxR, stat=iStat)

! If there are no descriptors, we create a new Paths_In with information for:
!   1) If dym type == 4: all basis atoms in the cell.
!   2) If dym type == 2: all unique atoms. NOTE IMPLEMENTED YET
!   For anything else:   ALL atoms in the system.
      select case ( dym_In%Type )

        case ( 4 )

          Paths_In%nDesc = size(dym_In%BasAtInd)
          allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
                   Paths_In%Desc(Paths_In%nDesc,0:0), &
                   Paths_In%Desc_mxR(Paths_In%nDesc)    )

          Paths_In%Desc_Len  = 1
          Paths_In%Desc(:,0) = dym_In%BasAtInd
          Paths_In%Desc_mxR  = 1.0_r08

! NOTE: This case is not implemented yet cause I don't have an easy way to
!       assign what the unique atom indices are based on the dym info. Will
!       figure it out later. For now just default to all atoms.

!       case ( 2 )

!         Paths_In%nDesc = size(dym_In%BasAtInd)
!         if ( .not. allocated(Paths_In%Desc_Len) ) then
!           allocate(Paths_In%Desc_Len(Paths_In%nDesc))
!         end if
!         if ( .not. allocated(Paths_In%Desc) ) then
!           allocate(Paths_In%Desc(Paths_In%nDesc,0:0))
!         end if
!         if ( .not. allocated(Paths_In%Desc_mxR) ) then
!           allocate(Paths_In%Desc_mxR(Paths_In%nDesc))
!         end if

!         Paths_In%Desc_Len  = 1
!         Paths_In%Desc(:,0) = dym_In%BasAtInd
!         Paths_In%Desc_mxR  = 1.0_r08

        case default

          Paths_In%nDesc = dym_In%nAt
          allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
                   Paths_In%Desc(Paths_In%nDesc,0:0), &
                   Paths_In%Desc_mxR(Paths_In%nDesc)    )

          Paths_In%Desc_Len  = 1
          Paths_In%Desc(:,0) = (/ (iAtom,iAtom=1,dym_In%nAt) /)
          Paths_In%Desc_mxR  = 1.0_r08

      end select

    else

! Check that everything is properly allocated and has the right dimensions
! This is a bit of overkill. We should just go with it.
      if ( (.not. allocated(Paths_In_%Desc_Len)) .or. &
           (.not. allocated(Paths_In_%Desc))     .or. &
           (.not. allocated(Paths_In_%Desc_mxR))        ) then
        write(IO_Err,fmt='(a)') &
          ' Error in RunTyp_PDOS: Unallocated Path_In arrays'
        stop
      end if

! Since only single atom "path" descriptors are allowed, here we remove all
! the other types.
      Paths_In%nDesc = count( Paths_In_%Desc_Len == 1 )
      allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
               Paths_In%Desc(Paths_In%nDesc,0:0), &
               Paths_In%Desc_mxR(Paths_In%nDesc)    )

      Paths_In%Desc_Len  = 1
      Paths_In%Desc(:,0) = pack(Paths_In_%Desc(:,0),Paths_In_%Desc_Len == 1)
      Paths_In%Desc_mxR  = 1.0_r08

    end if

! Debug
!   print *, Paths_In%nDesc
!   print *, allocated(Paths_In%Desc_Len), Paths_In%Desc_Len
!   print *, allocated(Paths_In%Desc), Paths_In%Desc
!   print *, allocated(Paths_In%Desc_mxR), Paths_In%Desc_mxR
!   stop

! Since we have no idea how many VFE components we will get for this Path_In,
! we do a dry run pass to get a count
    nTot_VFE = 0
    do iDesc=1,Paths_In%nDesc
      call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)
      nTot_VFE = nTot_VFE + Desc_Paths%N
      call Paths_DeInit(Desc_Paths)
    end do

! Take into account the three directional perturbations
    nTot_VFE = 3*nTot_VFE

! Allocate the array that will save all the outputs
    allocate(All_DW_Outs(nTot_VFE))

! Process all the VFE generated by each descriptor
! We need to call Calc_DW three times, for x, y, and z.
    iDW_Out = 1
    do iDesc=1,Paths_In%nDesc

! Initialize the paths list for this descriptor
      call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)

! Debug
!     write(IO_Out,*) iDesc
!     write(IO_Out,*) Desc_Paths%N
!     write(IO_Out,*) size(Desc_Paths%Ind,1), size(Desc_Paths%Ind,2)
!     write(IO_Out,*) Desc_Paths%Ind(:,:)
!     write(IO_Out,*) Desc_Paths%Len

! Process all paths for this descriptor
      do iPath=1,Desc_Paths%N

        write(IO_Out,fmt='(a)') &
              '=============================================================='
! Loop over all three perturbations (x, y, z, or 0, 1, 2)
        do iu2Pert=0,2

          call Calc_DW(Lanc_In, dym_In, &
                       Desc_Paths%Ind(iPath,:), Paths_In%Desc_Len(iDesc), &
                       0,iu2Pert,DW_Out)

          call Print_DW_Out(Lanc_In,DW_Out)

! Save this path's results for complete printout at the end
          All_DW_Outs(iDW_Out) = DW_Out
          iDW_Out = iDW_Out + 1

        end do

      end do

! Now we print the total VFE for these selection of paths/atoms
      Tot_VFE = 0.0_r08
      do iDW_Out=1,nTot_VFE
        Tot_VFE = Tot_VFE + All_DW_Outs(iDW_Out)%vfe
      end do

      write(IO_Out,fmt='(a)') &
            '=============================================================='
      write(IO_Out,fmt='(a)') ''
      write(IO_Out,fmt='(a)') ' Total VFE for the paths requested:'
      write(IO_Out,fmt='(a)') ''

!     write(IO_Out,fmt='(a)') ' Temp (K)        VFE (J/mol-c)'
      write(IO_Out,fmt='(a)') ' Temp (K)        VFE (eV)'
      do iT=1,Lanc_In%nT
        write(IO_Out,fmt='(f8.2,a5,f16.6)') Lanc_In%T(iT), ' ', &
              Tot_VFE(iT)/Jpmol2eV
      end do

! Deallocate stuff
      call Paths_DeInit(Desc_Paths)

    end do

  end subroutine RunTyp_VFE

  subroutine RunTyp_PDOS(Lanc_In,dym_In,Paths_In_)

    implicit none

    type(dym_Info),       intent(in) :: dym_In
    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(Paths_Info),     intent(in) :: Paths_In_

    type(DW_Out_Info)                :: DW_Out

    type(DW_Out_Info),   allocatable :: All_DW_Outs(:)
    type(DW_Out_Info)                :: DW_Out_Tot
    type(Paths_Info)   :: Paths_In
    type(Paths)        :: Desc_Paths
    integer(kind=i08)  :: iDesc
    integer(kind=i08)  :: iStat
    integer(kind=i08)  :: iPath, iAtom
    integer(kind=i08)  :: iPole, jPole
    integer(kind=i08)  :: iu2Pert
    integer(kind=i08)  :: iT
    integer(kind=i08)  :: nTot_pDOS, iDW_Out, nTot_Poles
    character          :: u2Pert_Label(0:2) = (/ 'x', 'y', 'z' /)

    integer(kind=i08), allocatable :: Swap(:)

! Debug
!   print *, Paths_In_%nDesc
!   print *, allocated(Paths_In_%Desc_Len), Paths_In_%Desc_Len
!   print *, allocated(Paths_In_%Desc), Paths_In_%Desc
!   print *, allocated(Paths_In_%Desc_mxR), Paths_In_%Desc_mxR
!   stop

! Initialize the Paths for the calculation of the PDOS:
    if ( Paths_In_%nDesc == 0 ) then

! Deallocate everything, just in case
      deallocate(Paths_In%Desc_Len, &
                 Paths_In%Desc,     &
                 Paths_In%Desc_mxR, stat=iStat)

! If there are no descriptors, we create a new Paths_In with information for:
!   1) If dym type == 4: all basis atoms in the cell.
!   2) If dym type == 2: all unique atoms. NOTE IMPLEMENTED YET
!   For anything else:   ALL atoms in the system.
      select case ( dym_In%Type )

        case ( 4 )

          Paths_In%nDesc = size(dym_In%BasAtInd)
          allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
                   Paths_In%Desc(Paths_In%nDesc,0:0), &
                   Paths_In%Desc_mxR(Paths_In%nDesc)    )

          Paths_In%Desc_Len  = 1
          Paths_In%Desc(:,0) = dym_In%BasAtInd
          Paths_In%Desc_mxR  = 1.0_r08

! NOTE: This case is not implemented yet cause I don't have an easy way to
!       assign what the unique atom indices are based on the dym info. Will
!       figure it out later. For now just default to all atoms.

!       case ( 2 )

!         Paths_In%nDesc = size(dym_In%BasAtInd)
!         if ( .not. allocated(Paths_In%Desc_Len) ) then
!           allocate(Paths_In%Desc_Len(Paths_In%nDesc))
!         end if
!         if ( .not. allocated(Paths_In%Desc) ) then
!           allocate(Paths_In%Desc(Paths_In%nDesc,0:0))
!         end if
!         if ( .not. allocated(Paths_In%Desc_mxR) ) then
!           allocate(Paths_In%Desc_mxR(Paths_In%nDesc))
!         end if

!         Paths_In%Desc_Len  = 1
!         Paths_In%Desc(:,0) = dym_In%BasAtInd
!         Paths_In%Desc_mxR  = 1.0_r08

        case default

          Paths_In%nDesc = dym_In%nAt
          allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
                   Paths_In%Desc(Paths_In%nDesc,0:0), &
                   Paths_In%Desc_mxR(Paths_In%nDesc)    )

          Paths_In%Desc_Len  = 1
          Paths_In%Desc(:,0) = (/ (iAtom,iAtom=1,dym_In%nAt) /)
          Paths_In%Desc_mxR  = 1.0_r08

      end select

    else

! Check that everything is properly allocated and has the right dimensions
! This is a bit of overkill. We should just go with it.
      if ( (.not. allocated(Paths_In_%Desc_Len)) .or. &
           (.not. allocated(Paths_In_%Desc))     .or. &
           (.not. allocated(Paths_In_%Desc_mxR))        ) then
        write(IO_Err,fmt='(a)') &
          ' Error in RunTyp_PDOS: Unallocated Path_In arrays'
        stop
      end if

! Since only single atom "path" descriptors are allowed, here we remove all
! the other types.
      Paths_In%nDesc = count( Paths_In_%Desc_Len == 1 )
      allocate(Paths_In%Desc_Len(Paths_In%nDesc), &
               Paths_In%Desc(Paths_In%nDesc,0:0), &
               Paths_In%Desc_mxR(Paths_In%nDesc)    )

      Paths_In%Desc_Len  = 1
      Paths_In%Desc(:,0) = pack(Paths_In_%Desc(:,0),Paths_In_%Desc_Len == 1)
      Paths_In%Desc_mxR  = 1.0_r08

    end if

! Debug
!   print *, Paths_In%nDesc
!   print *, allocated(Paths_In%Desc_Len), Paths_In%Desc_Len
!   print *, allocated(Paths_In%Desc), Paths_In%Desc
!   print *, allocated(Paths_In%Desc_mxR), Paths_In%Desc_mxR
!   stop

! Since we have no idea how many partial DOSs we will get for this Path_In,
! we do a dry run pass to get a count
    nTot_pDOS = 0
    do iDesc=1,Paths_In%nDesc
      call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)
      nTot_pDOS = nTot_pDOS + Desc_Paths%N
      call Paths_DeInit(Desc_Paths)
    end do

! Take into account the three directional perturbations
    nTot_pDOS = 3*nTot_pDOS

! Allocate the array that will save all the outputs
    allocate(All_DW_Outs(nTot_pDOS))

! Process all the pDOS generated by each descriptor
! We need to call Calc_DW three times, for x, y, and z.
    iDW_Out = 1
    do iDesc=1,Paths_In%nDesc

! Initialize the paths list for this descriptor
      call Paths_Init(dym_In,Paths_In,iDesc,Desc_Paths)

! Debug
!     write(IO_Out,*) iDesc
!     write(IO_Out,*) Desc_Paths%N
!     write(IO_Out,*) size(Desc_Paths%Ind,1), size(Desc_Paths%Ind,2)
!     write(IO_Out,*) Desc_Paths%Ind(:,:)
!     write(IO_Out,*) Desc_Paths%Len

! Process all paths for this descriptor
      do iPath=1,Desc_Paths%N

        write(IO_Out,fmt='(a)') &
              '=============================================================='
! Loop over all three perturbations (x, y, z, or 0, 1, 2)
        do iu2Pert=0,2

          call Calc_DW(Lanc_In, dym_In, &
                       Desc_Paths%Ind(iPath,:), Paths_In%Desc_Len(iDesc), &
                       0,iu2Pert,DW_Out)

          call Print_DW_Out(Lanc_In,DW_Out)

! Save this path's results for complete printout at the end
          All_DW_Outs(iDW_Out) = DW_Out
          iDW_Out = iDW_Out + 1

        end do

      end do

! Deallocate stuff
      call Paths_DeInit(Desc_Paths)

    end do
    write(IO_Out,fmt='(a)') &
          '=============================================================='

! If requested, print the individual components
    if ( Lanc_In%PDOS_Part ) then
      do iDW_Out=1,nTot_pDOS
        if ( Lanc_In%PDOS_Poles ) then
          call PrnPlot_PDOS_Poles(Lanc_In,All_DW_Outs(iDW_Out))
        end if
        if ( Lanc_In%PDOS_Rect ) then
          call PrnPlot_PDOS_Rect(Lanc_In,All_DW_Outs(iDW_Out))
        end if
        if ( Lanc_In%PDOS_Gauss ) then
          call PrnPlot_PDOS_Gauss(Lanc_In,All_DW_Outs(iDW_Out))
        end if
      end do
    end if

! Generate the total density of states and print
! To make the total density of states we collect all the poles and then
! renormalize their weights so they add up to 1 again.

! First count how many poles we have in total
! NOTE: This could be done in an easier way if we assume that all calculations
!       have the same number of poles, but I might want to generalize things
!       in the future, so I will do it the long way.
    nTot_Poles = 0
    do iDW_Out=1,nTot_pDOS
      nTot_Poles = nTot_Poles + size(All_DW_Outs(iDW_Out)%Poles_Frq)
    end do

! Now we set the contents of the total output that actually make sense
    DW_Out_Tot%Path_nAt  = 0
    DW_Out_Tot%Path_Pert = -1
    DW_Out_Tot%Path_Len  = 0.0_r08
    allocate(DW_Out_Tot%Poles_Frq(nTot_Poles), &
             DW_Out_Tot%Poles_Wgt(nTot_Poles)  )
    iPole = 1
    DW_Out_Tot%SPole_Frq = 0.0_r08
    do iDW_Out=1,nTot_pDOS
      do jPole=1,size(All_DW_Outs(iDW_Out)%Poles_Frq)
        DW_Out_Tot%Poles_Frq(iPole) = All_DW_Outs(iDW_Out)%Poles_Frq(jPole)
        DW_Out_Tot%Poles_Wgt(iPole) = All_DW_Outs(iDW_Out)%Poles_Wgt(jPole)
        iPole = iPole + 1
      end do
      DW_Out_Tot%SPole_Frq = DW_Out_Tot%SPole_Frq + &
                             All_DW_Outs(iDW_Out)%SPole_Frq
    end do
    DW_Out_Tot%Poles_Wgt = DW_Out_Tot%Poles_Wgt/sum(DW_Out_Tot%Poles_Wgt)
    DW_Out_Tot%SPole_Frq = DW_Out_Tot%SPole_Frq/nTot_pDOS

! For some of the PDOS plotting routines to work properly, the poles need to
! be sorted from lowest to highest freq.
    allocate(Swap(nTot_Poles))
    call Shell_Sort(DW_Out_Tot%Poles_Frq,Swap,DW_Out_Tot%Poles_Frq)
    DW_Out_Tot%Poles_Wgt = DW_Out_Tot%Poles_Wgt(Swap)

! Here we collect degenerate or near-degenerate poles to avoid problems with
! some of the PDOS printing formats (poles and rectangles)
    call Rmv_Deg_Poles(DW_Out_Tot)

    call Print_DW_Out(Lanc_In,DW_Out_Tot)
!   stop

    if ( Lanc_In%PDOS_Poles ) then
      call PrnPlot_PDOS_Poles(Lanc_In,DW_Out_Tot)
    end if
    if ( Lanc_In%PDOS_Rect ) then
      call PrnPlot_PDOS_Rect(Lanc_In,DW_Out_Tot)
    end if
    if ( Lanc_In%PDOS_Gauss ) then
      call PrnPlot_PDOS_Gauss(Lanc_In,DW_Out_Tot)
    end if

! Debug
! Test the different ways of printing the pDOS
!   call PrnPlot_PDOS_Rect(Lanc_In,All_DW_Outs(1))
!   call PrnPlot_PDOS_Poles(Lanc_In,All_DW_Outs(1))
!   call PrnPlot_PDOS_Gauss(Lanc_In,All_DW_Outs(1))

! Debug
!   do iDW_Out=1,nTot_pDOS
!     call Print_DW_Out(Lanc_In,All_DW_Outs(iDW_Out))
!   end do

  end subroutine RunTyp_PDOS

  subroutine Print_DW_Out(Lanc_In,DW_Out)

    implicit none

    type(Lanczos_Info),   intent(in) :: Lanc_In
    type(DW_Out_Info),    intent(in) :: DW_Out

    integer(kind=i08) :: iPole, nPoles, iAtom, iT
    integer(kind=i08) :: nMom
    real(kind=r08)    :: Mom, Corr

! Make sure this is DW_Out doesn't contain a total PDOS mashup.
    if ( DW_Out%Path_nAt > 0 ) then

      if ( Lanc_In%RunTyp == 0 ) then
        write(IO_Out,fmt='(a,$)') ' Path Indices: '
      end if
      if ( (Lanc_In%RunTyp == 1) .or. &
           (Lanc_In%RunTyp == 3) .or. &
           (Lanc_In%RunTyp == 5)      ) then
        write(IO_Out,fmt='(a,$)') ' Atom Index: '
      end if
      do iAtom=1,DW_Out%Path_nAt
        write(IO_Out,fmt='(i4,$)') DW_Out%Path(iAtom)
      end do
      write(IO_Out,fmt='(a)') ''
      if ( (Lanc_In%RunTyp == 1) .or. &
           (Lanc_In%RunTyp == 3) .or. &
           (Lanc_In%RunTyp == 5)      ) then
        write(IO_Out,fmt='(a,a,a)') '--------- Direction ', &
              u2Pert_Label(DW_Out%Path_Pert), &
              ' ---------------------------'
      end if

    else

      write(IO_Out,fmt='(a)') ' Total PDOS results: '

    end if

! Adding some extra useful output
    if ( Lanc_In%IOFlag .gt. 0 ) then

! Check that both arrays have the same size
      if ( size(DW_Out%Poles_Frq) /= size(DW_Out%Poles_Wgt) ) then
        write(IO_Err,fmt='(a)') &
          ' Error in Print_DW_Out: ', &
          ' Freq and weight arrays have different size'
        stop
      end if

! Number of poles
      nPoles = size(DW_Out%Poles_Frq)

! Print the pole positions and weights
      write(IO_Out,FMT='(A)') ' PDOS Poles:'
      write(IO_Out,FMT='(A)') '     Freq. (THz)    Weight'
      do iPole=1,nPoles
        write(IO_Out,FMT='(A5,F8.3,F18.9)') ' ', &
              DW_Out%Poles_Frq(iPole), DW_Out%Poles_Wgt(iPole)
      end do
      write(IO_Out,FMT='(A)') ''

! Print the single pole Einstein frequency
      if ( DW_Out%Path_nAt > 0 ) then
        write(IO_Out,FMT='(A)') &
          ' PDOS Einstein freq (single pole), associated temp and eff. force constant: '
      else
        write(IO_Out,FMT='(A)') &
          ' PDOS Einstein freq (avg. of single poles), associated temp and eff. force constant: '
      end if
      write(IO_Out,FMT='(A)') ' Freq (THz)   Temp (K)   Eff. FC (N/m)'
      write(IO_Out,FMT='(F8.3,A5,F8.2,A1,F12.4)') &
        DW_Out%SPole_Frq, ' ', DW_Out%SPole_Frq*47.9908741942_r08, &
        ' ', DW_Out%RedMass*(2.0_r08*pi*DW_Out%SPole_Frq)**2*0.00166053873000_r08
      write(IO_Out,FMT='(A)') ''

! Print print some of the important moments and their associated Einstein freqs
      write(IO_Out,FMT='(A)') &
        ' pDOS n Moments, associated Einstein freqs, temps and eff. force constants:'
      write(IO_Out,FMT='(A)') &
        '  n     Mom (THz^n)   Freq (THz)     Temp (K)    Eff. FC (N/m)'
      do nMom=-2,2
        Mom  = 0.0_r08
        Corr = 0.0_r08
! We ignore any imaginary frequencies in the calculation of the moments.
! This should make clear how bad the results are regarding the imag. frqs.
        do iPole=1,nPoles
          if ( DW_Out%Poles_Frq(iPole) > 0.0_r08 ) then
            Mom = Mom + DW_Out%Poles_Frq(iPole)**(nMom)*DW_Out%Poles_Wgt(iPole)
          else
            Corr = Corr + DW_Out%Poles_Wgt(iPole)
          end if
        end do
! Apply a correction so the moments are normalized to 1 after removing the
! imaginary freqs.
        Mom = Mom/(1-Corr)
        if ( nMom /= 0 ) then
          write(IO_Out,FMT='(I3,2F14.5,A5,F8.2,A1,F12.4)') nMom, &
                Mom, Mom**(1.0_r08/nMom), ' ', &
                Mom**(1.0_r08/nMom)*47.9908741942_r08, ' ', &
                DW_Out%RedMass*(2.0_r08*pi*Mom**(1.0_r08/nMom))**2*0.00166053873000_r08
        else
          write(IO_Out,FMT='(I3,F14.5,A14,A5,A8)') nMom, &
                Mom, '     ---------', ' ', &
                '--------'
        end if
      end do
      write(IO_Out,FMT='(A)') ''

    end if

    if ( Lanc_In%RunTyp == 0 ) then

! Add printout of path reduced mass
      write(IO_Out,fmt='(a,$)') ' Path Red. Mass (AMU):'
      write(IO_Out,fmt='(f12.6)') DW_Out%RedMass

! Print the results, depending on how many temperatures we have
      if ( Lanc_In%nT .eq. 1 ) then
        write(IO_Out,fmt='(a,$)') ' Path Length (Ang), s^2 (1e-3 Ang^2):'
        if ( allocated(DW_Out%s2) ) then
          write(IO_Out,fmt='(f8.4,f9.4)') &
            DW_Out%Path_Len/ang2au, DW_Out%s2*1000
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%s2'
          stop
        end if
      else
        write(IO_Out,fmt='(a,$)') ' Path Length (Ang):'
        write(IO_Out,fmt='(f8.4)') DW_Out%Path_Len/ang2au
        write(IO_Out,fmt='(a)') ' Temp (K)   s^2 (1e-3 Ang^2)'
        if ( allocated(DW_Out%s2) ) then
          do iT=1,Lanc_In%nT
            write(IO_Out,fmt='(f8.2,a5,f8.4)') Lanc_In%T(iT), ' ', &
                  DW_Out%s2(iT)*1000
          end do
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%s2'
          stop
        end if
      end if

    end if

    if ( Lanc_In%RunTyp == 1 ) then

! Print the results, depending on how many temperatures we have
      if ( Lanc_In%nT .eq. 1 ) then
!       write(IO_Out,fmt='(a,$)') ' VFE (J/mol-c):'
        write(IO_Out,fmt='(a,$)') ' VFE (eV):'
        if ( allocated(DW_Out%vfe) ) then
          write(IO_Out,fmt='(f16.6)') DW_Out%vfe/Jpmol2eV
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%vfe'
          stop
        end if
      else
!       write(IO_Out,fmt='(a)') ' Temp (K)        VFE (J/mol-c)'
        write(IO_Out,fmt='(a)') ' Temp (K)        VFE (eV)'
        if ( allocated(DW_Out%vfe) ) then
          do iT=1,Lanc_In%nT
            write(IO_Out,fmt='(f8.2,a5,f16.6)') Lanc_In%T(iT), ' ', &
                  DW_Out%vfe(iT)/Jpmol2eV
          end do
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%vfe'
          stop
        end if
      end if

    end if

    if ( Lanc_In%RunTyp == 3 ) then

! Print the results, depending on how many temperatures we have
      if ( Lanc_In%nT .eq. 1 ) then
        write(IO_Out,fmt='(a,$)') ' u^2 (1e-3 Ang^2):'
        if ( allocated(DW_Out%u2) ) then
          write(IO_Out,fmt='(f8.3)') DW_Out%u2*1000
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%u2'
          stop
        end if
      else
        write(IO_Out,fmt='(a)') ' Temp (K)   u^2 (1e-3 Ang^2)'
        if ( allocated(DW_Out%u2) ) then
          do iT=1,Lanc_In%nT
            write(IO_Out,fmt='(f8.2,a5,f7.3)') Lanc_In%T(iT), ' ', &
                  DW_Out%u2(iT)*1000
          end do
        else
          write(IO_Err,fmt='(a)') &
            'Error in Print_DW_Out: Unallocated DW_Out%u2'
          stop
        end if
      end if

    end if

    if ( Lanc_In%RunTyp == 5 ) then
      write(IO_Out,fmt='(a)') ' Projected DOS component computed.'
      write(IO_Out,fmt='(a)') ''
    end if

  end subroutine Print_DW_Out

  subroutine PrnPlot_PDOS_Rect(Lanc_In,DW_Out)

    implicit none

    type(Lanczos_Info), intent(in)           :: Lanc_In
    type(DW_Out_Info),  intent(in)           :: DW_Out

    integer(kind=i08) :: iPole, nPoles, iFrq, nFrq
    real(kind=r08), dimension(:), allocatable :: Frq, Wgt

    integer(kind=i08)  :: iAt
    character(len=300) :: Label, NumStr

! Number of digits for atom index padding.
! NOTE: This cold be set in a more intelligent way, using a fixed number for now
!       to simplify coding
    integer            :: nPad = 3

! NOTE: When using the rect. format, the resulting plot WILL NOTE BE CORRECT if
! poles are right on top of each other (degenerate). The degenerate
! poles should be eliminated from DW_Out before calling this routine.

! Initialize label with the format string
    Label = 'rect'

! If this is an atomic path, use a poles.nnnn.x format
    if ( DW_Out%Path_nAt < 1 ) then
      Label = trim(Label) // '.' // 'tot'
    elseif ( DW_Out%Path_nAt == 1 ) then
      NumStr = pad(nPad,DW_Out%Path(1))
      NumStr = adjustl(NumStr)
      Label = trim(Label) // '.' // &
              trim(NumStr) // '.' // u2Pert_Label(DW_Out%Path_Pert)
    else
! Otherwise use a poles.iiii.jjjj. ... .kkkk format
      do iAt=1,DW_Out%Path_nAt
        NumStr = pad(nPad,DW_Out%Path(iAt))
        NumStr = adjustl(NumStr)
        Label = trim(Label) // '.' // &
                trim(NumStr)
      end do
    end if
    Label = adjustl(Label)

    if ( Lanc_In%PDOS_DropL ) then
      Label = trim(Label) // '.wdl'
    end if

! Open the unit for printing the PDOS in this format
    call DMDW_Open_PDOS(trim(Label))

! Number of poles
    nPoles = size(DW_Out%Poles_Frq)

    do iPole=1,nPoles
      write(IO_PDOS,fmt='(a,2f12.6)') &
        '#', DW_Out%Poles_Frq(iPole), DW_Out%Poles_Wgt(iPole)
    end do
    write(IO_PDOS,fmt='(a)') ''

! This is gonna be ugly code, but I'll clean it up later.
    if ( DW_Out%Poles_Frq(1) > 0.0_r08 ) then

      nFrq = nPoles+2
      allocate(Frq(0:nFrq))

! Initialize, just in case
      Frq = 0.0_r08

      Frq(0) = 0.0_r08
      Frq(1) = DW_Out%Poles_Frq(1)/2.0_r08
      do iPole=2,nPoles
        Frq(iPole) = (DW_Out%Poles_Frq(iPole)+DW_Out%Poles_Frq(iPole-1))/2.0_r08
      end do
      Frq(nPoles+1) = 2.0_r08*DW_Out%Poles_Frq(nPoles) - Frq(nPoles)
      Frq(nPoles+2) = 3.0_r08*DW_Out%Poles_Frq(nPoles) - 2.0_r08*Frq(nPoles)

    else

      nFrq = nPoles+3
      allocate(Frq(0:nFrq))

! Initialize, just in case
      Frq = 0.0_r08

      Frq(0) = 2.0_r08*DW_Out%Poles_Frq(1) - DW_Out%Poles_Frq(2)
      Frq(1) = DW_Out%Poles_Frq(1) - &
               (DW_Out%Poles_Frq(2)-DW_Out%Poles_Frq(1))/2.0_r08
      iFrq = 1
      do iPole=1,nPoles-1
        if ( (DW_Out%Poles_Frq(iPole)   < 0.0_r08) .and. &
             (DW_Out%Poles_Frq(iPole+1) > 0.0_r08)       ) then
          iFrq = iFrq + 1
          Frq(iFrq) = &
            (DW_Out%Poles_Frq(iPole))/2.0_r08
          iFrq = iFrq + 1
          Frq(iFrq) = &
            (DW_Out%Poles_Frq(iPole+1))/2.0_r08
        else
          iFrq = iFrq + 1
          Frq(iFrq) = &
            (DW_Out%Poles_Frq(iPole)+DW_Out%Poles_Frq(iPole+1))/2.0_r08
        end if
      end do
      iFrq = iFrq + 1
      Frq(iFrq) = 2.0_r08*DW_Out%Poles_Frq(nPoles) - Frq(nPoles)
      iFrq = iFrq + 1
      Frq(iFrq) = 3.0_r08*DW_Out%Poles_Frq(nPoles) - 2.0_r08*Frq(nPoles)

    end if

    allocate(Wgt(0:nFrq))

! Initialize, just in case
    Wgt = 0.0_r08
    do iFrq=0,nFrq-1
      Wgt(iFrq) = 0.0_r08
      do iPole=1,nPoles
        if ( (Frq(iFrq) < DW_Out%Poles_Frq(iPole))   .and. &
             (DW_Out%Poles_Frq(iPole) < Frq(iFrq+1))       ) then
          Wgt(iFrq) = Wgt(iFrq) + &
                      DW_Out%Poles_Wgt(iPole)/(Frq(iFrq+1)-Frq(iFrq))
        end if
      end do
    end do

!   do iFrq=0,nFrq
!     write(IO_PDOS,fmt='(2f16.10)') Frq(iFrq), Wgt(iFrq)
!   end do
    do iFrq=0,nFrq-1
      if ( Lanc_In%PDOS_DropL ) then
        write(IO_PDOS,fmt='(2f16.10)') Frq(iFrq),   0.0_r08
      end if
      write(IO_PDOS,fmt='(2f16.10)') Frq(iFrq),   Wgt(iFrq)
      write(IO_PDOS,fmt='(2f16.10)') Frq(iFrq+1), Wgt(iFrq)
      if ( Lanc_In%PDOS_DropL ) then
        write(IO_PDOS,fmt='(2f16.10)') Frq(iFrq+1), 0.0_r08
      end if
    end do

    call DMDW_Close_PDOS

    deallocate(Frq,Wgt)

  end subroutine PrnPlot_PDOS_Rect

  subroutine PrnPlot_PDOS_Poles(Lanc_In,DW_Out)

    implicit none

    type(Lanczos_Info), intent(in) :: Lanc_In
    type(DW_Out_Info),  intent(in) :: DW_Out

    integer(kind=i08)  :: iPole, nPoles

    integer(kind=i08)  :: iAt
    character(len=300) :: Label, NumStr

! Number of digits for atom index padding.
! NOTE: This cold be set in a more intelligent way, using a fixed number for now
!       to simplify coding
    integer            :: nPad = 3

! NOTE: When using a pole format, the resulting plot is not very informative if
! certain poles are right on top of each other (degenerate). The degenerate
! poles should be eliminated from DW_Out before calling this routine.

! Initialize label with the format string
    Label = 'poles'

! If this is an atomic path, use a poles.nnnn.x format
    if ( DW_Out%Path_nAt < 1 ) then
      Label = trim(Label) // '.' // 'tot'
    elseif ( DW_Out%Path_nAt == 1 ) then
      NumStr = pad(nPad,DW_Out%Path(1))
      NumStr = adjustl(NumStr)
      Label = trim(Label) // '.' // &
              trim(NumStr) // '.' // u2Pert_Label(DW_Out%Path_Pert)
    else
! Otherwise use a poles.iiii.jjjj. ... .kkkk format
      do iAt=1,DW_Out%Path_nAt
        NumStr = pad(nPad,DW_Out%Path(iAt))
        NumStr = adjustl(NumStr)
        Label = trim(Label) // '.' // &
                trim(NumStr)
      end do
    end if
    Label = adjustl(Label)

! Open the unit for printing the PDOS in this format
    call DMDW_Open_PDOS(trim(Label))

! Number of poles
    nPoles = size(DW_Out%Poles_Frq)

    do iPole=1,nPoles
      write(IO_PDOS,fmt='(a,2f12.6)') &
        '#', DW_Out%Poles_Frq(iPole), DW_Out%Poles_Wgt(iPole)
    end do
    write(IO_PDOS,fmt='(a)') ''

! Print the low "anchor point"
    if ( DW_Out%Poles_Frq(1) > 0.0_r08 ) then
      write(IO_PDOS,fmt='(2f12.6)') 0.0_r08, 0.0_r08
    else
      write(IO_PDOS,fmt='(2f12.6)') &
        DW_Out%Poles_Frq(1) - &
          (DW_Out%Poles_Frq(2)-DW_Out%Poles_Frq(1))/2.0_r08, &
        0.0_r08
    end if

! Print all the pole points
    do iPole=1,nPoles
      write(IO_PDOS,fmt='(2f12.6)') &
        DW_Out%Poles_Frq(iPole), 0.0_r08
      write(IO_PDOS,fmt='(2f12.6)') &
        DW_Out%Poles_Frq(iPole), DW_Out%Poles_Wgt(iPole)
      write(IO_PDOS,fmt='(2f12.6)') &
        DW_Out%Poles_Frq(iPole), 0.0_r08
    end do

! Print the high "anchor point"
    write(IO_PDOS,fmt='(2f12.6)') &
      DW_Out%Poles_Frq(nPoles) + &
        (DW_Out%Poles_Frq(nPoles)-DW_Out%Poles_Frq(nPoles-1))/2.0_r08, &
      0.0_r08

    call DMDW_Close_PDOS

  end subroutine PrnPlot_PDOS_Poles

  subroutine PrnPlot_PDOS_Gauss(Lanc_In,DW_Out)

    implicit none

    type(Lanczos_Info), intent(in) :: Lanc_In
    type(DW_Out_Info),  intent(in) :: DW_Out

    integer(kind=i08) :: iPole, nPoles
    integer(kind=i08) :: iFrq, nFrq
    real(kind=r08)    :: Frq, Wgt
    real(kind=r08)    :: FrMin, FrMax
    real(kind=r08)    :: HArg, BArg

! This parameter determines the width of the plot edges, to ensure
! that all relevant DOS is plotted
    integer(kind=i08), parameter :: nBroad = 6

! sqrt(pi/ln 2)
    real(kind=r08), parameter :: spln2 = 2.12893403886245247140_r08

! Define ln 2
    real(kind=r08), parameter :: ln2   = 0.69314718055994528623_r08

    integer(kind=i08)  :: iAt
    character(len=300) :: Label, NumStr

! Number of digits for atom index padding.
! NOTE: This cold be set in a more intelligent way, using a fixed number for now
!       to simplify coding
    integer            :: nPad = 3

! Initialize label with the format string
    Label = 'gaussian'

! If this is an atomic path, use a poles.nnnn.x format
    if ( DW_Out%Path_nAt < 1 ) then
      Label = trim(Label) // '.' // 'tot'
    elseif ( DW_Out%Path_nAt == 1 ) then
      NumStr = pad(nPad,DW_Out%Path(1))
      NumStr = adjustl(NumStr)
      Label = trim(Label) // '.' // &
              trim(NumStr) // '.' // u2Pert_Label(DW_Out%Path_Pert)
    else
! Otherwise use a poles.iiii.jjjj. ... .kkkk format
      do iAt=1,DW_Out%Path_nAt
        NumStr = pad(nPad,DW_Out%Path(iAt))
        NumStr = adjustl(NumStr)
        Label = trim(Label) // '.' // &
                trim(NumStr)
      end do
    end if
    Label = adjustl(Label)

! Open the unit for printing the PDOS in this format
    call DMDW_Open_PDOS(trim(Label))

! Number of poles
    nPoles = size(DW_Out%Poles_Frq)

    do iPole=1,nPoles
      write(IO_PDOS,fmt='(a,2f12.6)') &
        '#', DW_Out%Poles_Frq(iPole), DW_Out%Poles_Wgt(iPole)
    end do
    write(IO_PDOS,fmt='(a)') ''

! First we find the limits of the frequency grid to be used
    FrMin = minval(DW_Out%Poles_Frq) - nBroad*Lanc_In%PDOS_Broad
    FrMax = maxval(DW_Out%Poles_Frq) + nBroad*Lanc_In%PDOS_Broad

! Unless we have negative freqs, we want it to start from 0
    FrMin = min(0.0_r08,FrMin)

! Compute arguments for the normal distribution
    BArg = ln2/(Lanc_In%PDOS_Broad/2)**2
    HArg = 1.0_r08/((Lanc_In%PDOS_Broad/2)*spln2)

! Number of frq points to print
    nFrq = int((FrMax-FrMin)/Lanc_In%PDOS_Res)
    do iFrq=0,nFrq+1
      Frq = FrMin + Lanc_In%PDOS_Res*iFrq
      Wgt = 0.0_r08
      do iPole=1,nPoles
        Wgt = Wgt + HArg*DW_Out%Poles_Wgt(iPole) * &
                    exp(-BArg*(Frq-DW_Out%Poles_Frq(iPole))**2)
      end do
      write(IO_PDOS,fmt='(2f12.6)') Frq, Wgt
    end do

    call DMDW_Close_PDOS

  end subroutine PrnPlot_PDOS_Gauss

  subroutine Rmv_Deg_Poles(DW_Out)

    implicit none

    type(DW_Out_Info),  intent(inout) :: DW_Out

    type(DW_Out_Info) :: DW_Out_Old
    real(kind=r08)    :: Freq_Thresh = 0.0001_r08
    integer(kind=i08) :: Types(size(DW_Out%Poles_Frq))

    integer(kind=i08), allocatable :: Idxs(:)

    integer(kind=i08) :: iPole, nPoles, iType, nPoles_Smp

! NOTE: This is temporary code, we will make something more intellingent later.

! Save the input into a working version
    DW_Out_Old = DW_Out

    nPoles = size(DW_Out%Poles_Frq)
    iType    = 1
    Types(1) = iType
    do iPole=2,nPoles
      if ( abs(DW_Out%Poles_Frq(iPole-1) - &
               DW_Out%Poles_Frq(iPole)) > Freq_Thresh ) then
        iType = iType + 1
      end if
      Types(iPole) = iType
! Debug
!     print *, DW_Out%Poles_Frq(iPole), Types(iPole)
    end do

    nPoles_Smp = maxval(Types)

! Debug
!   print *, nPoles_Smp

! Deallocate the old stuff in DW_Out and allocate the new sizes
    deallocate(DW_Out%Poles_Frq,DW_Out%Poles_Wgt)
    allocate(DW_Out%Poles_Frq(nPoles_Smp),DW_Out%Poles_Wgt(nPoles_Smp))

    do iType=1,nPoles_Smp
      allocate(Idxs(size(pack( (/ (iPole,iPole=1,nPoles) /), Types == iType))))
      Idxs = pack( (/ (iPole,iPole=1,nPoles) /), Types == iType)
!     print *, Idxs
      DW_Out%Poles_Wgt(iType) = sum(DW_Out_Old%Poles_Wgt(Idxs))
      DW_Out%Poles_Frq(iType) = &
        sum(DW_Out_Old%Poles_Frq(Idxs)*DW_Out_Old%Poles_Wgt(Idxs)) / &
        DW_Out%Poles_Wgt(iType)
      deallocate(Idxs)
    end do

  end subroutine Rmv_Deg_Poles

end module m_DMDW
