!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: reldos.f90,v $:
! $Revision: 1.11 $
! $Author: jorissen $
! $Date: 2012/01/30 22:04:08 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reldos    !KJ I put everything in modules 7-09

  use DimsMod, only: nheadx, nspx=>nspu
  use controls
  use struct,nphkevin=>nph
  use kklist
  use strfacs,only :  streta,strrmax,strgmax
  use constants
  use atoms_inp
  use ldos_inp
  use reciprocal_inp
  use potential_inp
  use hubbard_inp
  use xsph_inp
  implicit none

  ! Local stuff
  character*512 :: slog
  character*80  :: head(nheadx)
  integer :: lhead(nheadx)

  ! Added to satisfy implicit none
  integer :: i,iat,ios,iph,i1b,idum,j,nhead
  real*8  :: celvin

  call atoms_read    ! Read  geom.inp file
  call ldos_read     !     read ldos.inp
  call potential_read !KJ I want to know nohole
  call xsph_read
  call hubbard_read

  ! !KJ Next section added for k-space calculations
  call init_controls
  call reciprocal_read(celvin)  ! read reciprocal.inp file
  if (nohole.ge.0) corehole=.false. !KJ if NOHOLE 2, corehole needs to be used in fms but *not* in ldos !!
!  call init_struct(nph)
  if(ispace.eq.0) then
     !KJ next lines : initialize nsp in the struct module (old routines will use same value reinitialized by fmstot).
     nsp = 1
     if (abs(ispin).eq.1 ) nsp = nspx
     lpot(0:nph)=lmaxph(0:nph)
  endif

  makekmeshnow=.true.
  if(ispace.eq.0)  then !KJ go on and read the k-mesh!
     nphkevin=nph
     a1=a1/bohr  ! lattice constants in bohr
     a2=a2/bohr
     a3=a3/bohr
     celvin=celvin/(bohr**3)
     call crystalstructure(celvin)

     if(makekmeshnow) then
        call kmesh
     else
        !c  klist.inp  !KJ added 8/06
        open(3,file='klist.inp',form='formatted',status='old')
        read(3,*)
        read(3,*) nkp,usesym,ktype
        call init_kklist(nkp,nsym) !KJ 6-09
        read(3,*)
        do i=1,nkp
           read(3,*) bk(:,i),weight(i)
        enddo
        close(3)
        !   k-mesh in fractional units 0-1 for strfac
        call wlog('assuming a klist.inp in Angstrom')
        call wlog('update these instructions !!')
        do i=1,3
           bk(i,:)=bk(i,:) * alat(i)/(dble(2)*pi) *bohr  !remove *bohr if klist.inp is in a.u.
        enddo
     endif

  endif
  ! !KJ end my changes

  !     transform to code units (bohrs and hartrees - atomic units)
  rfms2 = rfms2 / bohr
  rdirec = rdirec / bohr
  emin  = emin  / hart
  emax  = emax  / hart
  eimag = eimag / hart
  rat = rat / bohr

  return
end subroutine reldos




subroutine changeklist(ktab,nktab,weight)

  implicit none
  integer nktab
  real*8 ktab(3,nktab),weight(nktab)
  real*8, dimension(3) :: brx,bry,brz,bgx,bgy,bgz
  real*8 a,zetax,zetay,zetaz,f1,f2,f3,alat,vol
  integer i,iq,il

  !KJ debugging : substitute the k-mesh for another one


  weight=dble(1)/DBLE(NKTAB)

  ! SIMPLE CUBIC :
  brx(:)=dble(0) ; bry(:)=dble(0) ; brz = dble(0)
  brx(1)=dble(1) ; bry(2)=dble(1) ; brz(3)=dble(1)
  alat=dble(4.5)

  ! Construct the reciprocal lattice - primitive vectors (BGX,BGY,BGZ) of reciprocal space
  DO I = 1,3
     Iq = 1 + MOD(I,3)
     Il = 1 + MOD(Iq,3)
     BGX(I) = BRY(Iq)*BRZ(Il) - BRZ(Iq)*BRY(Il)
     BGY(I) = BRZ(Iq)*BRX(Il) - BRX(Iq)*BRZ(Il)
     BGZ(I) = BRX(Iq)*BRY(Il) - BRY(Iq)*BRX(Il)
  END DO
  VOL = DABS(BRX(1)*BGX(1)+BRY(1)*BGY(1)+BRZ(1)*BGZ(1))
  DO I = 1,3
     BGX(I) = BGX(I)/VOL
     BGY(I) = BGY(I)/VOL
     BGZ(I) = BGZ(I)/VOL
  END DO

  ! Now construct the k-mesh
  ! ------------------------------------- initialize set up of weyl-k-mesh
  A = 0.359D0
  ZETAX = DSQRT(3.0D0)*A
  ZETAY = DSQRT(5.0D0)*A
  ZETAZ = 2.0D0*DSQRT(13.0D0)*A

  DO I = 1,NKTAB
     !
     ! --------------------------------------- create ->k in  units of  2PI/A
     ! -----------------------------------------------------  -0.5 < F  < 0.5
     F1 = I*ZETAX - IDINT(I*ZETAX) - 0.5D0
     F2 = I*ZETAY - IDINT(I*ZETAY) - 0.5D0
     F3 = I*ZETAZ - IDINT(I*ZETAZ) - 0.5D0

     !! For regular k-mesh :
     !         KTAB(1,I) = dble(0)
     !         KTAB(2,I) = dble(0)
     !         KTAB(3,I) = dble(-0.5)+(i-1)/dble(nktab-1)
     ! For Weyl-k-mesh :
     KTAB(1,I) = F1*BGX(1) + F2*BGX(2) + F3*BGX(3)
     KTAB(2,I) = F1*BGY(1) + F2*BGY(2) + F3*BGY(3)
     KTAB(3,I) = F1*BGZ(1) + F2*BGZ(2) + F3*BGZ(3)



  END DO



  return
end subroutine changeklist !subroutine changeklist
