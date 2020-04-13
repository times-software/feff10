!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: reafms.f90,v $:
! $Revision: 1.15 $
! $Author: jorissen $
! $Date: 2012/01/30 22:04:08 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reafms  !KJ 7-09 I put everything in modules.

  use DimsMod, only: nspx=>nspu
  use controls
  use struct, nphkevin=>nph
  use kklist
  use strfacs
  use constants,only: pi,bohr
  use global_inp
  use atoms_inp
  use fms_inp
  use reciprocal_inp
  use eels_inp
  use nrixs_inp
  use hubbard_inp
  implicit none

  !     Local stuff
  integer :: i
  real*8 :: celvin
  call atoms_read       ! read geom.dat
  call global_read      ! read global.inp
  call fms_read         ! read fms.inp
  call eels_read        ! read eels.inp 
  call hubbard_read
  if (do_nrixs.eq.1) call nrixs_init       ! initialize nrixs variables !KJ used to be in ffmod3jas.f90
  call init_controls
  call reciprocal_read(celvin)  ! read reciprocal.inp
!  call init_struct(nph)

 !write(*,*) 'reafms nph=',nph
  nphkevin=nph !KJ copy global value to struct module
  ! !KJ Next section added for k-space calculations
  if(ispace.eq.0) then
     !KJ next lines : initialize nsp in the struct module (old routines will use same value reinitialized  by fmstot).
     nsp = 1
     if (abs(ispin).eq.1 ) nsp = nspx
     lpot(0:nph)=lmaxph(0:nph)
  endif

  makekmeshnow=.true.
  if(ispace.eq.0)  then !KJ go on and read the k-mesh!
     a1=a1/bohr  ! lattice constants in bohr
     a2=a2/bohr
     a3=a3/bohr
     celvin=celvin/(bohr**3)
     call crystalstructure(celvin)

     if(makekmeshnow) then
        call kmesh
     else
        !c          klist.inp  !KJ added 8/06
        open(3,file='klist.inp',form='formatted',status='old')
        read(3,*)
        read(3,*) nkp,usesym,ktype
        call init_kklist(nkp,nsym) !KJ 6-09
        read(3,*)
        do i=1,nkp
           read(3,*) bk(:,i),weight(i)
        enddo
        close(3)
        !           k-mesh in fractional units 0-1 for strfac
        call wlog('assuming a klist.inp in Angstrom')
        call wlog('update these instructions !!')
        do i=1,3
           bk(i,:)=bk(i,:) * alat(i)/(dble(2)*pi) *bohr  !remove *bohr if klist.inp is in a.u.
        enddo
     endif

  endif

  !     transform to code units (bohrs and hartrees - atomic units)
  rfms2 = rfms2 / bohr
  rdirec = rdirec / bohr
  rat=rat / bohr
  return
end subroutine reafms
