subroutine importcif(cifname,cif_equivalence)

! This routine opens a .cif file (Crystallographic Information File).
! It then extracts the information needed to set up a k-space FEFF9 calculation.
! ciftbx v.4 is used for the manipulation of the .cif file.
! The second part of the routine reworks the info from the .cif file to the right format.
use dimsmod,only: nphx=>nphu
use struct,only:  m_a1=>a1,m_a2=>a2,m_a3=>a3,m_alat=>alat,m_alfalat=>alfalat, &
  m_sgroup=>sgroup,m_nats=>nats,m_nph=>nph,m_natom=>natom, &
  m_latticename=>latticename,m_ppos=>ppos,m_ppot=>ppot,m_absorber=>absorber, &
  m_lpot=>lpot,m_nsp=>nsp,m_celvol=>celvol,m_cryst_gr=>cryst_gr, &
  m_nsym=>nsym,m_sgroup_hm=>sgroup_hm,m_lattice=>lattice,m_label=>label, &
  m_izatom=>izatom,m_firstpos=>firstpos,m_init_struct=>init_struct, &
  m_R_spacegroup_is_hexagonal=>R_spacegroup_is_hexagonal 
use reciprocal_inp, only: ispace
     !contains the k-space variables specifying the structure
use par
use constants, only:pi
implicit none
include 'ciftbx.cmn' !to use the ciftbx library
character*120,intent(in) :: cifname  ! cif-file that will be opened to read the crystal structure
integer,intent(inout) :: cif_equivalence ! Different schemes for choosing unique potentials
! local variables:
logical f1,f2,f3,inside
character*80 line
real*8 a,b,c,sig,alpha,beta,gamma,vec(3),lattice(3,3),deg2rad,cellmin(3),cellmax(3),dummy
real*8,allocatable :: x(:),y(:),z(:),xn(:,:),dpos(:,:)
character*5,allocatable :: symbol(:),label(:),symbol_old(:),label_old(:)
character*5 zname
character*60 symcell,sym(200)
character*8 sgrhm
character*2 str2
character*10 sgrhm1 !dummy
integer syminttablesnr
integer nat,nsym,i,j,ki,kj,l,index,iorigin,i1,i2,i3
real*8,allocatable :: pos(:,:),cpos(:,:),dist(:)
integer,allocatable :: mult(:),indequiv(:),iz(:),pot(:),firstpos(:),iz_old(:),pot_old(:),firstpos_old(:)
integer nat_cif,nat_ineq,nat_all,nat_temp
integer,parameter :: nat_start=1000
real*8 latvect(3,27),d,dminimum,xtrans(3)
real*8 aa,bb,cc
integer iclosest
character*3 lattyp
character*52,parameter :: alphabet='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
character*512 slog
logical,parameter :: verbose=.false.
logical,parameter :: allow_negative_coordinates=.true.
logical R_spacegroup_is_hexagonal
logical test_sghrm
logical, external :: checksgrhm


! Assign the CIFtbx files
f1=init_(1,2,3,6)

! Request dictionary validation check
!if (.not. dict_('/Users/jorissen/science/feff90/src/RDINP/cif_core.dic','valid')) write(6,*) 'Requested Core dictionary not present'
!if (.not. dict_('cif_core.dic','valid')) write(6,*) 'Requested Core dictionary not present'


! Open the CIF
if (.not. ocif_(cifname)) stop 'CIF cannot be opened'
if (verbose) call wlog('cif file opened')

! Assign the data block to be accessed
120 continue
if (.not.data_(' ')) goto 200

! Read cell dimensions
f1=numd_('_cell_length_a',a,sig)
f2=numd_('_cell_length_b',b,sig)
f3=numd_('_cell_length_c',c,sig)
if(.not. (f1.and.f2.and.f3)) write(6,*) 'Cell lengths missing'
f1=numd_('_cell_angle_alpha',alpha,sig)
f2=numd_('_cell_angle_beta',beta,sig)
f3=numd_('_cell_angle_gamma',gamma,sig)
if(.not. (f1.and.f2.and.f3)) write(6,*) 'Cell angles missing'
deg2rad=acos(-1.0d0)/dble(180)
alpha=alpha*deg2rad
beta=beta*deg2rad
gamma=gamma*deg2rad



allocate(x(nat_start),y(nat_start),z(nat_start),symbol(nat_start),label(0:nat_start))
x=0. ; y=0. ; z=0. ; symbol='' ; label=''
nat_cif=0
! Extract atom site data
160 continue
nat_cif=nat_cif+1
f1=char_('_atom_site_type_symbol',symbol(nat_cif)) ! e.g. "Ge"
f2=char_('_atom_site_label',label(nat_cif))  ! e.g. "Ge1" or also "Ge"
if (f2 .and. (.not.f1)) symbol(nat_cif)=label(nat_cif)
if (.not.f1 .and. .not.f2) stop 'error type of atom not specified'
f2=numd_('_atom_site_fract_x',x(nat_cif),sig)
f2=numd_('_atom_site_fract_y',y(nat_cif),sig)
f2=numd_('_atom_site_fract_z',z(nat_cif),sig)
if (loop_) goto 160

!Correct fractions - cif's are often less accurate than theoretical symmetry routines
do i=1,nat_cif
index=12
do j=-2*index,2*index
   if (dabs(x(i)-dble(j)/dble(index)).lt.0.0002) x(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
   if (dabs(y(i)-dble(j)/dble(index)).lt.0.0002) y(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
   if (dabs(z(i)-dble(j)/dble(index)).lt.0.0002) z(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
enddo
index=8
do j=2*index,2*index
   if (dabs(x(i)-dble(j)/dble(index)).lt.0.0002) x(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
   if (dabs(y(i)-dble(j)/dble(index)).lt.0.0002) y(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
   if (dabs(z(i)-dble(j)/dble(index)).lt.0.0002) z(i)=dble(j)/dble(index)  ! so 0.6667 becomes 0.66666666666...
enddo
enddo   


! Space group
f1=char_('_symmetry_cell_setting',symcell)  ! e.g. "hexagonal" ; this variable is not currently used
f2=char_('_symmetry_space_group_name_H-M',sgrhm1)  ! e.g. "P 63/m m c"
f3=numd_('_symmetry_Int_Tables_number',dummy,sig) ! e.g. "194"
if (f3) syminttablesnr=nint(dummy)
if(.not.f2 .and. .not.f3) stop "Spacegroup not specified in CIF file."
!Now, some idiotic code to get rid of redundant whitespaces and whatnot:
j=0
sgrhm='        '
!write(*,*) 'sgrhm1',sgrhm1
do i=1,10
   if (sgrhm1(i:i).ne." ") then
      j=j+1
	  sgrhm(j:j)=sgrhm1(i:i)
   endif
   if (j.ge.8) exit
   if(i.lt.10) then
       if (sgrhm1(i+1:i+1).eq.'{'.or.sgrhm1(i+1:i+1).eq.':') exit  !added a condition to remove { ... } comments or things like " :1"
   endif
enddo
!write(*,*) 'sgrhm',sgrhm
if(f2 .and. .not.f3) call getsgnum(sgrhm,syminttablesnr)
if(f3 .and. .not.f2) call getsgrhm(syminttablesnr,sgrhm)
if(f2 .and. f3) then
   !We were given 2 equivalent names for the space group.
   !Let's use the number, which has no room for notational confusion, to enforce our clean naming of the H-M symbol:
   call getsgrhm(syminttablesnr,sgrhm1)
   test_sghrm = checksgrhm(syminttablesnr,sgrhm)
   if (.not. test_sghrm) then
      call wlog('Warning - H-M name of space group possibly sloppy or incorrect.')
      call wlog('Using provided International Tables number for space group, and proceeding.')
      write(*,*) 'International Tables space group number from cif file: ',syminttablesnr
      write(*,*) 'This may correspond to H-M symbol ',sgrhm1
      write(*,*) 'The CIF file provided ',sgrhm,'   (adjusted to',sgrhm1,')'
	  write(*,*) 'If this is not what you intended, please edit the cif file to avoid ambiguity.'
	  write(*,*) 'E.g., it is possible to #comment out the line given the International Tables number.'
      sgrhm = sgrhm1
   endif
endif
call getlattype(sgrhm,lattyp) ! lattyp is one of the  Bravais lattices, i.e., P, CXZ, CYZ, F, B, H, R
if(lattyp.eq."   ") then
   call wlog("Something unexpected happened in spacegroup input; setting lattice type to Primitive.")
   write(*,*) 'sgrhm',sgrhm
   write(*,*) 'lattyp',lattyp
   write(*,*) 'syminttablesnr',syminttablesnr
   call par_stop('Confusing space group in CIF file.')
   lattyp="P  "
endif
!write(*,*) 'lattyp',lattyp

! Symmetry operations
nsym=0
sym=""
170 continue
nsym=nsym+1
f1=char_('_symmetry_equiv_pos_as_xyz',sym(nsym))
if(.not.f1) sym(nsym)=""
if(.not.f1) nsym=nsym-1
if(loop_) goto 170

if(nsym.eq.0) then
   180 continue
   nsym=nsym+1
   f1=char_('_space_group_symop_operation_xyz',sym(nsym))
   if(.not.f1) nsym=nsym-1
   if(loop_) goto 180
endif

if(nsym.eq.0) then
   if(syminttablesnr .eq. 1) then
      !P1
   nsym=1
   sym(1)='x,y,z'
   else
      call wlog('Error - could not read symmetry operations from cif file.  Please include either :')
      call wlog('    _symmetry_equiv_pos_as_xyz')
      call wlog('or')
      call wlog('    _space_group_symop_operation_xyz')
      call par_stop('symmetry operations missing from CIF file')
   endif
endif


goto 120  !This is a loop over "data blocks", i.e. sections in the .cif file.  (Many .cif files will have only one section.)  We loop over all of them to make
! sure we find the info we need.  Once the whole file is read, we pass through the 200 continue statement and proceed with processing tasks.)
200 continue

call close_   ! Closes the .cif file; we are now done using the CIFTBX library.
if (verbose) call wlog ('Done reading .cif file.')


if(m_absorber .gt. nat_cif  .or. m_absorber .lt. 1) then
   call wlog('INPUT ERROR')
   call wlog('Invalid TARGET iabsorber ')
   write(*,*) 'target is ',m_absorber,' but nat_ineq is ',nat_cif
   stop
endif


! Now special treatment for R (rhombohedral) lattices.
! A spacegroup starting with "R" can be expressed in 2 ways:
! 1. hexagonal     a=b,c alpha=beta=90,gamma=120 x,y,z
! 2. rhombohedral  a=b=c alpha=beta=gamma        x,y,z
! The hexagonal expression is easier to work with in the future, so we choose that as default.
! If the rhombohedral expression is given, we convert it to hexagonal (while keeping the spacegroup the same).
if(lattyp.eq.'R') then
   R_spacegroup_is_hexagonal = ( (dabs(a-b).lt.0.0002) .and. (dabs((alpha/deg2rad)-90d0).lt.0.0002) &
         .and. (dabs((beta/deg2rad)-90.d0).lt.0.0002)  .and. (dabs((gamma/deg2rad)-120d0).lt.0.0002) )
   if(R_spacegroup_is_hexagonal) then
      call wlog('Found R spacegroup given in hexagonal coordinates.')
      lattyp='H  '
      call wlog('lattice type set to H.')
   else
      call wlog('Found R spacegroup given in rhombohedral coordinates.')
      call wlog('lattice type set to R.')
   endif
else
   R_spacegroup_is_hexagonal = .false.   !It doesn't matter - just to avoid errors about uninitialized variables
endif

! Now generate the FULL list of atom positions in the unit cell, from the basic list of positions in the CIF file.
! This is actually not so easy because the sym ops are given in a computer-unfriendly way ...
   nat_ineq=nat_cif
   nat_all=0
   allocate(pos(3,nsym*nat_ineq),pot(nsym*nat_ineq),firstpos(nat_ineq))
   allocate(xn(3,nsym),mult(nat_ineq),indequiv(nsym),iz(0:nat_ineq))
   do i=1,nat_ineq

      if(verbose) write(*,*) 'For atom type i= ',i

      vec(1)=x(i);vec(2)=y(i);vec(3)=z(i)
	  ! The cif-file may contain things like "O(1)" instead of "O" for atom type.  The following lines strip it down to something the findz routine will understand:
	  zname=symbol(i)
	  zname(3:5)='   '
	  f1=.false.
	  do j=1,52
	     if(zname(2:2).eq.alphabet(j:j)) f1=.true.
	  enddo
	  if(.not.f1) zname(2:2)=' '
	  label(i)=zname  !save stripped name
	  call findz(zname(1:3),iz(i))
      do j=1,nsym
         call apply_cifsymop(sym(j),vec,xn(1:3,j))         
      enddo
      if(verbose) then
         write(*,*) 'After apply_cifsymop:'
         do j=1,nsym
            write(*,'(i4,x,3(f12.5,x),a2,x,i3)') j,xn(:,j)
         enddo
         write(*,*)
      endif
      call gen_equiv(nsym,xn,lattyp,indequiv,mult(i))
      if(verbose) then
         write(*,*) 'After gen_equiv:'
         do j=1,mult(i)
            write(*,'(i4,x,3(f12.5,x),a2,x,i3)') j,xn(:,indequiv(j))
         enddo
         write(*,*)
      endif

	  firstpos(i)=nat_all+1
      do j=1,mult(i)
         nat_all=nat_all+1
         pos(1:3,nat_all)=xn(1:3,indequiv(j))       
		 pot(nat_all)=i
      enddo
   enddo
   if(verbose) then
!     Let's see what we have now ...
!    "index" is now the number of atoms in the unit cell
     do i=1,nat_all
        write(*,'(i4,x,3(f12.5,x),a2,x,i3)') i,pos(:,i),symbol(pot(i)),pot(i)
     enddo
     write(*,*)
   endif


if ((.not.R_spacegroup_is_hexagonal) .and. (lattyp.eq.'R')) then
    !call wlog('R spacegroup; converting from rhombohedral to hexagonal basis')
    ! Don't do the transformation until now -- after the symmetry operations have been applied :).

    ! This is probably not necessary for real-space calculations : the atoms list will be fine anyway.
    ! But I *think* the kmesh routine expects the hexagonal lattice and may produce incorrect
    ! k-list otherwise.  I'm just not sure.
    ! Also, the hexagonal cell is a little more regular and might "converge" better in some ways.
    call wlog('Proceeding in rhombohedral basis.  Fine for real space calculations.')
    call wlog('For k-space calculations the authors have not yet checked that the k-mesh will be correct.')
    if(ispace.eq.0) then
    call wlog('Send us a message if you need help with this.')
    call wlog('Quitting.')
    stop
    aa=a*2.0d0*dcos((pi-alpha)/2.0d0)
    bb=aa
    cc=3.0d0*dsqrt(a**2-(aa**2)/3.0d0)
    a=aa
    b=bb
    c=cc
    alpha=90.0d0 * deg2rad
    beta=90.0d0 * deg2rad
    gamma=120.0d0 * deg2rad
    !In WIEN2k the x,y,z would remain in rhomb basis.
    !Here, we undoubtedly need to convert them also.
    ! (x' y' z')rhomb  =  (x y z)hex  (-1  1   0)
    !                                 ( 1  0  -1)
    !                                 ( 1  1   1)
    ! The inverse of which is  1/3    (-1  1  1)
    !                                 ( 2  1  1)
    !                                 (-1 -2  1)
    ! xx = 1/3* ( -x + 2*y -  z)
    ! yy = 1/3* (  x +   y -2*z)
    ! zz = 1/3* (  x +   y +  z)
    do i=1,nat_all
    aa=pos(1,i)
    bb=pos(2,i)
    cc=pos(3,i)
    pos(1,i) = (-aa + 2*bb -   cc) / dble(3)
    pos(2,i) = ( aa +   bb - 2*cc) / dble(3)
    pos(3,i) = ( aa +   bb -   cc) / dble(3)
    enddo
    endif
endif


! Now it is time for a crucial step.  We must choose the unique potentials.
if (cif_equivalence .eq. 4) then
   !The practical approach.  Use "1" for small cells and "2" if there are a lot of atoms.
   cif_equivalence=1
   if (nat_ineq .gt. nphx) cif_equivalence=2
endif
if (cif_equivalence .eq. 1) then
   !Stick with what's crystallographically correct
   if(verbose) call wlog('using exact crystallographic equivalence to choose unique potentials.')
elseif (cif_equivalence .eq. 2) then
   !FEFF-style : only care about atomic number Z, not about crystallographic environment
   if(verbose) call wlog('using atomic number to choose unique potentials.')
   !Yuck, we must renumber the whole damn list.
   if (nat_ineq .gt. 1) then
   nat_temp=1  ! The new number of inequivalent atoms
   mult(:)=0
   allocate(iz_old(0:nat_ineq),symbol_old(nat_start),label_old(0:nat_start),pot_old(nsym*nat_ineq),firstpos_old(nat_ineq))
   iz_old=iz 
   symbol_old=symbol
   label_old=label
   pot_old=pot 
   firstpos_old=firstpos
   pot(1)=1
   mult(1)=1
   firstpos(1)=1
   do i=2,nat_all
      f1=.true.
      do j=1,i-1  !See if there's a previous atom with the same Z
	     if (iz(pot(j)) .eq. iz_old(pot(i)) .and. f1) then
		    pot(i)=pot(j)
			f1=.false.
			mult(pot(j))=mult(pot(j))+1
		 endif
	  enddo
	  if (f1) then  !First atom of this Z
	     nat_temp=nat_temp+1
		 iz(nat_temp)=iz_old(pot(i))
		 label(nat_temp)=label_old(pot(i))
		 symbol(nat_temp)=symbol_old(pot(i))
		 pot(i)=nat_temp
		 firstpos(nat_temp)=i
		 mult(nat_temp)=1
	  endif
   enddo

   nat_ineq=nat_temp
   endif
   
elseif (cif_equivalence .eq. 3) then
   !Hybrid method : use crystallographic information up to nearest neighbors only
   if(verbose) call wlog('using approximate crystallographic equivalence to choose unique potentials - first shell.')
   stop 'not yet implemented'
else 
   call wlog('The EQUIVALENCE parameter for choosing unique potentials must be between 1 and 4.  Quitting now.')
   stop 'not yet implemented'
endif


if (nat_ineq .gt. nphx) then
   write(slog,*) 'CIF file contains',nat_cif,' atom types, which is larger than the hardwired limit nphx=',nphx,'.'
   call wlog(slog)
   call wlog('You need to recompile FEFF or simplify the structure.  Exiting now.')
   stop
endif



! Now we have to be careful : at this point, we still have special lattice types (P,B,F,CXZ,CZY,CYZ,H).
! Meaning the basis vectors are given for the conventional cell, while the positions are given for the primitive cell.
! The way FEFF currently works, that's going to give problems.  Either have it all primitive, or all conventional.
! Primitive requires changing the basis (lattice) vectors.  Conventional requires creating extra atom positions.
! Typically, the primitive cell has fewer atoms (by definition) but the conventional cell has higher symmetry.
! This higher symmetry results in a smaller k-mesh; however, symmetry is currently not used to reduce the computational
! cost of calculating G(k).
! Overall I think that the primitive cell is the fastest way to go in the current version of FEFF.  (Verified in a few cases.)
! However we pass everything on as is at this point.
! Reciprocal.inp will contain stuff as it is here.
! Then in each module, during setup, the routine crystalstructure is called.
! This converts us to the primitive lattice.



! Set up the absorber:
! Be careful here - with CIF, "absorber" counts inequivalent positions; without CIF, "absorber" counts the positions in the ATOMS card.  This is altogether different ...
if (cif_equivalence.eq.2 .and. nat_cif.gt.1) then
   !Use the old information from the .cif
   iz(0)=iz_old(m_absorber)
   label(0)=label_old(m_absorber)		
else 
   iz(0)=iz(m_absorber)
   label(0)=label(m_absorber)
endif   


! So far, the lattice vectors have only been specified by length and angle.  Let's now make the vectors explicitly, in carthesian coordinates.
! choose c || e_z
lattice(:,:)=dble(0)
lattice(3,3)=c
! we have one more choice; let's make b _|_ e_x
lattice(2,3)=dcos(alpha)*b
lattice(2,2)=dsin(alpha)*b
! and then it follows for a that ...
lattice(1,3)=dcos(beta)*a
lattice(1,2)=((dcos(gamma)-dcos(beta)*dcos(alpha))/dsin(alpha))*a
lattice(1,1)=dsqrt(a**2-lattice(1,2)**2-lattice(1,3)**2)
if (verbose) then
   write(*,*) 'lattice vectors:'
   write(*,*) 'a1: ',lattice(1,:)
   write(*,*) 'a2: ',lattice(2,:)
   write(*,*) 'a3: ',lattice(3,:)
endif

! Now print out unit cell position in Carthesian coordinates
allocate(cpos(3,nat_all))
allocate(dist(nat_all)) ! for testing
cpos=dble(0)
do i=1,nat_all
do j=1,3
   cpos(j,i)=pos(1,i)*lattice(1,j)+pos(2,i)*lattice(2,j)+pos(3,i)*lattice(3,j)
enddo
enddo


! Shift origin
if (cif_equivalence.eq.2 .and. nat_cif.gt.1) then
   iorigin=firstpos_old(m_absorber)
else
   iorigin=firstpos(m_absorber)
endif
vec(:)=cpos(:,iorigin)
do i=1,nat_all
do j=1,3
   cpos(j,i)=cpos(j,i)-vec(j) !cpos(j,iorigin) - somehow that makes the compiler produce rubbish ...
enddo
enddo

if(.not.allow_negative_coordinates) then
! After this transformation pos -> cpos, some positions may not be in the unit cell anymore.  Bring them back in.
! The cell boundaries in carthesian coordinates :  ! [0,1]^3 in lattice coordinates
cellmin=dble(0)
cellmax=dble(0)
do i1=0,1 ; do i2=0,1 ; do i3=0,1
   do j=1,3
      vec(j)=i1*lattice(1,j)+i2*lattice(2,j)+i3*lattice(3,j)
	  if (vec(j).lt.cellmin(j)) cellmin(j)=vec(j)
	  if (vec(j).gt.cellmax(j)) cellmax(j)=vec(j)
   enddo
enddo ; enddo ; enddo

do i=1,nat_all
   inside=(cpos(1,i).ge.cellmin(1) .and. cpos(1,i).le.cellmax(1) &
      .and. cpos(2,i).ge.cellmin(2) .and. cpos(2,i).le.cellmax(2)  &
	  .and. cpos(3,i).ge.cellmin(3) .and. cpos(3,i).le.cellmax(3) )
   if (.not.inside) then
      do i1=-3,3 ; do i2=-3,3 ; do i3=-3,3
	     if(.not.inside) then
	     do j=1,3
	        vec(j)=cpos(j,i)+i1*lattice(1,j)+i2*lattice(2,j)+i3*lattice(3,j)
	     enddo
         inside=(vec(1).ge.cellmin(1) .and. vec(1).le.cellmax(1) &
		 .and. vec(2).ge.cellmin(2) .and. vec(2).le.cellmax(2)  &
		 .and. vec(3).ge.cellmin(3) .and. vec(3).le.cellmax(3) )
		 if (inside) cpos(:,i)=vec(:)
	     endif
	 enddo; enddo; enddo
   endif
   if (.not.inside) write(*,*) 'Insideing failed for atom ',i
enddo
endif ! enforce all coordinates inside [0,1].

! alternative strategy:
! assume that the input is a list of atoms in a unit cell.
! Because of periodicity, all coordinates are given up to factors of 1 in lattice vector units.
! This means that 0.9 and 0.1 are only a distance of "0.2" away, rather than a distance of "0.8".
! However, the code so far does not know this.
! For small unit cells, it is not a problem.  We will make many copies of the unit cell to fill a list of several hundred real-space positions:
! This will fill up the volume around the absorber.  However, for unit cells that are essentially a molecule in a big vacuum cell, this
! does not happen and the molecule is effectively cut in 2 or more pieces.  Therefore the list in atoms.dat is rubbish with distances that are
! completely wrong.  This can lead to algorithmic errors (e.g. NOVP) or just rubbish spectra.
! To solve, look at all the images of each atom created by simple lattice translations and keep only the one closest to the absorber, which is 
! by now at the origin.
! Create list of lattice translations:
i=0
do i1=-1,1 ; do i2=-1,1 ; do i3=-1,1
   i=i+1
   do j=1,3
      latvect(j,i)=i1*lattice(1,j)+i2*lattice(2,j)+i3*lattice(3,j)
   enddo
enddo ; enddo ; enddo

! Iterate over all atoms in the cluster (except for the absorber)
do i=1,nat_all
   if(i.ne.iorigin) then
   dminimum=100.
   do j=1,27
      xtrans=cpos(:,i)+latvect(:,j)
      d=dsqrt((xtrans(1)-cpos(1,iorigin))**2+(xtrans(2)-cpos(2,iorigin))**2+(xtrans(3)-cpos(3,iorigin))**2)
      if(d.lt.dminimum) then
         iclosest=j
         dminimum=d
      endif
   enddo
   xtrans=cpos(:,i)+latvect(:,iclosest)
   cpos(:,i)=xtrans(:)
   endif
enddo



! As a check, find distance from "central atom"
do i=1,nat_all
   dist(i)=dsqrt((cpos(1,i)-cpos(1,iorigin))**2+(cpos(2,i)-cpos(2,iorigin))**2+(cpos(3,i)-cpos(3,iorigin))**2)
enddo

if (verbose) then
   do i=1,nat_all
      write(*,'(i4,x,3(f12.5,x),a2,3x,f12.5)') i,cpos(:,i),symbol(pot(i)),dist(i)
   enddo
   write(*,*)
   write(*,*) 'lattyp ',lattyp
endif

! "Translate" the absorber to FEFF language
! i.e., instead of giving a potential index, give a position index
if (cif_equivalence.eq.2 .and. nat_cif.gt.1) then
   m_absorber=firstpos_old(m_absorber)
else
   m_absorber=firstpos(m_absorber)
endif

! Copy variables into the "struct" module so the rest of the program can use them:
m_nsym=nsym
m_nats=nat_all
m_nph=nat_ineq
m_a1=lattice(1,:)
m_a2=lattice(2,:)
m_a3=lattice(3,:)
!KJ Something really weird happens here ; using ifort12.5 on Mac 10.7.2, "cpos" gets corrupted, apparently during the call to m_init_struct.
!   I don't understand what I'm during wrong.
!   Creating "dpos" seems to avoid whatever the compiler otherwise trips up on.
!   If anyone understands why, please explain it to me ...
allocate(dpos(3,m_nats))
do i=1,m_nats
   dpos(:,i)=cpos(:,i)
enddo

m_latticename=lattyp
m_lattice=lattyp(1:1)
call m_init_struct(m_nats)  !Allocates the necessary arrays
m_natom(1:m_nph)=mult(1:nat_ineq)
m_ppot(1:m_nats)=pot(1:nat_all)
do i=0,nat_ineq ! copy 5-string to 2-string carefully:
   zname=label(i)
   str2=zname(1:2)
   m_label(i)=str2
enddo
!!!m_label(0:m_nph)=label(0:nat_ineq)
m_izatom(0:m_nph)=iz(0:nat_ineq)
!m_ppos(1:3,1:m_nats)=pos(1:3,1:nat_all) !for WIEN2k-like coordinates, i.e. "ICOORD=4" ; neither of below
m_ppos(1:3,1:m_nats)=dpos(1:3,1:nat_all) !for FEFF-like coordinates, i.e. "ICOORD=1" ; origin at absorber and within first unit cell 
m_ppos=m_ppos/a !Convert the FEFF-like coordinates to SPRKKR-like coordinates by dividing by length of first lattice vector
m_sgroup=syminttablesnr
m_sgroup_hm=sgrhm
m_firstpos(1:m_nph)=firstpos(1:nat_ineq)
m_R_spacegroup_is_hexagonal=R_spacegroup_is_hexagonal

if (cif_equivalence.eq.2 .and. nat_cif.gt.1) deallocate(iz_old,symbol_old,label_old,pot_old,firstpos_old)

return
end
