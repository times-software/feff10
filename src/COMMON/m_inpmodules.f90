!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_inpmodules.f90,v $:
! $Revision: 1.64 $
! $Author: bmattern $
! $Date: 2013/01/11 19:17:07 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This file written by Kevin Jorissen 7-09

!     THIS FILE REPLACES ALLINP.H AND MOST OF WRTALL.F90, INIALL.F90.


!     LAYOUT - PLEASE READ FIRST !!
!     Each module contains :
!      - a list of variables and parameters
!      - a module_write subroutine to write variables to a module.inp file
!      - a corresponding module_read subroutine
!      - a module_init subroutine to specify default values

!     Note to programmers : use of implicit none for each module MANDATORY
!     since a small mistake here could otherwise have very messy repercussions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     We're getting so many input parameters, some of them optional, that passing them through "rdmodxinp" subroutines is getting messy.
!     It also makes updating options difficult and messy.

!     FEFF is set up in the following way : each subprogram starts by reading a file containing all the settings that determines how it should run.
!     Hence, manipulation of this single file allows the user to tweak any subprogram.
!     However, obviously some variables control more than one subprogram.
!     Still they can belong to only one module!
!     To keep things manageable, such variables will either be put in the global_inp module ;
!     or in the module of the first subprogram (first w.r.t. normal program flow) that needs it.
!     All other modules will then have to call the first module with a "use module_x, only : var_x" statement.
!     The choice is sometimes a bit arbitrary.

!     Note that the only alternative would be to make the input files such that each variable occurs only once.
!     This would make the modules much cleaner, ie no use statements necessary.
!     However, it would then be less easy for the user to see which variables affect a particular subprogram he wants to run.
!     Also, changing the input file would then affect all subprograms needing that input file.
!     (This is already the case for quite a few input files, such as eels.inp, reciprocal.inp, global.inp, ...)


!     Ideally, in the future, defaults will be set as private parameters in each module's INIT subroutine.
!     Then, input files might not have to contain some optionals, and parameters could be defined free-format using XML.

!     Even more ideally, in a very distant future, each line will contain only one variable declaration
!     paired with a comment line describing its function ...



!     I REALLY THINK IT'S BEST TO KEEP ALL THESE MODULES IN ONE FILE ; or at least in a separate folder so they don't get mixed up
!     with everything else ...






!=======================================================================
!     GEOMETRY
!=======================================================================

      module geometry_inp
	    use dimsmod, only: nattx
		implicit none
	!c    atoms.dat
		integer  natt
		integer iphatx(nattx)
		double precision  ratx(3,nattx)

		contains

		subroutine geometry_write_atoms
			integer iat
            double precision distance
			!c    atoms.dat to be read by ffsort, which will write smaller geom.dat file
			open (file='atoms.dat', unit=3, status='unknown')
			write (3, 35) natt
		  35    format ('natx =  ', i7)
			write (3, 10) '    x       y        z       iph  '
			do iat = 1, natt
                   distance=dsqrt((ratx(1,iat)-ratx(1,1))**2+(ratx(2,iat)-ratx(2,1))**2+(ratx(3,iat)-ratx(3,1))**2) ! core hole should be in position 1 by now
				write(3,36) ratx(1,iat), ratx(2,iat), ratx(3,iat), iphatx(iat), distance
		  36      format( 3f13.5, i4, f13.5)
			  enddo
			close(3)
		  10  format(a)
		  20  format (20i4)
		  30  format (9f13.5)
		end subroutine geometry_write_atoms

		subroutine geometry_init
			natt = 0
			iphatx(:) = -1
			ratx(:,:) = 0.d0 !KJ added 7-09
		end subroutine geometry_init

      end module


!=======================================================================
!     ATOMS
!=======================================================================

	  module atoms_inp
	  ! The geom.dat file
	    use dimsmod,only: natx,nphx=>nphu,nheadx
        implicit none
	    integer nat, nph, iphat(natx), ibounc(natx)
        integer, allocatable :: iatph(:)  !, iatph(0:nphx)
		! ibounc is currently set to 1 for all atoms in ffsort.  Path uses it.  Probably discontinued variable but ah well. !KJ
	    double precision  rat(3,natx)
		character(*),parameter,private :: filename='geom.dat'
!		iphat(natx)  -  given specific atom, which unique pot?
!		rat(3,natx)  -  cartesian coords of specific atom
!		iatph(0:nphx)  - given unique pot, which atom is model?
!                      (0 if none specified for this unique pot)

		contains

        subroutine atoms_allocate
            if(.not.allocated(iatph)) allocate(iatph(0:nphx))
        end subroutine atoms_allocate

		subroutine atoms_read
	!		Read  geom.dat file
            implicit none
		    character*512 slog
			character*80 head(nheadx)
			integer lhead(nheadx),j,j1,nhead
            real*8 rdum1(3)
            integer idum1,idum2
            call atoms_allocate
			open (file=filename, unit=3, status='old')
!			read header
			nhead = nheadx
			call rdhead (3, nhead, head, lhead)
			nat = 0
		    nph = 0
			iatph(:)=0
  50		continue
!KJ I switched up statements below so that code doesn't falsely abort when nat=natx.
!KJ			nat = nat+1
			if (nat .gt. natx)  then
              write(slog,'(a, 2i10)') ' nat, natx ', nat, natx
              call wlog(slog)
              stop 'Bad input'
			endif
                        read(3,*,end=60) j1,rdum1(1:3),idum1,idum2
			nat = nat+1
                        rat(1:3,nat)=rdum1(1:3)
                        iphat(nat)=idum1
                        ibounc(nat)=idum2
!KJ			read(3,*,end=60)  j1, (rat(j,nat),j=1,3), iphat(nat), ibounc(nat) !KJ j2  !KJ put ibounc back in for program PATH
			if (iphat(nat).gt.nph) nph = iphat(nat)
			if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
			goto 50
  60		continue
!KJ			nat = nat-1
			close(3)
		end subroutine atoms_read

	  end module




!=======================================================================
!     GLOBAL
!=======================================================================

      module global_inp
        implicit none
	!	the variables evnorm, xivnorm, spvnorm and l2lp are exclusive to nrixs (feffq) calculations
	!	le2 has different meaning for feffq calculations
	!	xivec serves many different functions depending on spectroscopy : xas/eels/nrixs
		integer do_nrixs,lj,ldecmx !KJ 7-09 for feff8q
	!	configuration average data :
		integer nabs, iphabs
		real*8 rclabs
	!	global polarization data :
		integer ipol, ispin, le2, l2lp
		real*8 elpty, angks
		real*8 evec(3), xivec(3), spvec(3)
		complex*16 ptz(-1:1,-1:1)
		double precision evnorm, xivnorm, spvnorm
	!   moved here because I think it belongs here !KJ 7-09
		integer ispec
		character(*),parameter,private :: filename='global.inp'  !KJ used to be global.dat !!!
!       How many q-vectors:  (impulse transfer)
		integer nq
!       Are we doing direction averaged impulse transfer?  (Note: this means q || e_z, which is not really averaging!)
		logical qaverage
!       The list of q-vectors and their norm:
		real*8,allocatable :: qs(:,:),qn(:)
!       Weights of q-vectors in the cross-section (probably calculated by another code):
        complex*16,allocatable :: qw(:)
!       Are we doing q,q' crossterms in the NRIXS code?
        logical mixdff
!       If entering q, q' as length(q), length(q'), angle(q,q'), this is cosine(angle(q,q'))
        real*8,allocatable :: cosmdff(:,:)
!       and this is norm(q')
        real*8 qqmdff
!       A rotation matrix for each q-vector  (containing cos(theta),sin(theta),cos(fi),sin(fi) for each q)
        real*8,allocatable :: qtrig(:,:)  !compare to Adam's code qtrig(iq,1)=qcst(iq); 2)=qsnt; 3)=qcsf; 4)=qsnf
!       Should the mdff program run? !11-2010
        integer imdff

		contains

		subroutine init_feffq
	!	called to calculate some variables for nrixs
			integer i
			evnorm=0.0d0
			xivnorm=0.0d0
			spvnorm=0.0d0
			do i=1,3
			   evnorm=evnorm+evec(i)*evec(i)
			   xivnorm=xivnorm+xivec(i)*xivec(i)
			   spvnorm=spvnorm+spvec(i)*spvec(i)
			end do
			spvnorm=sqrt(spvnorm)
			xivnorm=sqrt(xivnorm)
			evnorm=sqrt(evnorm)
		end subroutine init_feffq

		subroutine global_write(iniq)
			integer i
			logical,intent(in) :: iniq
			if(iniq) call init_feffq
			open (file=filename, unit=3, status='unknown')
			write (3, 10) ' nabs, iphabs - CFAVERAGE data'
			write (3, 45) nabs, iphabs, rclabs
		  45	  format ( 2i8, f13.5)
			write (3,10) ' ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj' !KJ last 4 added for feff8q.  Note le2 new meaning in feff8q
			write (3, 50)  ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj !KJ
		  50	  format ( 3i5, 2f12.4, 10i5)  !KJ
			write (3, 10) 'evec		  xivec 	   spvec'
			do 60 i = 1,3
			write (3,30) evec(i), xivec(i), spvec(i)
		  60	  continue
			write (3, 10) ' polarization tensor '
			do 70 i = -1, 1
				write(3,30) dble(ptz(-1,i)), dimag(ptz(-1,i)), dble(ptz(0,i)), dimag(ptz(0,i)),  dble(ptz(1,i)), dimag(ptz(1,i))
		  70	  continue
		!KJ for feff8q - was in different place in file in feff8q:
			write(3,10) 'evnorm, xivnorm, spvnorm - only used for nrixs'
			write (3,30) evnorm, xivnorm, spvnorm !KJ
		!KJ for a list of q-vectors (NRIXS) and MDFF calculation (NRIXS) - only relevant for NRIXS calculations  12-2010
			write(3,10) "nq,    imdff,   qaverage,   mixdff"
			   write(3,*) nq,imdff,qaverage,mixdff
			write(3,*) 'q-vectors : qx, qy, qz, q(norm), weight, qcosth, qsinth, qcosfi, qsinfi'
			if(nq.gt.0) then    !note that this is redundant with xivec if nq=1, but ah well.
			   do i=1,nq
			      write(3,30) qs(i,:),qn(i),qw(i),qtrig(i,1:4)
			   enddo
	        endif
			if(mixdff) then
			   write(3,*) "   qqmdff,   cos<q,q'>"
			   write(3,*) qqmdff,cosmdff
		    endif
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (20f13.5)

		end subroutine global_write

		subroutine global_read
			real*8 aa1,bb1,aa2,bb2,aa3,bb3
			integer i
			open (file=filename, unit=3, status='old')
			read  (3,*)
			read  (3,*) nabs, iphabs, rclabs
			read  (3,*)
			read  (3,*)  ipol, ispin, le2, elpty, angks, l2lp, do_nrixs, ldecmx, lj
			read  (3,*)
			do i = 1,3
			  read  (3,*) evec(i), xivec(i), spvec(i)
			enddo
			read  (3,*)
			do i = -1, 1
			  read (3,*) aa1, bb1, aa2, bb2, aa3, bb3 !KJ changed names of dummies to avoid confusion with my WELL DEFINED arrays (f*cking "implicit" people !#$&%)
			  ptz(-1,i)= dcmplx(aa1,bb1)  !KJ changed cmplx to dcmplx to satisfy thorough compilers
			  ptz(0,i) = dcmplx(aa2,bb2)
			  ptz(1,i) = dcmplx(aa3,bb3)
			enddo
			read (3,*)
			read (3,*) evnorm, xivnorm, spvnorm !KJ
			if(do_nrixs .ne. 0) then !compatibility with (most) old files
		       read(3,*)
			   read(3,*) nq,imdff,qaverage,mixdff
			   read(3,*)
			   call make_qlist(nq)
			   if(nq.gt.0) then    !note that this is redundant with xivec if nq=1, but ah well.
			      do i=1,nq
			         read(3,30) qs(i,:),qn(i),qw(i),qtrig(i,1:4)
			      enddo
	           endif
			   if(mixdff) then
			      read(3,*)
				  read(3,*) qqmdff,cosmdff
			   endif
	  30  format (20f13.5)
			endif
			close(3)
		end subroutine global_read

		subroutine make_qlist(n)
		    implicit none
			integer,intent(in) :: n
			allocate(qs(n,3),qn(n),qw(n),qtrig(n,4),cosmdff(n,n))
			qs(:,:)=0.d0
			qn(:)=0.d0
			qw(:)=1.d0
			qtrig(:,:)=0.d0
			qtrig(:,1)=1.d0
			qtrig(:,3)=1.d0 !corresponding to not rotating at all
			cosmdff(:,:)=1.d0
		    return
		end subroutine make_qlist

        subroutine make_qtrig
!          simple routine to get rotation angles; copied from mkptz
           implicit none
           double precision rr,rsp
           integer iq
	       do iq=1,nq
              if (qn(iq).gt.0.0d0) then
                 rsp = qn(iq)
                 rr = qs(iq,1)**2 + qs(iq,2)**2
				 if (rr.eq. 0) then
					qtrig(iq,1) = - 1.d0
					qtrig(iq,2) = 0.d0
					qtrig(iq,3) = 1.d0
					qtrig(iq,4) = 0.d0
				elseif (qs(iq,3).lt.0) then !meaning forward scattering ??
!                  rotation is defined by angles theta and fi
				   rr = sqrt(rr)
				   qtrig(iq,1) = qs(iq,3) / rsp
				   qtrig(iq,2) = rr / rsp
				   qtrig(iq,3) = qs(iq,1) / rr
			       qtrig(iq,4) = qs(iq,2) / rr
				else
                   qtrig(iq,1)=1.0d0
                   qtrig(iq,2)=0.0d0
                   qtrig(iq,3)=1.0d0 !surely this is a bug??  Shouldn't this be 1? !KJ 12-2011 changed 0->1 because produces NaN in genfmt otherwise
                   qtrig(iq,4)=0.0d0
				end if
              else
                 call wlog(' FATAL error: one of the q-vectors is zero')
                 call par_stop(' ')
              endif
		   enddo !iq
           return
        end subroutine make_qtrig

		subroutine global_init
			ispec = 0
			ldecmx=-1 ! initialize the number of decomposition channels - KJ 7-09 for feff8q
			nabs = 1
			iphabs = 0
			rclabs = 0.d0
			ipol = 0
			ispin = 0
			le2 = 0
			l2lp = 0
			elpty = 0.d0
			angks = 0.d0
			evec(:) = 0.d0
			xivec(:) = 0.d0
			spvec(:) = 0.d0
			ptz(:,:) = cmplx(0.d0,0.d0)
			evnorm=0.0d0
			xivnorm=0.0d0
			spvnorm=0.0d0
			do_nrixs=0 ! no nrixs calculation
			lj = -1
			nq=0
			qaverage=.true.
			imdff=0 !no mdff
			mixdff=.false. !no mdff in NRIXS
		!	cosmdff=1.d0  ! q || q'  => cos(0)=1     !KJ 11-2011 this is now an allocatable array.  Instruction fails on gfortran.
			qqmdff=-1.d0 ! leads to q=q' (norm only)
		end subroutine global_init


	end module



!=======================================================================
!     RECIPROCAL
!=======================================================================

      module reciprocal_inp
	!     k-space variables :
		use controls  !KJ 8/06
		use struct, nphstr => nph
		use kklist,only: nkp,usesym,nkx,nky,nkz,ktype
        use strfacs,only: streta,strrmax,strgmax,init_strfacs
		implicit none
		integer icorehole
		real*8 streimag ! additional broadening for calculation KKR structure factors ; not recommended
		character(*),parameter,private :: filename='reciprocal.inp'

		contains

		subroutine reciprocal_write
		!KJ next file added 8/06
		    integer i
			open (file=filename, unit=3, status='unknown')
		!       in which space are we?
			write(3,10) 'ispace'
			write(3,20) ispace
			if(ispace.eq.0) then
			   write(3,10) 'lattice vectors  (in A, in Carthesian coordinates)'
			   write(3,30) a1
			   write(3,30) a2
			   write(3,30) a3
			   write(3,10) 'Volume scaling factor (A^3); eimag; core hole'
			   write(3,30) dble(-1),dble(0),dble(1)
			   write(3,10) 'lattice type  (P,I,F,R,B,CXY,CYZ,CXZ)'
			   write(3,11) latticename,sgroup_hm,sgroup
			   write(3,10) '#atoms in unit cell ; position absorber ; corehole?'
			   write(3,20) nats,absorber,icorehole
			   write(3,10) '# k-points total/x/y/z ; ktype; use symmetry?'
			   write(3,*) nkp,nkx,nky,nkz,ktype,usesym  ! format line 20 limits integer to 4 positions - not enough for nkp!
			   write(3,10) 'ppos'
			   do i=1,nats
				  write(3,30) ppos(:,i)
			   enddo
			   write(3,10) 'ppot'
			    !KJ bugfix 5/2012: It's important not to use formatting when there are more atoms than fit on one line!!
			   write(3,*) ppot
			   write(3,10) 'label'
			   write(3,'(100(a2,1x))') label
			   write(3,10) 'streta,strgmax,strrmax'
			   write(3,30) streta,strgmax,strrmax
				endif
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  11  format(a3,4x,a8,4x,i4)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine reciprocal_write

		subroutine reciprocal_read(celvin)
		use struct, nphstr => nph
		integer i
		real*8,intent(out) :: celvin
        open (3,file=filename,status='unknown',err=167)
        read(3,*,end=167,err=167)
        read(3,*,end=167,err=167) ispace
        if(ispace.eq.0) then
             read(3,*) ; read(3,*) a1(:)
        	 read(3,*) a2(:)
        	 read(3,*) a3(:)
        	 read(3,*) ; read(3,*) celvin,streimag,cholestrength
        	 read(3,*) ; read(3,11) latticename,sgroup_hm,sgroup
			 lattice=latticename(1:1)
             read(3,*) ; read(3,*) nats,absorber,icorehole
             read(3,*) ; read(3,*) nkp,nkx,nky,nkz,ktype,usesym
             read(3,*)
	  11     format(a3,4x,a8,4x,i4)
			 !Careful: the next statement used to be "if size(ppot).eq.0".  However, on ifort size(ppot)=0 but on gfortran it =1!!
			 !Hence the new instruction.
			 !I wish if(allocated(ppot)) would work here; I don't understand why it doesn't.
                 !if(size(ppot).lt.nats) call init_struct(nats) !KJ 7-09 bugfix call this only once ; I can't seem to use "allocated(ppos)" here?
                 ! JK - replace above here. Check allocation inside
                 ! init_struct.
             call init_struct(nats)
	     do i=1,nats
	         read(3,*) ppos(:,i)
	     enddo
             read(3,*) ; read(3,*) ppot
			 read(3,*) ; read(3,'(100(a2,1x))') label
             read(3,*) ; read(3,*) streta,strgmax,strrmax
			 if(icorehole.eq.1) then
				corehole=.true.
			 else
				corehole=.false.
			 endif
		endif
        return
167     ispace=1
        return
		end subroutine reciprocal_read

		subroutine reciprocal_init
			call init_controls
			call init_strfacs
			icorehole = 1  ! use core hole
			streimag = dble(0) ! no extra broadening for KKR struc factors
			cholestrength = dble(1) ! don't mess with core hole
		end subroutine reciprocal_init

	end module






!=======================================================================
!     POTENTIAL
!=======================================================================

      module potential_inp
		use dimsmod, only: nheadx, nphx=>nphu, novrx
		use global_inp, only: ispec
		use atoms_inp, only : nph
		implicit none
		character(*),parameter,private :: filename='pot.inp'

		character*80 title(nheadx)
		integer mpot, ntitle, ihole, ipr1, iafolp, iunf,             &
			nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, ixc

                integer iscfxc !-LC- Exchange correlation functional for the SCF cycle
                                  ! 11=vBH 12=PZ 21=PWD 22=KSDT

		integer, allocatable :: iz(:)
!		iz(0:nphx)    - atomic number, input
		integer, allocatable :: lmaxsc(:)
		real rfms1
		double precision gamach, rgrd, ca1, ecv, totvol
		double precision, allocatable :: xnatph(:), folp(:), spinph(:)
!		xnatph(0:nphx) - given unique pot, how many atoms are there
!                      of this type? (used for interstitial calc)
!		folp(0:nphx) -  overlap factor for rmt calculation
		double precision, allocatable ::  xion(:)
!		xion(0:nphx)  - ionicity, input
		logical ExternalPot, FiniteNucleus
	!     for OVERLAP option
	    logical StartFromFile
		! read potential from pot.bin file and start from there
		integer, allocatable :: novr(:), iphovr(:,:), nnovr(:,:)
		double precision, allocatable ::  rovr(:,:)
!		novr(0:nphx) -  number of overlap shells for unique pot
!		iphovr(novrx,0:nphx) -  unique pot for this overlap shell
!		nnovr(novrx,0:nphx) -   number of atoms in overlap shell
!		rovr(novrx,0:nphx)  -   r for overlap shell
		! Added by Fer
		! Used to correct the excitation energy for chemical shifts
		integer  ChSh_Type
		integer configtype !KJ 12-2010 : which method for choosing atomic configuration?
		double precision corval_emin  !KJ 12-2012 defines energy window for search for core-valence separation energy.


!       criteria for self-consistency
        real*8,parameter :: tolmu = 1.D-3  ! Fermi level (Ha)
        real*8,parameter :: tolq = 1.D-3   ! net charge on atom iph (e)
	    real*8,parameter :: tolqp = 2.D-4  ! partial charge (e.g. l=1) on atom iph (e)
	    real*8,parameter :: tolsum = 0.05  ! total valence charge in Norman sphere compared to formal valence charge

	    real*8 scf_temperature ! Electronic temperature in eV
	    integer scf_thermal_vxc

        real*8 xntol
        INTEGER iscfth ! TS 07/2020
	integer nmu ! cja 10/15/20
		contains

        subroutine potential_allocate
        !call wlog('in potential_allocate')
            if(.not.allocated(iz))  allocate(iz(0:nphx))
            if(.not.allocated(lmaxsc)) allocate(lmaxsc(0:nphx))
            if(.not.allocated(xnatph)) allocate(xnatph(0:nphx), folp(0:nphx), spinph(0:nphx), &
              xion(0:nphx),novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx), rovr(novrx,0:nphx))
              !call wlog('finished potential_allocate')
        end subroutine potential_allocate

		subroutine potential_write
			integer ititle,ip,iph,iovr
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mpot, nph, ntitle, ihole, ipr1, iafolp, ixc,ispec, iscfxc'
			  write(3,20) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec, iscfxc
			  write(3,10) 'nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf'
			  write(3,20)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf
			  do ititle = 1, ntitle
		         write(3,10) title(ititle)
			  enddo
			  write(3,10) 'gamach, rgrd, ca1, ecv, totvol, rfms1, corval_emin'
			  write(3,30)  gamach, rgrd, ca1, ecv, totvol, rfms1, corval_emin
			  write(3,10) ' iz, lmaxsc, xnatph, xion, folp'
		  120   format ( 2i5, 4f13.5)
			  do ip = 0, nph
		        write(3,120) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
			  enddo
			  write(3,10) 'ExternalPot switch, StartFromFile switch'
			  write(3,*) ExternalPot,StartFromFile
		!       for OVERLAP option
			  write(3,10) 'OVERLAP option: novr(iph)'
			  write(3,20) ( novr(iph), iph=0,nph)
			  write(3,10) ' iphovr  nnovr rovr '
		  140   format ( 2i5, f13.5)
			  do iph = 0, nph
			  do iovr = 1, novr(iph)
		         write(3,140) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
			  enddo
	          enddo
		! Added by Fer
		! Correction of the excitation energy for chemical shifts
			  write(3,10) 'ChSh_Type:'
			  write(3,20) ChSh_Type
			  write(3,10) 'ConfigType:'
			  write(3,20) configtype
        ! Electronic temperature for thermal scf routine
              write(3,10) 'Temperature (in eV):'
              write(3,*) scf_temperature, scf_thermal_vxc
              write(3,10) 'Scf_th,  xntol,  nmu'
              write(3,*) iscfth, xntol, nmu ! TS 07/2020
		close(3)
              write(3,10) 'FiniteNucleus'
              write(3,*) FiniteNucleus
              close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (21i4)
	  30  format (9f13.5)
		end subroutine potential_write

		subroutine potential_read
			integer ititle,ip,iph,iovr
            call potential_allocate
            call potential_init
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec, iscfxc
			  read(3,*) ; read(3,*)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf
			  do ititle = 1, ntitle
		         read(3,*) title(ititle)
			  enddo
			  read(3,*) ; read(3,*)  gamach, rgrd, ca1, ecv, totvol, rfms1, corval_emin
			  read(3,*)
			  do ip = 0, nph
		        read(3,*) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
			  enddo
			  read(3,*) ; read(3,*) ExternalPot, StartFromFile
			  read(3,*) ; read(3,*) (novr(iph), iph=0,nph)
			  read(3,*)
			  do iph = 0, nph
			  do iovr = 1, novr(iph)
		         read(3,*) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
			  enddo
	          enddo
			  read(3,*) ; read(3,*) ChSh_Type
			  read(3,*,end=55) ; read(3,*,end=55) configtype
			  read(3,*) ; read(3,*) scf_temperature, scf_thermal_vxc
              read(3,*) ; read(3,*) iscfth, xntol, nmu
                          read(3,*) ; read(3,*) FiniteNucleus
			  55 continue
			close(3)
		end subroutine potential_read

		subroutine potential_init
        !call wlog('in potential_init')
            call potential_allocate
			title(:) = ' '
			mpot = 1
			nph = 0
			ntitle = 0
			ihole = 1
			ipr1 = 0
			iafolp = 0
			iunf = 0
			nmix = 1
			nohole = -1
			jumprm = 0
			inters = 0
			nscmt = 0
			icoul = 0
			ixc = 0
                        iscfxc = 11 !-LC- default for SCF cycle XC functional
			lfms1 = 0
			iz(:) = -1
			lmaxsc(:) = 0
			rfms1 = -1 * 1.e0
			ca1 = 0.d0
			ecv = -40*1.d0
			rgrd = 0.05 * 1.d0
			totvol = 0.d0
			gamach = 0.d0 !initialized later by setgam
			xnatph(:) = 0.d0
			spinph(:) = -1.d10
			xion(:) = 0.d0
			folp(:) = 1.d0
		    ExternalPot = .false.
			StartFromFile = .false. !KJ added 12-10
			novr(:) = 0
			iphovr(:,:)=0 !KJ added 7-09
			nnovr(:,:)=0 !KJ
			rovr(:,:) = 0.d0 !KJ
			ChSh_Type = 0 !Fer : standard feff
			configtype=1 !KJ 12-2010 standard feff9
			corval_emin=-70.d0 ! eV
			scf_temperature = 0.0 ! eV
			scf_thermal_vxc = 1
            xntol = 1.D-4
            nmu = 100
            iscfth = 1
                        FiniteNucleus = .FALSE.
		end subroutine potential_init

	end module

!=======================================================================
!     LDOS
!=======================================================================

      module ldos_inp
	    use atoms_inp,only: nph
        use potential_inp,only: ixc, rgrd, iscfxc
		use global_inp,only: ispin
		use dimsmod,only : nphx=>nphu,nex
		implicit none
		character(*),parameter,private :: filename='ldos.inp'
		integer mldos, lfms2, minv
        integer, allocatable :: lmaxph(:)
		double precision emin, emax, eimag
		integer neldos
		real rdirec, toler1, toler2, rfms2

		contains

        subroutine ldos_allocate
        !call wlog('in ldos_allocate')
            if (.not.allocated(lmaxph)) allocate(lmaxph(0:nphx))
        end subroutine ldos_allocate

		subroutine ldos_write
			integer iph
            !Sanity check on number of points in energy mesh for LDOS:
            if(neldos.gt.nex) then
               neldos=nex
               call wlog('# energy points for LDOS reduced to compiled limit "nex".')
            endif
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mldos, lfms2, ixc, ispin, minv, neldos, iscfxc'
			  write(3,21)  mldos, lfms2, ixc, ispin, minv, neldos, iscfxc
			  write(3,10) 'rfms2, emin, emax, eimag, rgrd'
			  write(3,30)  rfms2, emin, emax, eimag, rgrd
			  write(3,10) 'rdirec, toler1, toler2'
			  write(3,30)  rdirec, toler1, toler2
			  write(3,10) ' lmaxph(0:nph)'
			  write(3,20)  (lmaxph(iph),iph=0,nph)
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
      20  format (20i4)
      21  format (5i4,1x,i7,1x,i4)
	  30  format (9f13.5)
		end subroutine ldos_write

		subroutine ldos_read
			integer iph
            call ldos_init
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*)  mldos, lfms2, ixc, ispin, minv, neldos, iscfxc
			  read(3,*) ; read(3,*)  rfms2, emin, emax, eimag, rgrd
			  read(3,*) ; read(3,*)  rdirec, toler1, toler2
			  read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
			close(3)
		end subroutine ldos_read

		subroutine ldos_init
        !call wlog('in ldos_init')
            call ldos_allocate
			mldos = 0
			lfms2 = 0
			minv = 0
			emax = 0.d0
			emin = 1000*1.d0
			eimag = -1*1.d0
			neldos = 101
			rfms2 = -1 * 1.e0
			rdirec = -1 * 1.e0
			toler1 = 1.d-3
			toler2 = 1.d-3
			lmaxph(:) = 0
		end subroutine ldos_init

	end module



!=======================================================================
!     SCREEN
!=======================================================================

      module screen_inp
	    use atoms_inp,only: nph
		implicit none

                TYPE ScreenInputVars
                   integer ner, nei, maxl, irrh, iend, lfxc, nrptx0
                   double precision emin, emax, eimax, ermin, rfms
                END TYPE ScreenInputVars

                character(*),parameter,private :: filename='screen.inp'
                TYPE(ScreenInputVars) ScreenI

		contains

		subroutine screen_write
		           open(unit=3,file=filename,status='unknown')
				   write(3,*) 'ner',ScreenI%ner
				   write(3,*) 'nei',ScreenI%nei
				   write(3,*) 'maxl',ScreenI%maxl
				   write(3,*) 'irrh',ScreenI%irrh
				   write(3,*) 'iend',ScreenI%iend
				   write(3,*) 'lfxc',ScreenI%lfxc
				   write(3,*) 'emin',ScreenI%emin
				   write(3,*) 'emax',ScreenI%emax
				   write(3,*) 'eimax',ScreenI%eimax
				   write(3,*) 'ermin',ScreenI%ermin
				   write(3,*) 'rfms',ScreenI%rfms
				   write(3,*) 'nrptx0',ScreenI%nrptx0
                   close(3)
				   return
		end subroutine screen_write

        subroutine screen_inp_parse(str,vars)
		   implicit none
		   character*3,intent(in) :: str
		   real*8,intent(in) ::  vars
				if (str .eq. 'ner') then
				   ScreenI%ner   = vars
				elseif (str .eq. 'nei') then
				   ScreenI%nei   = vars
				elseif (str .eq. 'max') then
				   ScreenI%maxl  = vars
				elseif (str .eq. 'irr') then
				   ScreenI%irrh  = vars
				elseif (str .eq. 'ien') then
				   ScreenI%iend  = vars
				elseif (str .eq. 'lfx') then
				   ScreenI%lfxc  = vars
				elseif (str .eq. 'emi') then
				   ScreenI%emin  = vars
				elseif (str .eq. 'ema') then
				   ScreenI%emax  = vars
				elseif (str .eq. 'eim') then
				   ScreenI%eimax = vars
				elseif (str .eq. 'erm') then
				   ScreenI%ermin = vars
				elseif (str .eq. 'rfm') then
				   ScreenI%rfms  = vars
				elseif (str .eq. 'nrp')then
				   ScreenI%nrptx0  = vars
				else
				   call wlog("Unrecognized keyword submitted to screen.inp in SCREEN_INP_PARSE ; aborting.")
				   stop
                endif
				return
		end subroutine screen_inp_parse

        subroutine screen_inp_parse_and_write(str,vars)
		!KJ No longer used (1-2012).  Used in a previous version of feff.
		   implicit none
		   character*3,intent(in) :: str
		   real*8,intent(in) ::  vars
		        open(unit=3,file=filename,status='unknown',access='append')
				if (str .eq. 'ner') then
				   ScreenI%ner   = vars
				   write(3,*) 'ner',ScreenI%ner
				elseif (str .eq. 'nei') then
				   ScreenI%nei   = vars
				   write(3,*) 'nei',ScreenI%nei
				elseif (str .eq. 'max') then
				   ScreenI%maxl  = vars
				   write(3,*) 'maxl',ScreenI%maxl
				elseif (str .eq. 'irr') then
				   ScreenI%irrh  = vars
				   write(3,*) 'irrh',ScreenI%irrh
				elseif (str .eq. 'ien') then
				   ScreenI%iend  = vars
				   write(3,*) 'iend',ScreenI%iend
				elseif (str .eq. 'lfx') then
				   ScreenI%lfxc  = vars
				   write(3,*) 'lfxc',ScreenI%lfxc
				elseif (str .eq. 'emi') then
				   ScreenI%emin  = vars
				   write(3,*) 'emin',ScreenI%emin
				elseif (str .eq. 'ema') then
				   ScreenI%emax  = vars
				   write(3,*) 'emax',ScreenI%emax
				elseif (str .eq. 'eim') then
				   ScreenI%eimax = vars
				   write(3,*) 'eimax',ScreenI%eimax
				elseif (str .eq. 'erm') then
				   ScreenI%ermin = vars
				   write(3,*) 'ermin',ScreenI%ermin
				elseif (str .eq. 'rfm') then
				   ScreenI%rfms  = vars
				   write(3,*) 'rfms',ScreenI%rfms
				elseif (str .eq. 'nrp')then
				   ScreenI%nrptx0  = vars
				   write(3,*) 'nrptx0',ScreenI%nrptx0
				else
				   call wlog("Unrecognized keyword submitted to screen.inp ; aborting.")
				   stop
                endif
				close(3)
				return
		end subroutine screen_inp_parse_and_write


		subroutine screen_read
		    ! Reads screen.inp.  This routine is set up a little different from its brothers in the other input modules.
			! This is to keep it compatible with situations where there either is no screen.inp file (in which case defaults are used for all variables),
			! and with situations where screen.inp contains only the variables for which non-default values are specified.
			! This is because I've only added mandatory screen.inp files being written by rdinp now 1-2012.  KJ
			integer i
			character*8 strs
			character*3 str
			double precision vars
			call screen_init  !KJ set defaults in case screen.inp doesn't exist!
			open (file=filename, unit=3, status='old', err=60)
!KJ			read (3,*)  !KJ 11-2011 removing header line from screen.inp because it is incompatible with screen_inp_and_parse above.
			do i = 1, 12
				read(3,*,end=60)  strs, vars
				str = strs(1:3)
				if (str .eq. 'ner') ScreenI%ner   = nint(vars)
				if (str .eq. 'nei') ScreenI%nei   = nint(vars)
				if (str .eq. 'max') ScreenI%maxl  = nint(vars)
				if (str .eq. 'irr') ScreenI%irrh  = nint(vars)
				if (str .eq. 'ien') ScreenI%iend  = nint(vars)
				if (str .eq. 'lfx') ScreenI%lfxc  = nint(vars)
				if (str .eq. 'emi') ScreenI%emin  = vars
				if (str .eq. 'ema') ScreenI%emax  = vars
				if (str .eq. 'eim') ScreenI%eimax = vars
				if (str .eq. 'erm') ScreenI%ermin = vars
				if (str .eq. 'rfm') ScreenI%rfms  = vars
				if (str .eq. 'nrp') ScreenI%nrptx0  = nint(vars)
			end do
		  60 continue
		  close(3)
                  return
		end subroutine screen_read

		subroutine screen_init
		  ScreenI%ner   = 40
		  ScreenI%nei   = 20
		  ScreenI%maxl  = 4
		  ScreenI%irrh  = 1
		  ScreenI%iend  = 0
		  ScreenI%emin  = -40.0d0 !KJ This and next 3 values are in eV ; converted to Ha at a later point in the code (screen/rdgeom.f90)
		  ScreenI%emax  = 0.0d0
		  ScreenI%eimax = 2.0d0
		  ScreenI%ermin = 0.001d0
		  ScreenI%lfxc  = 0
		  ScreenI%rfms  = 4.0d0
		  ScreenI%nrptx0 = 251
		end subroutine screen_init

	end module

!=======================================================================
!     CRPA
!=======================================================================

      module crpa_inp
		implicit none

        TYPE CRPAInputVars
           integer  l_crpa, do_CRPA
           double precision  rcut
        END TYPE CRPAInputVars

        character(*),parameter,private :: filename='crpa.inp'
        TYPE(CRPAInputVars) CRPAI

		contains

		subroutine crpa_write
		           open(unit=3,file=filename,status='unknown')
				   write(3,*) 'do_CRPA', CRPAI%do_CRPA
				   write(3,*) 'rcut', CRPAI%rcut
				   write(3,*) 'l_crpa', CRPAI%l_crpa
                   close(3)
				   return
		end subroutine crpa_write


		subroutine crpa_read
			integer i
			character*8 strs
			character*3 str
			double precision vars
			call crpa_init
			open (file=filename, unit=3, status='old', err=60)
			do i = 1, 13
				read(3,*,end=60)  strs, vars
				str = strs(1:3)
				if (strs .eq. 'do_CRPA') CRPAI%do_CRPA= nint(vars)
				if (strs .eq.  'rcut') CRPAI%rcut = vars
				if (strs .eq. 'l_crpa') CRPAI%l_crpa = nint(vars)

			end do
		  60 continue
		  close(3)
                  return
		end subroutine crpa_read

		subroutine crpa_init
		  CRPAI%do_CRPA= 0
		  CRPAI%rcut = 1.6
		  CRPAI%l_crpa= 3
		end subroutine crpa_init

	end module

!=======================================================================
!     OPCONS
!=======================================================================
      MODULE opcons_inp
        USE dimsmod, only: nphx=>nphu
        LOGICAL run_opcons, print_eps
        REAL(8), allocatable :: NumDens(:)
		character(*),parameter,private :: filename='opcons.inp'

        CONTAINS

           subroutine opcons_allocate
               if(.not.allocated(NumDens)) allocate(NumDens(0:nphx))
           end subroutine opcons_allocate

           SUBROUTINE opcons_init
              call opcons_allocate
              run_opcons = .FALSE.
              print_eps  = .FALSE.
              NumDens(:) = -1.d0
           END SUBROUTINE opcons_init

           SUBROUTINE opcons_write
              INTEGER iph

              OPEN(FILE=filename,UNIT=8,STATUS='REPLACE')

              WRITE(8,'(A)') 'run_opcons'
              WRITE(8,*) run_opcons
              WRITE(8,'(A)') 'print_eps'
              WRITE(8,*) print_eps
              WRITE(8,'(A)') 'NumDens(0:nphx)'
              WRITE(8,*) NumDens(0:nphx)

              CLOSE(8)
           END SUBROUTINE opcons_write

           SUBROUTINE opcons_read
              INTEGER iph

              OPEN(FILE=filename,UNIT=8,STATUS='OLD')
              call opcons_allocate
              READ(8,*)
              READ(8,*) run_opcons
              READ(8,*)
              READ(8,*) print_eps
              READ(8,*)
              READ(8,*) NumDens(0:nphx)
           END SUBROUTINE opcons_read

      END MODULE opcons_inp
!=======================================================================
!     XSPH
!=======================================================================

      module xsph_inp
        use dimsmod, only: nphx=>nphu
		use global_inp
		use potential_inp
		use ldos_inp
		implicit none
		character(*),parameter,private :: filename='xsph.inp'
		integer mphase, ipr2, ixc0, lreal, iPlsmn
		integer iGammaCH, iGrid, NPoles, iCoreState
		character*6, allocatable ::  potlbl(:)
!		potlbl(0:nphx)    -   label for user convienence
		double precision xkstep, xkmax, vixan, vr0, vi0, Eps0, EGap
		integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis

    real*8 electronic_temperature
!		!KJ for the energy grid card EGRID :
		integer iegrid,egrid3a
		real*8 egrid3b,egrid3c
		character*100 egridfile
        logical lopt, PrintRL

		contains

        subroutine xsph_allocate
           if(.not.allocated(potlbl)) allocate (potlbl(0:nphx))
        end subroutine xsph_allocate

		subroutine xsph_write
			integer iph
			open (file=filename, unit=3, status='unknown')
		!     Josh - added flag for PLASMON card (iPlsmn = 0, 1, or 2)
			  write(3,10) 'mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,iPlsmn,NPoles,iGammaCH,iGrid,iCoreState,iscfxc'
			  write(3,20)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,   &
			 &        iPlsmn, NPoles, iGammaCH, iGrid, iCoreState, iscfxc
			  write(3,10) 'vr0, vi0'
			  write(3,30)  vr0, vi0
			  write(3,10) ' lmaxph(0:nph)'
			  write(3,20)  (lmaxph(iph),iph=0,nph)
			  write(3,10) ' potlbl(iph)'
			  write(3,170)  (potlbl(iph),iph=0,nph)
		  170   format (13a6)
			  write(3,10) 'rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap'
			  write(3,30)  rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap
			  write(3,10) 'spinph(0:nph)'
			  write(3,30)  (spinph(iph),iph=0,nph)
			  write(3,10) 'izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis'
			  write(3,20)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
  		  write(3,10) 'electronic temperature'
  		  write(3,30) electronic_temperature
			! Commented out by Fer
			! The following lines are commented out because they are not being read
			! in rexsph (commented out by JK). This screws up anything that comes after
			! them in mod2.inp (for example, the ChSh parameters that I'm including.
			!!KJ next lines contain EGRID variables ; added 01-07
			!        write(3,10) 'iegrid,egrid3a,egrid3b,egrid3c'
			!          write(3,'(2i4,2f13.5)') iegrid,egrid3a,egrid3b,egrid3c !format statement is a mix of 20 and 30
			!          write(3,10) 'egridfile'
			!          write(3,10) egridfile
			!!KJ
			! Added by Fer
			! Correction of the excitation energy for chemical shifts
			  write(3,10) 'ChSh_Type:'
			  write(3,20) ChSh_Type
			!KJ 7-09 Next 2 lines for feff8q
			  write(3,'(a)') ' the number of decomposition channels ; only used for nrixs'
			  write(3,'(i5)') ldecmx
			  write(3,'(a)') 'lopt'
                          write(3,*) lopt
			  write(3,'(a)') 'PrintRL'
			  write(3,*) PrintRL
			  close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (21i4)
	  30  format (20f13.5)
		end subroutine xsph_write

		subroutine xsph_read
			integer iph
			open (file=filename, unit=3, status='old')
            call xsph_allocate
            call ldos_allocate  ! allocates lmaxph
            call potential_allocate ! allocates spinph
			  read(3,*) ; read(3,*)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,iPlsmn, NPoles, iGammaCH, iGrid, iCoreState, iscfxc
			  read(3,*) ; read(3,*)  vr0, vi0
			  read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
			  read(3,*) ; read(3,'(13a6)')  (potlbl(iph),iph=0,nph)
			  read(3,*) ; read(3,*)  rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap
			  read(3,*) ; read(3,*)  (spinph(iph),iph=0,nph)
			  read(3,*) ; read(3,*)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
			!!KJ next lines contain EGRID variables ; added 01-07
			!          read(3,*) ; read(3,'(2i4,2f13.5)') iegrid,egrid3a,egrid3b,egrid3c !format statement is a mix of 20 and 30
			!          read(3,*) ; read(3,10) egridfile
        read(3,*) ; read(3,*) electronic_temperature
			  read(3,*) ; read(3,*) ChSh_Type
			!KJ 7-09 Next 2 lines for feff8q
			  read(3,*) ; read(3,*) ldecmx
			  read(3,*) ; read(3,*) lopt
			  read(3,*) ; read(3,*) PrintRL
			close(3)
		end subroutine xsph_read

		subroutine xsph_init
            call ldos_allocate  ! allocates lmaxph
            call potential_allocate ! allocates spinph
            call xsph_allocate
            lopt = .false.
            PrintRl = .false. ! JK - this will be turned on by RIXS card.
            iCoreState = -1  ! JK - this is used to set core state for
                             ! Matrix elements
			izstd = 0
			ifxc = 0
			ipmbse = 0
			itdlda = 0
			nonlocal = 0
			ibasis = 0
			potlbl(0:nphx) = ' '
			mphase = 1
			ipr2 = 0
			ixc0 = -1
			lreal = 0
			iPlsmn = 0 ! Josh Kas
            NPoles = 100 ! JJK 3/9/2010
            EGap = 0.d0 ! JJK 4/2010
			iGammaCH = 0
			iGrid = 0
			vr0 = 0.d0
			vi0 = 0.d0
			xkmax = 20*1.d0
			xkstep = 0.07*1.d0
			vixan = 0.d0
			iegrid=0 !KJ for EGRID card 1-07
			egridfile=' '
			egrid3a=0
			egrid3b=dble(0)
			egrid3c=dble(0)
			electronic_temperature = dble(0)
		end subroutine xsph_init


      end module



!=======================================================================
!     FMS
!=======================================================================

      module fms_inp
        use ldos_inp, only: lmaxph, ldos_allocate, minv,rfms2,rdirec,toler1,toler2,nph
		use global_inp,only: ldecmx
		implicit none
		character(*),parameter,private :: filename='fms.inp'
		integer mfms, idwopt, ipr3 !ipr3 is currently dummy - not in fms.inp
		!KJ June 2013: do_fms added.  Note that program control over the fms module in the FEFF code is flawed.
		!The above "mfms" determines whether or not to execute the "fms" module, and this depends only on the
		!corresponding entry in the CONTROL card.  E.g. if the fms module has been executed earlier and its results
		!are still valid, the user saves much time by not recalculating fms.
		!However, a separate question is whether or not the fms algorithm is used to produce results rather than the path
		!summation algorithm.  This is governed by the choice between FMS and RPATH cards.  Even if the FMS card is not used,
		!one still needs to execute the fms and mkgtr modules in order to generate a fms.bin file.  This fms.bin file will then
		!simply contain a lot of 0's, but it is required for proper execution of FEFF.  Similarly, if the fms algorithm is
		!used for results, one still needs to execute the path expansion to obtain a zero result for PE.  The total result of the
		!FEFF calculation is always the sum of FMS+PE, one of which is normally zero.
		!This setup worked alright for r-space calculations, where in the case of PE (no FMS card in feff.inp) the fms module
		!wastes a second on xprep and then exits because of a control loop that checks that the cluster contains at least 1 atom.
		!Configuration routines in rdinp have set it to 0 in this case.
		!For k-space calculations however the cluster does not exist, and there is no obvious way to distinguish FMS/PE based on such
		!derivative variables.  Therefore the fms algorithm will always try to run.  Currently it seems mostly to crash during
		!kprep ("STOP IN STRAA").  This at least keeps it from overwriting old fms.bin files (e.g. obtained in real-space), but it is
		!not ideal.  Hence the new variable do_fms:
		!KJ do_fms =1 if the FMS card is in feff.inp, and =0 if the FMS card is not in feff.inp.
		!It is equivalent to (inclus.lt.1) for r-space calculations.
		integer do_fms
		real rprec
		!KJ rprec seems to be bogus input, i.e. not used anywhere in entire FEFF90.  Set to 0 and kept here for compatibility.
		double precision   tk, thetad, sig2g
		logical save_gg_slice

		contains

		subroutine fms_write
			implicit none
			integer iph
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mfms, idwopt, minv'
			  write(3,20)  mfms, idwopt, minv
			  write(3,10) 'rfms2, rdirec, toler1, toler2'
			  write(3,30)  rfms2, rdirec, toler1, toler2
			  write(3,10) 'tk, thetad, sig2g'
			  write(3,30)  tk, thetad, sig2g
			  write(3,10) ' lmaxph(0:nph)'
			  write(3,20)  (lmaxph(iph),iph=0,nph)
			  !KJ 7-09 Next 2 lines for feff8q
			  write(3,'(a24)') ' the number of decomposition channels'
			  write(3,'(i5)') ldecmx
			  write(3,10) ' save_gg_slice'
			  write(3,'(l)') save_gg_slice
			  write(3,10) 'do_fms'
			  write(3,20) do_fms
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine fms_write

		subroutine fms_read
		    integer iph
			open (file=filename, unit=3, status='unknown')
            call ldos_allocate
			  read(3,*) ; read(3,*)  mfms, idwopt, minv
			  read(3,*) ; read(3,*)  rfms2, rdirec, toler1, toler2
			  read(3,*) ; read(3,*)  tk, thetad, sig2g
			  read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
			  !KJ 7-09 Next line for feff8q
			  read(3,*) ; read(3,*) ldecmx
			  read(3,*) ; read(3,*) save_gg_slice
			  read(3,*) ; read(3,*) do_fms
			close(3)
		end subroutine fms_read

		subroutine fms_init
            call ldos_allocate
			mfms = 1
			idwopt = -1
			sig2g = 0.d0
			thetad = 0.d0
			tk = 0.d0
			ipr3 = 0
			rprec = 0.e0
			save_gg_slice = .false.
			do_fms = 0
		end subroutine fms_init

	end module




!=======================================================================
!     PATHS
!=======================================================================

      module paths_inp
        use ldos_inp
		implicit none
		character(*),parameter,private :: filename='paths.inp'
		integer  mpath, ms, nncrit, nlegxx, ipr4, ica  !KJ added ica 6-06
		!KJ nncrit seems to be bogus input, i.e. not set in rdinp at all ; fully internal to PATH.  Set to 0 and kept here for compatibility.
		real critpw, pcritk, pcrith,  rmax

		contains

		subroutine paths_write
			implicit none
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mpath, ms, nncrit, nlegxx, ipr4'
			  write(3,20)  mpath, ms, nncrit, nlegxx, ipr4
			  write(3,10) 'critpw, pcritk, pcrith,  rmax, rfms2'
			  write(3,30)  critpw, pcritk, pcrith,  rmax, rfms2
			  write(3,10) 'ica' !KJ 6-06
			  write(3,20)  ica  !KJ 6-06
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine paths_write

		subroutine paths_read
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*)  mpath, ms, nncrit, nlegxx, ipr4
			  read(3,*) ; read(3,*)  critpw, pcritk, pcrith,  rmax, rfms2
			  read(3,*) ; read(3,*)  ica  !KJ 6-06
			close(3)
		end subroutine paths_read

		subroutine paths_init
			mpath = 1
			ms = 0
			ipr4 = 0
			ica=-1 !KJ 6-06
			critpw = 2.5*1.e0
			pcritk = 0.e0
			pcrith = 0.e0
			rmax = -1 * 1.e0
			nlegxx = 7
			nncrit = 0
		end subroutine paths_init

	end module



!=======================================================================
!     GENFMT
!=======================================================================

      module genfmt_inp
		use global_inp
		implicit none
		character(*),parameter,private :: filename='genfmt.inp'
		integer  mfeff, ipr5, iorder
		logical  wnstar
		double precision critcw

		contains

		subroutine genfmt_write
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mfeff, ipr5, iorder, critcw, wnstar'
			  write(3,180)  mfeff, ipr5, iorder, critcw, wnstar
		      180   format ( 2i4, i8, f13.5, L5)
			  !KJ 7-09 Next 2 lines for feff8q
			  write(3,'(a24)') ' the number of decomposition channels'
			  write(3,'(i5)') ldecmx
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine genfmt_write

		subroutine genfmt_read
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*)  mfeff, ipr5, iorder, critcw, wnstar
			  !KJ 7-09 Next line for feff8q
			  read(3,*) ; read(3,*) ldecmx
			close(3)
		end subroutine genfmt_read

		subroutine genfmt_init
			mfeff = 1
			ipr5 = 0
			iorder = 2
			wnstar = .false.
			critcw = 4*1.d0
		end subroutine genfmt_init

	end module



!=======================================================================
!     FF2X
!=======================================================================

    module ff2x_inp
		use global_inp
		use xsph_inp
		use fms_inp
		use genfmt_inp
		implicit none
		character(*),parameter,private :: filename='ff2x.inp'
		integer  mchi, ipr6, mbconv, absolu !KJ added absolu 3-06
		double precision  vrcorr, vicorr, s02, alphat, thetae
		!KJ 1-2013 a global DW factor that will broaden chi(E) by exp^-(sig2gk k)^2  .  In Angstrom.  (Unlike sig2g, which is in Angstrom^2 .)
		double precision sig_gk


		contains

		subroutine ff2x_write
			integer i
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH' !KJ added absolu 3-06
			  write(3,20)  mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH !KJ added absolu 3-06
			  write(3,10) 'vrcorr, vicorr, s02, critcw'
			  write(3,30)  vrcorr, vicorr, s02, critcw
			  write(3,10) 'tk, thetad, alphat, thetae, sig2g, sig_gk'
			  write(3,30)  tk, thetad, alphat, thetae, sig2g, sig_gk
			  !KJ 7-09 next 4 lines for feff8q
			  write(3,10) 'momentum transfer'
			  write(3, '(3f13.5)') (xivec(i),i=1,3)
			  write(3,'(a24)') ' the number of decomposition channels'
			  write(3,'(i5)') ldecmx
			  write(3,10) 'electronic temperature'
			  write(3,30) electronic_temperature
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine ff2x_write

		subroutine ff2x_read
			integer i
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*)  mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH !KJ added absolu 3-06
			  read(3,*) ; read(3,*)  vrcorr, vicorr, s02, critcw
			  read(3,*) ; read(3,*)  tk, thetad, alphat, thetae, sig2g, sig_gk !KJ added sig_gk 1-13
			  read(3,*) ; read(3, *) (xivec(i),i=1,3)
			  read(3,*) ; read(3,*) ldecmx
			  read(3,*) ; read(3,*) electronic_temperature
			close(3)
		end subroutine ff2x_read

		subroutine ff2x_init
			absolu=0  !KJ 3-06 for ABSOLUTE card
			mchi = 1
			ipr6 = 0
			mbconv = 0
			vicorr = 0.d0
			vrcorr = 0.d0
			s02 = 1.d0
			alphat = 0.d0
			thetae = 0.d0
			sig_gk = 0.d0
		end subroutine ff2x_init

	end module


!=======================================================================
!     SFCONV
!=======================================================================

      module sfconv_inp
		use global_inp, only : ispec
		use ff2x_inp, only : ipr6
		implicit none
		character(*),parameter,private :: filename='sfconv.inp'
		integer  msfconv, ipse, ipsk
		double precision wsigk, cen
		character(12) cfname

		contains

		subroutine sfconv_write
		!c    sfconv.inp - Josh Kas
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'msfconv, ipse, ipsk'
			  write(3,20)  msfconv, ipse, ipsk
			  write(3,10) 'wsigk, cen'
			  write(3,30) wsigk, cen
			  write(3,10) 'ispec, ipr6'
			  write(3,20)  ispec, ipr6
			  write(3,10) 'cfname'
			  write(3,10) cfname
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine sfconv_write

		subroutine sfconv_read
			open (file=filename, unit=3, status='old')
			read (3,*) ; read (3,*)  msfconv, ipse, ipsk
			read (3,*) ; read (3,*)  wsigk, cen
			read (3,*) ; read (3,*)  ispec, ipr6
			read (3,*) ; read (3,*)  cfname
			close(3)
		end subroutine sfconv_read

		subroutine sfconv_init
			msfconv = 0 ! Josh Kas
			ipse = 0
			ipsk = 0
			wsigk = 0.d0 ! Josh Kas
			cen = 0.d0 ! Josh Kas
			cfname = 'NULL'
		end subroutine sfconv_init

	end module



!=======================================================================
!     EELS
!=======================================================================

      module eels_inp
!		Beam direction in crystal frame of feff.inp (a.u.)
		use global_inp,only: xivec
		implicit none
		character(*),parameter,private :: filename='eels.inp'
!		Beam energy in eV :
		real*8 ebeam
!		Convergence semiangle in rad :
		real*8 aconv
!		Collection semiangle in rad :
		real*8 acoll
!		Integration mesh for q-vectors (radial/angular mesh size)
        integer nqr,nqf
!		Detector position ; angles in rad w.r.t. x and y directions
		real*8 thetax,thetay
!       what kind of q-mesh : uniform (U), logarithmic (L), or one dimensional logarithmic (1)
!       not currently in eels.inp/feff.inp
        character*1      qmodus
!		Parameter for logarithmic mesh - not currently in eels.inp/feff.inp
		real*8 th0
!		Make magic angle plot if magic=1
		integer        magic
!		Evaluate magic angle at this energy point
		real*8        emagic
!		Orientation sensitive?
		integer        aver
!		Do we have cross-terms?
		integer        cross
!		Do we do anything at all?
		integer        eels
!		How many spectra to combine
		integer ipmin,ipmax,ipstep ,nip
!		Where do we take input from :
		integer iinput
!		Which column? - to be replaced by more advanced switch
		integer spcol
!       Relativistic calculation or not?  Converted into logical inside eels-module.
		integer relat

		contains

		subroutine eels_write
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'calculate ELNES?'
			  write(3,20) eels
			  write(3,10) 'average? relativistic? cross-terms? Which input?'
			  write(3,20) aver, relat, cross, iinput, spcol
			  write(3,10) 'polarizations to be used ; min step max'
			  write(3,20) ipmin,ipstep,ipmax
			  write(3,10) 'beam energy in eV'
			  write(3,30) ebeam
			  write(3,10) 'beam direction in arbitrary units'
			  write(3,30) xivec
			  write(3,10) 'collection and convergence semiangle in rad'
			  write(3,30) acoll,aconv
			  write(3,10) 'qmesh - radial and angular grid size'
			  write(3,20) nqr,nqf
			  write(3,10) 'detector positions - two angles in rad'
			  write(3,30) thetax,thetay
			  write(3,10) 'calculate magic angle if magic=1'
			  write(3,20) magic
			  write(3,10) 'energy for magic angle - eV above threshold'
			  write(3,30) emagic
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine eels_write

		subroutine eels_read
			open (file=filename, unit=3, status='old',err=100)
			read(3,*) ; read(3,*,end=100,err=100) eels
			read(3,*) ; read(3,*,err=209) aver, relat, cross, iinput,spcol ; goto 210
			209   iinput=1;spcol=4;relat=1;cross=1;aver=0  !restore defaults - this construction for older files.
			210 read(3,*) ; read(3,*) ipmin,ipstep,ipmax  !KJ this un-Kevin-like construction for older files ...
			nip=1+((ipmax-ipmin)/ipstep)
			read(3,*) ; read(3,*) ebeam
			read(3,*) ; read(3,*) xivec
			read(3,*) ; read(3,*) acoll,aconv
			read(3,*) ; read(3,*) nqr,nqf
			read(3,*) ; read(3,*) thetax,thetay
			read(3,*) ; read(3,*) magic
			read(3,*) ; read(3,*) emagic
			close(3)
			return
100			eels = 0  ; ipmin=1 ; ipmax=1 ; ipstep=1 ! no eels.inp -> don't do eels
            return
		end subroutine eels_read


        subroutine eels_init !default values for everything (except xivec)
		  ebeam=0.d0
		  aconv=0.d0
		  acoll=0.d0
		  nqr=0
		  nqf=0
		  magic=0
		  emagic=0.d0
		  eels=0
		  relat=1
		  cross=1
		  aver=0
		  thetax=0.d0
		  thetay=0.d0
		  ipmin=1
		  ipmax=1
		  nip=1
		  ipstep=1
          iinput=1  ! xmu.dat - files from ff2x
          spcol=4   ! xmu.dat - use spectrum mu(omega)
          qmodus='U'  !  U for uniform grid
          th0=0.d0
		end subroutine eels_init

	end module



!=======================================================================
!     COMPTON
!=======================================================================


    module compton_inp
      implicit none
      character(*),parameter,private :: filename='compton.inp'

      ! spatial and momentum grid parameters
      integer :: ns, nphi, nz, nzp, npq
      real :: smax, phimax, zmax, zpmax, pqmax

      ! flags
      logical :: do_compton, do_rhozzp
      logical :: force_jzzp
      logical :: do_rho_xy, do_rho_yz, do_rho_xz, do_rho_vol, do_rho_line
	  integer run_compton_module

      ! apodization function type
      integer :: window
      real :: window_cutoff

      real :: temperature
      logical ::  set_chemical_potential
      real :: chemical_potential

      real*8  :: qhat(3)

      integer, parameter :: WINDOW_STEP = 0, WINDOW_HANNING = 1
    contains
      subroutine compton_write
        if (do_compton .or. do_rhozzp .or. do_rho_xy .or. do_rho_yz .or. do_rho_xz .or. do_rho_vol .or. do_rho_line) then
		   run_compton_module=1
		else
		   run_compton_module=0
		endif
        open (file=filename, unit=3, status='unknown')
		  write(3,10) 'run compton module?'
		  write(3,*)  run_compton_module
          write(3,10) 'pqmax, npq'
          write(3,*) pqmax, npq
          write(3,10) 'ns, nphi, nz, nzp'
          write(3,20) ns, nphi, nz, nzp
          write(3,10) 'smax, phimax, zmax, zpmax'
          write(3,30) smax, phimax, zmax, zpmax
          write(3,10) 'jpq? rhozzp? force_recalc_jzzp?'
          write(3,*) do_compton, do_rhozzp, force_jzzp
          write(3,10) 'window_type (0=Step, 1=Hann), window_cutoff'
          write(3,*) window, window_cutoff
          write(3,10) 'temperature (in eV)'
          write(3,30) temperature
          write(3,10) 'set_chemical_potential? chemical_potential(eV)'
          write(3,*) set_chemical_potential, chemical_potential
          write(3,10) 'rho_xy? rho_yz? rho_xz? rho_vol? rho_line?'
          write(3,*) do_rho_xy, do_rho_yz, do_rho_xz, do_rho_vol, do_rho_line
          write(3,10) 'qhat_x qhat_y qhat_z'
          write(3,*) qhat(1), qhat(2), qhat(3)
        close(3)
		! standard formats for string, integers, real numbers
    10  format(a)
    20  format (20i4)
    30  format (9f13.5)
      end subroutine compton_write

      subroutine compton_read
        open (file=filename, unit=3, status='old',err=100)
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) run_compton_module
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) pqmax, npq
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) ns, nphi, nz, nzp
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) smax, phimax, zmax, zpmax
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) do_compton, do_rhozzp, force_jzzp
          read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) window, window_cutoff
          read(3,*,end=200,err=100) ; read(3,*,end=100,err=100) temperature
          read(3,*,end=200,err=100) ; read(3,*,end=100,err=100) set_chemical_potential, chemical_potential
          read(3,*,end=200,err=100) ; read(3,*,end=100,err=100) do_rho_xy, do_rho_yz, do_rho_xz, do_rho_vol, do_rho_line
          read(3,*,end=200,err=100) ; read(3,*,end=100,err=100) qhat(1), qhat(2), qhat(3)
    200 continue ! recoverable errors

        close(3)
        return

    100 run_compton_module=0  ! no (or invalid) compton.inp -> don't do compton
        return
      end subroutine compton_read

      subroutine compton_inp_init
        real, parameter :: pi  = 3.1415926535897932384626433832795
        ns   = 32
        nphi = 32
        nz   = 32
        nzp  = 144

        smax   = 0
        phimax = 2*pi
        zmax   = 0
        zpmax  = 10.0

        npq   = 1000
        pqmax = 5.0

        do_compton     = .false.
        do_rhozzp  = .false.
        force_jzzp = .false.
        do_rho_xy  = .false.
        do_rho_yz  = .false.
        do_rho_xz  = .false.
        do_rho_vol = .false.
        do_rho_line = .false.
		run_compton_module=0

        window = WINDOW_HANNING
        window_cutoff = 0

        temperature = 0.0
        set_chemical_potential = .false.
        chemical_potential = 0

        qhat(:) = 0.0; qhat(3) = 1.0

      end subroutine compton_inp_init
    end module compton_inp

!=======================================================================
!     RIXS
!=======================================================================

MODULE rixs_inp
  USE constants
  TYPE RixsInp
     INTEGER m_run, nEdges
     REAL(8) gam_ch, gam_exp(2), EMin(2), EMax(2), xmu, EdgeSplit(10)
     LOGICAL ReadPoles, SkipCalc, MBConv, ReadSigma
     CHARACTER(LEN=30) Edges(2)
  END TYPE RixsInp

  TYPE(RixsInp) RixsI
  CHARACTER(8), PARAMETER :: filename = 'rixs.inp'
  CONTAINS
    SUBROUTINE rixs_write
      INTEGER iEdge
      ! Change to code units.
      RixsI%gam_ch = RixsI%gam_ch/hart
      RixsI%gam_exp = RixsI%gam_exp/hart
      RixsI%EMin = RixsI%EMin/hart
      RixsI%EMax = RixsI%EMax/hart
      RixsI%xmu = RixsI%xmu/hart
      OPEN(3,FILE=filename,STATUS='UNKNOWN')
      WRITE(3,*) 'm_run'
      WRITE(3,*) RixsI%m_run
      WRITE(3,*) 'gam_ch, gam_exp(1), gam_exp(2)'
      WRITE(3,'(4f20.10)') RixsI%gam_ch, RixsI%gam_exp(1), RixsI%gam_exp(2)
      WRITE(3,*) 'EMinI, EMaxI, EMinF, EMaxF'
      WRITE(3,'(4f20.10)') RixsI%EMin(1), RixsI%EMax(1), RixsI%EMin(2), RixsI%EMax(2)
      WRITE(3,*) 'xmu'
      WRITE(3,*) RixsI%xmu
      WRITE(3,*) 'Readpoles, SkipCalc, MBConv, ReadSigma'
      WRITE(3,*) RixsI%Readpoles, RixsI%SkipCalc, RixsI%MBConv, RixsI%ReadSigma
      WRITE(3,*) 'nEdges'
      WRITE(3,*) RixsI%nEdges
      DO iEdge = 1, RixsI%nEdges
         WRITE(3,*) 'Edge', iEdge
         WRITE(3,*) TRIM(ADJUSTL(RixsI%Edges(iEdge)))
      END DO
      CLOSE(3)
    END SUBROUTINE rixs_write

    SUBROUTINE rixs_read
      INTEGER iEdge
      OPEN(3,FILE=filename,STATUS='OLD')
      READ(3,*)
      READ(3,*) RIXSI%m_run
      READ(3,*)
      READ(3,*) RixsI%gam_ch, RixsI%gam_exp(1), RixsI%gam_exp(2)
      READ(3,*)
      READ(3,*)RixsI%EMin(1), RixsI%EMax(1), RixsI%EMin(2), RixsI%EMax(2)
      READ(3,*)
      READ(3,*) RixsI%xmu
      READ(3,*)
      READ(3,*) RixsI%Readpoles, RixsI%SkipCalc, RixsI%MBConv, RixsI%ReadSigma
      READ(3,*)
      READ(3,*) RixsI%nEdges
      DO iEdge = 1, RixsI%nEdges
         READ(3,*)
         READ(3,*) RixsI%Edges(iEdge)
      END DO
      CLOSE(3)
    END SUBROUTINE rixs_read

    SUBROUTINE rixs_init
      RixsI%m_run = 0
      RixsI%gam_ch = 0.1d0/hart
      RixsI%gam_exp(:) = 0.1d0/hart
      RixsI%EMin = 0.d0
      RixsI%EMax = 0.d0
      RixsI%xmu = -1.d10
      RixsI%ReadPoles = .TRUE.
      RixsI%SkipCalc = .FALSE.
      RixsI%MBConv = .FALSE.
      RixsI%ReadSigma = .FALSE.
      RixsI%Edges = 'NULL'
      RixsI%nEdges = 0
    END SUBROUTINE rixs_init

    SUBROUTINE rixs_set(gam_ch, gam_exp, EMin, EMax, xmu, ReadPoles, SkipCalc, MBConv, ReadSigma)
      REAL(8) gam_ch, gam_exp(2), EMin(2), EMax(2), xmu, EdgeSplit(10)
      LOGICAL ReadPoles, SkipCalc, MBConv, ReadSigma
      gam_ch = RixsI%gam_ch
      gam_exp = RixsI%gam_exp
      EMin = RixsI%EMin
      EMax = RixsI%EMax
      xmu = RixsI%xmu
      ReadPoles = RixsI%ReadPoles
      SkipCalc = RixsI%SkipCalc
      MBConv = RixsI%MBConv
      ReadSigma = RixsI%ReadSigma
    END SUBROUTINE rixs_set
  END MODULE rixs_inp


!=======================================================================
!     BAND
!=======================================================================

      module band_inp
		implicit none
		character(*),parameter,private :: filename='band.inp'
		integer mband
		real*8 emin,emax,estep
		integer nkp,nep,ikpath
		logical freeprop

		contains

		subroutine band_write
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'mband : calculate bands if = 1'
			  write(3,20)  mband
			  write(3,10) 'emin, emax, estep : energy mesh'
			  write(3,30)  emin, emax, estep
			  write(3,10) 'nkp : # points in k-path'
			  write(3,20)  nkp
			  write(3,10) 'ikpath : type of k-path'
			  write(3,20)  ikpath
			  write(3,10) 'freeprop :  empty lattice if = T'
			  write(3,*)   freeprop
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine band_write

		subroutine band_read
			open (file=filename, unit=3, status='old')
			read (3,*) ; read (3,*)  mband
			read (3,*) ; read (3,*)  emin,emax,estep
			read (3,*) ; read (3,*)  nkp
			read (3,*) ; read (3,*)  ikpath
			read (3,*) ; read (3,*)  freeprop
			close(3)
		end subroutine band_read

		subroutine band_init
			emin = 0.d0
			emax = 0.d0
			estep = 0.d0
			nkp = 0
			nep = 0
			ikpath = -1
			mband = 0
			freeprop = .false.
		end subroutine band_init

	end module band_inp


!=======================================================================
!     HUBBARD
!=======================================================================

      module hubbard_inp
	    implicit none
		character(*),parameter,private :: filename='hubbard.inp'
		integer i_hubbard, mldos_hubb, l_hubbard
		real*8  U_hubbard,J_hubbard
        real*8  fermi_shift

		contains

		subroutine hubbard_write
			open (file=filename, unit=3, status='unknown')
			  write(3,10) 'i_hubbard mldos_hubb U_hubbard J_hubbard fermi_shift l_hubbard'
			  write(3,*)  i_hubbard, mldos_hubb, U_hubbard, J_hubbard, fermi_shift, l_hubbard
			close(3)
		! standard formats for string, integers and real numbers
	  10  format(a)
	  20  format (20i4)
	  30  format (9f13.5)
		end subroutine hubbard_write

		subroutine hubbard_read
			open (file=filename, unit=3, status='old')
			  read(3,*) ; read(3,*)   i_hubbard, mldos_hubb, U_hubbard, J_hubbard, fermi_shift, l_hubbard
			close(3)
		end subroutine hubbard_read

		subroutine hubbard_init
      i_hubbard=1
      mldos_hubb=1
      U_hubbard=0.0d0
      J_hubbard=0.0d0
      fermi_shift=0.0d0
      l_hubbard=0
		end subroutine hubbard_init

	end module


!=======================================================================
!     FULLSPECTRUM
!=======================================================================
    module fullspectrum_inp
     implicit none
     integer mFullSpectrum
     character(*),parameter,private :: filename='fullspectrum.inp'

     contains
       subroutine fullspectrum_init
          mFullSpectrum=0
       end subroutine
       subroutine fullspectrum_read
         open(file=filename,unit=3,status='old')
         read(3,*) ; read(3,*) mFullSpectrum
         close(3)
       end subroutine
       subroutine fullspectrum_write
        open(file=filename,unit=3,status='unknown')
        write(3,*) 'mFullSpectrum'
        write(3,*) mFullSpectrum
        close(3)
       end subroutine
    end module
