!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: crystalstructure.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:50 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CRYSTALSTRUCTURE(celvin)


! I believe this is the ONLY place in the entire FEFF code where the lattice of the crystal is changed (if change_to_primitive_lattice is true).
! The changes persist for the duration of the current FEFF module.  However they are not written to file.

        use struct,ncryst=>nsym
        use constants,only: pi2
        use controls,only:iprint
      implicit none
        real*8,intent(in)  :: celvin

!  Variables for the crystal structure :
      real*8, allocatable :: coorat(:,:,:)
      real*8 bvs(3,3),avs(3,3),b(3,3),tnp(48,3)
      integer maxnat
!  Variables for the symmetry information :
      real*8 point_gr(3,3,48)
      real*8,allocatable::car_rot(:,:,:),car_trans(:,:)
!  Unimportant locals :
      real*8 ascale,real_dum(3,3),xmult,ascinv,a1t(3),a2t(3),a3t(3)
	  real*8 alatt_old,alatt_new
      real*8,external :: dotthree
      integer i,j,k,ierr,np,mult(48,48),it,mtrx(48,3,3)
	  logical,parameter :: verbose=.false.
	  logical,parameter :: change_to_primitive_lattice=.true.
	  logical,parameter :: transform_coordinates=.false.


!********************* ORGANIZE CRYSTAL DATA *************************************

!  First, figure out if we are already in the primitive lattice.  If not, go there.
      a1t=a1;a2t=a2;a3t=a3
	  if (verbose) write(11,*) 'lattice type is ',lattice,'   ',latticename
      if(change_to_primitive_lattice) then
	      !The assumption is that the atom positions are given for the atoms belonging in the primitive cell but in the conventional cell coordinates.
		  !That is, to go to the primitive cell representation, one needs merely transform the coordinates, which is what is done here.
		  !Alternatively, one could find the "full" list of atom positions by applying "center" operations to the list of coordinates,
		  !E.g. for a CXY system where 40 positions are given, one can add 40 positions by doing "+(a,b,0)/2"; this fills the conventional cell.
		  
		  
		  !Correction!  It is necessary to transform the basis vectors (a1,a2,a3).
		  !However, the KKR routines expect the atom coordinates to be given in Carthesian coordinates, normalized to the length of the first lattice vector.
		  !So, the coordinates (ppos) shouldn't be changed, except for the normalization.
		  
		  
		  if(verbose) then
		     write(*,*) 'old basis'
		     write(*,*) 'a1',a1
		     write(*,*) 'a2',a2
		     write(*,*) 'a3',a3
             write(*,*) 'old positions'
		     do i=1,nats
		        write(*,*) ppos(:,i)
		     enddo
		  endif
		  alatt_old=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)  ! coordinates are expressed in units of |a1|.  However we are about to change |a1| ...
		  if(lattice.eq.'P') then   ! Primitive
		  elseif(lattice.eq.'F') then  ! Face centered
			 a1=(a2t+a3t)/dble(2)
			 a2=(a1t+a3t)/dble(2)
			 a3=(a1t+a2t)/dble(2)
			 if(transform_coordinates) then
			 do i=1,nats
				  a1t=ppos(:,i)
				  ppos(1,i)=-a1t(1)+a1t(2)+a1t(3)
				  ppos(2,i)= a1t(1)-a1t(2)+a1t(3)
				  ppos(3,i)= a1t(1)+a1t(2)-a1t(3)
			 enddo
			 endif
		  elseif(lattice.eq.'I'.or.lattice.eq.'B') then  ! Body centered
			   a1=(-a1t+a2t+a3t)/dble(2)
			   a2=( a1t-a2t+a3t)/dble(2)
			   a3=( a1t+a2t-a3t)/dble(2)
			   if(transform_coordinates) then
			   do i=1,nats
				  a1t=ppos(:,i)
				  ppos(1,i)=a1t(2)+a1t(3)
				  ppos(2,i)=a1t(1)+a1t(3)
				  ppos(3,i)=a1t(1)+a1t(2)
			   enddo
               endif
		  elseif(lattice.eq.'C') then ! xz Base centered (or yz?)
			   if(latticename.eq.'CXY') then
				  a1=( a1t+a2t)/dble(2)
				  a2=(-a1t+a2t)/dble(2)
				  a3=a3t 
			      if(transform_coordinates) then
				  do i=1,nats
					 a1t=ppos(:,i)
					 ppos(1,i)= a1t(1)+a1t(2)
					 ppos(2,i)= a1t(2)-a1t(1)
					 ppos(3,i)= a1t(3)
				  enddo
                  endif
			   elseif(latticename.eq.'CXZ') then
				  a1=( a1t+a3t)/dble(2)
				  a3=(-a1t+a3t)/dble(2)
				  a2=a2t 
			      if(transform_coordinates) then
				  do i=1,nats
					 a1t=ppos(:,i)
					 ppos(1,i)= a1t(1)+a1t(3)
					 ppos(3,i)= a1t(3)-a1t(1)
					 ppos(2,i)= a1t(2)
				  enddo
                  endif
			   elseif(latticename.eq.'CYZ') then
				  a3=( a3t+a2t)/dble(2)
				  a2=(-a3t+a2t)/dble(2)
				  a1=a1t 
			      if(transform_coordinates) then
				  do i=1,nats
					 a1t=ppos(:,i)
					 ppos(3,i)= a1t(3)+a1t(2)
					 ppos(2,i)= a1t(2)-a1t(3)
					 ppos(1,i)= a1t(1)
				  enddo		   
                  endif
			   endif
			   call wlog('crystalstructure - warning: C lattice not well tested')
		  elseif(lattice.eq.'R') then ! Rhombohedral ; the conventional cell is hexagonal
			   a1=(2.d0*a1t     -a2t  +a3t)/dble(3)
			   a2=(    -a1t+2.d0*a2t  +a3t)/dble(3)
			   a3=(    -a1t     -a2t  +a3t)/dble(3)
			   if(transform_coordinates) then
			   do i=1,nats
				  a1t=ppos(:,i)
				  ppos(1,i)=a1t(1)+a1t(3)
				  ppos(2,i)=a1t(2)+a1t(3)
				  ppos(3,i)=a1t(3)-a1t(1)-a1t(2)
			   enddo
               endif
			   call wlog('crystalstructure - warning: R lattice not well tested')
		  elseif(lattice.eq.'H') then ! Hexagonal
			 !Do nothing; a hexagonal lattice is primitive
		  else
			   stop 'unknown lattice type in crystalstructure'
		  endif
		  !At this point, I think we need to tell the rest of the program that we are in a primitive cell now:
		  lattice='P'
		  latticename='P  '
		  alatt_new=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
		  ppos=ppos*alatt_old/alatt_new  !Renormalize so coordinates are now in units of the new lattice vector
		  if(verbose) then
		     write(*,*) 'new basis'
		     write(*,*) 'a1',a1
		     write(*,*) 'a2',a2
		     write(*,*) 'a3',a3
             write(*,*) 'new positions'
		     do i=1,nats
		        write(*,*) ppos(:,i)
		     enddo		  
		  endif
	  endif

        a1t=dble(0);a2t=dble(0);a3t=dble(0)  ! destroy temporary variables(a1t is nonsense by this time anyway!)
!  We have assumed that only the primitive atoms were given!  Eg., only 2 atoms specified for diamond input with
!  cubic lattice and 'F' option.  If not, garbage will be produced.


!  input the crystal info
!  First, figure out how many atoms there are in the unit cell :
      allocate(coorat(nph,nats,3))
      coorat=dble(0)
      natom=0
      j=0 !test variable
      do it=1,nph
        do i=1,nats
         if(ppot(i).eq.it) then
              natom(it)=natom(it)+1
              coorat(it,natom(it),:)=ppos(:,i)
           endif
        enddo
        j=j+natom(it)
        if(natom(it).eq.0) stop 'error input crystalstructure'
      enddo
      if(j.ne.nats) stop 'error input crystalstructure'

      maxnat = 0
      do i = 1,nph
        if (natom(i).gt.maxnat) maxnat = natom(i)
      end do

      call cross(a2,a3,b1)
      call cross(a3,a1,b2)
      call cross(a1,a2,b3)
      celvol=dotthree(b1,a1)
      xmult=pi2/celvol
      do i=1,3
        b1(i)=b1(i)*xmult
        b2(i)=b2(i)*xmult
        b3(i)=b3(i)*xmult
      end do
      celvol=dabs(celvol)
      if (celvin.gt.0.d0) then
	    if(verbose) write(11,*) 'cell volume found ',celvol,' ; cell volume requested ',celvin
        ascale = (celvin/celvol)**(1.d0/3.d0)
        ascinv = 1.d0/ascale
        a1=a1*ascale
        a2=a2*ascale
        a3=a3*ascale
        b1=b1*ascinv
        b2=b2*ascinv
        b3=b3*ascinv
        celvol = celvin
      endif

        alat(1)=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2)
        alat(2)=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2)
        alat(3)=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2)
        alfalat(1)=dacos(dotthree(a2,a3)/(alat(2)*alat(3)))
        alfalat(2)=dacos(dotthree(a1,a3)/(alat(1)*alat(3)))
        alfalat(3)=dacos(dotthree(a1,a2)/(alat(1)*alat(2)))

!  Done reading input from file.

!    Now, summarize the structure to stdout :
      if(verbose)then
		  do i=1,nph
			write(11,'(/,1x,i2,5x,3x,1a12,/)')  natom(i),'in positions'
			do j=1,natom(i)
			  write(11,'(1x,3f10.5)') (coorat(i,j,k),k=1,3)
			end do
		  end do

		  write(11,'(/,1a23,/,3(/,1x,3f10.5))') ' lattice vectors (a.u.)',a1,a2,a3
		  write(11,'(/,1a23,/,3(/,1x,3f10.5))') ' in wave number space  ',b1,b2,b3

		  write(11,'(/,1a14,1f10.4)') ' cell volume =',celvol
      endif
	  
!    Save the individual (reciprocal) basis in a 3*3 matrix avs (bvs)
      do j=1,3
        avs(1,j)=a1(j)
        avs(2,j)=a2(j)
        avs(3,j)=a3(j)
        bvs(1,j)=b1(j)
        bvs(2,j)=b2(j)
        bvs(3,j)=b3(j)
      end do

!    Calculate the 'metric' of the reciprocal basis in b
      b=dble(0)
      do i=1,3
        do j=1,3
          do k=1,3
            b(i,j)=b(i,j)+bvs(i,k)*bvs(j,k)
          end do
        end do
      end do

!    Is the lattice orthogonal?
      ascale=dabs(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))   &! a1 . a2
     &   +dabs(a1(1)*a3(1)+a1(2)*a3(2)+a1(3)*a3(3))   &! a1 . a3
     &   +dabs(a3(1)*a2(1)+a3(2)*a2(2)+a3(3)*a2(3))   ! a2 . a3
      ortho=(dabs(ascale).lt.dble(0.00000001))
      if(verbose) write(11,*) 'lattice orthogonality (T/F) :',ortho

!******************************* DETERMINE SYMMETRY OF THE CRYSTAL *********************


!    Determine the point group of the lattice
      call pointgroup(bvs,b,48,point_gr,np)
      if(verbose) write(11,*) "There are ",np," operations in the point group"

!    Now figure out the space group of lattice AND basis
      call spacegroup(avs,bvs,nph,maxnat,coorat,natom,48,point_gr,np,cryst_gr(:,:,:,2),ncryst)
      if(verbose) write(11,*) "There are ",ncryst," sym. operations"

!******************************* GIVE OUTPUT *******************************************

      if(verbose) then
        open(71,file='avsbvs.txt')
        do i=1,3
        write(71,'(3f12.4,3x,3f12.4)') avs(i,:),bvs(i,:)
        enddo
        close(71)
      endif

!    Calculate the integer representation of the symmetry matrices in mtrx, tnp and 
!    write it to file 23 if iprint says so.  
      if (verbose) open(23,file='operations.txt',form='formatted')
      do i=1,ncryst
         do j=1,3
            do k=1,3
               mtrx(i,j,k)=nint(cryst_gr(j,k,i,2))
            end do
            tnp(i,j)=cryst_gr(j,4,i,1)
         end do
         if (verbose) then
            write(23,*) "Oper.: ",i 
            write(23,'(3(/,2x,3i10))') ((mtrx(i,j,k),k=1,3),j=1,3)
            write(23,'(7x,1a1,3f10.3,1a2)')  '(',(tnp(i,j),j=1,3),' )'
            write(23,*) ""
         endif
         do j=1,3
            tnp(i,j) = pi2*tnp(i,j)
         end do
      end do
      if (verbose) close(23)


      if(verbose) then
         write(11,*) 'mtrx :'
         do j=1,3
            write(11,*) mtrx(8,j,1:3)
         enddo

         write(11,*) 'avx:'
         do j=1,3
            write(11,*) avs(j,1:3)
         enddo
         write(11,*) 'bvx:'
         do j=1,3
            write(11,*) bvs(j,1:3)
         enddo
		endif


!    Calculate the Carthesian representation of the symmetry matrices in car_rot, car_trans and print it to stdout.      
      allocate(car_rot(1:3,1:3,1:ncryst),car_trans(1:3,1:ncryst))
      if(verbose) write(11,*) "Cartesian Operators"
      car_trans=dble(0)
      do i=1,ncryst
         call change_car(bvs,avs,mtrx,i,real_dum)
           car_rot(:,:,i)=real_dum/pi2
         do j=1,3
            do k=1,3
               car_trans(j,i)=car_trans(j,i)+avs(k,j)*tnp(i,k)
            end do
         end do
		 car_trans(:,i)=car_trans(:,i)/pi2
		 if(verbose) then
			 write(11,*) i
			 write(11,'(3(/,1x,3f10.4))') ((car_rot(j,k,i),k=1,3),j=1,3)
			 write(11,*)
			 write(11,'(3f10.4)') (car_trans(j,i),j=1,3)
			 write(11,*)
		 endif
      end do


!KJ  I believe that the rest of my program does not care about the translations ; only about the rotations.
!    However, I still return the whole thing in WIEN2k-compatible format
      do i=1,ncryst
           cryst_gr(:,1:3,i,1)=car_rot(:,:,i)
      enddo
!    I leave the translations in fractional units OF THE PRIMITIVE LATTICE.  (Eg., if input is given as F,
!    specifying the cubic (conventional) lattice vectors, the translation vectors will be in fractions of
!    the rhombohedric (primitive) lattice vectors.

      deallocate(car_rot,car_trans)


      call symmetrycheck(ierr,ncryst,mtrx,tnp,mult)


	if (ierr .ne. 0) then
		call wlog('FEFF cannot use symmetry.  The calculation proceeds with 1 symmetry operation.')
		np=1
		ncryst=1
		tnp=dble(0)
		mtrx=0
		mtrx(1,1,1)=1
		mtrx(1,2,2)=1
		mtrx(1,3,3)=1
		cryst_gr=dble(0)
		cryst_gr(:,1:3,1,1)=mtrx(1,:,:)
		call change_car(bvs,avs,mtrx,1,real_dum)
		cryst_gr(:,1:3,1,2)=real_dum/pi2
	endif

 
      return
        end

