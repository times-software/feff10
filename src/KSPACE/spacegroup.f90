!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: spacegroup.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     spacegroup attempts to find the improper rotations of the crystal and creates the symmetry group

!KJ Important note : I am not confident that the translations come out properly!!

!
!     input:
!     bvs  reciprocal lattice vectors
!     avs  real space lattice vectors
!     b    'metric' of the bs
!     types number of types of atoms
!     maxnat maximum number of types of atoms
!     coorat coordinates of the basis
!     natom(i) number of atoms of type 'i'
!     ntot max number of point groups
!     point_gr point group of the crystal
!     npoint number of point group oper
!
!     output:
!     cryst_gr crystal group (rot + rot*translation)


      subroutine  spacegroup(avs,bvs,types,maxnat,coorat_a,             &
     &     natom,ntot,point_gr,npoint,spacegr,nspacegr)

      use constants,only : pi2
      implicit none

      integer,intent(in) :: types,maxnat,ntot,npoint
      real*8,intent(in) :: bvs(3,3),avs(3,3),                           &
     &                coorat_a(types,maxnat,3),point_gr(3,3,ntot)
      integer,intent(in) :: natom(types)
      integer,intent(out) :: nspacegr
      real*8,intent(out) :: spacegr(3,4,ntot)

      real*8 c(types,maxnat,3)
      real*8 da,ident(3,3),vdum,v1(3),r_red(3),rr(3)
      real*8, allocatable:: difs(:,:,:),co(:,:,:),rotc(:,:,:)
      integer mxnatom,it,iatom,ipoint,ntrans,nats,ncells,cell1,nn
      integer i,j,k,l,m,p,ifirst,expansion,ix,iy,iz

      real*8,parameter :: eps=1.0E-5

      logical used(ntot),debug
        logical,allocatable :: mapped(:,:),istrans(:,:,:),mappedsmall(:)
        real*8 dist
        real*8, allocatable :: trans(:,:,:,:)



      ntrans=0
        ident=dble(0)
        ident(1,1)=dble(1)
        ident(2,2)=dble(1)
        ident(3,3)=dble(1)


!KJ The coordinates are already carthesian!  Therefore, I'm skipping the next section, and instead just doing :
      c=coorat_a*dsqrt(avs(1,1)**2+avs(1,2)**2+avs(1,3)**2)
!c Change the basis coordinates to cartesian
!      c=dble(0)
!      do it=1,types
!         do iatom=1,natom(it)
!            do i=1,3
!               do j=1,3
!                  c(it,iatom,i)=c(it,iatom,i)+
!     &                 coorat_a(it,iatom,j)*avs(j,i)
!               end do
!            end do
!         end do
!      end do
!	c=c/pi2




      mxnatom=0
        nats=0
      do iatom=1,types
           nats=nats+natom(iatom)
         if (natom(iatom).gt.mxnatom) mxnatom=natom(iatom)
      end do

      if(maxnat.ne.mxnatom) stop 'atoms wrong'


      expansion=0

        i=2*expansion+1
        cell1=expansion*(i**2+i+1)
        ncells=i**3
      p= ncells * mxnatom  ! upper estimate for total number of atoms in this routine
      allocate (co(3,p,types))
      co=dble(0)

      do it=1,types
        j=0
        do ix=-expansion,expansion
        do iy=-expansion,expansion
        do iz=-expansion,expansion
        do i=1,natom(it)
         j=j+1
           r_red=dble(0)
         do k=1,3
            do l=1,3
               r_red(k)=r_red(k)+bvs(k,l)*c(it,i,l)
            end do
            r_red(k)=r_red(k)/pi2
         end do
           r_red(1)=r_red(1)+ix
           r_red(2)=r_red(2)+iy
           r_red(3)=r_red(3)+iz
! back to normal coordinates
         rr=dble(0)
         do k=1,3
            do l=1,3
               rr(k)=rr(k)+avs(k,l)*r_red(l)
            end do
         end do
         co(:,j,it)=rr
      enddo
        enddo
        enddo
      enddo
        enddo

! Rotating the basis coordinates of each type of atom and then
! finding out the differences in the position of same type 
! of atoms in rotated and original basis. This done for each element
! of the point group.
      
        allocate (difs(3,p,p),rotc(types,p,3),mapped(p,p),mappedsmall(p))
      ntrans=natom(1)*ncells
        allocate (trans(3,ntrans,ntrans,ntot),istrans(ntot,ntrans,ntrans))
        difs=dble(0)
        rotc=dble(0)
        istrans=.true.
        trans=dble(0)

      do ipoint=1,npoint  ! loop over all point group sym ops


         do it=1,types
            do iatom=1,natom(it)*ncells
                 v1=dble(0)
               do i=1,3
                  do j=1,3
                     v1(i)=v1(i)+point_gr(i,j,ipoint)*co(j,iatom,it)
                  end do
               end do
                 rotc(it,iatom,:)=v1
            end do
         end do

         do it=1,types

            nn=ncells*natom(it)
            difs=dble(0)
              do i=1,nn
              do j=1,nn
                 difs(:,i,j)=rotc(it,i,:)-co(:,j,it)
                 do k=1,3
!	           if(dabs(difs(k,i,j)).lt.0.00000001) difs(k,i,j)=dble(0)
                 enddo
              enddo
              enddo
              if(it.eq.1) trans(:,:,:,ipoint)=-difs

      debug=.false.
        if(debug) then
      open(11,file='debug.txt')
        write(11,*) 'symmetry operation'
        do i=1,3
        write(11,17) ipoint,point_gr(i,:,ipoint)
        enddo
        write(11,*) 'translation vectors'
        do i=1,ntrans
        do j=1,ntrans
        write(11,17) i,trans(:,i,j,ipoint)
        enddo
        enddo
        write(11,*) 'original coordinates'
        do i=1,types
        do j=1,natom(i)
        write(11,17) i,c(i,j,:)
        enddo
        enddo
        write(11,*) 'full set of coordinates'
        do i=1,types
        do j=1,ncells*natom(it)
        write(11,17) i,co(:,j,i)
        enddo
        enddo
        write(11,*) 'rotated coordinates'
        do i=1,types
        do j=1,ncells*natom(it)
        write(11,17) i,rotc(i,j,:)
        enddo
        enddo
17    format(i3,3f14.4)
      write(11,*) 'difs'
      do j=1,ntrans
        do i=1,natom(1)*ncells
         write(11,17) j,difs(:,j,i)
        enddo
        enddo
!	close(11)
 !     stop
        endif


            do i=1,ntrans
              do j=1,ntrans
                 if(istrans(ipoint,i,j)) then   ! don't bother if the translation's already been proven useless
                    mapped=.false.
                    do l=1,nn
                    do m=1,nn
                       dist=dble(0)
!	               debug=(i.eq.2.and.j.eq.2.and.l.eq.1.and.m.eq.1)
                       call reduce(avs,bvs,difs(:,l,m)+trans(:,i,j,ipoint)&
     &                               ,v1,p)
                       do k=1,3
                          dist=dabs(v1(k))+dist
                       enddo
                       if(dist.lt.eps.and.mapped(l,m)) then
        			       stop 'error overlap'
        			   elseif(dist.lt.eps) then
        			       mapped(l,m)=.true.
                       endif
!	 write(11,*) 'i,j,l,m',i,j,l,m
!	write(11,*) 'difslm,transij',difs(:,l,m),trans(:,i,j,ipoint)
!	write(11,*) 'v1',v1
                    enddo
                    enddo
                    mappedsmall=.false.
                    do l=1,nn
                    do m=1,nn
                       mappedsmall(l)=mappedsmall(l).or.mapped(l,m)
                    enddo
                    enddo
                    do l=1,nn
                       istrans(ipoint,i,j)=(istrans(ipoint,i,j).and.    &
     &                                               mappedsmall(l))
                    enddo
!	            if(i.eq.2.and.j.eq.2) then
!	               write(11,*) 'mapped',mapped
!	               write(11,*) 'mappedsmall',mappedsmall
!	               write(11,*) 'istransij',istrans(ipoint,i,j)
!	            endif
                 endif
              enddo
              enddo
         enddo  ! atom types it=1,types
!      close(11)!;stop

        enddo  ! point group sym ops


! Now determine the space group ...

      used=.false.
      nspacegr=0
        spacegr=dble(0)
        do i=1,npoint
         do j=1,ntrans
           do it=1,ntrans
              if(istrans(i,j,it).and.(.not.used(i))) then
                 used(i)=.true.
                 nspacegr=nspacegr+1
               do p=1,3
               do k=1,3
                  do l=1,3
                  do m=1,3
                        spacegr(p,k,nspacegr)=spacegr(p,k,nspacegr)+    &
     &                    avs(p,l)*point_gr(l,m,i)*bvs(k,m)
                  end do
                  end do
                  spacegr(p,k,nspacegr)=spacegr(p,k,nspacegr)/pi2
                 enddo
                 enddo
                 call subtract_a(bvs,trans(:,j,it,i),                   &
     &                           spacegr(:,4,nspacegr),p)
!	         spacegr(:,4,nspacegr)=trans(:,j,it,i)
              endif
           enddo
           enddo
        enddo




      if (nspacegr.lt.1) stop "Did not find any symmetry operations!!!"

     
!     Finally shift the identity operator to the top.
!     This will be better for later purposes.
     
      ifirst=0
      do i=1,nspacegr
         da=0.0d0
         do j=1,3
            da=da+abs(spacegr(j,4,i))
         end do
         if(da.le.eps) then
            da=0.0d0
            do j=1,3
              do k=1,3
               da=da+abs(spacegr(k,j,i)-ident(k,j))
              enddo
            end do
            if (da.le.eps) then
               if (ifirst.gt.0) then
                  write(6,*) "Found more than one identity operator"
                  stop
               end if
               ifirst=i
            end if
         end if
      end do

      if (ifirst.lt.1) stop "Did not find the identity operator"

      do i=1,3
        do j=1,4
         vdum=spacegr(i,j,1)
         spacegr(i,j,1)=spacegr(i,j,ifirst)
         spacegr(i,j,ifirst)=vdum
        enddo
      enddo

      return 
      end
      







