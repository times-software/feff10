!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: pointgroup.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/01/30 06:01:58 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     pointG_find attempts to find the point group of the crystal
!     input:
!     bvs  reciprocal lattice vectors
!     b    'metric' of the bs
!     ntot max number of point groups
!
!     output:
!     point_gr point group of the crystal
!     noper    number of operations in the point group

      subroutine pointgroup(bvs,b,ntot,point_gr,nopr)
      implicit none
      integer,intent(in)  :: ntot
      integer,intent(out) :: nopr
      real*8,intent(in)   :: bvs(3,3),b(3,3)
      real*8,intent(out)  :: point_gr(3,3,ntot)
!     Internal parameters
      real*8,parameter :: eps=1.0E-8
      real*8 tbvs(3,3),inv_tbvs(3,3),dum_tbvs(3,3)
      real*8,allocatable::Gdbl(:,:),Gnorm(:)
      real*8 bhelp,dum,diag(3),g1(3),g2(3),dot1,dot2,maxB
      integer i,j,k,imax,jmax,kmax,ig1,ig2,ig3,icar,nGmax

!KJ     calculate inv_tbvs := bvs^-1     :
      do i=1,3
         do j=1,3
            tbvs(i,j)=bvs(j,i)
         end do
      end do
      call invertmatrix(3,3,tbvs,dum_tbvs,inv_tbvs)

!KJ  Prepare some variables for the upcoming test :
      maxB=dble(0)
!      write(6,*) "b"
      do i=1,3
         if (b(i,i).gt.maxB) maxB=b(i,i)
         diag(i)=dble(1)/b(i,i)
!         write(6,*) (b(i,j),j=1,3)
      end do
!      write(6,*)
      if (maxB.le.eps) stop "b is zero in pointG_find"
      bhelp=b(1,2)*b(1,3)*b(2,3)


!KJ  To determine the point group, a number of alternative lattices will be set up and
!KJ  compared to the input lattice  (all in reciprocal space).  Which alternative lattices
!KJ  make sense?  Basically, all those consisting of vectors small enough that they could be
!KJ  rotations of the input lattice.  Here we determine the corresponding cutoff in three steps.
!KJ  (I'm not sure what exactly happens here.)

      dum=b(1,1)-b(1,2)*b(1,2)*diag(2)- b(3,1)*b(3,1)*diag(3)+ 2.0d0*bhelp*diag(2)*diag(3)
      if (dum.le.eps) stop "dum(1) is zero in pointG_find"
      imax=int(sqrt(maxB/dum)+1.0d0)
!      write(6,*) maxB,dum,dble(imax)*dble(imax)*dum

      dum=b(2,2)-b(1,2)*b(1,2)*diag(1)- b(3,2)*b(3,2)*diag(3)+ 2.0d0*bhelp*diag(1)*diag(3)
      if (dum.le.eps)  stop "dum(2) is zero in pointG_find"
      jmax=int(sqrt(maxB/dum)+1.0d0)
!      write(6,*) maxB,dum,dble(jmax)*dble(jmax)*dum

      dum=b(3,3)-b(1,3)*b(1,3)*diag(1)- b(3,2)*b(3,2)*diag(2)+ 2.0d0*bhelp*diag(1)*diag(2)
      if (dum.le. eps) stop "dum(3) is zero in pointG_find"
      kmax=int(sqrt(maxB/dum)+1.0d0)
!      write(6,*) maxB,dum,dble(kmax)*dble(kmax)*dum

      nGmax=(2*imax+1)*(2*jmax+1)*(2*kmax+1)


!KJ In the following block, a number (nGmax) of reciprocal lattice vectors are saved in
!KJ Gdbl, along with their norm^2 in Gnorm.
      allocate(Gdbl(3,nGmax),Gnorm(nGmax))
      ig1=0
      do i=-imax,imax
         do icar=1,3
           g1(icar)=dble(i)*tbvs(icar,1) 
         end do
         do j=-jmax,jmax
            do icar=1,3
               g2(icar)=g1(icar)+ dble(j)*tbvs(icar,2) 
            end do
            do k=-kmax,kmax
               ig1=ig1+1

               do icar=1,3
                  Gdbl(icar,ig1)= g2(icar)+dble(k)*tbvs(icar,3)
               end do
               dum=0.0d0
               do icar=1,3
                  dum=dum+Gdbl(icar,ig1)*Gdbl(icar,ig1)
               end do
               Gnorm(ig1)=dum
            end do
         end do
      end do

      nopr=0

!     Triple loop over the G vectors. Testing for each set 
!     of three DIFFERENT G vectors if they are some rotation 
!     of the b's (bvs). This done by comparing the 'metric' of
!     of the two sets of three vectors. If match is found the 
!     rotation matrix associated with group is stored.
 

      do ig1=1,nGmax

         if (ABS((Gnorm(ig1)-b(1,1))).le.eps) then  
         do ig2=1,nGmax

            if (ig2.ne.ig1 .and. ABS(Gnorm(ig2)-b(2,2)).le.eps) then
               dot1=dble(0)
               do icar=1,3
                  dot1=dot1+Gdbl(icar,ig1)*Gdbl(icar,ig2)
               end do
               if (ABS(dot1-b(2,1)).le.eps) then

               do ig3=1,nGmax

                  if (ig3.ne.ig2 .and. ig3.ne.ig1 .and.  ABS(Gnorm(ig3)-b(3,3)).le.eps) then
                     dot1=dble(0)
                     dot2=dble(0)
                     do icar=1,3
                        dot1=dot1+Gdbl(icar,ig1)*Gdbl(icar,ig3)
                        dot2=dot2+Gdbl(icar,ig2)*Gdbl(icar,ig3)
                     end do

                     
                     if (ABS(dot1-b(1,3)).le.eps  .AND. ABS(dot2-b(2,3)).le.eps) then
                     nopr=nopr+1
                     if (nopr.gt.ntot) then
                        write(6,*) "There are over ",ntot," rotations in the point group. Something must be wrong."
                        STOP
                     end if
                     call add_oper(inv_tbvs,Gdbl,nGmax,ig1,ig2,ig3,nopr,point_gr,ntot)
                     end if
                  end if
               end do
               end if
            end if
         end do
         end if
      end do
!     By this time the point group of the lattice should be in hand.
!     If not error is raised !!!
      if (nopr.lt.1) then
         write(6,*) "There is no point group for this system.  Something must have gone wrong!!"
         STOP
      end if
      deallocate(Gnorm,Gdbl)
      return
      end

!  Subroutine add_oper adds a rotation operation matrix 
!  to the collection point_gr (total size ntot, current size npoint-1).
!  This done by using the inverse of the b's (inv_tbvs).
!  The Gvectors index ig1,ig2,ig3 (from the matrix Gdbl) are used in  
!  process.

!KJ  We seek the matrix m that transforms (b1 b2 b3) into (Gdble(ig1) Gdble(ig2) Gdble(ig3)) =: g
!KJ  This is done by multiplying    m := point_gr(:,npoint) = g * (b1 b2 b3)^-1   (=: inv_tbvs)
!KJ  Finally, for some reason, m is saved as a 9 element array instead of a 3*3 matrix.  (Silly people ...)

      subroutine add_oper(inv_tbvs,Gdbl,nGmax,ig1,ig2,ig3,npoint,point_gr,ntot)
      implicit none
      integer ntot,npoint,ig1,ig2,ig3,nGmax
      real*8 inv_tbvs(3,3)
      real*8 Gdbl(3,nGmax)
      real*8 point_gr(3,3,ntot)
      
      real*8 g(1:3,1:3)
      integer i,j,icar
      
      do icar=1,3
         g(icar,1)=Gdbl(icar,ig1)
         g(icar,2)=Gdbl(icar,ig2)
         g(icar,3)=Gdbl(icar,ig3)
      end do
      
        point_gr(:,:,npoint)=dble(0)
      do i=1,3
         do j=1,3
            do icar=1,3
               point_gr(i,j,npoint)=point_gr(i,j,npoint)+ g(i,icar)*inv_tbvs(icar,j)
            end do
         end do
      end do

      return 
      end 


