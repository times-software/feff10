!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmskspace.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2012/03/16 22:55:02 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fmsband(ispin,xphase,em,eigenvalues,bk,freeprop)  ! KJ
  use DimsMod, only: lx, nphx=>nphu, nspx=>nspu
  use struct
  use strfacs, only : eimag  !,gkk=>gk
  use boundaries
  use wigner3j
  use trafo
  use controls, only : corehole,cholestrength,sprkkrpot,fullpot
  USE IOMod
  use constants,only : bohr 

  implicit none

 
  !  INPUT
  integer,intent(in) ::    ispin
  complex,intent(in) ::    xphase(nspx, -lx:lx, 0:nphx) ! array of partial wave phase shifts for each unique potential
  real*8,intent(in)  ::    bk(3)  ! k-point ( single point)
  logical,intent(in) ::    freeprop
  complex*16,intent(in) :: em !energy at which we work - one point, not a mesh
 !  OUTPUT
  complex,intent(out) ::   eigenvalues(msize)
   !********************************************************************
  !  LOCAL
  real*8,parameter :: one=1_8
  real*8,parameter :: pi = 3.1415926535897932384626433d0
  complex,parameter :: coni=(0.0,1.0)
  complex,parameter :: c0=(0.0,0.0),c1=(1.0,0.0)
  complex Gfms(msize,msize),Gtemp(msize,msize)  ! the fms-matrix of the Green's function
  complex tmatrx(nsp,msize+mls)   ! scattering matrix; defined otherwise than in real space calculation
  complex Tmat(msize,msize)   ! scattering matrix of the perfect lattice (without corehole)
  integer ist1,ist2,iat1,l1,m1,iat2,m2,isp1,iph,i,j,a,b,isp2,n,ik,k
  integer iatom,nt,lun,j1,jj1,j2,jj2
  complex*16 eryd,p
  complex prefac,aa
  COMPLEX*16, parameter :: CI2PI = (0.d0,2.d0)*PI
  integer irels,nts,zs(nats)
  integer ia,ib,na,nb,niphabs,info
  complex,allocatable :: work(:)
  real,allocatable :: rwork(:)
  logical,allocatable :: lwork(:)
  logical,parameter :: debug=.false.

  character*7 fname
  character*75 messg

  INTEGER istate





  if(freeprop) goto 12345
  !********************* TAKE CARE OF THE T-MATRIX :
  ! Make tmatrx
  ist1=0
  IAT_LOOP: do iat1=1,nats+1  ! last field is for core hole potential
     if(iat1.le.nats) then
        iph=ppot(iat1)
     else
        iph=0
     endif
     L1_LOOP: do l1=0,maxl
   !all potentials have same l now - replace by lpot(iph) - fix later
        M1_LOOP: do m1=-l1,l1
           ISP_LOOP: do isp1=1,nsp
              ist1=ist1+1
              if(nsp.eq.1.and.ispin.eq.0) then
                 tmatrx(1,ist1)= (exp(2*coni*xphase(isp1,l1,iph))-one)  &
                      &			                             /(2*coni)
              elseif(nsp.eq.1.and.ispin.gt.0) then
                 tmatrx(1,ist1)=  ( exp(2*coni*xphase(isp1,l1,iph)) - one)&
                      &  / (2*coni) * t3jm (l1, m1, 2)**2  +           &
                      &  ( exp(2*coni*xphase(isp1,-l1,iph)) - one )    &
                      &  / (2*coni) * t3jp (l1, m1, 2)**2 
              elseif(nsp.eq.1.and.ispin.lt.0) then
                 tmatrx(1,ist1)=  ( exp(2*coni*xphase(isp1,l1,iph)) - one)&
                      &  / (2*coni) * t3jm (l1, m1, 1)**2  +           &
                      &  ( exp(2*coni*xphase(isp1,-l1,iph)) - one )    &
                      &  / (2*coni) * t3jp (l1, m1, 1)**2 
              else  ! nsp=2

                 tmatrx(1,ist1)=  ( exp(2*coni*xphase(isp1,l1,iph)) - one) &! for l1=l2,m1=m2,isp1=isp2 diagonal
                      &  / (2*coni) * t3jm (l1, m1, isp1)**2  +        &
                      &  ( exp(2*coni*xphase(isp1,-l1,iph)) - one )    &
                      &  / (2*coni) * t3jp (l1, m1, isp1)**2 
                 if(isp1.eq.1) then
                    isp2=2
                 else
                    isp2=1
                 endif
                 m2=m1+isp1-isp2

                 tmatrx(2,ist1)= ( exp(2*coni*xphase(isp1, l1,iph)) - one  &! for l1=l2,m1+isp1=m2+isp2 off-diagonal
                      &  + exp(2*coni*xphase(isp2,l1,iph)) - one ) / (4*coni) &
                      &  * t3jm (l1, m1, isp1) * t3jm (l1, m2, isp2)  +       &
                      &  ( exp(2*coni*xphase(isp1,-l1,iph)) - one +           &
                      &   exp(2*coni*xphase(isp2,-l1,iph)) - one ) / (4*coni) &
                      &   * t3jp (l1, m1, isp1) * t3jp (l1, m2, isp2)

              endif  ! which case
           enddo ISP_LOOP
        enddo M1_LOOP
     enddo L1_LOOP
  enddo IAT_LOOP

  !  Blow up t-matrix to T-matrix.
  Tmat=cmplx(0,0)
  ist1=0
  do iat1=1,nats
     do l1=0,maxl
        do m1=-l1,l1
           do isp1=1,nsp
              ist1=ist1+1
              Tmat(ist1,ist1)=tmatrx(1,ist1)
              if(nsp.eq.2) then
                 if(isp1.eq.1) then
                    isp2=2
                    ist2=ist1-1
                 else
                    isp2=1
                    ist2=ist1+1
                 endif
                 if(ist2.le.msize.and.ist2.gt.0) Tmat(ist1,ist2)=tmatrx(2,ist1)
              endif  ! nsp=2
           enddo
        enddo
     enddo
  enddo
  
  

12345 continue
  !********************************** TAKE CARE OF THE STRUCTURE FACTORS :

  !    Calculate structure factors G(k), returned in matrix gk

  eryd = dble(2)*em
  eryd = dble(2)*dble(em)
  eimag=0.03d0
  !if (dabs(dimag(eryd)).lt.2*eimag) eryd=dcmplx(dble(eryd),dimag(eryd)+eimag)
  p = sqrt(eryd)

  ! ------------------ calculate energy - dependent terms of str.constants
  call strcc(eryd,alat(1),.false.)

  ! for SPRKKR units, multiply by p rather than dimensionless:
  !Tmat=Tmat*p

  !************************  NOW CALCULATE GG USING STRUCTURE FACTORS AND T-MATRIX

  call kkrband(Gfms,Tmat,p,bk,freeprop,eryd)


!        Diagonalize Gfms by finding its eigenvalues : 
         ! Multiply by p to go to SPRKKR units
		 ! This matters at negative energies ...
         Gtemp=Gfms * cmplx(p) 
	     allocate(work(10*msize),rwork(msize),lwork(msize))
         call cgees('N','N',.false.,msize,Gtemp,msize,i,eigenvalues,Tmat,msize,work,10*msize,rwork,lwork,info)
!        i,Tmat  : dummy arguments for options used - will not be overwritten ??
         if(info.ne.0)  stop 'info crash in cgees yikes'

         ! order eigenvalues large to small by real part
         do i=1,msize
	     do j=i,msize
	        if(dble(eigenvalues(j)).gt.dble(eigenvalues(i)))then
	           aa=eigenvalues(j)
	           eigenvalues(j)=eigenvalues(i)
	           eigenvalues(i)=aa
	        endif
	     enddo
	     enddo

		           if(debug) write(45,440) eryd,dble(eigenvalues(1:msize))
440      format(100(f11.5,2x))
  return
end subroutine fmsband   

