!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmskspace.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2012/03/16 22:55:02 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fmskspace(ispin,ck,xphase,ie,em,gg,iverb,sig2g)  ! KJ
  use DimsMod, only: lx, nspx=>nspu, nphx=>nphu
  use struct
  use kklist
  use strfacs, only : eimag  !,gkk=>gk
  use boundaries
  use wigner3j
  use trafo
  use kgenwork, only : bkf,wf,nkf,lsymf,linkf,nki,nka
  use controls, only : corehole,cholestrength,sprkkrpot,fullpot
  USE IOMod
  use constants,only : bohr 

  implicit none

  !  OUTPUT
  ! the 'scattering path' matrix  (cf. tau_ij in Beeby)
  complex, intent(out) ::  gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)       

  !  INPUT
  integer,intent(in) ::    ispin, ie, iverb
  complex,intent(in) ::    xphase(nspx, -lx:lx, 0:nphx) ! array of partial wave phase shifts for each unique potential
  complex,intent(in) ::    ck(nspx) ! momentum of current energy point
  real*8,intent(in)  ::    sig2g  ! constant Debye-Waller factor
  complex*16 :: em !energy at which we work - one point, not a mesh
  !********************************************************************
  !  LOCAL
  real*8,parameter :: one=1_8
  real*8,parameter :: pi = 3.1415926535897932384626433d0
  complex,parameter :: coni=(0.0,1.0)
  complex,parameter :: c0=(0.0,0.0),c1=(1.0,0.0)
  complex Gfms(msize,msize)  ! the fms-matrix of the Green's function
  complex tmatrx(nsp,msize+mls)   ! scattering matrix; defined otherwise than in real space calculation
  complex Tmat(msize,msize)   ! scattering matrix of the perfect lattice (without corehole)
  complex Cmat(mls,mls)       ! core hole matrix (Cmat = t_gs - t_ch )
  complex r1(mls,mls),rr1(msize,msize) !rr1(mls,mls,nats,nats)
  real r2(mls,mls),rr2(msize,msize) !rr2(mls,mls,nats,nats)
  complex v1(mls,mls),v2(mls,mls),v3(mls,mls),v4(mls,mls),v5(mls,mls),expfact
  integer ist1,ist2,iat1,l1,m1,iat2,m2,isp1,iph,i,j,a,b,isp2,n,ik,k
  integer iatom,nt,lun,netab,j1,jj1,j2,jj2
  complex*16 eryd,p
  complex*16,allocatable :: tsst(:,:,:),w1(:,:),w2(:,:)
  character*14 fefffile
  real*8 riajb(3),bvec(3),bvecl(3),bvecl2(3),dif,p2a(3)
  real*8 metric(3,3),brav(3,3),pospos(3),brav2(3,3),posvec(3)
  real norm,dotprod
  complex gk(msize,msize)
  complex prefac
  COMPLEX*16, parameter :: CI2PI = (0.d0,2.d0)*PI
  character*6 norms
  integer irels,nts,zs(nats),msizes
  integer ia,ib,na,nb,niphabs
  logical copyhelp(0:nph) !controls copying Gfms into gg
  logical justone !same deal

  !    Set program control :
  logical,parameter ::  teststrfacs=.false.
!  logical usesprkkrtmat

  logical,parameter ::  calcg0r=.false.
  logical,parameter ::  calcgr=.false.
  logical offdiagonal
  logical blabla
  character*7 fname
  character*75 messg

  LOGICAL,PARAMETER :: FullG = .FALSE.  !KJ 7-09 I think someone else put this here??  Switching it off just like in gglu.f90
  INTEGER istate

  save netab,msizes !KJ bah!


  ! --- runtime message if requested
  if (iverb.gt.0) then
     write(messg, 4010) ie, msize 
4010 format('  FMS matrix at point ', i3,' , number of state kets =', i4)
     call wlog(messg)
  endif


  do i=1,3
     p2a(i)=dble(2)*pi/alat(i)  ! 2 pi / a_j
  enddo

!  usesprkkrtmat=(fullpot.and.(sprkkrpot.eq.1))

  brav(1:3,1)=a1
  brav(1:3,2)=a2
  brav(1:3,3)=a3
  brav2(1:3,1)=b1
  brav2(1:3,2)=b2
  brav2(1:3,3)=b3
  metric=dble(0)
  do i=1,3
     do j=1,3
        do k=1,3
           metric(i,j)=metric(i,j)+brav(k,i)*brav2(k,j)
        enddo
     enddo
  enddo


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
  Cmat=cmplx(0,0)
  if(corehole) then
     !  Make core hole matrix Cmat
     i=ppot(absorber)
     do j=1,nats
        if(ppot(j).eq.i) exit
     enddo
     i=(j-1)*mls
     ist2=0

     do l1=0,maxl
        do m1=-l1,l1
           do isp1=1,nsp
              ist1=ist1+1
              ist2=ist2+1
              Cmat(ist2,ist2)=tmatrx(1,ist1)
              if(nsp.eq.2) then
                 if(isp1.eq.1) then
                    isp2=2
                    ist2=ist1-1
                 else
                    isp2=1
                    ist2=ist1+1
                 endif
                 if(ist2.le.msize.and.ist2.gt.0) Cmat(ist1,ist2)=tmatrx(2,ist1)
              endif  ! nsp=2
           enddo
        enddo
     enddo
     Cmat=Tmat(1+i:mls+i,1+i:mls+i)-Cmat
     Cmat=Cmat*cholestrength
  endif


  !      write(*,*) 'corehole,cholestrength,absorber,ppot(absorber),i,ist1,msize,j'
  !      write(*,*) corehole,cholestrength,absorber,ppot(absorber),i,ist1,msize,j
  !      call writematrix(Tmat,msize,'Tmat.txt')
  !      call writematrix(Cmat,mls,'Cmat.txt')
  !      write(72,*) tmatrx
  !      stop


!   if (usesprkkrtmat) then
!      !         Tmat=dcmplx(0,0)
!      call wlog('Using SPRKKR T-matrix')
!      if(corehole) call wlog('Still using FEFF core hole')
!      DO IATOM=1,nats !nph
!         LUN=70+IATOM
! 2019    FORMAT(3x,5i6,2x,a6)	  
! 2020    FORMAT(10000e14.6)	
! 2021    FORMAT(a2,1x,5i6,2x,a6)	  
!         IF (IE .EQ. 1) then
!            if(iatom.lt.10) then
!               FEFFFILE(1:14)='TFORFEFFX.DAT '
!               FEFFFILE(9:9)=CHAR(48+iatom)
!            elseif(iatom.lt.100) then
!               FEFFFILE(1:14)='TFORFEFFXX.DAT'
!               FEFFFILE(9:9)=char(48+iatom/10)
!               FEFFFILE(10:10)=char(48+iatom-10*(iatom/10))
!            else
!               stop 'FEFF cannot handle more than 99 atoms'
!            endif
!            OPEN(LUN,FILE=FEFFFILE,FORM='FORMATTED',STATUS='old')
!            READ(LUN,2019) MSIZES,NETAB,zs(iatom),irels,nts,norms
!            FEFFFILE(2:4)='VIA'
!            OPEN(LUN+20,FILE=FEFFFILE,FORM='FORMATTED',STATUS='unknown')
!            WRITE(LUN+20,2021) '# ',MSIZES,NETAB,zs(iatom),irels,nts,norms
!            if(msizes.ne.mls) write(*,*) 'msizes',msizes,'ne mls',mls
!            FEFFFILE(2:4)='BBB'
!            OPEN(LUN+40,FILE=FEFFFILE,FORM='FORMATTED',STATUS='unknown')
!            WRITE(LUN+40,2021) '# ',MSIZES,NETAB,zs(iatom),irels,nts,norms
!            if(nts.ne.nph) write(*,*) 'nts',nts,'ne nph',nph
!         ENDIF
!         if(iatom.eq.1) then
!            allocate(tsst(msizes,msizes,nats),w1(msizes,msizes),w2(msizes,msizes))
!            tsst=dcmplx(0,0);w1=dcmplx(0,0);w2=dcmplx(0,0)
!         endif
!         read(LUN,2020) ERYD,P,((w2(I,J),J=1,msizes),I=1,msizes)
!         write(LUN+20,2020) ERYD,P,((Tmat((iatom-1)*mls+I,(iatom-1)*mls+J),J=1,mls),I=1,mls)
!         IF (IE .EQ. NETAB.and.iatom.eq.nph) then
!            write(*,*) 'closing files at energy point',netab
!            CLOSE(LUN)
!         endif
!         !	    call writematrixdblefree(w2,msizes,'w2.txt')
!         !	    call writematrixdblefree(Tmat,msize,'Tmat.txt')
!         IF (IE .EQ. NETAB.and.iatom.eq.nph) CLOSE(LUN+20)
!         !           Convert from SPRKKR basis to FEFF basis
!         !           For nsp=1 : from real to spherical harmonics
!         call changerep(w2,'RLM>CLM',w1,mls,msizes,rc,crel,rrel,'          ',0)
!         !           To get the normalization of the spherical harmonics right, add factor i^l :
!         !            call writematrixdblefree(w1,msizes,'w1.txt')
!         ia=0;ib=0
!         do a=0,maxl;do na=-a,a;ia=ia+1;ib=0
!            do b=0,maxl;do nb=-b,b;ib=ib+1
!               tsst(ia,ib,iatom)=w1(ia,ib)*coni**(a-b)
!            enddo;enddo
!         enddo;enddo	    
!         !	    call writematrixdblefree(tsst(:,:,1),msizes,'tsst.txt')
!         !	    if(msizes.eq.mls) then 
!         Tmat((iatom-1)*mls+1:iatom*mls,(iatom-1)*mls+1:iatom*mls)=tsst(1:mls,1:mls,iatom)
!         !            call writematrixdblefree(Tmat,msize,'Tmat2.txt')
!         !	    stop
!         write(LUN+40,2020) ERYD,P,((Tmat((iatom-1)*mls+I,(iatom-1)*mls+J),J=1,mls),I=1,mls)
!         !	    elseif(msizes.eq.(mls*2)) then
!         !	       if (nsp.eq.2) stop 'whoops 2 !'
!         !	       do l1=0,maxl
!         !	       do m1=-l1,l1
!         !Tmat((iatom-1)*mls+1+l1*(l1+1)+m1,(iatom-1)*mls+1:iatom*mls)=tsst(,1:mls,iatom)
!         !	       do 
!         !	    else
!         !	       stop 'whoops!'
!         !	    endif
!      ENDDO
!      deallocate(tsst)
!   endif

  !********************************** TAKE CARE OF THE STRUCTURE FACTORS :

  !    Calculate structure factors G(k), returned in matrix gk

  eryd = dble(2)*em
  if (dabs(dimag(eryd)).lt.2*eimag) eryd=dcmplx(dble(eryd),dimag(eryd)+eimag)
  p = sqrt(eryd)

  ! ------------------ calculate energy - dependent terms of str.constants
  call strcc(eryd,alat(1),.false.)


!  !KJ
!  goto 1587
!  ! some other test ...
!  riajb=dble(0.3)
!  call structurefactor(p,riajb,gk)
!  open (99,file='g011kk.txt',position='append')
!  write(99,167) ck,(gk(ist1,1:mls),  ist1=1,mls)
!  close(99)
!  open (99,file='g022kk.txt',position='append')
!  write(99,167) ck,(gk(ist1,mls+1:2*mls),  ist1=mls+1,2*mls)
!  close(99)
!  open (99,file='g012kk.txt',position='append')
!  write(99,167) ck,(gk(ist1,mls+1:2*mls), ist1=1,mls)
!  close(99)
!  open (99,file='g021kk.txt',position='append')
!  write(99,167) ck,(gk(ist1,1:mls), ist1=mls+1,2*mls)
!  close(99)
!  return
!1587 continue
!  !KJ



!  if(calcg0r) then
!     ! another debugging test : try to reconstruct g0(r)
!     !   riajb = ria - rjb
!
!     !!   For sc Po, r_12 with (1)=(0,0,0) and (2)=(-a,0,0)
!     riajb=dble(0);  riajb(1)=dble(1) ; iat1=1 ; iat2=1 ; ! sc Po r_12
!     !! Following lines for sc diamond (8 atoms) :
!     !      riajb=dble(0) ; riajb(1)=dble(1) ; iat1=1 ; iat2=1 ! diamond r_1,30  = a lattice vector
!     !      riajb=dble(0); iat1=1 ; iat2=3  ! diamond r_1,16 = in-cell vector
!     !      riajb=dble(0); iat1=1 ; iat2=4  ! diamond r_1,17 = in-cell vector
!     !      riajb=dble(0); iat1=1 ; iat2=5  ! diamond r_1,84 = in-cell vector (1)=(0,0,0) and (2)=(0.75,0.75,0.75)
!     !      riajb=dble(0) ; iat1=1 ; iat2=6 ! diamond r_1,27 in-cell vector (1)=(0,0,0) and (2)=(0.75,0.25,0.25)
!
!     !! Following lines for fcc diamond (2 atoms) :
!     riajb=dble(0);iat1=1;iat2=2  ! 1,2   (1)=(0,0,0) and (2)=-(0.89,0.89,0.89) carth. = -(0.25,0.25,0.25) rhomb.
!     !	riajb=dble(0);riajb(3)=-dble(1);iat1=1;iat2=1  ! 1,15  (1)=(0,0,0) and (2)=(1.78,1.78,0.00) (carth.)=(1,0,0) rhomb - a lattice vector
!     !	riajb=-dble(1);iat1=1;iat2=1  ! 1,167  (1)=(0,0,0) and (2)=(3.56,3.56,3.56) (carth.)=(1,1,1) rhomb - a lattice vector
!     !      riajb=dble(0);iat1=2;iat2=1  ! 1,2   (1)=(0,0,0) and (2)=-(0.89,0.89,0.89) carth. = -(0.25,0.25,0.25) rhomb.- in-cell vector
!     riajb=dble(-1);riajb(1)=dble(1);iat1=1;iat2=1  ! 1-35 propagator  (1)=(0,0,0); (35)=(3.56,0,0) carth
!     iat1=1 ; iat2=2 ; riajb=dble(0) ; riajb(1)=dble(1)  ! 1-36 propagator  (1)=(0,0,0);(2)=(-0.89,-2.65,-2.65) carth
!
!     n=mls !shorthand
!
!     !    riajb is the lattice vector between atom 1 and 2 in lattice coordinates (relative units).
!     !    posvec is the same vector in the carthesian basis (a.u.).
!     !    pospos is posvec, plus the sublattice vector (carthesian) between atom 1 and 2.
!     pospos=dble(0);posvec=dble(0)
!     do i=1,3
!        do j=1,3
!           posvec(i)=posvec(i)+brav(i,j)*riajb(j)
!           pospos(i)=pospos(i)+brav(i,j)*(riajb(j)+ppos(j,iat1)-ppos(j,iat2))
!        enddo
!     enddo
!     write(*,*) 'distance between atoms in a.u.',pospos
!
!     blabla=.false.
!     fname='gxx.txt'
!
!     !      do isp2=1,nsym
!     !	if(isp2.lt.10) then
!     !	   write(fname(2:2),'(a1)') 'x' !'0'
!     !	   write(fname(3:3),'(a1)') char(48+isp2)
!     !	else
!     !	   write(fname(2:2),'(a1)') char(48+(isp2/10))
!     !	   write(fname(3:3),'(a1)') char(48+(isp2-10*(isp2/10)))
!     !	endif
!
!     if(blabla) open(55,file='testsymmetryno.txt')
!     r1=cmplx(0,0)
!     !    KKKKKKKKKKKKKKKKKKKKKK  big loop over k-points
!     !      do ik=4,4 !
!     do ik=1,nkp
!
!        call structurefactor(p,bk(:,ik),gk)
!        bvec=bk(:,ik)
!
!        jj1=(iat1-1)*mls
!        jj2=(iat2-1)*mls
!        v2(1:mls,1:mls)=gk(jj1+1:jj1+mls,jj2+1:jj2+mls)
!        if(usesym.eq.1) then
!           if(blabla) then
!              do j=1,mls
!                 write(55,'(1000f12.4)') v2(j,:)
!              enddo
!           endif
!           !	         do isp1=1,nsym
!           v5=cmplx(0,0)
!           !               do i = 2,2 !n(ik)  !2,2 !1,intn(ik)
!           do i = 1,intn(ik)  !2,2 !1,intn(ik)
!              if(blabla) write(*,*) 'inti',inti(ik,i,2),lsymf(inti(ik,i,1))
!              bvecl=bkf(:,inti(ik,i,1)) !  / p2a(1)  !from lattice units to atomic units
!              v4=cmplx(mrot(:,:,isp2)) !inti(ik,i,2))) !isp1))
!              call cgemm('C','N',n,n,n,c1,v4,n,v2,n,c0,v1,n)
!              call cgemm('N','N',n,n,n,c1,v1,n,v4,n,c0,v3,n)
!              ! sprsym(:,:,k)=wiensym(:,:,symid(1,k))
!              ! wiensym(i)=sprsym(symid(2,i))
!              dotprod=real(0)  ! calculate the product in the carthesian basis
!              do j=1,3
!                 dotprod=dotprod+real(bvecl(j)*posvec(j))
!              enddo
!              v5=v5+v3*real(intw(ik,i))*exp(-coni*dotprod)
!              if(blabla) then
!                 write(55,*) ik,bvec,i,isp1,bvecl,dotprod
!                 do j=1,mls
!                    write(55,'(1000f12.4)') v3(j,:)
!                 enddo
!              endif
!           enddo  ! i
!           !	         enddo  ! isp1
!        else
!           dotprod=real(0)  ! calculate the product in the carthesian basis
!           do i=1,3
!              dotprod=dotprod+real(bvec(i)*posvec(i))
!           enddo
!           if(blabla) then
!              write(55,*) ik,bvec,dotprod
!              do j=1,mls
!                 write(55,'(1000f12.4)') v2(j,:)
!              enddo
!           endif
!           v5=v2*real(weight(ik))  *exp(-coni*dotprod)
!        endif !usesym
!
!        r1=r1+v5
!
!     enddo
!     r1=r1/real(sumweights)
!
!
!     open(99,file='g012k.txt',position='append')
!     !      open(99,file=fname,position='append')
!     do a=1,mls
!        do b=1,mls
!           if (cabs(r1(a,b)).lt.0.00001) r1(a,b)=cmplx(0,0)
!        enddo
!     enddo
!     write(99,167) ck,(r1(a,:),a=1,mls)
!     close(99)
!167  format(5000(e12.4,1x,e12.4,3x))
!
!     !      enddo !isp2
!
!     if(.not.(calcgr.or.teststrfacs)) return 
!  endif



!  if(teststrfacs) then
!     ! Here, we test the structure factors.  The following requirement may be derived theoretically :
!     ! Integral {Brillouin Zone } Gk(k)  dk  = 0.  Only for iat1=iat2 !!
!
!
!
!     rr1=cmplx(0,0)
!     rr2=real(0)
!     n = mls
!     do ik=1,nkp
!        call structurefactor(p,bk(:,ik),gk)
!        do j1=1,nats
!           do j2=1,nats
!              jj1=(j1-1)*mls
!              jj2=(j2-1)*mls
!              !            expfact=EXP(cmplx(-CI2PI*(bk(1,ik)*(ppos(1,j1)-ppos(1,j2))
!              !     1        +bk(2,ik)*(ppos(2,j1)-ppos(2,j2))
!              !     1	    +bk(3,ik)*(ppos(3,j1)-ppos(3,j2)))))
!
!              rr1(jj1+1:jj1+mls,jj2+1:jj2+mls)=rr1(jj1+1:jj1+mls,         &
!                   &        jj2+1:jj2+mls)+gk(jj1+1:jj1+mls,jj2+1:jj2+mls)            &
!                   &        *real(weight(ik)) !*expfact
!           enddo
!        enddo
!        rr2=rr2+cabs(gk)*real(weight(ik))
!
!     enddo
!     norm=real(1)/real(sumweights)
!     rr1=cabs(rr1)*norm
!     rr2=rr2*norm
!
!
!     !         do j1=1,nats
!     !	   do j2=1,nats
!     !	      jj1=(j1-1)*mls
!     !	      jj2=(j2-1)*mls
!     !            r1=cmplx(0,0)
!     !	      r2=real(0)
!     !	      v2(1:mls,1:mls)=gk(jj1+1:jj1+mls,jj2+1:jj2+mls)  !,ik)
!     !            if(usesym.eq.1) then
!     !               do i = 1,intn(ik)
!     !	            v3=cmplx(mrot(1:mls,1:mls,inti(ik,i,2)))
!     !                  call cgemm('C','N',n,n,n,c1,v3,n,v2,n,c0,v1,n)
!     !                  call cgemm('N','N',n,n,n,c1,v1,n,v3,n,c0,v2,n)
!     !	            r1=r1+v2*real(intw(ik,i))
!     !                  r2=r2+cabs(v2)*real(intw(ik,i))
!     !               enddo
!     !	      else
!     !	         r1=r1+v2*real(weight(ik))
!     !	         r2=r2+cabs(v2)*real(weight(ik))
!     !            endif
!     !	      rr1(:,:,j1,j2)=cabs(r1)
!     !	      rr2(:,:,j1,j2)=r2
!     !         enddo
!     !	   enddo
!     !      enddo
!     !	norm=real(1)/real(sumweights)
!     !      rr1=rr1*norm
!     !	rr2=rr2*norm
!
!
!     write(43,'(1000e12.3)') em,                                     &
!          &    ((((abs(rr1(j1*mls+i,j2*mls+i))),            &!plot just residue
!          &     i=1,mls),j1=0,nats-1),j2=0,nats-1)
!     write(44,'(1000e12.3)') em,                                     &
!          &    ((((abs(rr2(j1*mls+i,j2*mls+i))),            &!plot just norm
!          &     i=1,mls),j1=0,nats-1),j2=0,nats-1)
!     write(42,'(1000e12.3)') em,                                     &
!          &    ((((abs(rr1(j1*mls+i,j2*mls+i))/abs(rr2(j1*mls+i,j2*mls+i))), &!plot ratio
!          &     i=1,mls),j1=0,nats-1),j2=0,nats-1)
!     !     1    (((abs(rr1(i,i,j1,1))/abs(rr2(i,i,j1,1))),i=1,mls),j1=1,nats)
!
!     ! The test is finished
!     if (.not.calcgr) return 
!  endif
!146 continue


  !************************  NOW CALCULATE GG USING STRUCTURE FACTORS AND T-MATRIX
  riajb=dble(0)
  offdiagonal=.false.
  posvec=dble(0)

!  if(calcgr) then
!     !       For sc Po, (1) = (0,0,0) and (2) = (-a,0,0)
!     iat1=1;iat2=1;riajb=dble(0);riajb(1)=dble(1)
!     offdiagonal=.true.  ! for a lattice vector
!
!     !       sc diamond, 1-1 propagator  (1)=(0,0,0)
!     iat1=1 ; iat2=1 ; riajb=dble(0) ; offdiagonal=.false.
!     !       sc diamond, 1-2 propagator  (1)=(0,0,0);(2)=(-a/4,-a/4,-a/4)
!     iat1=1 ; iat2=1 ; riajb=dble(0) ; offdiagonal=.false.
!     !       sc diamond, 1-36 propagator  (1)=(0,0,0);(2)=(-a/4,-3a/4,-3a/4)
!     iat1=1 ; iat2=7 ; riajb=dble(0) ; offdiagonal=.false.
!     !       sc diamond, 1-35 propagator  (1)=(0,0,0); (35)=(a,0,0)
!     iat1=1 ; iat2=1 ; riajb=dble(0) ; riajb(1)=-dble(1)
!     offdiagonal=.true.
!     !       sc diamond, 1-58 propagator  (1)=(0,0,0); (35)=(-a/2,a,-a/2)
!     iat1=1 ; iat2=3 ; riajb=dble(0) ; riajb(2)=-dble(1)
!     offdiagonal=.true.
!
!     !       r diamond, 1-1 propagator  (1)=(0,0,0)
!     iat1=1 ; iat2=1 ; riajb=dble(0) ; offdiagonal=.false.
!
!     pospos=dble(0);posvec=dble(0)
!     do i=1,3
!        do j=1,3
!           posvec(i)=posvec(i)+brav(i,j)*riajb(j)
!           pospos(i)=pospos(i)+brav(i,j)*(riajb(j)+ppos(j,iat1)-ppos(j,iat2))
!        enddo
!     enddo
!     write(*,*) 'distance between atoms in a.u.',pospos
!  endif

  ! Determine if there is more than one atom of type ppot(absorber)  !KJ this section 3-2012
  niphabs=0
  do j=1,nats
     if(ppot(j).eq.ppot(absorber)) niphabs=niphabs+1
  enddo
  justone=(niphabs.lt.2)
  if (justone) then
     posvec=0.d0
	 riajb=0.d0
	 riajb(1)=1.d0
	 do i=1,3
	 do j=1,3
	    posvec(i)=posvec(i)+brav(i,j)*riajb(j)
	 enddo
	 enddo
	 ! posvec is now simply the first lattice vector - guess we could've just copied that instead :-P
  endif
  if (offdiagonal.and.justone) then
     call wlog('offdiagonal and justone in fmskspace - you probably do not want this. Stopping.')
	 stop 'contact developers'
  endif

  call kkrintegral(Gfms,Tmat,Cmat,offdiagonal,posvec,p,justone)

  if (dabs(sig2g).gt.0.001) then
     ! Use correlated debye model, sigsqr is in AA^2
     prefac = exp(-1 * sig2g * ck(1)**2 / bohr**2)
!	 write(*,*) 'applying SIG2 factor of ',prefac,' and sig2g is ',sig2g
     Gfms = Gfms * prefac
  endif

!  if(calcgr) then
!     jj1=(iat1-1)*mls
!     jj2=(iat2-1)*mls
!     open(99,file='g12k.txt',position='append')
!     do a=1,msize
!        do b=1,msize
!           if (cabs(Gfms(a,b)).lt.0.00001) Gfms(a,b)=cmplx(0,0)
!        enddo
!     enddo
!     write(99,167) ck,(Gfms(a,jj2+1:jj2+mls),a=jj1+1,jj1+mls)
!     close(99)
!  endif



  !     Now select the parts of G_fms that are needed for the calculation of the spectrum
  !     and put them in the array gg.
  !     Just copy the diagonal of G_fms into gg.
  !     This code is expected to fail if there is only one atom in the unit cell of the type
  !     of the absorber.  fix later
  
  
 copyhelp=.false.
 do j=1,nats
     jj1=(j-1)*mls
     if(j.eq.absorber) then
        do n=1,mls
           do ik=1,mls
              gg(ik,n,0)=Gfms(jj1+ik,jj1+n)
           enddo
        enddo
        copyhelp(0)=.true.
     elseif(.not.copyhelp(ppot(j))) then
        do n=1,mls
           do ik=1,mls
              gg(ik,n,ppot(j))=Gfms(jj1+ik,jj1+n)
           enddo
        enddo
        copyhelp(ppot(j))=.true.              
     end if
  enddo
  istate = nats*(lx+1)**2
  if(FullG) then
     CALL WriteData('g0.bin',Int1 = ie)
     CALL Write2D('g0.bin',Gfms(:istate,:istate))
  end if
  
  if(justone) then
    ! There's only one atom in the unit cell of the type of the absorber.
	! Therefore, the corresponding block of G_fms was used for the information of the core hole.
	! Now we need to retrieve the information for a "neutral" atom of the same type.
	! It will be passed from kkrintegral in the Cmat variable, which earlier contained the t-matrix of the core hole.
	if(copyhelp(absorber)) stop 'Something suspicious in justone section of fmskspace - contact developers'
	gg(1:mls,1:mls,ppot(absorber))=Cmat(:,:)
	copyhelp(absorber)=.true.
  endif

  do j=1,nph
     copyhelp(0)=copyhelp(0).and.copyhelp(j)
  enddo
  if(.not.copyhelp(0)) stop 'copying failed in fmskspace'


  !  gg is ready to be used by fmstot

  return
end subroutine fmskspace   !  subroutine fms2

