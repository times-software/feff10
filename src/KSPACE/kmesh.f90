!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: kmesh.f90,v $:
! $Revision: 1.16 $
! $Author: jorissen $
! $Date: 2012/04/03 22:39:49 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine kmesh
! given the lattice of type lattic in lat, construct the reciprocal lattice and sample the first brillouin
! zone using nk points, gathered in bk.  If usesym and use_external_symmetry_input, use the nsym symmetry operations given in file symfil to reduce
! this mesh of nk points to the irreducible part of the BZ (nk then becomes the number of irreducible points).  
        use struct,ncrystsym=>nsym
        use kklist
        use kgenwork
        use tetrahedra
        use controlkgen
      implicit none
! this routine is to a great extent identical to the program kgen of the wien2k package,
! which is an adapted version of P. Bloechl's programs.

! Some comments about the list of k-vectors :
! The bk are in carthesian units and are limited to the first BZ.
! if eg the lattice is simple cubic of size a, then each of the coordinates
! of bk(.,i) is within [0,2pi/a[.
! for non simple lattices, coordinates may appear weird :
! eg F (diamond) takes 'cubic' lattice parameter input, but can actually
! be reduced to a rhombohedric simple lattice, which has a non cubic BZ.
! In this case, negative coordinates may occur.


!  INPUT
      real*8 lat(6) ! real space lattice specified as 3 lengths and 3 angles
      character*3 lattic    ! defines bravais lattice, ie H(exagonal),R(hombohedric),F(ace centered),
                            ! I(body centered),P(rimitive),CXZ/CYZ/CXY (one face centered)
!      logical usesym     ! true : reduce mesh by symmetry ; false : don't.
        character*20 symfil   ! (optional) file containing symmetry operations (Note: no longer input; now set below; not used anyway ; KJ 1-2012)
!	integer nk  ! at calling time, requested number of k-vectors for total mesh
!  OUTPUT
!      real*8 bk(3,:),volbz,weight(:) ! list of irr. k-vectors; volume of Brillouin zone; integration weight of each vector
!     integer nk : now contains the number of irreducible k-vectors (already declared in 'INPUT' section)
!  LOCALS
      integer iz(3,3,48), iio(3,3,48)
      integer nsym !local copy
        integer i,iarb(3),ndiv(3),ir,ikp,k1,k2,k3,idiv,nkj,j
      integer,allocatable :: klist(:,:)
      real*8 sym4(3,48),pi,afact,sumwgt,ak1,ak2,ak3,am(3,3),bm(3,3),cm(3,3),dm(3,3)
        real*8,allocatable :: wei(:)
      character*512 slog
	  logical ran_already
	  logical, parameter :: write_kmesh_dat=.true.
        logical,parameter :: writewienklist=.false.
        logical,parameter :: writesymops=.false.
        logical,parameter :: use_external_symmetry_input=.false.
		logical,parameter :: debug=.false.

      pi=dacos(dble(-1))                                           
      ran_already=(allocated(bk))
      symfil='di.wiensymmetry     '	     

!****************** GET K-LIST FROM FILE
!  First, we check if the file klist.in exists; if it does, we take k-points from this file instead of calculating them.
!  That means the input from feff.inp is ignored completely.  CAREFUL : the symmetry information in module kklist should be consistent with whatever's in the file ...
      open(18,file='klist.in',form='formatted',status='old',err=1010)
      read(18,*,err=1010,end=1010) nkj,volbz ! we call this nkj instead of nk so that we can fall back to input nk in case of trouble with the file
      call init_kklist(nkj,ncrystsym)  !KJ 6-09
      do i=1,nkj
        read(18,*,err=1009,end=1009) bk(1,i),bk(2,i),bk(3,i),wei(i)
      enddo
      close(18)
        nkp=nkj
        write(slog,25) nkp
25    format(i6,'k-vectors taken from klist.in ; feff.inp ignored')
        call wlog(slog)
      return
1009  stop 'error while reading klist.in in kgen'
1010  continue

!****************** CREATE NEW K-LIST     
      sumwgt=dble(0)
      if(debug) open(66,file='outputkgen.txt',form='formatted',status='unknown')
!	  open(15,file='file.kgen',form='formatted',status='unknown')



      if(use_external_symmetry_input) then
		  iz=0
		  if (usesym .eq. 1) then
			 call reasym(symfil,nsym,iz,sym4)
		  else
			 nsym = 1
			   do k1=1,3
				  iz(k1,k1,1)=1
			   enddo
		  endif
		  write(slog,26) nsym
	26    format('Using',i4,' symmetry operation(s) to reduce the k-mesh.')
		  call wlog(slog)
	  endif


	  if(usesym.eq.1) then
		 nsym=ncrystsym
	  else
		 nsym=1
	  endif


!     copy array
      lat(1:3)=alat
	  lat(4:6)=alfalat !/dble(180)*pi
	  lattic=latticename
!      if (nsym.eq.1) then
!		 lattic='P  ' !KJ this line added KJ - don't use symmetry, primitive lattice
!		 if(.not.ran_already) call wlog ('Lattice type set to Primitive in kgen.')
!	  endif
!KJ I think this would be a big mistake ??

!*********************** PREPARE BRAVAIS AND SYMMETRY MATRICES :

      rbas(1,:)=a1
      rbas(2,:)=a2
      rbas(3,:)=a3
      gbas(1,:)=b1
      gbas(2,:)=b2
      gbas(3,:)=b3
      volbz=(2*pi)**3 / celvol
	  if(debug) then
        open(77,file='kdebug2.txt')
        write(77,*) lattic,lat,ortho,volbz,celvol
        do i=1,3
        write(77,'(3f12.4,7x,3f12.4)') rbas(i,:),gbas(i,:)
        enddo 
        write(77,*) nsym
        do i=1,nsym
        do j=1,3
        write(77,'(i3,3x,3f12.4)') i,iio(j,:,i)
        enddo
        enddo
        close(77)
      endif

      call bravais(lattic,lat(1),lat(2),lat(3),rbas,gbas,afact,iarb,lat(4),lat(5),lat(6),ortho,volbz)          
      call gbass(rbas,gbas)   
      call sdef(iz,nsym,lattic)
      call sdefl(rbas,gbas,iio,nsym,iz,lattic,ortho)

      if(debug)then
        open(77,file='kdebug.txt')
        write(77,*) lattic,lat,ortho,volbz
        do i=1,3
        write(77,'(3f12.4,7x,3f12.4)') rbas(i,:),gbas(i,:)
        enddo 
        write(77,*) nsym
        do i=1,nsym
        do j=1,3
        write(77,'(i3,3x,3f12.4)') i,iio(j,:,i)
        enddo
        enddo
        close(77)
      endif

      iz=iio !backup
      rbas(1,:)=a1
        rbas(2,:)=a2
        rbas(3,:)=a3
        gbas(1,:)=b1
        gbas(2,:)=b2
        gbas(3,:)=b3
        volbz=celvol/(2*pi)**3
        do i=1,48
           sym4(:,i)=cryst_gr(:,4,i,1)
           do j=1,3
!	   iio(j,1:3,i)=nint(cryst_gr(j,1:3,i,1))
           iio(j,1:3,i)=nint(cryst_gr(j,1:3,i,2))
           enddo
        enddo
        
        if(nsym.eq.1) then
!          enforce that the identity is the chosen one
           iio(:,:,1)=0
           sym4(:,1)=dble(0)
           do i=1,3
              iio(i,i,1)=1
           enddo
        endif	

      if(debug) then        
        open(76,file='symoptest.txt')
        do i=1,48
        do j=1,3
        write(76,'(i3,4x,3i4,4x,3i4)') i,iz(j,:,i),iio(j,1:3,i)
        enddo
        enddo
        close(76)
      endif


!************************ FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA
	  call destroy_kklist  !KJ 1-2012 : if running kmesh for the second time, delete allocatable arrays first
	  call destroy_meshes !id.                                               
      nka=nkp      
      call arbmsh(gbas,nsym,iio,iarb,ndiv,sumwgt)       
!     Now nka is the number of points asked for; nkf is the nofp on the
!     actual mesh; nki of these points are in the IBZ.
!     The irreducible points are stored in bk.
      write(slog, 4020) nki,nkf,nka
 4020 format(i6,' k-points in IBZ (',i6,' in full BZ; ',i6,' requested)')
      if(.not.ran_already) call wlog(slog)

     
!********************* COPY RESULTS TO FEFF ARRAYS
      call init_kklist(nki,nsym) !KJ 6-09
	  nkp=nki
	  weight=wi
	  bk=bki


!      do i=1,48
!	do j=1,3
!	do k1=1,3
!	afact=dble(0)
!	do k2=1,3
!      afact=afact+cryst_gr(j,k2,i,1)*cryst_gr(k1,k2,i,1)
!      enddo
!	if((j.eq.k1.and.(dabs(afact-dble(1)).gt.dble(0.0001)))
!     1   .or.(j.ne.k1.and.(dabs(afact).gt.dble(0.0001))))
!     2   write(*,*) 'no inverse',i
!	enddo
!	enddo
!	enddo
!	write(*,*) 'all cleared!'
!	stop

!      call wlog('For convenience, inverting all symmetry matrices.')
        do i=1,nsym
         cm=cryst_gr(:,1:3,i,1)
           do j=1,3
           do k1=1,3
              cryst_gr(j,k1,i,1)=cm(k1,j)
           enddo
           enddo
        enddo
        afact=dble(0)

      if(writesymops) then
           open(55,file='symops_kgen.txt',form='formatted')
           do i=1,3
           do j=1,3
              cm(i,j)=gbas(j,i)/dble(6.28318)
              dm(i,j)=rbas(j,i)
           enddo
           enddo
!           write(*,*) 'possibly need transpose of cm ...'
           write(55,*) 'RBAS ------------ GBAS ------------ CM'
           do i=1,3
              write(55,1018) rbas(i,:),gbas(i,:),cm(i,:)
           enddo
1018     format(3f12.6,5x,3f12.6,5x,3f12.6)
           do i=1,nsym
               write(55,*) 'SYMMETRY OPERATION ',i
               do j=1,3
               do k1=1,3
        	      am(j,k1)=dble(iio(k1,j,i)) !the transpose of iio is its inverse
               enddo
               enddo
!           calculate am = gbas' iio(i)  rbas'  /(2pi)
               call matmm(bm,cm,am)
               call matmm(am,bm,dm)
               do j=1,3
!	          write(55,1019) iio(j,:,i),am(j,:),cryst_gr(j,1:3,i,1),
                  write(55,1019) iio(:,j,i),am(j,:),cryst_gr(j,1:3,i,1),cryst_gr(j,1:3,i,2)
               enddo
             cryst_gr(:,1:3,i,2)=am
           enddo
           close(55)
1019     format(3i6,3(5x,3f12.6))
        endif
                 
		if (allocated(wei)) deallocate(wei,klist)
		allocate(wei(nki),klist(nki,3))
		do ikp=1,nki                                                
            wei(ikp)=wi(ikp)*sumwgt/2.  ! weight=>wwi
		enddo


	  if(writewienklist) then

         IDIV=10000  
         if(iarb(1).eq.1.and.iarb(2).eq.1.and.iarb(3).eq.1) IDIV=NDIV(1)*2  
         if(iarb(1).eq.1.and.iarb(3).eq.0) IDIV=NDIV(1)*ndiv(3)*2          
         if(iarb(3).eq.1.and.iarb(1).eq.0) IDIV=NDIV(1)*ndiv(2)*2          
         if(iarb(2).eq.1.and.iarb(1).eq.0) IDIV=NDIV(2)*ndiv(3)*2          
         if(iarb(1).eq.0.and.iarb(2).eq.0.and.iarb(3).eq.0) IDIV=NDIV(1)*ndiv(2)*ndiv(3)*2                                                
         if (debug) write(66,*) 'nki,NDIV,afact ',nki,NDIV,afact                      
           open(8,file='file.klist',form='formatted',status='unknown')
         do ikp=1,nki !nk                                                  
            ak1=bki(1,ikp)/2./pi*lat(1)  !bk=>bki
            ak2=bki(2,ikp)/2./pi*lat(2)  !bk=>bki
            ak3=bki(3,ikp)/2./pi*lat(3)  !bk=>bki
            if (debug) write(66,6221) AK1,AK2,AK3,ak1*idiv,ak2*idiv,ak3*idiv
            K1=NINT(AK1*IDIV)                                                 
            K2=NINT(AK2*IDIV)                                                 
            K3=NINT(AK3*IDIV)
            if (.not.ortho) then
               klist(ikp,1)=NINT(bKi2(1,ikp)*IDIV)  !bki=>bki2
               klist(ikp,2)=NINT(bKi2(2,ikp)*IDIV)  !bki=>bki2
               klist(ikp,3)=NINT(bKi2(3,ikp)*IDIV)  !bki=>bki2
            else
               klist(ikp,1)=K1
               klist(ikp,2)=K2
               klist(ikp,3)=K3
            endif
         enddo
        
         call divisi(nki,nki,idiv,klist)  !nk,nk

         do ikp=1,nki  !nk
            if(ikp.eq.1) then 
               write(8,1523) IKP,(klist(ikp,ir),ir=1,3),IDIV,wei(ikp),-7.,1.5,nka,ndiv
            else
               write(8,1520) IKP, (klist(ikp,ir),ir=1,3),idiv,wei(ikp)
            endif                               
         enddo                                                          
         write(8,1521)
           close(8)
      endif
	  
	  if(write_kmesh_dat) then
	  open(44,file='kmesh.dat')
	  do ikp=1,nki
            if(ikp.eq.1) then 
               write(44,1525) IKP,(bki(ir,ikp),ir=1,3),wei(ikp),nka,nki,ndiv
            else
               write(44,1525) IKP, (bki(ir,ikp),ir=1,3),wei(ikp)
            endif                               
         enddo                                                          
	  close(44)
	  1525 format(i10,4f9.4,5i7)
	  endif
	  

 1520 format(I10,4I5,f5.1)                                              
 1521 format('end',/)
 1522 format(4f10.7)                                              
 1523 format(I10,4I5,3f5.1,4x,i6,' k, div: (',3i3,')')                  
 6221 format(3f12.5,10x,3f13.5)
      close(66);close(15)

!*********************** CLEAN UP UNWANTED ARRAYS
!      call destroy_meshes 
        if (dotets) call destroy_tetrahedra     

      return                                              
      end     ! subroutine kgen                                                            





!*******************************************************************************************
      subroutine ARBMSH(gbas,nsym,iio,iarb,n,sumwgt)       
!     **  CALCULATE IRREDUCIBLE K-POINTS AND FINDS INEQUIVALENT TETRAHEDRA  **
!     **  FOR BRILLOUIN ZONE INTEGRATION                              **
!     **                                                              **
!     **  INPUT :                                                     **
!     **    nsym        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IIO         SYMMETRY OPERATIONS                           **
!     **    iarb        DEPENDENCIES FOR DIVISIONS                    **
!     **                OF RECIPROCAL LATTICE VECTORS                 **
!     **                ( if iarb(1)=1 then 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  if iarb(3)=1 then SAME FOR 2ND AND 3RD;     **
!     **                  if iarb(2)=1 then SAME FOR 3RD AND 1ST)     **
!     **                                                              **
!     **   AUTHOR: PETER E. BLOECHL                                   **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **     GBASS,BASDIV,REDUZ,TETDIV,TETCNT,ORD1         **
      use kgenwork
        use controlkgen
      implicit none
      integer,intent(in) :: nsym,iio(3,3,nsym),iarb(3)
        integer,intent(inout) :: n(3)
        real*8,intent(in) :: gbas(3,3)
        real*8,intent(out) :: sumwgt
        integer ishift(3),tet0(3,4,6)
        integer nmshp
		logical,parameter :: debug=.false.         
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      nmshp=nka                                                         
      call basdiv(n,nmshp,gbas,iarb)
        nkw=nmshp
        nkf=n(1)*n(2)*n(3)
        call init_workmesh(nkw)
        call init_fullmesh(nkf)
!     ------------------------------------------------------------------
!     -- FIND IRREDUCIBLE K-POINTS                                    --
!     ------------------------------------------------------------------
      call REDUZ(n,ishift,nsym,IIO,sumwgt,gbas)             
      if(debug) write(66,*) ' NO. OF INEQUIVALENT K-POINTS ', nki                                                 
        if(.not.dotets) return
!     ----------------------------------------------------------------- 
!     --  CHOOSE TETRAHEDRA                                          -- 
!     ----------------------------------------------------------------- 
      call TETDIV(n,gbas,tet0)                                      
!     ----------------------------------------------------------------- 
!     --  FIND INEQUIVALENT TETRAHEDRA                               -- 
!     ----------------------------------------------------------------- 
      call TETCNT(TET0,n)               
      return                                                            
      end                                                               
  



!     ..................................................................
      subroutine basdiv(n,nmshp,gbas,iarb)                              
!     **                                                              **
!     **  BASDIV DETERMINES DIVISION OF BASEVECTORS OF REC. LATT.     **
!     **  SO THAT THE NUMBER OF MESHPOINTS IS JUST BELOW nmshp AND    **
!     **  TAKES INTO ACCOUNT THE DEPENDENCY BY POINT-SYMMETRY         **
!     **  INPUT :                                                     **
!     **    gbas        RECIPROCAL LATTICE VECTORS                    **
!     **    iarb        DEPENDENCIES FOR DIVISIONS                    **
!     **                OF RECIPROCAL LATTICE VECTORS                 **
!     **                ( if iarb(1)=1 then 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  if iarb(2)=1 then SAME FOR 1ST AND 3RD;     **
!     **                  if iarb(3)=1 then SAME FOR 3RD AND 2ND)     **
!     **    nmshp       TOTAL NUMBER OF GRID POINTS                   **
!     **  OUTPUT:                                                     **
!     **    n           NUMBER OF DIVISIONS OF RECIPROCAL LATTICE     **
!     **                VECTORS FOR DEFINITION OF SUBLATTICE          **
!     **    nmshp       NUMBER OF SUBLATTICE POINTS in A REC. UNIT CELL*
!     **                AND ON ALL ITS FACES                          **
!     **                                                              **
      implicit none
        integer,intent(in) ::   iarb(3)
        integer,intent(out) ::  n(3)
        integer,intent(inout) ::nmshp                                       
      real*8,intent(in) ::    gbas(3,3)
        integer i                                                 
      real*8,parameter :: opar=1.d-6 
      real*8 betr(3),rn(3),svar  
	  logical,parameter :: debug=.false.                                  
!     ==================================================================
!     == FIND UPPER LIMIT FOR THE LENGTH OF SUBLATTICE VECTORS        ==
!     ==================================================================
      do i=1,3                                                       
         betr(i)=sqrt(gbas(i,1)**2+gbas(i,2)**2+gbas(i,3)**2)            
      enddo                                                         
      svar=(dble(nmshp)/(betr(1)*betr(2)*betr(3)))**(1.D0/3.D0)         
      do i=1,3                                                       
         rn(i)=betr(i)*svar                                                
      enddo                                                          
!     ==================================================================
!     == FIND DIVISIONS OF LATTICE VECTORS                            ==
!     ==================================================================
      if(iarb(1).eq.1.and.iarb(2).eq.1)then                             
        n(1)=int((rn(1)*rn(2)*rn(3))**(1.D0/3.D0)+opar)                 
        n(2)=n(1)                                                       
        n(3)=n(1)                                                       
      else if(iarb(1).eq.1) then                                        
        n(1)=int(sqrt(rn(1)*rn(2))+opar)                               
        n(2)=n(1)                                                       
        n(3)=int(rn(3))                                                 
      else if(iarb(2).eq.1) then                                        
        n(1)=int(sqrt(rn(1)*rn(3))+opar)                               
        n(2)=int(rn(2))                                                 
        n(3)=n(1)                                                       
      else if (iarb(3).eq.1) then                                       
        n(1)=int(rn(1))                                                 
        n(2)=int(sqrt(rn(2)*rn(3))+opar)                               
        n(3)=n(2)                                                       
      else                                                              
        n(1)=int(rn(1)+opar)                                            
        n(2)=int(rn(2)+opar)                                            
        n(3)=int(rn(3)+opar)                                            
      end if                                                            
      n(1)=max0(1,n(1))                                                 
      n(2)=max0(1,n(2))                                                 
      n(3)=max0(1,n(3))                                                 
      if(debug) write(66,*) ' length of reciprocal lattice vectors:',rn
      nmshp=(n(1)+1)*(n(2)+1)*(n(3)+1)                                  
      return                                                            
      end                                                               




!*******************************************************************************

      subroutine REDUZ(n,ishift,nsym,IO,sumwgt,gbas)
!KJ  Constructs the following arrays :
!KJ  bki,bkf,bkw  ;  wi,wf,ww  ;  linkf,linkw  ;  lsymf,lsymw
!     **                                                              **
!     **  REDUZ CREATES THE RELATION BETWEEN THE MESHPOINTS AND       **
!     **  THE POINTS in THE "IRREDUCIBLE ZONE"                        **
!     **  INPUT :                                                     **
!     **    n           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    nmshp       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    nsym        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)*
!     **    in          WORK ARRAY                                    **
!     **  OUTPUT :                                                    **
!     **    num(i)      MAPPING FROM A GENERAL POINT (i) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (num)         **
!     **  REMARKS :                                                   **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :     **
!     **   (X,Y,Z)=rbas*(i,J,K)                                       **
!     **   (i,J,K)  <->  num = i*(n(2)+1)*(n(3)+1)+J*(n(3)+1)+K+1     **
!     **                                                              **
      use kgenwork
        use controlkgen
      implicit none
        integer,intent(in) ::  n(3),nsym,IO(3,3,nsym)
        integer,intent(out) :: ishift(3)
      real*8,intent(out) ::  sumwgt
        real*8,intent(in) ::   gbas(3,3)
        real*8 wsum(nkw)
        integer i,j,i1,i2,i3,j1,j2,j3,k,l,nirr,linkwf(nkw)
        real*8 wgt,summ(3)
      integer nf
        real*8 rinda,rindb,rindc
		logical,parameter :: debug=.false.
!     ------------------------------------------------------------------
!     -- TEST WHETHER SHIFT OF THE SUBLATTICE BY (1/2,1/2,1/2) IS     --
!     --  COMPATIBLE WITH THE SYMMETRYGROUP                           --
!     ------------------------------------------------------------------
!            write(11,*) 'SYMMETRY MATRICES :'  !KJ changed lun 68 to 11
!            do i=1,nsym           
!            write(11,6003) I,((IO(K,L,I),K=1,3),L=1,3)    !KJ changed lun 68 to 11               
!            enddo
      ishift=1                                                       
      do i=1,8                                                       
         bkwi(1,i)=(i-1)/4                                                   
         bkwi(2,i)=(i-bkwi(1,i)*4-1)/2                                         
         bkwi(3,i)=i-bkwi(1,i)*4-bkwi(2,i)*2-1                                   
      enddo                                                          
      do i=1,nsym                                                    
      do J=1,8                                                       
         I1=2*bkwi(1,J)+ishift(1)                                            
         I2=2*bkwi(2,J)+ishift(2)                                            
         I3=2*bkwi(3,J)+ishift(3)                                            
         J1=IO(1,1,i)*I1+IO(1,2,i)*I2+IO(1,3,i)*I3                         
         J2=IO(2,1,i)*I1+IO(2,2,i)*I2+IO(2,3,i)*I3                         
         J3=IO(3,1,i)*I1+IO(3,2,i)*I2+IO(3,3,i)*I3                         
         if(DMOD(dble(J1-ishift(1)),2.D0).ne.0.D0.or.                   &
     &     DMOD(dble(J2-ishift(2)),2.D0).ne.0.D0.or.                    &
     &     DMOD(dble(J3-ishift(3)),2.D0).ne.0.D0) then                    
            ishift=0                                                     
            if(debug) write(66,*) 'SUBMESH SHIFT CONFLICTS WITH POINTGROUP'           
            if(debug) write(66,6003) I,((IO(K,L,I),K=1,3),L=1,3)                 
6003        format(1H ,'SYMMETRYMATRIX NR. : ',I5/3(1H ,3I10/))             
            goto 30                                                         
         end if                                                            
      enddo
        enddo
                                                          
30    continue
!!       ishift=0                                                   
      if(ishift(1).eq.1.or.ishift(2).eq.1.or.ishift(3).eq.1) then       
         ishift=1
         if(debug) write(66,*) ' SUBMESH SHIFTED; SHIFT: ',ishift                  
      else                                                              
         if(debug) write(66,*)' SUBMESH NOT SHIFTED; SHIFT: ',ishift               
      end if        

!     ==================================================================
!     ==  INITIALIZE                                                  ==
!     ==================================================================
      wsum=dble(0)
      do i=1,nkw
         linkw(i)=i                                                          
         bkwi(1,i)=(i-1)/((n(3)+1)*(n(2)+1))                                 
         bkwi(2,i)=(i-bkwi(1,i)*(n(2)+1)*(n(3)+1)-1)/(n(3)+1)                  
         bkwi(3,i)=i-bkwi(1,i)*(n(2)+1)*(n(3)+1)-bkwi(2,i)*(n(3)+1)-1            
      enddo                                                          
!     ==================================================================
!     ==  REDUCTION OF REDUCIBLE K-POINTS                             ==
!     ==================================================================
      l=0
        lsymw=0
      do i=1,nsym                                                    
      do J=1,nkw                                                   
         I1=2*bkwi(1,J)+ishift(1)                                            
         I2=2*bkwi(2,J)+ishift(2)                                            
         I3=2*bkwi(3,J)+ishift(3)                                            
         J1=MOD(IO(1,1,i)*I1+IO(1,2,i)*I2+IO(1,3,i)*I3,2*n(1))             
         J2=MOD(IO(2,1,i)*I1+IO(2,2,i)*I2+IO(2,3,i)*I3,2*n(2))             
         J3=MOD(IO(3,1,i)*I1+IO(3,2,i)*I2+IO(3,3,i)*I3,2*n(3))             
         J1=J1+(1-ISIGN(1,J1))*n(1)                                        
         J2=J2+(1-ISIGN(1,J2))*n(2)                                        
         J3=J3+(1-ISIGN(1,J3))*n(3)                                        
         J1=(J1-ishift(1))/2                                               
         J2=(J2-ishift(2))/2                                               
         J3=(J3-ishift(3))/2
           k=J1*(n(2)+1)*(n(3)+1)+J2*(n(3)+1)+J3+1
           if (k.lt.linkw(j).or.(k.eq.linkw(j).and.lsymw(j).eq.0)) then
              lsymw(j)=i
        	  linkw(J)=min0(linkw(J),k)
           endif
           if (k.eq.j.and.i.eq.1) then
              l=l+1
              linkwf(j)=l
           elseif(i.eq.1.and.k.lt.j) then
              linkwf(j)=min(linkwf(k),k)
           elseif(i.eq.1.and.k.gt.j) then
              l=l+1
              linkwf(j)=l
              write(*,*) 'warning - retest now!'
         endif
      enddo
        enddo


      if(debug) write(66,6000) nkf,(n(i),i=1,3)                         
 6000 format(1H ,' NO. OF MESH POINTS in THE BRILLOUIN ZONE =',I6/      &
     & '  DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)=',3I5)      
      if(verbose) write(66,*) ' point    coordinates     relation'
      sumwgt=dble(0)
        wf=dble(0) 
      nirr=0 
      do j=1,nkw
         if(linkw(j).eq.j) nirr=nirr+1  !calculate the number of irreducible points 
         wgt=1.
         if(mod(bkwi(1,j),n(1)).eq.0) wgt=wgt/2.     ! in center : weight=1 ;  in face : weight=0.5                        
         if(mod(bkwi(2,j),n(2)).eq.0) wgt=wgt/2.     ! on side   : weight=0.25; on corner : weight=0.125
         if(mod(bkwi(3,j),n(3)).eq.0) wgt=wgt/2.                            
           ww(j)=wgt
           wf(linkwf(j))=wf(linkwf(j))+wgt
         wsum(linkw(j))=wsum(linkw(j))+wgt  !the precursor of wi
         sumwgt=sumwgt+wgt                                                
         if(verbose .or. debug) then 
               write(66,100) j,bkwi(1,j),bkwi(2,j),bkwi(3,j),linkw(j),wgt,&
     &                           wsum(linkw(j))
           endif
100      format(i5,5x,3i4,i8,2f10.5)
      enddo 
        nki=nirr
      call init_irrmesh(nki)

      do j=1,nkf
           if(wf(j).ne.1.) write(*,*) j,wf(j)
        enddo

      nirr=0
      do j=1,nkw
         if(linkw(j).eq.j) then
            nirr=nirr+1 
            wi(nirr)=wsum(linkw(j))  !wien2k has a factor *2 here, which I have removed.
         endif       
      enddo
        ww=ww/sumwgt
        wf=wf/sumwgt
        wi=wi/sumwgt
        summ=dble(0)
        do i=1,nkw;summ(1)=summ(1)+ww(i);enddo
        do i=1,nkf;summ(2)=summ(2)+wf(i);enddo
        do i=1,nki;summ(3)=summ(3)+wi(i);enddo
      if(verbose.or.debug) write(66,153) summ
         

      nirr=0
        nf=0                                                            
      i=0    ! point i in the worklist is point nf in the full list and/or point nirr in the irr list                                                           
      if(verbose.or.debug) write(66,*)' internal and cartesian k-vectors:'
      do I1=1,n(1)+1                                                   
      do I2=1,n(2)+1                                                 
      do I3=1,n(3)+1                                                 
         i=i+1                                                             
         RINDA=(dble(I1-1)+dble(ishift(1))/2.D0)/dble(n(1))              
         RINDB=(dble(I2-1)+dble(ishift(2))/2.D0)/dble(n(2))              
         RINDC=(dble(I3-1)+dble(ishift(3))/2.D0)/dble(n(3))              
         bkw(1,i)=gbas(1,1)*RINDA+gbas(2,1)*RINDB+gbas(3,1)*RINDC      
         bkw(2,i)=gbas(1,2)*RINDA+gbas(2,2)*RINDB+gbas(3,2)*RINDC      
         bkw(3,i)=gbas(1,3)*RINDA+gbas(2,3)*RINDB+gbas(3,3)*RINDC      
         if(i.eq.linkw(i))then                                               
            nirr=nirr+1                                                     
            linkw(i)=nirr
        	  bki(:,nirr)=bkw(:,i)                                                     
            bki2(1,nirr)=rinda
            bki2(2,nirr)=rindb
            bki2(3,nirr)=rindc
            if(verbose.or.debug) write(66,103) rinda,rindb,rindc,                &
     &                        bki(1,nirr),bki(2,nirr),bki(3,nirr)
 103        format(3f10.5,10x,3f10.5)
         else                                                              
            linkw(i)=linkw(linkw(i))                                              
         end if
           if(linkwf(i).eq.(nf+1)) then
              nf=nf+1
        	  linkf(nf)=linkw(i)
              lsymf(nf)=lsymw(i)
        	  bkf(:,nf)=bkw(:,i)
           endif                                                            
      enddo
        enddo
        enddo                                                          
      if(nirr.ne.nki) stop 'problem in reduc'


      if(nkf.gt.5000) then
         write(*,*) 'skipping sanity check in reduz - too expensive'
        else
! sanity check :
      nirr=0
      do i=1,nkf
        do j=1,nkf
         if(i.ne.j.and.linkf(i).eq.linkf(j).and.lsymf(i).eq.lsymf(j)) then
            nirr=nirr+1
            write(*,*) 'points ',i,' and ',j,' problematic in reduz'
         endif
        enddo
        enddo
        if(nirr.gt.0) stop
!
      endif


      return                                                            
  150 format('  weights of k-points:')
  152 format(2i8,2f12.6)  
  153 format('  sum of weights: ',3f14.6)
      end                                                               


!*****************************************************************************************
      subroutine TETDIV(n,gbas,TET0)                                    
!     **                                                             ** 
!     **  TETDIV DETERMINES THE DIVISION OF THE PARALLELEPIPEDS      ** 
!     **  in TETRAHEDRONS ACCORDING TO THE SHORTEST DIAGONAL         ** 
!     **                                                             ** 
!     **  INPUT:                                                      **
!     **    n           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    gbas        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT:                                                    ** 
!     **    TET0(i,J,K) COORDINATES (i) in THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) in A SUBLATTICE UNIT CELL    ** 
!     **                                                             ** 
      implicit none                               
        real*8,intent(in) :: gbas(3,3)
        integer,intent(in) :: n(3)
        integer,intent(out) :: tet0(3,4,6)
      real*8 P(8,3),DIAG(4)                                
      integer IACHT(8),TET(4,6),i,j,k,l,isvar,mndg
!     ------------------------------------------------------------------
!     -- SEARCH FOR THE SHORTEST DIAGONAL                             --
!     ------------------------------------------------------------------
      do i=0,1                                                       
      do J=0,1                                                       
      do K=0,1                                                       
      ISVAR=4*i+2*J+K+1                                                 
      do L=1,3                                                       
         P(ISVAR,L)=gbas(1,L)*dble(i)/dble(n(1))+gbas(2,L)*             &
     &           dble(J)/dble(n(2))+gbas(3,L)*dble(K)/dble(n(3))
      enddo
        enddo
        enddo
        enddo
                                                                      
      do i=1,4                                                      
         DIAG(i)=0.D0                                                     
         do J=1,3                                                      
            DIAG(i)=DIAG(i)+(P(i,J)-P(9-i,J))**2                             
           enddo
        enddo

      MNDG=1                                                           
      do i=2,4                                                      
         if(DIAG(i).LT.DIAG(MNDG)) then                                   
            MNDG=i                                                         
         end if                                                           
      enddo                                                         
!     ------------------------------------------------------------------
!     -- ROTATE PARALLELEPIPED                                        --
!     ------------------------------------------------------------------
      if(MNDG.eq.1)then                                                 
        do i=1,8                                                     
           IACHT(i)=i                                                      
          enddo
      else if(MNDG.eq.2) then                                           
        do i=1,4                                                     
           IACHT(2*i-1)=2*i                                                
           IACHT(2*i)=2*i-1                                                
          enddo
      else if(MNDG.eq.3) then                                           
        do i=0,1                                                     
        do J=1,2                                                     
           IACHT(4*i+J)=4*i+J+2                                            
           IACHT(4*i+J+2)=4*i+J                                            
          enddo
          enddo
      else if(MNDG.eq.4) then                                           
        do i=1,4                                                     
           IACHT(i)=i+4                                                    
           IACHT(i+4)=i                                                    
          enddo

      end if                                                            
!      **  CREATION OF TETRAHEDRA  **                                   
      do i=1,6                                                      
         TET(1,i)=IACHT(1)                                                
         TET(4,i)=IACHT(8)                                                
        enddo

      TET(2,1)=IACHT(2)                                                
      TET(3,1)=IACHT(4)                                                
      TET(2,2)=IACHT(4)                                                
      TET(3,2)=IACHT(3)                                                
      TET(2,3)=IACHT(3)                                                
      TET(3,3)=IACHT(7)                                                
      TET(2,4)=IACHT(7)                                                
      TET(3,4)=IACHT(5)                                                
      TET(2,5)=IACHT(5)                                                
      TET(3,5)=IACHT(6)                                                
      TET(2,6)=IACHT(6)                                                
      TET(3,6)=IACHT(2)                                                
                                                                      
      do i=1,4                                                       
      do J=1,6                                                       
         TET0(1,i,J)=(TET(i,J)-1)/4                                        
         TET0(2,i,J)=(TET(i,J)-TET0(1,i,J)*4-1)/2                          
         TET0(3,i,J)=TET(i,J)-TET0(1,i,J)*4-TET0(2,i,J)*2-1                
        enddo
        enddo

      return                                                            
      end                                                               



!******************************************************************************


      subroutine TETCNT(TET0,n)
!     **  TETCNT CALCULATES ALL DIFFERENT TETRAHEDRA AND COUNTS THEM  **
!     **  INPUT :                                                     **
!     **    nmshp       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    num(i)      MAPPING FROM A GENERAL POINT (i) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (num)         **
!     **    TET0(i,J,K) COORDINATES (i) in THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) in A SUBLATTICE UNIT CELL    ** 
!     **    n           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    MWRIT       INFORMATION FOR MWRIT TETRAHEDRA ARE WRITTEN  **
!     **                AT ONE TIME.                                  **
      use kgenwork
        use controlkgen
        use tetrahedra
      implicit none
        integer,intent(in) :: n(3)   
      integer TET0(3,4,6)
      integer ntt,nrec,mwrittest
        integer i,ind,ip,ipp,isvar1,isvar2,ixx,j,jmax,k,k1,k2,k3,l,m,ni,&
     &            nkitest,ntmax
      real*8 sum,v
      integer,parameter :: ICHK=0
        
      call init_tetrahedra(n(3),nkf)               
        
                                                        
      NTMAX=n(1)*n(2)*n(3)*6                                            
      IPP=0                                                             
      do K1=1,n(1)                                                   
      do K2=1,n(2)                                                   
         IND=0                                                             

         do K3=1,n(3)     
            IP=K3+(n(3)+1)*((K2-1)+(n(2)+1)*(K1-1))                           
                                                                        
            do i=1,6                                                       
            IND=IND+1                                                         
            do J=1,4                                                       
            IXX=TET0(1,J,i)*(n(2)+1)*(n(3)+1)+TET0(2,J,i)*(n(3)+1)      &
     &                 +TET0(3,J,i)                                                   
            ITET(J,IND)=IP+IXX                                                
            enddo
              enddo
           enddo                                                          
                                                                        
!     --  TRANSFORM THE EDGEPOINTS ONTO THE IRREDUCIBLE POINTS          
         do M=1,4                                                       
         do J=1,n(3)*6                                                  
             ITET(M,J)=linkw(ITET(M,J))                                          
         enddo
           enddo
                                                                        
!     --  ORDER THE POINTS OF EACH TETR. ACC. TO INREASING NUMBER       
         do K=1,3                                                       
         do J=K+1,4                                                     
         do L=1,n(3)*6                                                  
            ISVAR1=ITET(K,L)                                                  
            ISVAR2=ITET(J,L)                                                  
            ITET(K,L)=min0(ISVAR1,ISVAR2)                                     
            ITET(J,L)=max0(ISVAR1,ISVAR2)                                     
         enddo
           enddo
           enddo
                                                               
!     --  IDENTIFY THE TETRAHEDRA WITH INTEGERS                         
         do M=1,n(3)*6
            IPP=IPP+1                                                         
            IY(1,IPP)=ITET(1,M)
            IY(2,IPP)=ITET(2,M)
            IY(3,IPP)=ITET(3,M)
            IY(4,IPP)=ITET(4,M)
         enddo
         if(IPP.GE.NTMAX) then                                             
            IPP=NTMAX                                                       
            goto 100                                                        
         end if                                                            
                                                                      
      enddo
        enddo

      print*,'UNNORMAL end OF LOOP.......................STOP in TETCNT'
      STOP                                                              
                                                                        
100   continue                                                          
!     ------------------------------------------------------------------
!     --  ORDER TETRAHEDRA                                            --
!     ------------------------------------------------------------------
      NTET=IPP                                                          
      call ORD1(NTET,IY)

!     ==  CHECK ORDERING AND CALCULATE NUMBER OF INEQUIVALENT TETRAHEDRA
      NTT=1                                                             
      do i=1,NTET-1                                                 
         NTT=NTT+1                                                         
         if(IY(4,i+1).eq.IY(4,i)) then                                     
         if(IY(3,i+1).eq.IY(3,i)) then                                     
         if(IY(2,i+1).eq.IY(2,i)) then                                     
         if(IY(1,i+1).eq.IY(1,i)) then                                     
            ntt=ntt-1
         end if
         end if
         end if
         end if
      enddo
!      write(66,1018) NTT                                                 
1018  format(1H ,'NUMBER OF DIFFERENT TETRAHEDRA :',I5)                 
!     ------------------------------------------------------------------
!     --  write ON FILE                                               --
!     ------------------------------------------------------------------
      if((ntt/mwrit)*mwrit.eq.ntt) then
         NREC=NTT/MWRIT  
      else
         NREC=NTT/MWRIT+1
      end if                                                            
      V=1.D0/dble(6*n(1)*n(2)*n(3))                                     
      SUM=V*dble(NTET)                                          
      if(DABS(SUM-1.D0).GT.1.d-5) then                                  
        print*,'SUMRULE NOT FULLFILLED...................STOP in TETCNT'
        print*,' SUM ',SUM,' SHOUD BE EQUAL TO 1'                       
        STOP                                                            
      end if
                                                                    
      REWIND 15                                                         
      write(15,1234) nki,NTT,V,MWRIT,NREC                                
      NTT=0                                                             
      NREC=0
        ittfl=0                                                            
      do i=1,NTET                                                   
         if(i.eq.1) then
             NTT=NTT+1                                                       
                                                                        
             if(NTT.GT.MWRIT*(NREC+1)) then                                  
                NREC=NREC+1                                                  
                write(15,1235)(ITTFL(J),J=1,5*MWRIT)  
        		  ittfl=0                       
             end if                                                          
                                                                        
             NI=5*(NTT-1-NREC*MWRIT)                                         
             ITTFL(NI+1)=1                                               
             ITTFL(NI+2)=IY(1,i)
             ITTFL(NI+3)=IY(2,i)
             ITTFL(NI+4)=IY(3,i)
             ITTFL(NI+5)=IY(4,i)
           else if(IY(1,i).eq.IY(1,i-1).and.IY(2,i).eq.IY(2,i-1)        &
     &   .and.IY(3,i).eq.IY(3,i-1).and.IY(4,i).eq.IY(4,i-1)) then         
            NI=5*(NTT-1-NREC*MWRIT)+1                                       
            ITTFL(NI)=ITTFL(NI)+1                                       
         else 
            NTT=NTT+1                                                       
                                                                        
            if(NTT.GT.MWRIT*(NREC+1)) then                                  
               NREC=NREC+1                                                  
               write(15,1235)(ITTFL(J),J=1,5*MWRIT)                         
               ittfl=0
            end if                                                          
                                                                        
            NI=5*(NTT-1-NREC*MWRIT)                                         
            ITTFL(NI+1)=1                                               
            ITTFL(NI+2)=IY(1,i)
            ITTFL(NI+3)=IY(2,i)
            ITTFL(NI+4)=IY(3,i)
            ITTFL(NI+5)=IY(4,i)
                                                                     
         end if                                                            
      enddo                                                          
      NREC=NREC+1                                                       
      write(15,1235)ITTFL                                               
                                                                        
      if(ICHK.eq.0) return                                              
!     ------------------------------------------------------------------
!     --  CHECK OF i/O                                                --
!     ------------------------------------------------------------------
      REWIND(15)  
      read(15,1234) nkitest,NTT,V,MWRITtest,NREC     
        if(mwrittest.ne.mwrit) stop 'mwrit is wrong in tetcnt'              
        if(nkitest.ne.nki) stop 'nki is wrong in tetcnt'
 1234 format(2i10,e20.12,2i10) 
      SUM=dble(0)
                                                                  
      do i=1,NREC                                                   
         ittfl=0
         read(15,1235)(ITTFL(J),J=1,5*MWRIT)                               
 1235    format(6i10)
         JMAX=5*min0(MWRIT,NTT-(i-1)*MWRIT)                                
         do J=1,JMAX/5                                                 
            SUM=SUM+dble(ITTFL(5*(J-1)+1))*V                                  
         enddo
        enddo
                                                             
      if(DABS(SUM-1.D0).GT.1.d-5) then                                  
        print*,'SUMRULE NOT FULLFILLED...................STOP in TETCNT'
        print*,' SUM ',SUM,' SHOUD BE EQUAL TO 1'                       
        STOP                                                            
      end if                                                            
      return                                                            
      end                                                               




!************************************************************************************





      subroutine bravais(latti,AX,BX,CX,rbas,gbas,afact,iarb,           &
     &  alpha,beta,gamma,ortho,v)
!      implicit real*8 (A-H,O-Z)
!     **                                                              **
!     **  CONSTRUCTION OF TRANSLATION VECTORS : rbas (AS COLUMNS)     **
!     **  AND RECIPROCAL LATTICE VECTORS      : gbas (AS ROWS   )     **
!     **   V: VOLUME OF THE BRILLOUIN ZONE                            **
!     **  INPUT IS THE NAME OF THE BRAVAIS LATTICE ACCORDING TO       **
!     **  TABLE 3.3 OF BRADLEY AND CRACKNELL : THE MATHEMATICAL       **
!     **  THEORY OF SYMMETRY in SOLIDS (OXFORD)                       **
!     ..................................................................
!:UB  The direct lattice vectors     : ai = rbas(i,*)
!:UB  The reciprocal lattice vectors : bi = gbas(*,i)
!:UB           redefined by GBASS as : bi = gbas(i,*) !
!:UB  if( ndiv1 = ndiv2 ) iarb(1) = 1
!:UB  if( ndiv1 = ndiv3 ) iarb(2) = 1
!:UB  if( ndiv2 = ndiv3 ) iarb(3) = 1
!:UB  ..................................................................
      implicit none
      real*8,intent(out) ::     gbas(3,3),rbas(3,3),v
        real*8,intent(inout) ::   ax,bx,cx
        real*8,intent(in) ::      alpha,beta,gamma
        integer,intent(out) ::    iarb(3)
      character*3,intent(in) :: latti
      logical,intent(out) :: ortho
      real*8 eps(3,3,3),det,ay,az,by,bz,cy,cz,pi,a1,a2,cosg1,afact,     &
     &   gamma0
        integer i,j,k

      eps=dble(0)
      afact=dble(1)
        iarb=1
        gbas=dble(0)
        rbas=dble(0)
      det=dble(0)
      ay=dble(0)
      az=dble(0)
      by=dble(0)
      bz=dble(0)
      cy=dble(0)
      cz=dble(0)

      
        if (latti(1:1).eq.'H') then

!        HEXAGONAL : GH
         rbas(1,1)=AX*sqrt(.75E0)
         rbas(1,2)=-AX/2.
         rbas(2,2)=AX
         rbas(3,3)=CX
         iarb(2)=0
         iarb(3)=0
         ortho=.FALSE.

        elseif (latti(1:1).eq.'F') then

!        ORTHOROMBIC : GOF
         AX=AX*0.5E0
         BX=BX*0.5E0
         CX=CX*0.5E0
         rbas(1,2)=BX
         rbas(1,3)=CX
         rbas(2,1)=AX
         rbas(2,3)=CX
         rbas(3,1)=AX
         rbas(3,2)=BX
         afact=0.5
         ortho=.TRUE.

        elseif (latti(1:1).eq.'B') then

!        ORTHORHOMBIC : GOV
         AX=AX*0.5E0
         BX=BX*0.5E0
         CX=CX*0.5E0
         rbas(1,1)=-AX
!        ai = rbas(i,*) >>
         rbas(1,2)=BX
         rbas(1,3)=CX
         rbas(2,1)=+AX
         rbas(2,2)=-BX
         rbas(2,3)=+CX
         rbas(3,1)=AX
         rbas(3,2)=+BX
         rbas(3,3)=-CX
         afact=0.5
         ortho=.TRUE.

        elseif ((latti(1:1).eq.'P'.and.abs(gamma-1.570796d0).gt.0.0001) &
     &    .or.(latti(1:1).eq.'P'.and.abs(beta-1.570796d0).gt.0.0001)    &
     &    .or.(latti(1:1).eq.'P'.and.abs(alpha-1.570796d0).gt.0.0001)) then

!     TRICLINIC : GT
         cosg1=(cos(gamma)-cos(alpha)*cos(beta))/sin(alpha)/sin(beta)
         gamma0=acos(cosg1)
         rbas(1,1)=AX*sin(gamma0)*sin(beta)
         rbas(1,2)=AX*cos(gamma0)*sin(beta) 
         rbas(2,2)=BX*sin(alpha)
         rbas(1,3)=AX*cos(beta)
         rbas(2,3)=BX*cos(alpha)
         rbas(3,3)=CX
         iarb=0
         ortho=.FALSE.

      elseif ((latti(1:1).eq.'C'.and.abs(gamma-1.570796d0).gt.0.0001)   &
     &        .or.       (latti(1:3).eq.'MXZ')) then

!        MONOCLINIC : GMB
!        AX = A * SIN(GAMMA) / 2
!        AY = A * COS(GAMMA) / 2
!        CX = C / 2
!        ay=ax(orig)*cos(gamma)/2 thus ay must be evaluated first
         ay=ax*cos(gamma)/2.
         ax=ax*sin(gamma)/2.
         cx=cx/2.
         rbas(1,1)=ax
         rbas(1,2)=ay
         rbas(1,3)=-cx
         rbas(2,2)=bx
         rbas(3,1)=aX
         rbas(3,2)=ay
         rbas(3,3)=CX
         iarb(1)=0
!:UB[    |b1| = |b3|
         iarb(3)=0
         ortho=.FALSE.

      elseif ((latti(1:1).eq.'S').or.(latti(1:1).eq.'P')) then

!        ORTHORHOMBIC : GO
         rbas(1,1)=AX
         rbas(2,2)=BX
         rbas(3,3)=CX
         iarb=0
         ortho=.TRUE.

        elseif (latti(1:1).eq.'C') then

!        ORTHORHOMBIC : GOB
         if(latti(2:3).eq.'XZ') then
  
            rbas(1,1)=AX*0.5E0
!           ai = rbas(i,*) >>
            rbas(1,3)=-CX*0.5E0
            rbas(3,1)=AX*0.5E0
            rbas(3,3)=CX*0.5E0
            rbas(2,2)=BX
            iarb(1)=0
            iarb(3)=0
            ortho=.TRUE.
 
         elseif(latti(2:3).eq.'YZ') then

            rbas(2,2)=BX*0.5E0
!           ai = rbas(i,*) >>
            rbas(2,3)=-CX*0.5E0
            rbas(3,2)=BX*0.5E0
            rbas(3,3)=CX*0.5E0
            rbas(1,1)=AX
            iarb(1)=0
            iarb(2)=0
            ortho=.TRUE.

         else

            rbas(1,1)=AX*0.5E0
!           ai = rbas(i,*) >>
            rbas(1,2)=-BX*0.5E0
            rbas(2,1)=AX*0.5E0
            rbas(2,2)=BX*0.5E0
            rbas(3,3)=CX
            iarb(2)=0
            iarb(3)=0
            ortho=.TRUE.

         endif

        elseif(latti(1:3).eq.'M  ') then

!        MONOCLINIC : GM
!        AX = A * SIN(GAMMA)
!        AY = A * COS(GAMMA)
         a1=ax*sin(gamma)
         a2=ax*cos(gamma)
         rbas(1,1)=a1
         rbas(1,2)=a2
         rbas(2,2)=bx
         rbas(3,3)=CX
         iarb=0
         ortho=.FALSE.

        elseif(latti(1:1).eq.'R') then

!        TRIGONAL : GRH
         rbas(1,1)=ax/2.d0/sqrt(3.d0)
         rbas(1,2)=-AX/2.d0
         rbas(1,3)=CX/3.d0
         rbas(2,1)=AX/2.d0/sqrt(3.d0)
         rbas(2,2)=AX*0.5E0
         rbas(2,3)=CX/3.d0
         rbas(3,1)=-AX/sqrt(3.d0)
         rbas(3,2)=0.d0
         rbas(3,3)=CX/3.d0
         ortho=.FALSE.

        else

!        CUBIC : GCV
         AX=AX*0.5E0
         rbas(1,1)=-AX
         rbas(2,1)=AX
         rbas(3,1)=AX
         rbas(1,2)=AX
         rbas(2,2)=-AX
         rbas(3,2)=AX
         rbas(1,3)=AX
         rbas(2,3)=AX
         rbas(3,3)=-AX
         afact=0.5
         ortho=.TRUE.

        endif

      do J=1,3
!:UB[ << ai = rbas(i,*) >>
!         write (66,240) J,(rbas(J,i),i=1,3)
      enddo

! Careful : below, gbas is calculated.  However, in wien, this gbas is not used
! (at least not in the kgen program).  So I am not 100 % sure about consistency
! with other data, although comments in the first lines of this file claim it is
! just the transpose of what it 'should' be.
  240 format (1H ,' R',I1,' = ',3F10.6)
      pi=4.*ATAN(1.)
      EPS(1,2,3)=1.E0
      EPS(2,3,1)=1.E0
      EPS(3,1,2)=1.E0
      EPS(1,3,2)=-1.E0
      EPS(3,2,1)=-1.E0
      EPS(2,1,3)=-1.E0
      do i=1,3
      do J=1,3
      do K=1,3
         det=det+EPS(i,J,K)*rbas(1,i)*rbas(2,J)*rbas(3,K)
         gbas(i,1)=gbas(i,1)+EPS(i,J,K)*rbas(2,J)*rbas(3,K)
         gbas(i,2)=gbas(i,2)+EPS(i,J,K)*rbas(3,J)*rbas(1,K)
         gbas(i,3)=gbas(i,3)+EPS(i,J,K)*rbas(1,J)*rbas(2,K)
      enddo
        enddo
        enddo
      do i=1,3
      do J=1,3
         gbas(J,i)=2*pi*gbas(J,i)/det
      enddo
        enddo
      V=(2*pi)**3/det
290   format (1H ,' G',I1,' = ',3F10.6)
!      write (66,300) iarb
300   format (1H ,' DEPENDENCE OF DIVISION OF TRANSLATION VECTORS iarb=' ,3I3)
      return
      end



!********************************************************************************


      subroutine divisi(idkp,nkp,idiv,klist)
      implicit none
        integer,intent(in) ::    idkp,nkp
        integer,intent(out) ::   klist(idkp,3)
        integer,intent(inout) :: idiv
      integer,parameter ::     nprim=16
      integer,parameter ::     niter=10
        integer iprim(nprim),idivi,ip,ik,ir,idummy,itest

      idivi=1
        iprim(1)=2
        iprim(2)=3
        iprim(3)=5
        iprim(4)=7
      iprim(5)=11      
      iprim(6)=13     
      iprim(7)=17     
      iprim(8)=19     
      iprim(9)=23     
      iprim(10)=29     
      iprim(11)=31     
      iprim(12)=37     
      iprim(13)=41     
      iprim(14)=43     
      iprim(15)=47     
      iprim(16)=53 
!    
        do ip=1,nprim
         do idummy=1,niter
            do ik=1,nkp 
              do ir=1,3
                 itest=mod(klist(ik,ir),iprim(ip))
                 if (itest.ne.0) goto 1
            enddo
              enddo

            idivi=idivi*iprim(ip)
            do ik=1,nkp 
              do ir=1,3
               klist(ik,ir)=klist(ik,ir)/iprim(ip)
            enddo
        	  enddo
         enddo
      enddo
1     continue
      idiv=idiv/idivi
      if(idiv.eq.0) idiv=1
      return
        end


!********************************************************************************


      subroutine GBASS(rbas,gbas)                                       
!     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM real SPACE        **
!     **  LATTICE VECTORS OR VICE VERSA                               **

      implicit none
      real*8,intent(in) ::  rbas(3,3)
        real*8,intent(out) :: gbas(3,3)
        real*8 pi,det
        integer i
      pi=4.D0*DATAN(1.D0)                                               
      gbas(1,1)=rbas(2,2)*rbas(3,3)-rbas(3,2)*rbas(2,3)                 
      gbas(2,1)=rbas(3,2)*rbas(1,3)-rbas(1,2)*rbas(3,3)                 
      gbas(3,1)=rbas(1,2)*rbas(2,3)-rbas(2,2)*rbas(1,3)                 
      gbas(1,2)=rbas(2,3)*rbas(3,1)-rbas(3,3)*rbas(2,1)                 
      gbas(2,2)=rbas(3,3)*rbas(1,1)-rbas(1,3)*rbas(3,1)                 
      gbas(3,2)=rbas(1,3)*rbas(2,1)-rbas(2,3)*rbas(1,1)                 
      gbas(1,3)=rbas(2,1)*rbas(3,2)-rbas(3,1)*rbas(2,2)                 
      gbas(2,3)=rbas(3,1)*rbas(1,2)-rbas(1,1)*rbas(3,2)                 
      gbas(3,3)=rbas(1,1)*rbas(2,2)-rbas(2,1)*rbas(1,2)                 
      det=0.D0                                                          
      do i=1,3                                                      
         det=det+gbas(i,1)*rbas(i,1)                                       
      enddo                                                          
      gbas=gbas*2.D0*pi/det                                   
      return                                                            
      end                                                               


!********************************************************************************

      subroutine ORD1(nmax,IX)                                     
!     **  ORD1 ORDERS THE ARRAY IX WITH SIZE nmax ACCORDING TO  INCREASING NUMBER   **
!     **  SUBROUTINES USED:            indexx (from Numerical recipes)               **
      implicit none
        integer,intent(in) :: nmax
        integer,intent(inout) :: ix(4,nmax)
      integer WORK(nmax),index(nmax)
        integer i,i2,i3,i4,ichang
      
!     sort first index of ix
      work=ix(1,:)
      call indexx(work,index,nmax)
      do i=1,nmax
         ix(1,i)=work(index(i))
      enddo
      work=ix(2,:)
      do i=1,nmax
         ix(2,i)=work(index(i))
        enddo
      work=ix(3,:)
      do i=1,nmax
         ix(3,i)=work(index(i))
      enddo
      work=ix(4,:)
      do i=1,nmax
         ix(4,i)=work(index(i))
      enddo
! 
!     sort higher indices with trivial procedure
 501  ichang=0
      do i=2,nmax
      if(ix(1,i).eq.ix(1,i-1)) then
          if(ix(2,i).lt.ix(2,i-1)) then
          i2=ix(2,i-1)
          i3=ix(3,i-1)
          i4=ix(4,i-1)
          ix(2,i-1)=ix(2,i)
          ix(3,i-1)=ix(3,i)
          ix(4,i-1)=ix(4,i)
          ix(2,i)=i2
          ix(3,i)=i3
          ix(4,i)=i4
          ichang=1
          else if(ix(2,i).eq.ix(2,i-1)) then
            if(ix(3,i).lt.ix(3,i-1)) then
            i3=ix(3,i-1)
            i4=ix(4,i-1)
            ix(3,i-1)=ix(3,i)
            ix(4,i-1)=ix(4,i)
            ix(3,i)=i3
            ix(4,i)=i4
            ichang=1
            else if(ix(3,i).eq.ix(3,i-1)) then
              if(ix(4,i).lt.ix(4,i-1)) then
              i4=ix(4,i-1)
              ix(4,i-1)=ix(4,i)
              ix(4,i)=i4
              ichang=1
              end if
            end if
          end if
      end if
      enddo
      if(ichang.eq.1) goto 501
      return                                                            
      end                                                               


!********************************************************************************

      subroutine indexx(arrin,indx,n)
!     sort with heapspot (from Numerical Recipes, p.233)
      implicit none
        integer,intent(in) :: n,arrin(n)
        integer,intent(out) :: indx(n)
        integer q,j,l,ir,indxt,i

      do j=1,n
         indx(j)=j
      enddo
      L=n/2+1
      ir=n

  10  continue
      if(L.gt.1) then
         l=l-1
         indxt=indx(l)
         q=arrin(indxt)            
      else
         indxt=indx(ir)
         q=arrin(indxt)
         indx(ir)=indx(1)
         ir=ir-1
         if(ir.eq.1) then
            indx(1)=indxt
            return
         endif
      endif
      i=l
      j=l+l
  20  if(j.le.ir) then
         if(j.lt.ir) then
            if(arrin(indx(j)).lt.arrin(indx(j+1))) j=j+1
         end if
         if(q.lt.arrin(indx(j))) then
            indx(i)=indx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         goto 20
      end if
      indx(i)=indxt
      goto 10
      end    


!********************************************************************************

      subroutine sdef(iio,nsym,lattic)
!    redefines symmetry operation to include cxz and cyz
!    bravais lattices (from cxy)
!
      implicit none
      character*3,intent(in) :: lattic
        integer,intent(in) ::     nsym
      integer,intent(inout) ::  iio(3,3,48)
      integer ihelp,j

      if(LATTIc(1:3).eq.'CXZ'.or.LATTIc(1:3).eq.'BO ') then
         do j=1,nsym
            ihelp=iio(2,2,j)                                   
            iio(2,2,j)=iio(3,3,j)
            iio(3,3,j)=ihelp
            ihelp=iio(1,2,j)                                   
            iio(1,2,j)=iio(1,3,j)
            iio(1,3,j)=ihelp
            ihelp=iio(2,3,j)                                   
            iio(2,3,j)=iio(3,2,j)
            iio(3,2,j)=ihelp
            ihelp=iio(2,1,j)                                   
            iio(2,1,j)=iio(3,1,j)
            iio(3,1,j)=ihelp
         enddo
      else if(LATTIc(1:3).eq.'CYZ'.or.LATTIc(1:3).eq.'AO ') then
         do j=1,nsym
            ihelp=iio(1,1,j)                                   
            iio(1,1,j)=iio(3,3,j)
            iio(3,3,j)=ihelp
            ihelp=iio(1,2,j)                                   
            iio(1,2,j)=iio(3,2,j)
            iio(3,2,j)=ihelp
            ihelp=iio(1,3,j)                                   
            iio(1,3,j)=iio(3,1,j)
            iio(3,1,j)=ihelp
            ihelp=iio(2,3,j)                                   
            iio(2,3,j)=iio(2,1,j)
            iio(2,1,j)=ihelp
         enddo
      end if
      return
      end


!******************************************************************************

      subroutine sdefl(rbas,gbas,iio,nsym,iz,lattic,ortho)
!    redefines symmetry operation from lapw-struct with unitary transformation  u(-1) * S * u

      implicit none 
      logical,intent(in) :: ortho
      character*3,intent(in) :: lattic                
      integer,intent(out) :: iio(3,3,48)
        integer,intent(in) :: iz(3,3,48),nsym
        real*8,intent(in) :: rbas(3,3),gbas(3,3)
      real*8 gbas1(3,3),pi,a(3,3),b(3,3)
        integer i,j,k,ind,i1,i2,i3,i4
      pi=acos(-1.d0)
    
      do i=1,3
      do j=1,3          
         gbas1(j,i)=gbas(i,j)/2.d0/pi
        enddo
        enddo
 
      do ind=1,nsym 
         a=iz(:,:,ind)
         if(ortho.or.(.not.ortho.and.lattic(1:3).eq.'CXZ')) then
             call matmm(b,rbas,a)
             call matmm(a,b,gbas1)
         end if
         do i=1,3;do j=1,3
             iio(i,j,ind)=nint(a(i,j))
         enddo;enddo
      enddo
!     Write symm matrices iio to file 66 :
      do i=1,int(FLOAT(nsym)/4.+.9)                                 
         I1=4*i-3                                                          
         I2=4*i-2                                                          
         I3=4*i-1                                                          
         I4=4*i                                                            
 !        write (66,120) I1,I2,I3,I4                                        
         do J=1,3                                                      
 !           write (66,140) (Iio(J,K,I1),K=1,3),(IiO(J,K,I2),K=1,3),     &
 !    &         (IiO(J,K,I3),K=1,3),(IiO(J,K,I4),K=1,3)                          
         enddo
        enddo                                                          

  120 format (T5,'SYMMETRY MATRIX NR.',I3,T30,'SYMMETRY MATRIX NR.'     &
     &  ,I3,T55,'SYMMETRY MATRIX NR.',I3,T80,'SYMMETRY MATRIX NR.',I3)   
  140 format (T5,3I5,T30,3I5,T55,3I5,T80,3I5)                           
      return
      end


!******************************************************************************

      subroutine  MATMM (C,A,B)                                        
!     FORM C = A B, WHERE C MAY OVERLAP WITH EITHER A OR B, OR BOTH,    
!     SINCE THE PRODUCT IS DEVELOPED in A TEMPORARY MATRIX.             
      real*8  A(3,3),AB(3,3),B(3,3),C(3,3)    
      AB(1,1) = A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      AB(1,2) = A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      AB(1,3) = A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      AB(2,1) = A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      AB(2,2) = A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      AB(2,3) = A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      AB(3,1) = A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      AB(3,2) = A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      AB(3,3) = A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      C(1,1) = AB(1,1)                                                  
      C(2,1) = AB(2,1)                                                  
      C(3,1) = AB(3,1)                                                  
      C(1,2) = AB(1,2)                                                  
      C(2,2) = AB(2,2)                                                  
      C(3,2) = AB(3,2)                                                  
      C(1,3) = AB(1,3)                                                  
      C(2,3) = AB(2,3)                                                  
      C(3,3) = AB(3,3)                                                  
      return                                                            
      end                                                               


!****************************************************************************** 
 
      subroutine reasym(fn,nsym,sym,sym4)
!   read symmetry operations from wien2k case.struct formatted file

      implicit none
!  INPUT
      character*20,intent(in) :: fn  ! file to use
!  OUTPUT
      integer,intent(out) :: nsym   ! number of symmetry operations
      integer,intent(out) :: sym(3,3,48)  !  rotation matrices
        real*8,intent(out)  :: sym4(3,48)   !  translation vectors
!  LOCALS
        integer j,j1,j2

      open(20,file=fn,form='formatted',status='old')
      read(20,1151) nsym                                                
      do j=1,nsym                                                     
         read(20,1101) ( (sym(J1,J2,J),J1=1,3),sym4(J2,J),J2=1,3 )
      enddo
 
      close(20)
        return
1101  format(3(3I2,F10.7,/)) 
1151  format(I4) 
        end



    
 
 
