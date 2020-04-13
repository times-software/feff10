!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: kprep.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/02/17 07:39:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine kprep(em,ne,prepare_structure_factors)


  use struct,lat=>alat
  use kklist
  use strfacs
  use workstrfacs
  use workstrfacs2
  use workstrfacssimple,only: init_workstrfacssimple
  use boundaries
  use wigner3j
  use trafo
  use controls,only : allocated,irel,singleprec
  use kgenwork
  use energygrid

  implicit none
  integer,intent(in) :: ne ! J Kas - got rid of nex, as it wasn't used.
  complex*16,intent(in) :: em(ne)
  logical,intent(in) :: prepare_structure_factors   !Allows us to exit after allocation of energy grid etc. to save time when we don't really need to venture out into k-space.

  ! LOCALS
  integer l1,mm,isp1,j1,j2,j3p,j3m,m1,m2
  real*8,external :: cwig3j
  logical,parameter :: debug=.false.
  ! For structure factors :
  integer i,iq,il,a,na,ia,b,nb,ib
  integer NK,NKM
  integer,allocatable :: NLQ(:),LTAB(:),KAPTAB(:),NMUETAB(:)
  real*8 ALAT,BOA,COA,VOL,FACT(0:100),etop,gmulti,rmulti,rrr(6)
  real*8,allocatable :: CGC(:,:)
  complex*16,allocatable :: crellocal(:,:),rclocal(:,:),rrellocal(:,:)
 ! complex*16 w1(32,32),w2(32,32)
  complex*16 w1(2*(maxl+1)**2,2*(maxl+1)**2),w2(2*(maxl+1)**2,2*(maxl+1)**2)
  integer savenph !KJ don't like this but am stumped ...

  !     debugging
  integer j,k,ik,l
  real*8 wiensym(3,3,48),dif
  logical addi,rtoc,fromrel

  !      nsym=1
  !	if(usesym.eq.1) nsym=48


  !  debugging
  if (debug) open(101,file='instrfeff.txt')


  !#######################################################################
  !       FIRST:    INITIALIZE DIMENSIONS 
  !#######################################################################
  maxl=0
  do i=0,nph
     maxl=max(maxl,lpot(i))
  enddo

  nl=maxl+1
  msize=nsp*nats*(maxl+1)**2
  mls=nsp*(maxl+1)**2

!Note that a fair amount of allocations was moved to strvecgen (workstrfacs,workstrfacs2 - 2/2012) KJ
  if (.not.allocated) then
     call init_strfacs  !KJ 11-2011 : this just zeroes streta,strrmax,strgmax,eimag - after they've been read from file! || update: changed init_strfacs
     call init_boundaries(maxl,nats)
     !	   call init_gk(msize,nkp)
     call init_workstrfacs(nats)
     call init_workstrfacs2
     call init_wigner3j(maxl)
     call init_trafo(mls)
     call init_energygrid(ne)
     allocated=.true.
  endif


  !* PROCESS energy grid
  egrid=em(1:ne)
  emin=dcmplx(500,0)
  emax=-emin
  do i=1,ne
     if(dble(emin-em(i)).gt.0) emin=em(i)
     if(dble(em(i)-emax).gt.0) emax=em(i)
  enddo
  if(debug) write(6,*) 'emin,emax',emin,emax


if (.not.prepare_structure_factors) return  !KJ 2-2012


  !* STRUCTURE
  ALAT = lat(1)
  BOA = lat(2)/lat(1)
  COA = lat(3)/lat(1)

  brx=bramat(:,1)
  bry=bramat(:,2)
  brz=bramat(:,3)

  brx(1)=a1(1)
  brx(2)=a2(1)
  brx(3)=a3(1)
  bry(1)=a1(2)
  bry(2)=a2(2)
  bry(3)=a3(2)
  brz(1)=a1(3)
  brz(2)=a2(3)
  brz(3)=a3(3)
  brx=brx /alat
  bry=bry /alat
  brz=brz /alat


  qx(1:nats)=ppos(1,1:nats)  ! I believe these need to be in CARTHESIAN fractional coordinates
  qy(1:nats)=ppos(2,1:nats)
  qz(1:nats)=ppos(3,1:nats)
  !* CALCULATION
  ! In SPRKKR, irel takes values from 0 to 4.  For the routines used here, only > or < 2 matters.
  if (nsp.le.1) then
     irel=1
  elseif (nsp.gt.1) then
     irel=3
  endif
  !  Next variables allow to expand the clusters used for the calculation of the structure constants
  !  in real (rmulti) and reciprocal (gmulti) space.
  rmulti = dble(1)
  gmulti = dble(1)



  ! allocate some locals
  allocate(NLQ(NQMAX),LTAB(NMUEMAX),KAPTAB(NMUEMAX),NMUETAB(NMUEMAX),CGC(NKMPMAX,2) )
  allocate(crellocal(nkmmax,NKMMAX),rclocal(NKMMAX,NKMMAX),rrellocal(NKMMAX,NKMMAX))


  !*************************************************************************************
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

  ! Calculate the sum of the k-mesh integration weights to normalize integrals :
  sumweights=dble(0)
  do i=1,nkp
     sumweights=sumweights+weight(i)
  enddo



  if (debug) then
     write(101,*) 'bgx',bgx
     write(101,*) 'bgy',bgy
     write(101,*) 'bgz',bgz
     write(101,*) 'vol',vol
  endif


  FACT(0) = 1.0D0
  DO I=1,100
     FACT(I) = FACT(I-1) * DBLE(I)
  END DO
  NLQ(1:nqmax) = 0



  !************************** For FEFF's t-matrix (in fms2) :
  !     Calculate Clebsch-Gordon coefficients <LS|J>
  do l1 = 0, maxl  !KJ used to be lx instead of maxl 11-06
     do mm = -l1, l1
        do isp1 = 1, 2
           j1 = 2 * l1
           j2 = 1
           j3p = j1 + 1
           j3m = j1 - 1
           m1 = 2*mm
           m2 = 2*isp1 - 3
           !  j = l+1/2
           t3jp( l1, mm, isp1) = sqrt( j3p + 1.0e0 ) *    real( cwig3j( j1, j2, j3p, m1, m2, 2) )
           if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)     t3jp( l1, mm, isp1) = - t3jp( l1, mm, isp1)
           !  j = l-1/2
           t3jm( l1, mm, isp1) = sqrt( j3m + 1.0e0 ) *   real( cwig3j( j1, j2, j3m, m1, m2, 2) )
           if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)    t3jm( l1, mm, isp1) = - t3jm( l1, mm, isp1)
        enddo
     enddo
  enddo

  !*************************



  !   ********************************************************************
  !   *                                                                  *
  !   *   rel. quantum numbers    up to                                  *
  !   *                                                                  *
  !   ********************************************************************
  !                     s   p   p   d   d   f   f   g   g   h
  !     DATA LTAB    /  0,  1,  1,  2,  2,  3,  3,  4,  4,  5 /
  !     DATA LBTAB   /  1,  0,  2,  1,  3,  2,  4,  3,  5,  4 /
  !     DATA KAPTAB  / -1,  1, -2,  2, -3,  3, -4,  4, -5,  5 /
  !     DATA NMUETAB /  2,  2,  4,  4,  6,  6,  8,  8, 10, 10 /

  DO I = 1,NMUEMAX
     LTAB(I) = I/2
     IF( 2*LTAB(I) .EQ. I ) THEN
        KAPTAB(I) = LTAB(I)
     ELSE
        KAPTAB(I) = - LTAB(I) - 1
     END IF
     NMUETAB(I) = 2*ABS( KAPTAB(I) )
  END DO

  if (debug) then
     write(101,*) 'ltab',ltab
     write(101,*) 'kaptab',kaptab
     write(101,*) 'nmuetab',nmuetab
  endif


  ! calculate some gaunt symbols before ltab and nmuetab get overwritten for non-sprel calcul.
  CALL CALCCGC( LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX )

  !=======================================================================
  !                initialize tables of quantum numbers
  !=======================================================================
  NLM  = NL**2
  ! Fix site specific momentum cutoffs
  DO IQ=1,nats
     NLQ(IQ) = NL   !KJ for now, here all sites get the same l-cutoff.  We could change this using lipotx or lpot.
     NKMQ(IQ) = 2*NLQ(IQ)**2
  ENDDO

  if(debug) then
     write(101,*) 'nlm',nlm
     write(101,*) 'nlq',nlq
     write(101,*) 'nkmq',nkmq
  endif

  IF( IREL .LT. 2 ) THEN
     !
     !-------------------------------------- non-relativistic kkr-calculation

     DO 220 IQ=1,nats
        IF( IQ .EQ. 1 ) THEN
           IND0Q(IQ) = 0
        ELSE
           IND0Q(IQ) = IND0Q(IQ-1) + NLQ(IQ-1)**2
        END IF
        NKMQ(IQ)  = NLQ(IQ)**2
        DO 220 IL=1,NLMAX
           LTAB(IL)    = IL - 1
           NMUETAB(IL) = 2*LTAB(IL)+1
220        CONTINUE

   ELSE
           !------------------------------------------ relativistic kkr-calculation
           NK   = 2*NL-1
           NKM  = 2*NLM
           WRITE(6,108)
           I = 0
           DO 230 IQ=1,nats
              IF( IQ .EQ. 1 ) THEN
                 IND0Q(IQ) = 0
              ELSE
                 IND0Q(IQ) = IND0Q(IQ-1) + 2*NLQ(IQ-1)**2
              END IF
230           CONTINUE
   END IF

           !=======================================================================
!           WRITE(6,8)  ALAT,nats,NL,NLM,NK,NKM
8          FORMAT(//,10X,'ALAT =',F10.6,//,                                  &
                &  10X,'nats =',I3,5X,'NL =',I3,5X,'NLM =',I4,5X,/,                &
                & 22X,'NK =',I3,5X,'NKM =',I4,5X,//)


           !sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
           !                         structure constants
           !sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
		   !next parameters communicated to strinit via module workstrfacs(2)
           ETA= 0D0
           RMAX=0D0 !rmax and gmax will be set in strinit 
           GMAX=0D0
           !If eta=0, an algorithm in strinit will find a value automatically.  However, if the user set a value in feff.inp, use it here: 
		   ETA=streta

           etop=dble(0)
           do i=1,ne
              if(dabs(dble(em(i))).gt.etop)  etop=dabs(dble(em(i)))
           enddo
           etop=etop*2  ! SPRKKR works in Ry, FEFF in Ha

           CALL STRINIT(ETOP,nats,ALAT,FACT,CGC,gmulti,rmulti ) !KJ eliminating arguments that I think are unnecessary
!           CALL STRINIT( ETA, RMAX, GMAX,ETOP,                               &
!                &                 BRX,BRY,BRZ,R1,R2,R3,NRDLTAB,QX,QY,QZ,           &
!                &                 BGX,BGY,BGZ,G1,G2,G3,NGRLTAB,                    &
!                &                 nats,ALAT,FACT,CGC,gmulti,rmulti )

           !     Reset icall - strcc won't actually calculate anything
           call STRCC(dcmplx(0),ALAT,.true.)  ! 1st argument etop chosen arbitrarily - doesn't matter if 3rd arg. = true

           if(debug) then
              write(101,*) 'eta ',eta
              write(101,*) 'rmax ',rmax
              write(101,*) 'gmax ',gmax
              write(101,*) 'etop ',etop
              write(101,*) 'qx',qx
              write(101,*) 'qy',qy
              write(101,*) 'qz',qz
              write(101,*) 'nats',nats
              write(101,*) 'cgc',cgc

           endif

           call bastrmat(maxl,CGC,rclocal,crellocal,rrellocal,nkmmax,nkmpmax)
           rc=rclocal(1:mls,1:mls)
           crel=crellocal(1:mls,1:mls)
           rrel=rrellocal(1:mls,1:mls)



           drot=dcmplx(0,0)
           savenph=nph
           call makerotations
           nph=savenph
           addi=.false.
           rtoc=.false.
           fromrel=.false.
           do i = 1,48
              w1=dcmplx(0,0)
              if(fromrel) then
                 w2=drot(:,:,i,1)
                 CALL CHANGEREP(w2,'REL>CLM',W1,nkmmax,nkmmax,RClocal,CRELlocal,RRELlocal,'          ',0)
                 w2=w1
              else
                 w2=drot(:,:,i,2)
              endif
              !        Convert from SPRKKR basis to FEFF basis
              !        For nsp=1 : from real to spherical harmonics
              if(rtoc) then
                 call changerep(w2,'RLM>CLM',w1,nkmmax,nkmmax,rclocal,crellocal,rrellocal,'          ',0)	   
              else
                 w1=w2
              endif
              if(addi) then
                 !        To get the normalization of the spherical harmonics right, add factor i^l :
                 ia=0;ib=0
                 do a=0,maxl;do na=-a,a;ia=ia+1;ib=0
                    do b=0,maxl;do nb=-b,b;ib=ib+1
                       w1(ia,ib)=w1(ia,ib)*dcmplx(0,1)**(a-b)
                    enddo;enddo
                 enddo;enddo
              endif



              drot(:,:,i,2)=w1
           enddo

           !!KJ debugging
           !      open(78,file='symops_working.txt')
           !	do i=1,48
           !	write(78,*) i
           !	do j=1,nkmmax
           !	write(78,'(100f10.4)') drot(j,:,i,2)
           !      enddo
           !	enddo
           !	close(78)
           !!KJ



!           write(*,*) 'mls= ',mls
           allocate(mrot(mls,mls,48))
           if(nsp.le.1) then
              mrot=drot(1:mls,1:mls,1:48,2)
           else
              mrot=drot(1:mls,1:mls,1:48,1)
           endif




           ! kill some locals
           deallocate(NLQ,LTAB,KAPTAB,NMUETAB,CGC,crellocal,rrellocal,rclocal)

           !     following section only works if you've already run kgen
           if(usesym.eq.1) then
              open(55,file='symops_kgen.txt',form='formatted')
!              write(*,*) 'possibly need transpose of cm ...'
              do i=1,4
                 READ(55,*)
              enddo
!              write(*,*) 'nsym in kprep',nsym
              do i=1,nsym
                 read(55,*)
                 do j=1,3
                    read(55,*) j1,j2,j3p,rrr(1:6)
                    wiensym(j,1,i)=rrr(4)
                    wiensym(j,2,i)=rrr(5)
                    wiensym(j,3,i)=rrr(6)
                    !	         read(55,1019) wiensym(j,:,i)
                 enddo
                 do k=1,48
                    dif=dble(0)
                    do j=1,3
                       do l=1,3
                          dif=dif+dabs(wiensym(j,l,i)-mrotr(j,l,k)) !mrotr(l,j,k))
                          !	            dif=dif+dabs(wiensym(j,l,i)-mrotr(l,j,k)) !mrotr(l,j,k))
                       enddo
                    enddo
                    if(dif.lt.0.0001) then
                       symid(1,k)=i  ! sprsym(:,:,k)=wiensym(:,:,symid(1,k))
                       symid(2,i)=k  ! wiensym(i)=sprsym(symid(2,i))
                    endif
                 enddo
              enddo
              close(55)
1019          format(18x,5x,3f12.6)

              if(debug) then
              open(55,file='symidcheck.txt')
              do i=1,nsym !48
                 write(55,*) i,symid(2,i),'  wiensym - mrotr'
                 do j=1,3
                    write(55,1020) wiensym(j,:,i),mrotr(j,:,symid(2,i))
                 enddo
              enddo
              close(55)
1020          format (3f12.6,5x,3f12.6)
              endif
           endif !usesym=1


           !    Process the k-mesh :
           if (usesym.eq.1) then
!              write(*,*) 'nkp,nka,nki,nkf',nkp,nka,nki,nkf
              !    Prepare a table that will help carry out the integral
              inti=0
              intn=0
              intw=dble(0)
              symact=0
              do ik=1,nkp
                 do j=1,nkf
                    if(linkf(j).eq.ik)then
                       intn(ik)=intn(ik)+1
                       inti(ik,intn(ik),1)=j
                       inti(ik,intn(ik),2)=symid(2,lsymf(j))
                       symact(symid(2,lsymf(j)))=1
                       intw(ik,intn(ik))=wf(j)
                    endif
                 enddo
              enddo
              if(debug) then
              open(55,file='matchk.txt')
              do ik=1,nkp
                 write(55,*) 'K-POINT ',ik
                 do j=1,intn(ik)
                    write(55,*) inti(ik,j,1),inti(ik,j,2),intw(ik,j)
                 enddo
              enddo
              close(55)
              endif
           endif



        if(debug) then
           open(44,file='bkf.txt')
           do i=1,nkf
              write(44,1818) i,bkf(:,i),wf(i),linkf(i),lsymf(i)
           enddo
           close(44)
           open(44,file='bkw.txt')
           do i=1,nkw
              write(44,1818) i,bkw(:,i),ww(i),linkw(i),lsymw(i)
           enddo
           close(44)
           open(44,file='bki.txt')
           do i=1,nki
              write(44,1818) i,bki(:,i),wi(i)
           enddo
           close(44)
1818       format(i5,4f14.4,2i5)
         endif


           ! kill some variables of the k-mesh that are no longer needed
           deallocate(bki,bki2,wi,wf,linkf) !lsymf
           if(usesym.ne.1) deallocate(bkf)
           !     if(debug) stop


           if(singleprec) call init_workstrfacssimple

           return
108        FORMAT(1H ,/,                                                     &
                &10X,'RELATIVISTIC QUANTUM NUMBERS',/,                             &
                &10X,'============================',//,                            &
                &10X,' I  L  KAP  J   MUE    I',                                   &
                &    '  GDIA   GOFF   GMDIA  GMOFF',                               &
                &    '  FDIA   FOFF   FMDIA  FMOFF ')

             end subroutine kprep

