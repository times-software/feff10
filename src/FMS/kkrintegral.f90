!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: kkrintegral.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/03/16 22:55:02 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine kkrintegral(G,Tmat,Cmat,offdiagonal,rij,p,justone)

  !  * create Gamma(k) * T (matrices of dimension (lx+1)**1 * natlat
  !  * invert (1-Gamma(k)T) for all k in the array bk
  !  * take an integral over the Brillouin Zone of    exp(-i k . rij)  (1-Gamma(k)T)^-1  Gamma(k)
  !  * if there is a core hole atom in the cell, make the difference t-matrix and construct the new G_corehole

  use struct,only: nats,alat,absorber,ppot
  use kklist
  use boundaries,only: maxl,msize,mls
  use strfacs,only: eimag
  use kgenwork ,only : bkf
  use trafo
  use controls,only : corehole

  implicit none

  ! INPUT
  complex,intent(in) ::  Tmat(msize,msize)   ! scattering matrix
  complex,intent(inout) ::  Cmat(mls,mls)       ! Core hole scattering matrix  (note : difference between core hole and regular potential)
                                                ! If (justone), information will be passed back in this matrix
  logical,intent(in) ::  offdiagonal         ! include factor exp(-ik.r) or not
  real*8,intent(in)  ::  rij(3)              ! the r in exp(-ik.r)
  complex*16,intent(in) :: p !energy at which we work - one point, not a mesh
  logical,intent(in) ::  justone             ! If =true, means there is only 1 position corresponding to the type of the excited atom.  This requires extra care.
  ! OUTPUT
  complex,intent(out) :: G(msize,msize)  ! The fms-matrix G_fms - also used as a work array in this routine
  ! LOCAL
  complex Mmat(msize,msize),Gfull(msize,msize),TT(msize,msize),SS(msize,msize)
  integer info,i,j,ik,ipiv(msize)
  character*13 trans
  complex expfact
  logical bandstructure,transfo,useusesym
  integer,parameter :: method=1  ! 1 to calculate G_fms ; 0 to calculate tau_fms
  integer n,m
  integer k,j1,j2,jj1,jj2,i1
  real*8 bvec(3),bvecl(3),bvecl2(3),dif,p2a(3)
  complex v1(mls,mls),v2(mls,mls),v3(mls,mls),v4(mls,mls),vjustone(2,mls,mls)
  real*8,parameter :: pi = 3.1415926535897932384626433d0
  logical,parameter :: useuseusesym=.true. !false.  ! T to get the whole G_fms right ; F to only get Tr G_fms right (much faster)
  complex,parameter :: ie=(0.0,1.0),c0=(0.0,0.0),c1=(1.0,0.0)  ! careful : these are defined as complex*16 in other routines!
  integer,parameter :: coreholemethod=1  ! must be 1,2,3 or 4
  logical,parameter :: debug=.false.
  logical :: usevjustone


  if (justone .and. offdiagonal) call wlog('Fishy in kkrintegral : justone and offdiagonal both use rij!')
  usevjustone= justone .and. corehole .and. (coreholemethod.eq.1 .or. coreholemethod.eq.4)  ! In these cases, it is necessary to construct the matrix vjustone
  if(debug) call writematrix(Tmat,msize,'Tmat.txt')
  if(debug) call writematrix(Cmat,mls,'Cmat.txt')

  !***************** PREPARATION **************************************
  !   Process input and initialize variables.
  Gfull=cmplx(0,0)
  G=cmplx(0,0)
  if(usevjustone) vjustone=cmplx(0,0)
  trans = 'NotTransposed'
  bandstructure=.false.
  transfo=.false.
  useusesym=useuseusesym.or.offdiagonal
  n=mls
  m=msize
  do i=1,3
     p2a(i)=dble(2)*pi/alat(i)  ! 2 pi / a
  enddo
  expfact=c1


  if(method.eq.1) then  !FEFF calculation




     !******************************** START BZ INTEGRATION ******************************************

     if(debug) open(55,file='dededebug.txt')

     do ik=1,nkp  !summing over all irreducible k-vectors
!	 write(*,*) 'starting k-point #',ik  !for debugging
!write(*,*) 'p',p
!write(*,*) 'bk(:,ik)',bk(:,ik)
!call writematrix(G,msize,'G_1.txt')

        !        prepare Gamma(k)
        call structurefactor(p,bk(:,ik),G) !note : G contains (-1) * structure factor
        if(debug) call writematrix(G,msize,'G.txt')

        call unity(TT,msize)
        call cgemm('N','N',m,m,m,-C1,G,m,Tmat,m,C1,TT,m)   ! TT = 1 - (-G) (-t)  !!!KJ here I'm messing with FEFF

        !        matrix inversion of TT=(1-Gamma(k)*T)  into SS=((1-Gamma(k)*T)^-1
        call cgetrf(msize,msize,TT,msize,ipiv,info)
        if (info.ne.0) then
           call wlog('*** Error in cgetrf when processing 1-Gamma(k)T')
           call par_stop('Error in cgetrf when processing 1-Gamma(k)T')
        endif
        call unity(SS,msize)
        call cgetrs(trans,msize,msize,TT,msize,ipiv,SS,msize,info)
        if (info.lt.0) then
           call par_stop('*** Error in cgetrs when computing (1-Gamma(k)T)^-1')
        endif

        call cgemm('N','N',m,m,m,C1,SS,m,G,m,C0,TT,m)   ! TT = SS * (-G)  !!!c1 used to be -c1 - messing with FEFF again!

        bvec=bk(:,ik)

        if(usesym.eq.1.and.useusesym) then
           do i = 1,intn(ik)
              if(offdiagonal) then
                 bvecl2=dble(0)
                 do j=1,3
                    do k=1,3
                       bvecl2(j)=bvecl2(j)+mrotr(k,j,inti(ik,i,2))*bvec(k)
                       !    The transpose of a unitary matrix is its inverse ; this is the transformation from irreducible point to its image,
                       !    whereas the mrotr matrices store the opposite transformation that brings the image to its representative point.
                    enddo
                    bvecl2(j)=mod(bvecl2(j),p2a(j))
                    if(bvecl2(j).lt.dble(0)) bvecl2(j)=bvecl2(j)+p2a(j)
                 enddo
                 do j=1,3
                    bvecl(j)=bkf(j,inti(ik,i,1))  !* p2a(j)  !from fractional to atomic units
                 enddo
                 dif=dble(0)
                 do j=1,3
                    dif=dif+dabs(bvecl(j)-bvecl2(j))
                 enddo
                 if(dif.gt.0.00001) then
                    write(*,1021) ik,i,bvecl,inti(ik,i,2),bvecl2
                    write(*,*) (mrotr(j,:,inti(ik,i,2)),j=1,3)
1021                format(2i5,3f12.4,i5,3f12.4)
                    stop
                 endif
                 expfact= exp(-ie*real((bvecl(1)*rij(1)+ bvecl(2)*rij(2)+bvecl(3)*rij(3))))
              endif
              do j1=1,nats
                 do j2=1,nats
                    jj1=(j1-1)*mls
                    jj2=(j2-1)*mls
                    v2=TT(1+jj1:mls+jj1,1+jj2:mls+jj2)
                    !                 sprsym(:,:,k)=wiensym(:,:,symid(1,k))
                    !                 wiensym(i)=sprsym(symid(2,i))
                    v4=cmplx(mrot(1:mls,1:mls,inti(ik,i,2)))
                    call cgemm('C','N',n,n,n,C1,v4,n,v2,n,C0,v1,n)
                    call cgemm('N','N',n,n,n,C1,v1,n,v4,n,C0,v3,n)
                    SS(1+jj1:mls+jj1,1+jj2:mls+jj2)=v3
                 enddo
              enddo
              Gfull=Gfull+SS*real(intw(ik,i))*expfact
           enddo
        else   ! No symmetry :
           if(offdiagonal) expfact=exp(-ie*real((bvec(1)*rij(1)+  bvec(2)*rij(2)+bvec(3)*rij(3))))

           Gfull=Gfull+TT*real(weight(ik))*expfact
		   !KJ 03-2012 next line calculates G_c,c' where c is the absorber position and c' is the same but removed by 1x the first lattice vector:
		   if (usevjustone) then
		      vjustone(1,:,:)=vjustone(1,:,:)+ &
                        TT(1+((absorber-1)*mls):(absorber*mls),1+ &
                        ((absorber-1)*mls):(absorber*mls)) *real(weight(ik)) &
                         * exp(-ie*real((bvec(1)*rij(1)+  bvec(2)*rij(2)+bvec(3)*rij(3))))
		      vjustone(2,:,:)=vjustone(2,:,:)+ &
                         TT(1+((absorber-1)*mls):(absorber*mls),1+ &
                         ((absorber-1)*mls):(absorber*mls)) *real(weight(ik)) &
                         * exp(ie*real((bvec(1)*rij(1)+  bvec(2)*rij(2)+bvec(3)*rij(3))))
           endif
        endif !usesym
     enddo
     Gfull=Gfull/sumweights

     if(debug) call writematrix(Gfull,msize,'Gfull.txt')
     !****************************** BZ INTEGRATION IS FINISHED *******************************

     !    Deal with the core hole :
     if (.not.corehole) then
        G=Gfull
		if (justone) then
		   !just copy the "core hole" to the "neutral atom":
		   jj1=(absorber-1)*mls
		   Cmat=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)
		endif 
     else

        if(coreholemethod.eq.1) then
           !             Correct the whole G-matrix for the presence of the core hole	   	   

           !           Prepare   v1 := (1+ t_ch G_fms_ch,ch)^-1 t_ch
           i=(absorber-1)*mls
           v1=Gfull(1+i:mls+i,1+i:mls+i)
           call unity(v2,mls)
           call cgemm('N','N',n,n,n,c1,Cmat,n,v1,n,c0,v3,n)
           v1=v2+v3
           call cgetrf(n,n,v1,n,ipiv,info)
           if (info.ne.0) then
              call par_stop('Error in cgetrf when processing v1')
           endif
           call cgetrs(trans,n,n,v1,n,ipiv,v2,n,info)
           if (info.lt.0) then
              call par_stop('*** Error in cgetrs for (v1)^-1')
           endif
           call cgemm('N','N',n,n,n,c1,v2,n,Cmat,n,c0,v1,n)

           !       Reconstruct the fms-matrix with core hole impurity :
           do j1=1,nats
              do j2=1,nats
                 jj1=(j1-1)*mls
                 jj2=(j2-1)*mls
                 v3=Gfull(1+i:mls+i,1+jj2:mls+jj2)      ! v3 = G_ch,j2
                 call cgemm('N','N',n,n,n,c1,v1,n,v3,n,c0,v4,n)  ! v4 = v1 G_ch,j2
                 v3=Gfull(1+jj1:mls+jj1,i+1:i+mls)      ! v3 = G_j1,ch
                 call cgemm('N','N',n,n,n,c1,v3,n,v4,n,c0,v2,n)  ! v2 = G_j1,ch  v1  G_ch,j2
                 v3=Gfull(1+jj1:mls+jj1,1+jj2:mls+jj2)  ! v3 = G_j1,j2
                 G(1+jj1:mls+jj1,1+jj2:mls+jj2)=v3-v2
              enddo
           enddo
		   
		   if (justone) then
           ! G_c',c' =   G_c',c'   +   G_c',c   (1 - t_c G_c,c)^-1   t_c   G_c,c'
                 v3=vjustone(1,:,:)      ! v3 = G_ch,c'
                 call cgemm('N','N',n,n,n,c1,v1,n,v3,n,c0,v4,n)  ! v4 = v1 G_ch,c'
                 v3=vjustone(2,:,:)      ! v3 = G_c',ch
                 call cgemm('N','N',n,n,n,c1,v3,n,v4,n,c0,v2,n)  ! v2 = G_c',ch  v1  G_ch,c'
		         jj1=(absorber-1)*mls
		         v3=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)  ! v3 = G_c',c'
                 Cmat=v3-v2		   
		   endif


        elseif(coreholemethod.eq.2) then
           !           ALTERNATIVE : "lazy code" : fix only ch,ch block!
		      
           G=Gfull
           !           Prepare   v1 := (1+ t_ch G_fms_ch,ch)^-1
           i=(absorber-1)*mls
           v1=Gfull(1+i:mls+i,1+i:mls+i)
           call unity(v2,mls)
           call cgemm('N','N',n,n,n,c1,Cmat,n,v1,n,c0,v3,n)
           v4=v2+v3
           call cgetrf(n,n,v4,n,ipiv,info)
           if (info.ne.0) then
              call par_stop('Error in cgetrf when processing v1')
           endif
           call cgetrs(trans,n,n,v4,n,ipiv,v2,n,info)
           if (info.lt.0) then
              call par_stop('*** Error in cgetrs for (v4)^-1')
           endif
           call cgemm('N','N',n,n,n,c1,v1,n,v2,n,c0,v4,n)
           G(1+i:mls+i,1+i:mls+i)=v4
           if (justone) then
		      !just copy the "core hole" to the "neutral atom":
		      jj1=(absorber-1)*mls
		      Cmat=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)
		   endif 
           !           NOTE THAT this code returns with only the ch,ch block correct - all the rest is ground-state!!	   
		   

        elseif(coreholemethod.eq.3) then    
           !           ALTERNATIVE : "efficient code" : sum over scatterers at this level ...
		      
           G=Gfull
           do i1=1,nats
              if(ppot(i1).eq.ppot(absorber)) then
                 i=(i1-1)*mls
                 v1=Gfull(1+i:mls+i,1+i:mls+i)
                 call unity(v2,mls)
                 call cgemm('N','N',n,n,n,c1,Cmat,n,v1,n,c0,v3,n)
                 v4=v2+v3
                 call cgetrf(n,n,v4,n,ipiv,info)
                 if (info.ne.0) then
                    call par_stop('Error in cgetrf when processing v1')
                 endif
                 call cgetrs(trans,n,n,v4,n,ipiv,v2,n,info)
                 if (info.lt.0) then
                    call par_stop('*** Error in cgetrs for (v4)^-1')
                 endif
                 call cgemm('N','N',n,n,n,c1,v1,n,v2,n,c0,v4,n)
                 G(1+i:mls+i,1+i:mls+i)=v4
              endif
           enddo
           if (justone) then
		      !In this case, coreholemethod=3 is the same as coreholemethod=2
		      !just copy the "core hole" to the "neutral atom":
		      jj1=(absorber-1)*mls
		      Cmat=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)
		   endif 
           !           NOTE THAT this code returns with only the diagonal blocks of the atoms with the same potential type
           !                  as the absorber correct - all the rest is ground-state!!	   

        elseif(coreholemethod.eq.4) then
           !             Correct the diagonal of the G-matrix for the presence of the core hole	   	   

           !           Prepare   v1 := (1+ t_ch G_fms_ch,ch)^-1 t_ch
           i=(absorber-1)*mls
           v1=Gfull(1+i:mls+i,1+i:mls+i)
           call unity(v2,mls)
           call cgemm('N','N',n,n,n,c1,Cmat,n,v1,n,c0,v3,n)
           v1=v2+v3
           call cgetrf(n,n,v1,n,ipiv,info)
           if (info.ne.0) then
              call par_stop('Error in cgetrf when processing v1')
           endif
           call cgetrs(trans,n,n,v1,n,ipiv,v2,n,info)
           if (info.lt.0) then
              call par_stop('*** Error in cgetrs for (v1)^-1')
           endif
           call cgemm('N','N',n,n,n,c1,v2,n,Cmat,n,c0,v1,n)

           !       Reconstruct the diagonal of the fms-matrix with core hole impurity :
           do j1=1,nats
              jj1=(j1-1)*mls
              v3=Gfull(1+i:mls+i,1+jj1:mls+jj1)      ! v3 = G_ch,j1
              call cgemm('N','N',n,n,n,c1,v1,n,v3,n,c0,v4,n)  ! v4 = v1 G_ch,j1
              v3=Gfull(1+jj1:mls+jj1,i+1:i+mls)      ! v3 = G_j1,ch
              call cgemm('N','N',n,n,n,c1,v3,n,v4,n,c0,v2,n)  ! v2 = G_j1,ch  v1  G_ch,j1
              v3=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)  ! v3 = G_j1,j1
              G(1+jj1:mls+jj1,1+jj1:mls+jj1)=v3-v2
           enddo

		   if (justone) then
           ! G_c',c' =   G_c',c'   +   G_c',c   (1 - t_c G_c,c)^-1   t_c   G_c,c'
                 v3=vjustone(1,:,:)      ! v3 = G_ch,c'
                 call cgemm('N','N',n,n,n,c1,v1,n,v3,n,c0,v4,n)  ! v4 = v1 G_ch,c'
                 v3=vjustone(2,:,:)      ! v3 = G_c',ch
                 call cgemm('N','N',n,n,n,c1,v3,n,v4,n,c0,v2,n)  ! v2 = G_c',ch  v1  G_ch,c'
		         jj1=(absorber-1)*mls
		         v3=Gfull(1+jj1:mls+jj1,1+jj1:mls+jj1)  ! v3 = G_c',c'
                 Cmat=v3-v2		   
		   endif

        else
           stop 'invalid core hole method in recinv2'
        endif  ! which core hole method	          


     endif  ! is there a core hole?




  elseif (method.eq.0) then  !SPRKKR calculation :

     !********************************** SET UP MATRICES *******************************
     if(nats.gt.1) stop 'this part not yet adapted'
     !       matrix inversion of Tmat  into Mmat
     call cgetrf(msize,msize,Tmat,msize,ipiv,info)
     if (info.ne.0) then
        call wlog('*** Error in cgetrf when processing T')
        call par_stop('Error in cgetrf when processing T')
     endif
     call unity(Mmat,msize)
     call cgetrs(trans,msize,msize,Tmat,msize,ipiv,Mmat,msize,info)
     if (info.lt.0) then
        call par_stop('*** Error in cgetrs when computing T^-1')
     endif


     !******************************** START BZ INTEGRATION ******************************************

     do ik=1,nkp  !summing over all irreducible k-vectors


        !        prepare Gamma(k)
        !!         G=gk(:,:,ik)  !note : G contains (-1) * structure factor
        call structurefactor(p,bk(:,ik),G)

        !        now add Mmat        
        TT = G + Mmat

        !        matrix inversion of TT=((T^-1)-Gamma(k))  into SS=((T^-1)-Gamma(k))^-1
        call cgetrf(msize,msize,TT,msize,ipiv,info)
        if (info.ne.0) then
           call wlog('*** Error in cgetrf when processing 1-Gamma(k)T')
           call par_stop('Error in cgetrf when processing 1-Gamma(k)T')
        endif
        call unity(SS,msize)
        call cgetrs(trans,msize,msize,TT,msize,ipiv,SS,msize,info)
        if (info.lt.0) then
           call par_stop('*** Error in cgetrs when computing (1-Gamma(k)T)^-1')
        endif

        !        now simply add to perform the integral.
        Gfull=Gfull+real(weight(ik))*SS
        if(bandstructure.and.(bk(3,ik).ne.0)) write(77,'(7e13.6)')  bk(3,ik),imag(SS(1,1))

     enddo    !sum over all irreducible points - integral of BZ


     !     Now restore the full BZ from the IBZ by rotating the results :
     if(usesym.eq.1) then
        SS=Gfull  !use for temporary storage
        Gfull=cmplx(0,0)
        do i = 1,48
           n=mls
           v4(:,:) = cmplx(mrot(:,:,i))
           ! JK - Changed mrot to v4 below because of type discrepancy. 09/09
           call cgemm('N','N',N,N,N,C1,v4,n,SS,n,C0,TT,n) !DROT(1,1,i) !cgemm for single precision, zgemm for double
           !call cgemm('N','N',N,N,N,C1,mrot(1,1,i),n,SS,n,C0,TT,n) !DROT(1,1,i) !cgemm for single precision, zgemm for double
           call cgemm('N','C',N,N,N,C1,TT,n,v4,n,C1,Gfull,n)
           !call cgemm('N','C',N,N,N,C1,TT,n,mrot(1,1,i),n,C1,Gfull,n)
        enddo
        Gfull=Gfull/dble(48)
     endif

     Gfull=Gfull/sumweights

     !****************************** BZ INTEGRATION IS FINISHED *******************************


     if (transfo) then
        !        transform from tau to G_fms formalism
        !        below, Gfull = Tmat^-1 * Gfull * Tmat^-1
        !        the Gfull on the LHS is G_fms ; the one on the RHS is tau.
        call cgemm('n','n',msize,msize,msize,c1,Gfull,msize, Mmat,msize,c0,SS,msize)
        call cgemm('n','n',msize,msize,msize,c1,Mmat,msize,SS,msize,c0,Gfull,msize)
     endif


     !    Deal with the core hole :
     if (.not.corehole) then
        G=Gfull
     else
        stop 'core hole not yet implemented in recinv2'
     endif







  else
     stop 'unknown method in recinv2'
  endif


  if(debug) stop

  return
end subroutine kkrintegral  ! subroutine recinv





!**********************************************************************************************************
      subroutine unity(a,n)
! Initialize a to the unit matrix of dimension n.
        implicit none
        integer n,i
        complex a(n,n)
        complex,parameter :: nul=(0.0,0.0),een=(1.0,0.0)

      a=nul
        do i=1,n
           a(i,i)=een
        enddo

        return
        end   ! unity


