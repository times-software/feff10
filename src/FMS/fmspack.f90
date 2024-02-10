!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmspack.f90,v $:
! $Revision: 1.19 $
! $Author: jorissen $
! $Date: 2011/12/11 00:37:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fms(lfms, nsp, ispin, inclus, npot, ck, lipotx, xphase,&
     &   ik, iverb, minvin, rdirec, toler1, toler2, lcalc, gg)
  use DimsMod, only: lx, nspx=>nspu, nphx=>nphu, nclusx, istatx, nphasx
  use controls,only : sprkkrpot,fullpot

  use stkets
  use rotx
  use lnlm
  use xstruc
  use t3j
  use par
  use fms_mod, only: gg_full
    
  implicit none

  !--------------------------------------------------------------------
  !  Compute full multiple scattering within some cluster at some energy 
  !
  !  This uses the LU decomposition package from LAPACK.
  !  Driver routines: cgetrf (Decomposition), cgetrs: (Backsubstitution)
  !  coded by b.ravel
  !  modified by a.l.ankudinov to include spin and SO interactions
  !  feb 2000
  !
  !  Most of the information needed by this package is set into modules
  !  In package xprep, the lists of
  !  atomic coordinates and potential indeces are organized so that the
  !  first npot+1 entries are examples of each of the unique potentials.
  !  Consequently, only the upper left hand corner of the FMS matrix
  !  need be recomposed to get the set of submatrices necessary to
  !  compute chi for each type of atom in the cluster.
  !
  !  See subroutine fmstot.f for an example of decoding the output of this
  !  subroutine. The third index of gg refers to the unique potential with
  !  element 0 being the absorbing atom.  
  !  The first two indeces are related to the |lms> state by the
  !  formula:
  !       nsp=1, no spin indeces
  !       lm  = ( l**2 + 1 ) + ( l + m )
  !            thus {1..(lx+1)^2} ==>
  !            {0,0 1,-1 1,0 1,1 2,-2 2,-1 2,0 2,1 2,2 ...
  !                   lx,lx-1 lx,lx}
  !       nsp=2, with spin indeces
  !       lms  = 2*( l**2 + 1 ) + 2*( l + m ) + (s-1/2)
  !            thus {1...2*(lx+1)^2} ==>
  !            {0, 0,-1/2  0. 0,1/2
  !             1,-1,-1/2  1,-1,1/2  1,0,-1/2  1,0,1/2  1,1,-1/2 1,1,1/2
  !             2,-2,-1/2  2,-2,1/2  2,-1,-1/2 2,-1,1/2 ...    lx,lx,1/2}
  !
  !  The calling protocol for xpreppack and fmspack is;
  !          ...
  !          call xprep(nat, inclus, npot, iphat, rmax, rat,
  !     $            xnrm, izx, temper, thetad)
  !          energy loop {
  !             ...
  !             call fms(nsp, inclus, npot, ck, lipotx, xphase,
  !                      ik, iverb, gg)
  !             ... }
  !
  !  fmspack contains the following routines:
  !    fms.f:     main routine of fmspack
  !    kets.f:    compute all state kets for current energy
  !    xclmz.f:   compute hankle-like polynomials for current energy
  !    xgllm.f:   compute z-axis propagators for current energy
  !    cgetrf.f:  LU decomposition of matrix
  !    cgetrs.f:  backsubstitution of LU decomposed matrix
  !    lu_misc.f: various routines called by LU package
  !
  !---------------------------------------------------------------------
  !  input
  !    nsp:    1) no spin indeces 2) with spin indeces
  !    inclus: number of atoms in cluster
  !    npot:   number of unique potentials in cluster
  !    ck:     complex momentum of current energy point
  !    lipotx: (0:nphasx) max l for each unique potential
  !    xphase: (0:lx, 0:nphasx) single complex array of partial wave
  !            phase shifts for each unique potential
  !    ik:     current energy grid index, used for run-time messages
  !    iverb:  do nothing when iverb <= 0
  !            1  => write a message about grid point and matrix size
  !
  !  passed in common from xprep package (xstruc.h)
  !    xrat:   (3,nclusx) array of coordinates with first npot+1
  !            elements each a unique potential
  !    xphi:   (nclusx, nclusx) angles between z axis and vectors
  !            connecting the atoms in the cluster
  !    iphx:   (nclusx) potential index of each atom in the cluster
  !    drix:   huge matrix containing all rotation matrix elements
  !            needed for computation of free electron propagators
  !    xnlm:   matrix of legendre polynomial normalization factors
  !    xpsile: matrix containing wave functions for hybridization
  !            calculation
  !    sigsqr: (nclusx,nclusx) matrix of pair-wise mean square
  !            displacements about interatomic distances.  Currently only
  !            calculated by the correlated debye model.
  !
  !  output
  !    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
  !          angular momentum basis for each unique potential

  ! Inputs:
  integer, intent(in) :: lfms, nsp,ispin,inclus,npot
  integer, intent(in) :: ik,iverb,minvin
  complex, intent(in) :: xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
  real,    intent(in) :: rdirec,toler1,toler2
  logical, intent(in) :: lcalc(0:lx)
  
  integer, intent(inout) :: lipotx(0:nphx)

  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)

  ! Local Variables:
  complex   term, prefac, gllmz
  complex, allocatable :: clm(:,:), xclm(:,:,:,:,:)
  complex, allocatable :: xrho(:,:,:)

  complex,allocatable ::   tmatrx(:,:) !KJ made allocatable 5-07 (nspx, istatx)

  !     big work matrices
  complex, allocatable :: g0(:,:), g0t(:,:)

  integer minv
  integer i0(0:nphx)
  character*3  cerr, dec
  character*13 trans
  character*75 messg

  !KJ next variables are all mine :
  complex*16 jl(0:3),hl1(0:3),xxx(3),yy(16)  !KJ debugging
  logical, parameter :: makeg0=.false.,makeg=.false.,makegk0=.false.

!  logical usesprkkrtmat

  complex v1(9,9),vv(9,9,2,2)
  real*8 kx,ky,kz,x(3),m(3,3),y(3)
  complex,allocatable :: tsst(:,:,:)
  complex eryd,p
  integer nts,msizes,irels,lun,iatom,zs(npot)
  character*6 norms
  character*13 fefffile

  ! These are in the constants module But kept for now
  ! so that precision matches earlier verions for tests

   real,    parameter :: pi = 3.1415926535897932384626433e0
   real,    parameter :: bohr = 0.529177249e0
   real,    parameter :: one = 1, zero = 0
   complex, parameter :: coni = (0,1)

   ! Added to satisfy implicit none:
   integer :: i,j,k,isp,isp1,ist1,ist2,iatom1,iatom2 ! Loop indecies
   integer :: ix,ip,ll,mm,mu ! Loop indecies
   integer :: iat1,iat2,l1,l2,m1,m2,muabs,nstates,netab,mls,isp2
   integer :: iph,is,jj1,jj2,inlat,msord,lplus1,mplus1,ipi,ipf
   real    :: r,rr,rdir2

   allocate(g0(istatx,istatx), g0t(istatx,istatx))
   allocate(clm(lx+2, 2*lx+3), xclm(0:lx, 0:lx, nclusx, nclusx,nspx))
   allocate(xrho(nclusx, nclusx, nspx))
   g0t = 0.d0 ! JK - added initialization of g0t to get rid of valgrind errors.

   minv = minvin
   if (allocated(gg_full) .and. minv .ne. 0) then
! Modified by FDV
! Split line to compile in Solaris Studio
     call wlog("Calculating full scattering matrix (required by " // &
               "DENSITY and COMPTON cards) is currently only supported " // &
               "for LU inversion. Changing minv to 0 here and using LU" // &
               " (this may take longer than your chosen method).")
     minv = 0
   end if

400 format(i4)
  do i=0,nphx
     if (lipotx(i).lt.0)  lipotx(i) = lx
     if (lipotx(i).gt.lx) lipotx(i) = lx
     i0(i) = -1
  enddo

  !     initialize gg to zero
  do i = 0, nphasx
     do j = 1, nspx*(lx+1)**2
        do k = 1, nspx*(lx+1)**2
           gg( k, j, i) = cmplx( zero, zero)
        enddo
     enddo
  enddo

!JPR Comment usesprkkrtmat out
  !KJ next section initializes the t-matrix
!  usesprkkrtmat=(fullpot.and.(sprkkrpot.eq.1)) !KJ
!  if(.not.usesprkkrtmat) then

     allocate(tmatrx(nspx, istatx))
!  else
!     allocate(tmatrx(istatx,istatx))
!     if(minv.ne.0) stop 'fms cannot do full potential for minv > 0 yet'
!  endif
  tmatrx=cmplx(0,0)
  !KJ

  if (lfms.eq.0) then
     ipi = iphx(1)
     ipf = iphx(1)
  else
     ipi = 0
     ipf = npot
  endif

  ! --- get basis kets; output array 'lrstat' passed through module
  call getkts(nsp, inclus, npot, iphx, lipotx, i0)

  ! --- sanity check for i0(ip)
  do ip = ipi, ipf
     if (i0(ip) .lt. 0) then
        call wlog (' Cannot find all representative atoms')
        call wlog (' Increase FMS radius and rerun.')
        call par_stop(' In subroutine FMS')
     endif
  end do

  ! --- runtime message if requested
!KJ  if (iverb.gt.0 .and. minv.eq.0) then
  if (iverb.gt.0) then  !KJ let's show this for any inversion method (minv)
     dec = 'LUD'
     write(messg, 4010)this_process,dec, ik, istate
4010 format('  ',i3,'   FMS matrix (', a, ') at point ', i3, ', number of state kets =', i4)
     !call wlog(messg)
  endif

  ! --- get all c_lm(z) values for this energy, i,j sum over all atom
  !     pairs xrho and xclm are symmetric in ij
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  !  nota bene, in the code for setting the clmz, the indexing starts
  !  at 1 rather than 0.  To my mind, that is confusing, so here I
  !  reindex when I copy from clm to xclm.  See the note about this in
  !  subroutine xclmz
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  lplus1 = lx+1
  mplus1 = lx+1
  do i=1,inclus
     do j=1,i

        ! ------- get and store rho for this pair of atoms   bohr units on
        !         r and ck
        r   = zero
        do ix=1,3
           r = r + (xrat(ix,i) - xrat(ix,j))**2
        enddo
        r   = sqrt(r)

        do isp = 1,nsp
           xrho(i,j,isp) = ck(isp) * r
           xrho(j,i,isp) = xrho(i,j,isp)

           ! ------- store the c_lm(z) for all the rhos at this energy
           !            xclm(i,j) = xclm(j,i) by symmetry
           if (i.ne.j) call xclmz(lx,lplus1,mplus1,xrho(i,j,isp),clm)
           do ll = 0,lx
              do mm = 0,lx
                 if (i.eq.j) then
                    xclm(mm,ll,j,i,isp) = cmplx(zero,zero)
                 else
                    xclm(mm,ll,j,i,isp) = clm(ll+1,mm+1)
                    xclm(mm,ll,i,j,isp) = clm(ll+1,mm+1)
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo



  !KJ next section reads a t-matrix from file, output of (FP-)SPRKKR
!   if (usesprkkrtmat) then
!      call wlog('Using SPRKKR T-matrix')
!      call wlog('Still using FEFF core hole')
!      if(npot.gt.9) stop 'cannot read more than 9 files'
!      DO IATOM=1,npot
!         LUN=70+IATOM
! 2019    FORMAT(3x,5i6,2x,a6)	  
! 2020    FORMAT(1000e14.6)	
!         IF (ik .EQ. 1) then
!            FEFFFILE(1:13)='TDIAFEFFX.DAT'
!            FEFFFILE(9:9)=CHAR(48+IATOM)
!            OPEN(LUN,FILE=FEFFFILE,FORM='FORMATTED',STATUS='UNKNOWN')
!            READ(LUN,2019) MSIZES,NETAB,zs(iatom),irels,nts,norms
!            if(iatom.eq.1) then
!               allocate(tsst(msizes,msizes,nts))
!               tsst=dcmplx(0,0)
!            endif
!            if(msizes.ne.mls) write(*,*) 'msizes',msizes,'ne mls',mls
!            if(nts.ne.npot) write(*,*) 'nts',nts,'ne npot',npot
!         ENDIF
!         read(LUN,2020) ERYD,P,(TSST(I,I,IATOM),I=1,MSIZEs)
!         IF (ik .EQ. NETAB) CLOSE(LUN)
!      ENDDO
!      deallocate(tsst)
!   endif


  ! --- fill the G0 and T matrices for this energy
  rdir2 = rdirec**2
  !	call wlog('KJ modified big loop in fms!!!')
  STATE1_LOOP: do ist1=1,istate
     iat1 = lrstat(1, ist1)
     l1   = lrstat(2, ist1)
     m1   = lrstat(3, ist1)
     isp1 = lrstat(4, ist1)

     STATE2_LOOP: do ist2=1,istate
        iat2 = lrstat(1, ist2)
        l2   = lrstat(2, ist2)
        m2   = lrstat(3, ist2)
        isp2 = lrstat(4, ist2)

        rr =   (xrat(1,iat1)-xrat(1,iat2))**2 +        &
             & (xrat(2,iat1)-xrat(2,iat2))**2 +     &
             & (xrat(3,iat1)-xrat(3,iat2))**2

        ! Equation 9 in Rehr, Albers
        ! <LR| G |L'R'>

        if (iat1.eq.iat2) then
           ! Same atom: G=0, calculate T-matrix 
           g0(ist1,ist2)     = cmplx(zero,zero)

          ! if(usesprkkrtmat) then
              ! CODE HERE???

          ! else  ! construct FEFF t-matrix

              ! Notice that T is tri-diagonal, due to conservation of
              ! total momentum.(will be broken by nonspherical potential)
              !  --- potential index for this atom
              iph = iphx(iat1)
              if (nsp.eq.1.and.ispin.eq.0) then
                 if (ist1.eq.ist2) then
                    tmatrx(1, ist1) =  &
                         & (exp(2*coni*xphase(isp1,l1,iph)) - one ) &
                         &  / (2*coni) 
                 end if
              else
                 if (ist1.eq.ist2) then
                    ! Set spin index for t3jm and t3jp
                    is = isp1
                    if (nsp.eq.1) then
                       ! Special case
                       is = 1
                       if (ispin.gt.0) is = 2
                    endif

                    ! Diagonal matrix element
                    tmatrx(1, ist1) =                &
                         &  ( exp(2*coni*xphase(isp1,l1,iph)) - one )     &
                         &  / (2*coni) * t3jm (l1, m1, is)**2  +          &
                         &  ( exp(2*coni*xphase(isp1,-l1,iph)) - one )    &
                         &  / (2*coni) * t3jp (l1, m1, is)**2 
                 elseif (nsp.eq.2.and.l1.eq.l2.and.m1+isp1.eq.m2+isp2) then
                    ! Same orb. mom. and total momentum projections conserved
                    ! Calculate off-diagonal T-matrix element
                    ! tmatrx(2, ist1) = here only if nspx equal to 2
                    tmatrx(nsp, ist1) =                    &
                         & (exp(2*coni*xphase(isp1, l1,iph)) - one +     &
                         & exp(2*coni*xphase(isp2, l1,iph)) - one ) / (4*coni) &
                         & * t3jm (l1, m1, isp1) * t3jm (l1, m2, isp2)  + &
                         & (exp(2*coni*xphase(isp1,-l1,iph)) - one +     &
                         &  exp(2*coni*xphase(isp2,-l1,iph)) - one ) / (4*coni) &
                         & * t3jp (l1, m1, isp1) * t3jp (l1, m2, isp2)
                 endif
              endif
           !endif  ! FEFF or SPRKKR t-matrix


        elseif (isp1.eq.isp2 .and. rr.le.rdir2) then
           ! Different atoms, same spin: T=0, calculate G
           g0(ist1,ist2) = cmplx(zero,zero)
           do mu=-l1,l1
              ! --- third arg in drix: 0==>beta, 1==>-beta
              muabs = abs(mu)
              call xgllm(muabs, ist1, ist2, lrstat,                     &
                   &                   xclm(0,0,1,1,isp1), gllmz )
              g0(ist1,ist2) = g0(ist1,ist2) +                           &
                   & drix(mu,m1,l1,1,iat2,iat1) *  gllmz *                &
                   & drix(m2,mu,l2,0,iat2,iat1)
           enddo
           prefac = exp(coni*xrho(iat1,iat2,isp1)) /                   &
                &   xrho(iat1,iat2,isp1)
           ! Use correlated debye model, sigsqr is in AA^2
           prefac = prefac * exp(-1 * sigsqr(iat1,iat2) *              &
                &   ck(isp1)**2 / bohr**2)
           g0(ist1,ist2) = prefac * g0(ist1,ist2)
        else
           ! Different atoms, different spins:T=G=0
           g0(ist1,ist2) = cmplx(zero,zero)
        endif

     enddo STATE2_LOOP
  enddo STATE1_LOOP



  if(makeg.or.makeg0.or.makegk0) write(*,*) 'ck',ck

  if(makegk0) then
     if(nsp.eq.1) nstates=9 !9 diamond !16 Po
     if(nsp.eq.2) nstates=18 !18 diamond !32 Po

     kx=dble(0.3) ; ky=dble(0.3) ; kz=dble(0.3)

     m=dble(1)
     do i=1,3 
        m(i,i)=-dble(1)
     enddo

     m=m/dble(1.7835)*dble(0.529177)/dble(2)
     vv=cmplx(0,0)

     IATOM1_LOOP: do iatom1=1,2
        IATOM2_LOOP: do iatom2=1,2
           jj1=(iatom1-1)*nstates ; jj2=(iatom2-1)*nstates

           v1=cmplx(0,0)
           inlat=0

           INCLUS_LOOP: do i=1,inclus
              jj2=(i-1)*nstates
              x=xrat(:,i)-xrat(:,iatom2) !-xrat(:,iatom1)
              !        Check that x is a lattice vector :
              y=dble(0)
              do j=1,3
                 do k=1,3
                    y(j)=y(j)+m(j,k)*x(k)
                 enddo
              enddo
              do j=1,3
                 y(j)=dabs(y(j)-nint(y(j)))
              enddo
              if(dabs(y(1))+dabs(y(2))+dabs(y(3)).gt.dble(0.001)) then
                 cycle  !not a lattice vector
              else
                 inlat=inlat+1
              endif

              v1=v1+g0(jj1+1:jj1+nstates,jj2+1:jj2+nstates)*exp(-coni*    &
                   &   real((kx*x(1)+ky*x(2)+kz*x(3))))
           enddo INCLUS_LOOP

           vv(:,:,iatom1,iatom2)=v1
        enddo IATOM2_LOOP
     enddo IATOM1_LOOP

     open (99,file='g011rk.txt',position='append')
     write(99,167) ck,(vv(ist1,:,1,1),ist1=1,nstates)
     close(99)
     open (99,file='g022rk.txt',position='append')
     write(99,167) ck,(vv(ist1,:,2,2),ist1=1,nstates)
     close(99)
     open (99,file='g012rk.txt',position='append')
     write(99,167) ck,(vv(ist1,:,1,2),ist1=1,nstates)
     close(99)
     open (99,file='g021rk.txt',position='append')
     write(99,167) ck,(vv(ist1,:,2,1),ist1=1,nstates)
     close(99)

     return

  end if


  if(makeg0) then
     !KJ debugging
     open (99,file='g0136r.txt',position='append')
     if(nsp.eq.1) nstates=9 !9 diamond !16 Po
     if(nsp.eq.2) nstates=18 !18 diamond !32 Po
     iatom1=1
     iatom2=36
     write(*,*) 'doing atoms',iatom1,iatom2
     jj1=(iatom1-1)*nstates ; jj2=(iatom2-1)*nstates
     write(100,167) xrat(:,iatom1),xrat(:,iatom2),ck(1)
     do ist1=jj1+1,jj1+nstates
        do ist2=jj2+1,jj2+nstates
           if (cabs(g0(ist1,ist2)).lt.0.00001) g0(ist1,ist2)=cmplx(0,0)
        enddo
     enddo
     write(99,167) ck,(g0(ist1,jj2+1:jj2+nstates),                   &
          &               ist1=jj1+1,jj1+nstates)

     close(99)
251  format(5000(F14.7,1x,F14.7,3x))       


167  format(5000(e12.4,1x,e12.4,3x))

     open(99,file='hankel.txt')
     xxx(1)=ck(1)*(xrat(1,1)-xrat(1,2))
     xxx(2)=ck(1)*(xrat(2,1)-xrat(2,2))
     xxx(3)=ck(1)*(xrat(3,1)-xrat(3,2))

     call besjh(dcmplx(xrho(1,2,1)),3,jl,hl1)
     call ylm(dble(xxx),3,yy)
     yy=yy*sqrt(4*3.1415)  !feff's Ylm's are not normalized ...
     ist1=0
     do j=0,3
        do i=-j,j
           ist1=ist1+1
           yy(ist1)=yy(ist1)*dcmplx(0,1)**j * hl1(j)
        enddo
     enddo
     write(99,167) yy
     close(99)


     !KJ
     if (.not.makeg) return
  endif

  if (minv.eq.0) then
     !if(usesprkkrtmat) then
     !   call gglufullpot(nsp, i0, ipi, ipf, lipotx, g0,tmatrx,g0t,gg,ck(1))
     !else
        call gglu(nsp, i0, ipi, ipf, lipotx, g0,tmatrx,g0t,gg,ck(1),ik)
     !endif
  elseif (minv.eq.1) then
     dec = 'VdV'
     call ggbi ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg, toler1, toler2, lcalc, msord)
  elseif (minv.eq.2) then
     dec = 'LLU'
     call ggrm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg, toler1, toler2, lcalc, msord)
  elseif (minv.eq.3) then
     dec = 'GMS'
     call gggm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg, toler1, toler2, lcalc, msord)
  else
     dec = 'TF'
     call ggtf ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg, toler1, toler2, lcalc, msord)
  endif

!KJ commented out bc way too annoying ...
!KJ (years later) commenting out again, it still annoys me ...
!  if (minv.ne.0) then
!     write(messg, 410)this_process,dec, ik, istate, msord
410  format('  ',i3,'. Iterative FMS (', a, ') at point ', i3,      &
          & '; matrix size =', i4,'; MS order =',i5)
!     call wlog(messg)
!  endif


  deallocate(tmatrx)
  deallocate(g0, g0t)
  deallocate(clm, xclm,xrho)

  return
end subroutine fms

!--------------------------------------------------------------------

!    ----------------------------------------------------------------
subroutine xclmz(lx,lmaxp1,mmaxp1,rho,clm)

  use constants, only: one,zero,coni

  implicit none

  
  !     calculates energy dependent factors needed in subroutine gllm
  !     c(il,im) = c_l^(m)z**m/m!=c_lm             by recursion
  !     c_l+1,m  = c_l-1,m-(2l+1)z(c_l,m-c_l,m-1)  l ne m
  !     c_m,m    = (-z)**m (2m)!/(2**m m!)         with z=1/i rho
  !
  !  input:
  !    lmaxp1, mmaxp1:  largest angular momentum under consideration + 1
  !    rho:  distance between atoms * complex momentum at this energy
  !          point
  !  output:
  !    clm(lx+1,lx+1):  Hankle-like polynomials from RA

  integer :: ltotb,mtotb,ntotb,mntot

  ! Input:
  integer, intent(in) :: lmaxp1,mmaxp1
  complex, intent(in) :: rho
  integer, intent(in) :: lx
  ! Output:
  complex, intent(out) :: clm(lx+2,2*lx+3)
  ! jpr: Original dimensions of clm were: clm(ltotb,mntot) => clm(lx+1,2*lx+2)
  !      These were changed to match the size of the input array.

  ! Local variables
  complex z, cmm

  ! Added to satisfy implicit none
  integer :: i,im,il  ! Loop indecies
  integer :: imp1,lmax,m,mmxp1

  ltotb=lx+1
  mtotb=ltotb
  ntotb=ltotb
  mntot=mtotb+ntotb

  clm(:,:)=cmplx(0.0d0,0.0d0) !KJ 7-09 taken from JAS feffq - don't think it's ever necessary unless array used inappropriately (eg. queried for higher element than calculated).

  cmm  = cmplx(one, zero)
  z    = (-coni)/rho

  clm(1,1) = cmplx(one,zero)
  clm(2,1) = clm(1,1) - z

  lmax = lmaxp1-1
  do il=2,lmax
     clm(il+1,1) = clm(il-1,1) - (z * (2*il-1) * clm(il,1))
  enddo
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  !  nota bene:  the 2l-1 factor above is correct, even though in Rehr,
  !  Albers equation 4 appears with a 2l+1.  The reason has to do with
  !  the indexing.  in RA the subscripts on the c's start at 0.  In this
  !  piece of code, the subscripts start at 1.  If you sub l-1 for
  !      l, 2l+1 --> 2l-1
  !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


  mmxp1 = min(lmaxp1, mmaxp1)
  IMLOOP: do im=2,mmxp1
     m    = im-1
     imp1 = im+1
     cmm  = (-cmm) * (2*m-1) * z
     clm(im,im)   = cmm
     clm(imp1,im) = cmm * (2*m+1) * (1-im*z)
     ILLOOP: do il=imp1,lmax
        clm(il+1,im) = clm(il-1,im) - (2*il-1) * z *                  &
             &                               ( clm(il,im)+clm(il,m) )
        !           l = il-1
        !           clm(il+1,im) = clm(l,im) - (2*l+1) * z *
        !      $                               ( clm(il,im)+clm(il,m) )
     enddo ILLOOP
  enddo IMLOOP

  return
end subroutine xclmz

subroutine xgllm(mu, ist1, ist2, lrstat, xclm, gllmz)
  !--------------------------------------------------------------------
  !  this calculates equations 11,12 from Rehr, Albers PRB v.41,#12,
  !  p.8139,  the output is the G term in equation 9 from that paper
  !
  !  input:
  !    mu:         abs val of magnetic state in sum in eqn 11 RA, mu>=0
  !    ist1, ist2: state indices of mat. elem., first index of lrstat
  !    lrstat:     (4,istatx,nkmin:nex) array of LR states
  !    xclm:       (0:lx,0:lx,nclusx,nclusx) array of c_lm(z) for
  !                present energy value
  !  output:
  !    gllmz:      g_ll'^|m|(z), for present state & energy, eqn 11 RA
  !--------------------------------------------------------------------
  !  this requires that N_lm normalization factors and c_lm(z)
  !  polynomials have already been calculated.
  !--------------------------------------------------------------------
  use DimsMod, only: istatx, nclusx, lx
  use rotx
  use lnlm
  use xstruc

  !  use constants

  implicit none

  integer, intent(in) :: mu,ist1,ist2
  integer, intent(in) :: lrstat(4, istatx)
  complex, intent(in) :: xclm(0:lx, 0:lx, nclusx, nclusx)

  complex, intent(out) :: gllmz

  real, parameter :: zero=0.e0

  ! Local variables
  complex sum
  complex gam, gamtl

  ! Added to satisfy implicit none:
  integer :: iat1,iat2,l1,l2,nu,numax,mn

  iat1     = lrstat(1, ist1)
  l1       = lrstat(2, ist1)
  iat2     = lrstat(1, ist2)
  l2       = lrstat(2, ist2)
  numax    = min(l1, l2-mu)

  sum      = cmplx(zero, zero)
  do nu=0,numax
     mn    = mu+nu

     !       bug for xnlm with nspx=2
     gamtl = (2*l1+1) * xclm(nu, l1, iat2, iat1) / xnlm(mu, l1)
     gam   = (-1)**mu * xclm(mn, l2, iat2, iat1) * xnlm(mu, l2)

     sum   = sum + gamtl * gam
  enddo

  gllmz = sum

  return
end subroutine xgllm
