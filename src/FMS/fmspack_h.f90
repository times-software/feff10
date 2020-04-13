!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmspack.f90,v $:
! $Revision: 1.19 $
! $Author: jorissen $
! $Date: 2011/12/11 00:37:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fms_h(lfms, nsp, ispin, inclus, npot, ck, lipotx, xphase_m, &
        ik, iverb, rdirec, toler1, toler2, lcalc, gg_m, TFrm,  &
        TFrmInv, UseTFrm)

!KJ removed "minv" argument -- since only minv=0 was implemented anyway
  use DimsMod, only: lx, nphx=>nphu, nspx=>nspu, nphasx, istatx, nclusx
  use stkets
  use rotx
  use lnlm
  use xstruc
  use t3j
  use par
  use hubbard_inp, only: l_hubbard
    
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
  integer, intent(in) :: ik,iverb
  complex, intent(in) :: xphase_m(nspx,-lx:lx,(lx+1)**2,0:nphasx),ck(nspx)
  real,    intent(in) :: rdirec,toler1,toler2
  logical, intent(in) :: lcalc(0:lx)
  logical, intent(in) ::   UseTFrm(0:lx,0:nphx)
  complex, intent(in) ::   TFrm(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx)
  complex, intent(in) ::   TFrmInv(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx)
  integer, intent(inout) :: lipotx(0:nphx)
  ! Return matrix containing info about each unique potential
  complex, intent(out) :: gg_m(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

  ! Local Variables:
  complex   term, prefac, gllmz
  complex, allocatable :: clm(:,:), xclm(:,:,:,:,:)
  complex, allocatable :: xrho(:,:,:)

  complex,allocatable ::   tmatrx_m(:,:) 
  complex,allocatable ::   tmatrxfull(:,:) 
  complex,allocatable ::   tmatrxtmp(:,:)
  integer, allocatable::   iStartTFrm(:,:) 
  integer, parameter :: iTestTFrm = 0
  ! iTestTFrm = {0 - No test, 1, identity, 2 small off diagonal on all atoms.}

  !     big work matrices
  complex, allocatable :: g0(:,:), g0t(:,:)

  integer i0(0:nphx)
  character*3  cerr, dec
  character*13 trans
  character*75 messg

  ! These are in the constants module But kept for now
  ! so that precision matches earlier versions for tests

   real,    parameter :: pi = 3.1415926535897932384626433e0
   real,    parameter :: bohr = 0.529177249e0
   real,    parameter :: one = 1, zero = 0
   complex, parameter :: coni = (0,1)

   ! Added to satisfy implicit none:
   integer :: i,j,k,isp,isp1,ist1,ist2,iatom1,iatom2,im1,im2 ! Loop indices
   integer :: ix,ip,ll,mm,mu ! Loop indices
   integer :: iat1,iat2,l1,l2,m1,m2,muabs,nstates,netab,mls,isp2
   integer :: iph,is,jj1,jj2,inlat,msord,lplus1,mplus1,ipi,ipf
   real    :: r,rr,rdir2

   ! Again, adding to satisfy implicit none:
   integer irow, icol, iStart, k1, k2, il, ik1
   real tmp

   allocate(g0(istatx,istatx), g0t(istatx,istatx))
   allocate(clm(lx+2, 2*lx+3), xclm(0:lx, 0:lx, nclusx, nclusx,nspx))
   allocate(xrho(nclusx, nclusx, nspx))
   allocate(tmatrx_m(nspx, istatx),tmatrxtmp(2*lx+1,2*lx+1))
   allocate(tmatrxfull(istatx, istatx))
   allocate(iStartTFrm(nclusx,0:lx))

400 format(i4)

      tmatrxfull(:,:) = cmplx(0.0,0.0)
!!     Testing transformation of t-matrices.
!      IF(iTestTFrm.EQ.1) THEN
!!        First test - identity matrix
!         TFrm(:,:,:,:) = cmplx(0.0,0.0)
!         TFrmInv(:,:,:,:) = cmplx(0.0,0.0)
!         UseTFrm(:,:) = .FALSE.
!         UseTFrm(2,:) = .TRUE.
!         DO il = 0, lx
!            DO iph = 0, nphx
!               DO ik1 = 1, 2*(l_hubbard+1)
!                  TFrm(ik1,ik1,il,iph) = cmplx(1.0,0.0)
!                  TFrmInv(ik1,ik1,il,iph) = cmplx(1.0,0.0)
!               END DO
!            END DO
!         END DO         
!      ELSEIF(iTestTFrm.EQ.2) THEN
!!        for cubic system, all elements are the same.
!!        Test will off diagonal elements x on the 1,-1 and -1, 1.
!!        This gives the following transformation matrix.
!         UseTFrm(:,:) = .FALSE.
!         IF(iphx(1).EQ.0) UseTFrm(1,0) = .TRUE.
!
!         DO il = 1, lx
!            DO iph = 1, nphx
!               DO ik1 = 1, 2*(l_hubbard+1)
!                  TFrm(ik1,ik1,il,iph) = cmplx(1.0,0.0)
!               END DO
!            END DO
!         END DO
!         TFrm(1,1,1,0) = -1.0/SQRT(2.0)
!         TFrm(2,1,1,0) = 0.0
!         TFrm(3,1,1,0) = 1.0/SQRT(2.0)
!         TFrm(1,2,1,0) = 0.0
!         TFrm(2,2,1,0) = 1.0
!         TFrm(3,2,1,0) = 0.0
!         TFrm(1,3,1,0) = 1.0/SQRT(2.0)
!         TFrm(2,3,1,0) = 0.0
!         TFrm(3,3,1,0) = 1.0/SQRT(2.0)
!         TFrmInv(1,1,1,0) = -1.0/SQRT(2.0)
!         TFrmInv(2,1,1,0) = 0.0
!         TFrmInv(3,1,1,0) = 1.0/SQRT(2.0)
!         TFrmInv(1,2,1,0) = 0.0
!         TFrmInv(2,2,1,0) = 1.0
!         TFrmInv(3,2,1,0) = 0.0
!         TFrmInv(1,3,1,0) = 1.0/SQRT(2.0)
!         TFrmInv(2,3,1,0) = 0.0
!         TFrmInv(3,3,1,0) = 1.0/SQRT(2.0)
!  END IF
  do i=0,nphx
     if (lipotx(i).le.0)  lipotx(i) = lx
     if (lipotx(i).gt.lx) lipotx(i) = lx
  enddo
  !     initialize gg to zero
  gg_m( :, :, :) = cmplx( zero, zero)

  if (lfms.eq.0) then
     ipi = iphx(1)
     ipf = iphx(1)
  else
     ipi = 0
     ipf = npot
  endif

  ! --- get basis kets; output array 'lrstat' passed through module
  call getkts(nsp, inclus, npot, iphx, lipotx, i0)
!  write(*,*) 'for getkts,',nsp, inclus, npot
!  write(*,*)'iphx', iphx
!  write(*,*)'lipotx', lipotx
!  write(*,*) 'ipi,ipf,i0',ipi,ipf,i0

  ! --- sanity check for i0(ip)
  do ip = ipi, ipf
     if (i0(ip) .lt. 0) then
        call wlog (' Cannot find all representative atoms. Increase FMS radius and rerun.')
        call par_stop(' In subroutine FMS')
     endif
  end do

  ! --- runtime message if requested
!KJ  if (iverb.gt.0 .and. minv.eq.0) then
  if (iverb.gt.0) then  !KJ let's show this for any inversion method (minv)
     dec = 'LUD'
     write(messg, 4010)this_process,dec, ik, istate
4010 format('  ',i3,'   FMS matrix (', a, ') at point ', i3, ', number of state kets =', i4)
     call wlog(messg)
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
        ! ------- get and store rho for this pair of atoms   bohr units on  r and ck
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
           if (i.ne.j) call xclmz_h(lx,lplus1,mplus1,xrho(i,j,isp),clm)
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


  ! --- fill the G0 and T matrices for this energy
  rdir2 = rdirec**2
  STATE1_LOOP: do ist1=1,istate
     iat1 = lrstat(1, ist1)
     l1   = lrstat(2, ist1)
     m1   = lrstat(3, ist1)
     isp1 = lrstat(4, ist1)
     if(m1.eq.-l1) iStartTFrm(iat1,l1)=ist1
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

              ! Notice that T is tri-diagonal, due to conservation of
              ! total momentum.(will be broken by nonspherical potential)
              !  --- potential index for this atom
              iph = iphx(iat1)
              if (nsp.eq.1.and.ispin.eq.0) then
                 if (ist1.eq.ist2) then
                    tmatrx_m(1, ist1) =  &
                    & (exp(2*coni*xphase_m(isp1,l1,l1**2+l1+m1+1,iph))-one) / (2*coni)
                     tmatrxfull(ist1,ist1) = tmatrx_m(1,ist1)
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
                    tmatrx_m(1, ist1) =                &
                          ( exp(2*coni*xphase_m(isp1,l1,l1**2+l1+m1+1,iph)) - one )     &
                           / (2*coni) * t3jm (l1, m1, is)**2  +          &
                           ( exp(2*coni*xphase_m(isp1,-l1,l1**2+l1+m1+1,iph)) - one )    &
                           / (2*coni) * t3jp (l1, m1, is)**2
                        tmatrxfull(ist1,ist1) = tmatrx_m(1,ist1)
                 elseif (nsp.eq.2.and.l1.eq.l2.and.m1+isp1.eq.m2+isp2) then
                    ! Same orb. mom. and total momentum projections conserved
                    ! Calculate off-diagonal T-matrix element
                    ! tmatrx(2, ist1) = here only if nspx equal to 2
                    tmatrx_m(nsp, ist1) =                    &
                          (exp(2*coni*xphase_m(isp1,l1,l1**2+l1+m1+1,iph)) - one +     &
                           exp(2*coni*xphase_m(isp2,l1,l1**2+l1+m1+1,iph)) - one ) / (4*coni) &
                          * t3jm (l1, m1, isp1) * t3jm (l1, m2, isp2)  + &
                          (exp(2*coni*xphase_m(isp1,-l1,l1**2+l1+m1+1,iph))-one +     &
                           exp(2*coni*xphase_m(isp2,-l1,l1**2+l1+m1+1,iph))-one)/(4*coni) &
                          * t3jp (l1, m1, isp1) * t3jp (l1, m2, isp2)
                       tmatrxfull(ist1+(-1)**isp1,ist1) = tmatrx_m(nsp,ist1)
                 endif
              endif

        elseif (isp1.eq.isp2 .and. rr.le.rdir2) then
           ! Different atoms, same spin: T=0, calculate G
           g0(ist1,ist2) = cmplx(zero,zero)
           do mu=-l1,l1
              ! --- third arg in drix: 0==>beta, 1==>-beta
              muabs = abs(mu)
              call xgllm_h(muabs, ist1, ist2, lrstat, xclm(0,0,1,1,isp1), gllmz )
              g0(ist1,ist2) = g0(ist1,ist2) +   drix(mu,m1,l1,1,iat2,iat1) *  gllmz *  drix(m2,mu,l2,0,iat2,iat1)
           enddo
           prefac = exp(coni*xrho(iat1,iat2,isp1)) /   xrho(iat1,iat2,isp1)
           ! Use correlated debye model, sigsqr is in AA^2
           prefac = prefac * exp(-1 * sigsqr(iat1,iat2) *   ck(isp1)**2 / bohr**2)
           g0(ist1,ist2) = prefac * g0(ist1,ist2)
        else
           ! Different atoms, different spins:T=G=0
           g0(ist1,ist2) = cmplx(zero,zero)
        endif

     enddo STATE2_LOOP
  enddo STATE1_LOOP

!!     Testing transformation of t-matrix
!      IF(iTestTFrm.EQ.1) THEN
!!        First test - identity matrix
!      ELSEIF(iTestTFrm.EQ.2) THEN
!!        Second test - alpha is same on diagonal, with off
!!        diagonal elements in the 1,-1 and -1,1 corners only.
!         tmp = 1.0
!         tmatrxfull(2,2) = tmatrxfull(2,2)*(1.0 - tmp) 
!         tmatrxfull(3,3) = tmatrxfull(3,3)
!         tmatrxfull(4,4) = tmatrxfull(4,4)*(1.0 + tmp)
!      END IF
      DO iat1 = 1, inclus
         iph = iphx(iat1)
         DO il = 0, lipotx(iph)
            IF(UseTFrm(il,iph)) THEN
!              Transform t-matrices, i.e., TFrmInv.tmatrx.TFrm
               iStart = iStartTFrm(iat1,il) - 1
               DO icol = 1, 2*il + 1
                  DO irow = 1, 2*il + 1
                     tmatrxtmp(irow,icol) = 0.0
                     tmp = 0.0
                     DO k1 = 1, 2*il + 1
                        DO k2 = 1, 2*il + 1
!                          This transforms from the diagonal
!                          basis to mm'.
                            tmatrxtmp(irow,icol) =  tmatrxtmp(irow,icol) +  TFrm(irow,k1,il,iph) *  &
                                     tmatrxfull(k1+iStart,k2+iStart)*  TFrmInv(k2,icol,il,iph)
                        END DO
!                       Check transform matrices.
                        tmp = tmp + TFrm(irow,k1,il,iph)*  TFrmInv(k1,icol,il,iph)
                     END DO
                     IF(irow.EQ.icol) THEN
                        IF(ABS(tmp-1.0).GT.1.E-4) PRINT*, tmp,il, iph, irow,  'Possible problem with transformation matrices'
                     ELSE
                        IF(ABS(tmp).GT.1.E-5) PRINT*, tmp, irow,icol,  'Possible problem with transformation matrices'
                     END IF
                  END DO
               END DO
               DO k1 = 1, 2*il + 1
                  DO k2 = 1, 2*il + 1
                     tmatrxfull(k1+iStart,k2+iStart) =  tmatrxtmp(k1,k2)
                  END DO
               END DO
            END IF
         END DO
      END DO
      


      call gglu_h(nsp, i0, ipi, ipf, lipotx, g0,tmatrx_m, tmatrxfull, g0t,gg_m,ck(1),ik)

      DO iph = 0, nphasx
         DO il = 0, lipotx(iph)
            IF(UseTFrm(il,iph)) THEN
!              back transform gg, i.e., TFrmInv.gg.TFrm
               iStart = il**2
               DO icol = 1, 2*il + 1
                  DO irow = 1, 2*il + 1
                     tmatrxtmp(irow,icol) = 0.0
                     DO k1 = 1, 2*il + 1
                        DO k2 = 1, 2*il + 1
                           tmatrxtmp(irow,icol) =  tmatrxtmp(irow,icol) + TFrmInv(irow,k1,il,iph)* &
                               gg_m(k1+iStart,k2+iStart,iph) * TFrm(k2,icol,il,iph)
                        END DO
                     END DO
                  END DO
               END DO
               DO k1 = 1, 2*il + 1
                  DO k2 = 1, 2*il + 1
                     gg_m(k1+iStart,k2+iStart,iph) =  tmatrxtmp(k1,k2)
                  END DO
               END DO
            END IF
         END DO
      END DO
!KJ commented out bc way too annoying ...
!KJ (years later) commenting out again, it still annoys me ...
!  if (minv.ne.0) then
!     write(messg, 410)this_process,dec, ik, istate, msord
410  format('  ',i3,'. Iterative FMS (', a, ') at point ', i3,      &
          & '; matrix size =', i4,'; MS order =',i5)
!     call wlog(messg)
!  endif

  deallocate(g0, g0t)
  deallocate(clm, xclm,xrho)
  deallocate(iStartTFrm) 
  deallocate(tmatrxfull)
  deallocate(tmatrxtmp)
  deallocate(tmatrx_m)

  return
end subroutine fms_h

!--------------------------------------------------------------------

!    ----------------------------------------------------------------
subroutine xclmz_h(lx,lmaxp1,mmaxp1,rho,clm)

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
end subroutine xclmz_h 

subroutine xgllm_h(mu, ist1, ist2, lrstat, xclm, gllmz)
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
end subroutine xgllm_h  

