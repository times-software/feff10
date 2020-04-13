!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmstot.f90,v $:
! $Revision: 1.32 $
! $Author: jorissen $
! $Date: 2012/05/15 22:57:34 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fmstot

  ! Calculates FMS contribution to absorption
  ! uses Bruce Ravel subroutine to do FMS in self-consistency loop
  ! notice that it can do FMS with polarization dependence and
  ! always include l-->l-1 channel.
  ! written by alexei ankudinov 06.1997

!KJ cleaned up a bunch of unneeded (big!) arrays that were left over from old version fms. 7-09
! Kevin Jorissen 7-09 included feff-q (NRIXS).
!  2 differences :  1/ call rdxsphjas / rdxsph
!					2/ setting lcalc
!  I hope there's no additional junk like hidden differences in basis set or so, but I don't think so.

  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx, natx, ltot, dump_dimensions, nex
  use controls,only : ispace,sprkkrpot,sprkkrklist,fullpot
  use kklist !KJ for debugging only
  use kgenwork,only: nka
  use boundaries,only : maxl
  use IOMod
  use ErrorMod
  use constants
  use par
  use stkets
  use rotx
  use fms_inp,rclust=>rfms2,sigma2=>sig2g
  use global_inp,only: ipol,ispin,le2,angks,ptz,do_nrixs,nq
  use atoms_inp,only: nat,iphat,ratdbl=>rat
  use nrixs_inp,only: jmax,kfinmax,jinit,indmax,lgind
  use errorfile
  use fms_mod, only: gg_full
  use hubbard_inp

  implicit none
  real*8 wall_commend, wall_commst
  integer, parameter :: npadx=8, iblock=1
  real :: rat(3,natx)
  real :: rpart,aipart, rnrmax, thetax, temper, sig2
  integer :: ne, ne1, ne3, ihole ! ,nph  !KJ now in modules
  integer :: iz(0:nphx)
  character*512 :: slog
  ! Work space
  complex*16, allocatable :: ph(:,:,:,:)
  integer, allocatable    :: lmax(:,:)
  ! Complex energy grid emg is decomposed into em and eref to have
  !  the same structure in phase.bin
  complex*16, allocatable  :: em(:), eref(:,:)
  character*6, allocatable  :: potlbl(:)
  ! FMS staff
  integer, allocatable :: map(:)
  complex :: ck(nspx)
  complex, allocatable, dimension(:,:,:)   :: xphase
  complex, allocatable, dimension(:,:,:,:) :: gg,ggloc,ggremember 
  complex, allocatable, dimension(:,:,:)   :: gg_slice, gg_slice_loc
  complex, allocatable, dimension(:,:,:,:)   :: gg_diag, gg_diag_loc

  complex*16 :: dck
  integer :: lind(8)
  logical, allocatable :: lcalc(:)
  complex*16, allocatable :: rkk(:,:,:),rkk_nrixs(:,:,:,:)

  integer :: npot, nsp, ie, iverb, lfms
  integer :: NEPts, iggSize, igg_slice_size, igg_diag_size
  integer :: nkp_max,nkp_min,nkp_total,nkp_fixed

  ! Added to satisfy implicit none:
  integer :: i,j,iph,iph0,ik,ik0,nk,ipot,imj,inclus,lmaxp1,iat, imm
  integer :: idwopx,isp,ipp,ill, k1,nphtmp,indx
  integer :: kinit,linit
  real*8  :: rnrmav,xmu,edge,e_relative,e_largest
  integer :: istart,isize,ii,is
  integer :: i_for_next_report
  integer, parameter :: i_report_granularity=10
  logical, parameter :: write_dims=.false.

  ! Hubbard Definitions
  complex ::   TFrm(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx), TFrmInv(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx)
  complex ::   Trans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx),InvTrans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
  complex    ::   xphase_m(nspx, -lx:lx, (lx+1)**2, 0:nphx)
  complex, allocatable ::   gg_m(:,:,:)
  logical    ::   UseTFrm(0:lx,0:nphx)
  real*8     ::   max_gap_up, max_gap_down, aa, xmuO
  complex*16, allocatable ::   aph(:,:,:,:,:)
  real*8     ::   Vnlm(0:lx,(lx+1)**2,2,0:nphx)
  integer    ::   ios, ip, il, im1, im2, iss, lll, mmm, ispinp, nspp
  character*30 :: fname

!  !KJ following variables for reading SPRKKR phase shifts :
!  real*8, allocatable :: pstab(:,:,:,:)
!  complex*16, allocatable :: etab(:),ptab(:)
!  integer    :: sprne,sprnk,sprnmj
!  character  :: adsun
!  allocate(pstab(nex,2*lx+1,3,0:nphx))
!  allocate(etab(nex),ptab(nex))

  integer :: ldim
  ldim = nspx * (lx+1)**2
  ! Allocate local variables
  allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(lmax(nex, 0:nphx))
  allocate(em(nex), eref(nex, nspx))
  allocate(potlbl(0:nphx))
  allocate(map(nex))
  allocate(rkk(nex,8,nspx),rkk_nrixs(nex,nq,kfinmax,nspx))
  allocate(xphase(nspx, -lx:lx, 0:nphx), lcalc(0:lx))
  allocate( aph(nex,lx+1,(lx+1)**2,2,0:nphx))

  NEPts = 0
  iggSize = 0
  xphase = cmplx(0,0)
  iph0 = 0

  if(write_dims) call dump_dimensions

  ! Need less data than rphbin.f provides, also dimensions of ph array are different.
  !KJ Here nph could be changed - not too keen on the whole arrangement ... 7-09
  if(do_nrixs .eq. 1) then !NRIXS calculation
      call rdxsphjas (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk_nrixs, lmax, lmaxp1)
  elseif(do_nrixs .ne. 1 .and. i_hubbard .eq. 1) then  !regular FEFF calculation
      call rdxsph (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
  elseif(do_nrixs .ne. 1 .and. i_hubbard .eq. 2) then
      call rdxsph_h (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1,   &
      &                Trans, InvTrans, UseTFrm, aph)
  endif

  ! J Kas - moved allocation of the following to below rdxsph to allow use of ne/nph intstead of nex
  allocate(ggremember(ldim, ldim, 0:nphx, ne))
  allocate(     ggloc(ldim, ldim, 0:nphx, ne))
  allocate(        gg(ldim, ldim, 0:nphx, ne))
  allocate(gg_m(nspx*(lx+1)**2,nspx*(lx+1)**2, 0:nphx))
  call setkap (ihole, kinit, linit)
  npot = nph
  deallocate(rkk,rkk_nrixs)  !These are not used in "fms"; only here for compatibility with rdxsph as used elsewhere

  ! allocate work array for saving central-atom slice of fms matrix (used by COMPTON)
  if (save_gg_slice) then
    allocate(gg_full(istatx, istatx))
    gg_full(:,:) = 0.0
  end if

!! *KJ To import SPRKKR phase shifts, grab commented code block 1 at end of file and paste here.


  if(master) write(slog,'(a,i5,a)') 'Using ', ne,' energy points.'
  call wlog(slog)

  ! Transform to single precision
  rat(:,1:nat) = real (ratdbl(:,1:nat))
  rnrmax = real(rnrmav)
  temper = real(tk)
  thetax = real(thetad)
  sig2 = real(sigma2)
  idwopx = idwopt

  ! Prepare major work arrays:
  if (ispace.eq.0) then 
     !KJ 8/06 added kprep, which is a reduced version of xprep
     call kprep(em,ne,(do_fms.eq.1))  !KJ only initialize structure factors if we're actually doing work here.
     inclus=0
	 if (ktype.eq.3) then
	    nkp_max=nka 
		nkp_total=0
		nkp_fixed=nkp 
		nkp_min=max(10,nkp_max/50)
		e_largest=-1000.d0
		do ie=1,ne
		   if(dble(em(ie)).gt.e_largest) e_largest=dble(em(ie))
		enddo
	 endif
  else
     call xprep(iph0, idwopx, nat, inclus, npot, iphat, rclust, rat, iz, rnrmax, temper, thetax, sig2, minv, rdirec)
  endif
  call wlog('xprep done')

  lfms = 0
! JK - Making lcalc true for now. This way we can use Lanczos and eels, reuse Lanczos calculated
! Greens function for different edges, etc. 9/2009
  lcalc(:) = .TRUE.


  if (inclus.gt.1 .or. (ispace.eq.0 .and. do_fms.eq.1) ) then  !KJ 8/06 added "or ispace" !KJ 6/13 added do_fms
     !      call fms for a cluster around central atom
     if (ispace.ne.0) write (slog,25) inclus
25   format('FMS for a cluster of ',i4,' atoms')
     if (master .and. ispace.ne.0) call wlog(slog)

     nsp = 1
     if (abs(ispin).eq.1 ) nsp = nspx
     !KJ new code 10-2014:
     !nsp = nspu
     !KJ Note that now nspx=nspu redirected in "use" statement
     !nsp = nspx
     if (nsp.gt.nspx) then
        call wlog('ERROR - FEFF wants to do a spin-polarized calculation in fmstot but nspx=1 (arrays are too small).')
        stop
     endif

!KJ  feffq here has the following constraint :   7-09
!KJ        if (nsp.gt.2) stop "nsp too large!"
!KJ  I'm assuming that this isn't a problem for fmstot ; will check again for getgtr!

     istart = this_process*iblock + 1   ! =1 for proc1, =2 for proc2, etc.
     isize = 0 ! relic - no longer used
     i_for_next_report=1  !for writing progress to stdout below

     do ii=istart,ne,numprocs*iblock   ! e.g. for 4 processors, proc1 does energy points 1,5,9,... , proc2 does 2,6,10, ...  etc.
        do ie = ii, MIN(ii+iblock-1,ne) ! ie==ii since currently iblock=1; so this is a "false" loop that can be ignored
           if (worker) par_type = 3
           do  isp = 1, nsp
              dck=sqrt(2*(em(ie)-eref(ie,isp)))
              rpart  = real( dble(dck))
              aipart = real(dimag(dck))
              ck(isp) = cmplx(rpart, aipart)
           enddo
           do ipp = 0,nph
              do isp = 1, nsp
                 do ill = -lmaxph(ipp), lmaxph(ipp)
                    rpart  = dble( ph(ie, ill, isp, ipp))
                    aipart = dimag(ph(ie, ill, isp, ipp)) 
                    xphase(isp, ill, ipp) = cmplx(rpart, aipart)
                 enddo
              enddo
           enddo
           ! Report progress to stdout:
           iverb=0
           ! Complicated test because in || calculation we don't know which node will get e.g. point "10".
           ! I only let master report because I don't want several nodes reporting at once (e.g. if numprocs large).
           if(ie.ge.i_for_next_report) then
              iverb=1 ! report to stdout if master
              if(i_for_next_report.eq.1) then
                 i_for_next_report=10
              else
                 i_for_next_report=i_for_next_report+i_report_granularity * max(int(numprocs/i_report_granularity),1)
              endif
           endif
           write(slog,'(a,i4,a,i4)') 'Energy point ',ie,'/',ne
           if(master .and. iverb.eq.1) call wlog(slog)

           ! Here at last real space and reciprocal space calculations separate :   KJ 8/06
           if (ispace.eq.0) then 
		      if(ktype.eq.3) then
		         !KJ 1-2012 adaptive grid that uses fewer points as we go away from the threshold.
				 e_relative = (dble(em(ie)-em(1))/dble(e_largest-eref(ie,1)))  ! goes from 0 -> 1 if use em(ie)-em(1) ; could also use em(ie)-eref(ie,1)
		         nkp=nkp_min + nint( (nkp_max - nkp_min) * ( 1.d0 - e_relative**1.8)**5.0 )  ! just a random formula - saves 35% time in graphite CK, essentially exact
		         !nkp=nkp_min + nint( (nkp_max - nkp_min) * ( 1.d0 - e_relative**1.1)**6.0 )  ! just a random formula - saves 50% time in graphite CK, but some loss accuracy
			     nkx=nkp
			     nky=0
			     nkz=0
                 call kmesh ! Generate a new k-mesh. 
				 write(14,'(i4,2(f13.6),2i6)') ie,dble(em(ie)),e_relative,nkp,nkx
				 nkp_total=nkp_total+nkp
                 ! Recalculate the sum of the k-mesh integration weights to normalize integrals :
                 sumweights=dble(0)
                 do i=1,nkp
                    sumweights=sumweights+weight(i)
                 enddo
		      endif
              call fmskspace(ispin,ck,xphase,ie,em(ie)-eref(ie,1),gg(:,:,:,ie),iverb,sigma2)
           else
              if (i_hubbard.eq.1) then 
                 call fms(lfms,nsp,ispin,inclus,npot,ck,lmaxph,xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc,gg(:,:,:,ie))
              else if (i_hubbard.eq.2) then
               do isp=1,nsp
                  TFrm(:,:,:,:)=Trans(:,:,isp,:,:)
                  TFrmInv(:,:,:,:)=InvTrans(:,:,isp,:,:)
                do ip=0,nph
!                 do il=0,lx
!                  do im1=1,2*l_hubbard+1
!                   do im2=1,2*l_hubbard+1
!                    TFrm(im1,im2,il,ip)=Trans(im1,im2,isp,il,ip)
!                    TFrmInv(im1,im2,il,ip)=InvTrans(im1,im2,isp,il,ip)
!                   enddo
!                  enddo
!                 enddo
                 do il=-lmaxph(ip), lmaxph(ip)
                  do  imm = (abs(il)**2+1),(abs(il)+1)**2
                   rpart  = dble(aph(ie, abs(il)+1, imm, isp, ip))
                   aipart = dimag(aph(ie, abs(il)+1, imm, isp, ip))                
                   xphase_m(isp, il, imm, ip) = cmplx(rpart, aipart)
                  enddo
                 enddo
                enddo 
                enddo
                
                call fms_h(lfms,nsp,ispinp,inclus,npot,ck,lmaxph,xphase_m, ie, iverb, rdirec, toler1, toler2, lcalc, &
                          gg_m, TFrm,TFrmInv,UseTFrm)
               gg(:,:,:,ie)=gg_m(:,:,:) 
              
              end if
              if (save_gg_slice) then
                ! delay allocation of these arrays to here since istate hasn't been set prior to calling fms
                if (.not. allocated(gg_slice)) then
                  allocate(gg_slice(ldim, istate, ne), gg_slice_loc(ldim, istate, ne))
                  allocate(gg_diag(ldim, ldim, inclus, ne), gg_diag_loc(ldim, ldim, inclus, ne))
                end if
                gg_slice(:,:,ie) = gg_full(:ldim,:istate)
                do iat=1,inclus
                  i = (iat-1)*ldim + 1
                  j = iat * ldim
                  !j = i + ldim!  Brian version EK I think it is not right since we need
                  ! slice up to ldim not ldim + 1 for example
                  gg_diag(:,:,iat,ie) = gg_full(i:j,i:j)
                end do
              end if
           endif
           if (worker) then
              NEPts = NEPts + 1
              ggloc(:,:,:,NEPts) = gg(:,:,:,ie)
              if (save_gg_slice) then
                gg_slice_loc(:,:,NEPts) = gg_slice(:,:,ie)
                gg_diag_loc(:,:,:,NEPts) = gg_diag(:,:,:,ie)
              end if
           end if

        enddo  ! ie
     enddo  ! ii
     if (worker) par_type = 2
  endif

  if (ktype.eq.3) then
     write(slog,*) 'time used ',(100*nkp_total)/(ne*nkp_fixed)	,'% of fixed grid'
     call wlog(slog)
  endif
  if (numprocs .gt. 1) then
     call seconds(wall_commst)
     ! Collect gg
     if (worker) then
        !-- Send pointers for gg buffer to master 
        iggSize = NEPts*ldim**2*(nphx+1)
        call par_send_int_scalar(iggSize,1,0,this_process)  !KJ 11-2011 added _scalar to enable interface checking at compile time
        if (iggSize.gt.0) call par_send_cmplx(ggloc,iggSize,0,this_process)

        if (save_gg_slice) then
          ! send central-atom slice of gg back
          igg_slice_size = NEPts*(ldim)*istate
          call par_send_int_scalar(igg_slice_size,1,0,this_process)
          if (igg_slice_size.gt.0) call par_send_cmplx(gg_slice_loc,igg_slice_size,0,this_process)

          ! send site-diagonal piece of gg back
          igg_diag_size = NEPts*ldim**2*inclus
          call par_send_int_scalar(igg_diag_size,1,0,this_process)
          if (igg_diag_size.gt.0) call par_send_cmplx(gg_diag_loc,igg_diag_size,0,this_process)
        end if
     else
        ! Set up map of (!KJ gg?) to processor 
        do is = 1,numprocs-1
           istart = is*iblock + 1 
           do ii=istart,ne,numprocs*iblock
              do ie = ii, MIN(ii+iblock-1,ne) 
                 map(ie) = is
              enddo
           enddo
        enddo
        do i = 1,numprocs-1
           !-- Receive pointers for gg buffer from i
           call par_recv_int_scalar(iggSize,1,i,i)  !KJ 11-2011 added _scalar to enable interface checking at compile time
           !-- Receive buffer from i
           ! J. Kas - added gg saving.
           if (iggSize .gt. 0) then
              call par_recv_cmplx(ggloc,iggSize,i,i)
              indx = 1
              do j = 1,ne
                 if (map(j) .eq. i) then
                    gg(:,:,:,j) = ggloc(:,:,:,indx)
                    indx = indx + 1
                 endif
              enddo
           end if

           !-- Now repeat for gg_slice and gg_diag --
           if (save_gg_slice) then
             call par_recv_int_scalar(igg_slice_size,1,i,i)
             if (igg_slice_size .gt. 0) then
                 call par_recv_cmplx(gg_slice_loc,igg_slice_size,i,i)
                 indx = 1
                 do j = 1,ne
                   if (map(j) .eq. i) then
                       gg_slice(:,:,j) = gg_slice_loc(:,:,indx)
                       indx = indx + 1
                   endif
                 enddo
             end if

             call par_recv_int_scalar(igg_diag_size,1,i,i)
             if (igg_diag_size .gt. 0) then
                 call par_recv_cmplx(gg_diag_loc,igg_diag_size,i,i)
                 indx = 1
                 do j = 1,ne
                   if (map(j) .eq. i) then
                       gg_diag(:,:,:,j) = gg_diag_loc(:,:,:,indx)
                       indx = indx + 1
                   endif
                 enddo
             end if
           end if
        enddo
     endif
     call seconds(wall_commend)
     wall_comm = wall_comm + wall_commend - wall_commst      
  end if

  ! Write gg.bin nspx*(lx+1)**2
563 continue
  IF(master) THEN
!     PRINT*, 'Writing gg.bin ...'
     DO ie = 1, ne
        IF(.FALSE.) THEN
           nphtmp = nph
        ELSE
           nphtmp = 0
        END IF
        DO iph = 0, nphtmp
           CALL Write2D('gg.bin', gg(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,iph,ie))
        END DO
     END DO
     CALL CloseFl('gg.bin')
!     PRINT*, 'Done writing gg.bin.'


    if (save_gg_slice) then
      call wlog("Saving gg_slice.bin")
      open(unit=8, file="gg_slice.bin", form="unformatted")
      write(8) ldim, istate, ne
      write(8) gg_slice
      close(8)

      open(unit=8, file="gg_diag.bin", form="unformatted")
      write(8) ldim, ldim, inclus, ne
      write(8) gg_diag
      close(8)
    end if

  END IF
  
  ! Deallocate local variables
  deallocate(ph,lmax,em, eref,potlbl)
  deallocate(map)
  deallocate(gg,gg_m)
  deallocate(ggloc)
  deallocate(ggremember)
 ! deallocate(pstab,etab,ptab)
  deallocate(xphase, lcalc)
  deallocate(aph)
  RETURN

!! * KJ for using SPRKKKR phase shifts, grab code block 2 at end of file and paste here.

end subroutine fmstot



! ** KJ code block 1 for using SPRKRR phase shifts:
!!********** KJ this whole section is my experimental stuff and can be ignored for normal runs ******
!  if (sprkkrpot.eq.1) then   !KJ 8-06 read phase shifts from SPRKKR code
!     open(77,file='phasefeff.dat',form='formatted')
!     do iph=1,nph
!        do ie=1,ne
!           write(77,2078) em(ie),cdsqrt(2*(em(ie)-eref(ie,1))),ph(ie,0:ltot,1,iph)
!        enddo
!     enddo
!     close(77)
!2078 format(100e14.5)
!     if(.not.fullpot) then
!        open(78,file='phaseforfeff.dat',form='formatted',status='old',err=1475)
!        do iph=1,nph
!           read(78,*) adsun,ipot,sprne,sprnk,sprnmj
!           nk=sprnk
!           if(sprne.ne.ne) write(*,*) 'sprne and ne differ',sprne,ne
!           if(ipot.ne.iph) write(*,*) 'ipot and iph differ',ipot,iph
!           do ie=1,sprne
!              read(78,2078) etab(ie),ptab(ie),((pstab(ie,ik,imj,iph),ik=1,nk),imj=1,3)
!           enddo
!        enddo
!        close(78)
!        ! Now copy them into the ph array :
!        call wlog('copying of phaseshifts is probably wrong.')
!        do ik=1,nk,2
!           ph(1:ne,(ik-1)/2,1,0:nph)=dcmplx(pstab(1:ne,ik,1,0:nph),dble(0))
!        enddo
!        ph(:,:,:,0)=ph(:,:,:,1)
!        call wlog('Just copying GS to core hole phase shifts.')
!
!        open(77,file='phasefeff2.dat',form='formatted')
!        do iph=1,nph
!           do ie=1,ne
!              write(77,2078) em(ie),eref(ie,1),cdsqrt(2*(em(ie)-eref(ie,1))),ph(ie,0:ltot,1,iph)
!           enddo
!        enddo
!        close(77)
!     else
!        open(78,file='TFORFEFF1.DAT',status='old',err=1475)
!        close(78)
!     endif
!  endif
!!********** KJ end my stuff that can be ignored for normal runs ******

! ** KJ code block 2 for using SPRKRR phase shifts:
!1475 continue
!  if(master) then
!  call wlog ('SPR-KKR must calculate phase shifts first.')
!  open(76,file='energygridfms.dat',form='formatted')
!  write(76,*) ne
!  do ie=1,ne
!     write(76,2078) ((em(ie)-eref(ie,isp)),isp=1,nspx),em(ie),eref(ie,1)
!  enddo
!  close(76)
!  call wlog ('Energy grid written to energygridfms.dat .')
!  open(77,file='SPRKKR_FEFF.inp',status='unknown')
!  write(77,*) 'CONTROL  DATASET = SPRKKR'
!  write(77,*) '         ADSI = DOS'
!  write(77,*) '         POTFIL = SPRKKR.pot'
!  write(77,*) '         PRINT = 0'
!  write(77,*) '         PARA'
!  write(77,*) 'MODE  NREL'
!  write(77,*) 'SITES  NL = {',maxl,'}'
!  write(77,*) 'TAU  BZINT= POINTS  NKTAB= 250'
!  write(77,*) 'ENERGY  GRID={8}   NE={',ne,'}'
!  write(77,*) '         EMIN=0.2  EMAX=1.0  ImE=0.01 Ry'
!  write(77,*) 'TASK  FEFF'
!  close(77)
!  call wlog('Default SPRKKKR input file SPRKKR_FEFF.inp has been written.')
!  call wlog('Please customize, run SPRKKR, and then rerun FEFF.')
!  endif
!
!  ! If this is a stop call, is the memory automatically released?
!  ! Deallocate local variables
!  deallocate(ph,lmax,em, eref,potlbl)
!  deallocate(map,gg,ggloc,ggremember)
!  deallocate(pstab,etab,ptab)
!  deallocate(xphase, lcalc)
!  if (save_gg_slice) then
!    deallocate(gg_slice, gg_slice_loc, gg_full)
!  end if
!  call WipeErrorfileAtFinish  !KJ considering this a valid stop
!  stop


