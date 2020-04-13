
  subroutine fmsdos_h_step2 (rclust, lfms,iph0, &
        lmaxph, nat, iphat, ratdbl, inclus,    &
        ne, nph, em, eref_sp, aph_sp,        &
        rdirec, toler1, toler2, Trans, InvTrans, gtr_m)

! This is a FMS routine for the second step of the Hubbard LDOS calculation.
! Looking at the program flow in ldossub_h (which calls fmsdos_h),
! only gtr_m is needed as output.
! Therefore, it is not necessary to track/keep/MPI-copy gtr or gtr_off.
! As a result, this routine may look a bit off-canon compared to other FMS drivers.

  use DimsMod, only: natx, nex, nphx=>nphu, lx, nspx=>nspu, ltot
  use controls,only : ispace
  use constants
  use par
  use errorfile
  use hubbard_inp
  implicit none

  real*8 wall_commend, wall_commst
  integer, parameter :: iblock=1

  ! Input
  integer, intent(in) :: ne,lfms
  integer, intent(in) :: iph0,nat,nph
  integer, intent(out) :: inclus
  real,    intent(in) :: rdirec,toler1,toler2,rclust
  integer, intent(in) :: iphat(natx)
  real*8,  intent(in) :: ratdbl(3,natx)
  complex*16, intent(in) :: em(nex)
  integer, intent(in) :: lmaxph(0:nphx)
  complex*16, intent(in) :: aph_sp(nex,lx+1,(lx+1)**2,2,0:nphx)
  complex*16, intent(in) ::  eref_sp(nex,2)
  complex, intent(in) :: Trans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
  complex, intent(in) :: InvTrans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx)
  complex, intent(out) :: gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex)

  real :: rat(3,natx)
  real :: rpart,aipart
  integer :: ihole,iph,lmaxphpass(0:nphx) !KJ the last one to avoid sloppy interfaces
  character*512 slog
  complex, parameter :: conis=(0,1)
  
  ! Added to satisfy implicit none
  integer :: iat
  integer :: ipot,imj,nk,ik,istart,isize,iverb
  integer :: nsp,ispin,j,maxl,ip,ios,isp
  integer :: is,i,indx,length,maxlen,ixl,iyl
  integer :: il,ix,im,ii,ie,ipp,ill
  real*8  :: wall_prep,wall_yprep
  integer :: i_for_next_report
  integer, parameter :: i_report_granularity=20

  integer mapMPI(nex)
  logical, allocatable :: lcalc(:)

  complex, allocatable, dimension(:,:,:) :: gg_m
  complex, allocatable, dimension(:,:,:,:,:) :: gtr_mloc
  complex, allocatable :: xphase_m(:,:,:,:), ck(:)

  logical   UseTFrm(0:lx,0:nphx)
  complex, allocatable, dimension(:,:,:,:) :: TFrm, TFrmInv
  complex*16 dck
  integer im1,im2,imm,length_off,maxlen_off,ixl_off,iyl_off,length_m,maxlen_m,ixl_m,iyl_m

 ! Allocate local variables
  allocate( lcalc(0:lx), ck(nspx) )
  allocate( gg_m(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx), &
          xphase_m(nspx,-lx:lx,(lx+1)**2,0:nphx))
  allocate( TFrm(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx), &
          TFrmInv(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx) )
  allocate ( gtr_mloc(0:lx,nspx*(lx+1)**2,2,0:nphx,nex))
  UseTFrm(:,:) = .false.
  !KJ Again here we hardcode that the HUBBARD term is on ATOM 1:
  UseTFrm(l_hubbard,1) = .true.

  inclus=0
  gtr_m = cmplx(0,0)
  gtr_mloc=cmplx(0,0)

  write(slog,'(a,i5,a)') 'Using ', ne,' energy points in ldos_hubbard step 2'
  if(master .and. iph0.eq.0) call wlog(slog)  !write this message only on first pass

  do iat=1,nat
    do j=1,3
       rat(j,iat) = real (ratdbl(j,iat))
    enddo
  enddo

  call seconds(wall_prep)
  if (ispace.eq.0) call par_stop('ERROR - KSPACE calculations not allowed with HUBBARD card.')
  call yprep(iph0, nat, inclus, nph, iphat, rclust, rat )
  call seconds(wall_yprep)
  wall_yprep = wall_yprep - wall_prep
  i_for_next_report=1
  if (inclus.le.1 ) goto 900

     ! call fms for a cluster around central atom
     write (slog,'("FMS for a cluster of ",i3," atoms around atom type ",i2)') inclus, iph0
     if(master) call wlog(slog)
     lcalc(0:lx) = .true.


     do is=1,2  ! spin up / spin down

        do ip=0,nphx
            do il=0,lx
               do im1=1,(2*l_hubbard+1)
                  do im2=1,(2*l_hubbard+1)
                     TFrm(im1,im2,il,ip)=Trans(im1,im2,is,il,ip)
                     TFrmInv(im1,im2,il,ip)=InvTrans(im1,im2,is,il,ip)
                  enddo
               enddo
            enddo
        enddo


        istart = this_process*iblock + 1
	    isize = 0
        do ii=istart,ne,numprocs*iblock
	    do ie = ii, min(ii+iblock-1,ne) ! ineffective loop, as iblock=1 => ie=ii always
	      if (worker) par_type = 3
          dck=sqrt(2*(em(ie)-eref_sp(ie,is)))
          rpart  = real( dble(dck))
          aipart = real(dimag(dck))
          ck(1) = cmplx(rpart, aipart)
          do ipp = 0,nph
              do ill = -lmaxph(ipp), lmaxph(ipp)
                    do imm = (abs(ill))**2+1,(abs(ill)+1)**2
                       rpart  = dble( aph_sp(ie, abs(ill)+1,imm,is,ipp))
                       aipart = dimag(aph_sp(ie, abs(ill)+1,imm,is,ipp))
                       xphase_m(1,ill, imm, ipp) = cmplx(rpart, aipart)
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

          ! calculate gg_m
          call fms_h(1, 1,0,inclus,nph,ck,lmaxphpass,xphase_m,ie,iverb,rdirec,toler1,toler2,lcalc,gg_m,  &
                 TFrm, TFrmInv, UseTFrm)
!          call fms_m(lfms, 1,0,inclus,nph,ck,lmaxph,xphase_m,ie,iverb,rdirec,toler1,toler2,lcalc,gg_m,  &
!                 TFrm, TFrmInv, UseTFrm)

          if (worker) isize = isize + 1
          !make gtr_m from gg_m_sp
          do ip=0,nph
              do il=0,lmaxph(ip)
                 ix = il**2
                 do im=1,2*il+1
                    gtr_m(il,ix+im,is,ip,ie)=  gg_m(ix+im,ix+im,ip)
                 enddo
                 do im=1,2*il+1
                       gtr_m(il,ix+im,is,ip,ie)= gtr_m(il,ix+im,is,ip,ie) * exp(2*conis*xphase_m(1,il,ix+im,ip))/(2*il+1)
                       if (worker) gtr_mloc(il,ix+im,is,ip,isize) =   gtr_m(il,ix+im,is,ip,ie)
                 enddo
              enddo ! il
          enddo ! ip

        enddo  !ii
        enddo  !ie double energy loop

     enddo  ! is=1,2  spin up / spin down
	 if (worker) par_type = 2




  if (numprocs .gt. 1) then
     call seconds(wall_commst)

     length_m=(lx+1)*(nphx+1)*((lx+1)**2)*2
     maxlen_m = (lx+1)*((lx+1)**2)*(nphx + 1)*2
     ixl_m = length_m * isize
     iyl_m = maxlen_m * ne
     if (worker) then
        !-- Send pointers for gtr_m buffer to master
        !-- Send buffer
        call par_send_int_scalar(ixl_m,1,0,this_process)
        if (ixl_m .ne. 0) call par_send_cmplx(gtr_mloc,ixl_m,0,this_process)
     else ! master
        ! Set up map of gtr_m to processor
        do is = 1,numprocs-1
            istart = is*iblock + 1
            do ii=istart,ne,numprocs*iblock
               do ie = ii, MIN(ii+iblock-1,ne)
                  mapMPI(ie) = is
               enddo
            enddo
        enddo
        do i = 1,numprocs-1
            !-- Receive pointers for gtr_m buffer from i
            !-- Receive buffer from i
            call par_recv_int_scalar(ixl_m,1,i,i)
            if (ixl_m .ne. 0) then
               call par_recv_cmplx(gtr_mloc,ixl_m,i,i)
               indx = 1
               do j = 1,ne
                  if (mapMPI(j) .eq. i) then
                     do ip=0,nph
                        do il=0,lmaxph(ip)
                           do im=(il**2+1),(il+1)**2
                              do is=1,2
                                 gtr_m(il,im,is,ip,j)=gtr_mloc(il,im,is,ip,indx)
                              enddo
                           enddo
                        enddo
                     enddo
                     indx = indx + 1
                  endif
               enddo
            endif
        enddo
     endif  ! master/slave?

     call seconds(wall_commend)
     wall_comm = wall_comm + wall_commend - wall_commst
  endif


  900  continue
  call par_barrier

  ! Deallocate local variables
  deallocate(gtr_mloc)
  deallocate(gg_m)
  deallocate(lcalc,xphase_m,ck)
  deallocate(TFrm, TFrmInv)


  return
  end

