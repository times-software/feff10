
  subroutine fmsdos_h_step1 (rclust, lfms,iph0, lmaxph, nat, iphat, ratdbl, inclus,    &
        ne, nph, em, eref_sp, iz, ph, rdirec, toler1, toler2, gtr, gtr_m, gtr_off)

  use DimsMod, only: natx, nex, nphx=>nphu, lx, nspx=>nspu, ltot
  use controls,only : ispace
  use constants
  use par
  use errorfile
  use hubbard_inp, only: l_hubbard
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
  integer, intent(in) :: iz(0:nphx)

  complex*16, intent(inout) :: ph(nex,ltot+1,2,0:nphx)
  complex*16, intent(in) ::  eref_sp(nex,2)
  complex, intent(out) :: gtr(0:lx, 2, 0:nphx, nex)
  complex, intent(out) :: gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex)
  complex, intent(out) :: gtr_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,0:nphx,nex)

  real :: rat(3,natx)
  real :: rpart,aipart
  integer :: ihole,iph,lmaxphpass(0:nphx) !KJ the last one to avoid sloppy interfaces

  character*30  fname
  character*512 slog
  complex, parameter :: conis=(0,1)
  
  ! Added to satisfy implicit none
  integer :: iat,ipot,imj,nk,ik,istart,isize,iverb
  integer :: nsp,ispin,j,maxl,ip,ios,isp
  integer :: is,i,indx,length,maxlen,ixl,iyl
  integer :: il,ix,im,ii,ie,ipp,ill
  real*8  :: wall_prep,wall_yprep
  integer :: i_for_next_report
  integer, parameter :: i_report_granularity=20

  !     fms staff
  integer, allocatable :: mapMPI(:)
  logical, allocatable :: lcalc(:)

  complex, allocatable, dimension(:,:,:) :: gg
  complex, allocatable, dimension(:,:,:,:) :: gtrloc
  complex, allocatable, dimension(:,:,:,:,:) :: gtr_mloc
  complex, allocatable, dimension(:,:,:,:,:,:) :: gtr_offloc
  complex, allocatable :: xphase(:,:,:), ck(:)

  complex*16 dck
  complex bmat(-1:1,-2:2,-4:4,-4:4)
  integer lkap(-1:1), lkapq(-2:2)
  integer im1,im2,imm,length_off,maxlen_off,ixl_off,iyl_off,length_m,maxlen_m,ixl_m,iyl_m

 ! Allocate local variables
  allocate(lcalc(0:lx),  &
          gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx),           &
          xphase(nspx, -lx:lx, 0:nphx), ck(nspx),               &
          mapMPI(nex)     )
  allocate(gtrloc(0:lx, 2, 0:nphx, nex))
  allocate(gtr_offloc(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2, 0:nphx,nex), &
           gtr_mloc(0:lx,nspx*(lx+1)**2,2,0:nphx,nex) )

  gtr = cmplx(0,0)
  gtrloc = cmplx(0,0)

  gtr_m = cmplx(0,0)
  gtr_mloc=cmplx(0,0)
  gtr_off = cmplx(0,0)
  gtr_offloc = cmplx(0,0)

  if (rclust.le.0.0) goto 900

  write(slog,'(a,i5,a)') 'Using ', ne,' energy points in ldos_hubbard step 1'
  if(master .and. iph0.eq.0) call wlog(slog)  !write this message only on first pass

  do iat=1,nat
    do j=1,3
       rat(j,iat) = real (ratdbl(j,iat))
    enddo
  enddo

  call seconds(wall_prep)
  if(ispace.eq.0) call par_stop('ERROR - KSPACE calculations not allowed with HUBBARD card.')
  inclus=0
  call yprep(iph0, nat, inclus, nph, iphat, rclust, rat )
  call seconds(wall_yprep)
  wall_yprep = wall_yprep - wall_prep
  i_for_next_report=1
  if (inclus.le.1) goto 900

     ! call fms for a cluster around central atom
     write (slog,'("FMS for a cluster of ",i3," atoms around atom type ",i2)') inclus, iph0
     if(master) call wlog(slog)


     do is=1,2

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
                 rpart  = dble( ph(ie, abs(ill)+1, is, ipp))
                 aipart = dimag(ph(ie, abs(ill)+1, is, ipp))
                 xphase(1, ill, ipp) = cmplx(rpart, aipart)
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

          lcalc(0:lx) = .true.
!          call fms(lfms, 1, 0, inclus, nph, ck, lmaxph, xphase, ie, iverb, 0, rdirec, toler1, toler2, lcalc, gg)
          call fms(1, 1, 0, inclus, nph, ck, lmaxphpass, xphase, ie, iverb, 0, rdirec, toler1, toler2, lcalc, gg)

          if (worker) isize = isize + 1

          do ip=0,nph
              do il=0,lmaxph(ip)
                 ix = il**2
                 do im=1,2*il+1
                    gtr(il,is,ip,ie)=gtr(il,is,ip,ie)+gg(ix+im,ix+im,ip)
                    gtr_m(il,ix+im,is,ip,ie)=gg(ix+im,ix+im,ip)
                 enddo

                 if (il.eq.l_hubbard) then
                    ix = il**2
                    do im1=1,2*il+1
                       do im2=1,2*il+1
                            gtr_off(il,ix+im1,ix+im2,is,ip,ie)= gg(ix+im1,ix+im2,ip)
                       enddo
                    enddo
                 endif

                 gtr(il,is,ip,ie)= gtr(il,is,ip,ie) * exp(2*conis*xphase(1, il,ip))/(2*il+1)
                 if (worker) gtrloc(il,is,ip,isize) = gtr(il,is,ip,ie)
                 do im=1,2*il+1
                       gtr_m(il,ix+im,is,ip,ie)= gtr_m(il,ix+im,is,ip,ie) * exp(2*conis*xphase(1,il,ip))/(2*il+1)
                       if (worker) gtr_mloc(il,ix+im,is,ip,isize) =  gtr_m(il,ix+im,is,ip,ie)
                 enddo

                 if(il.eq.l_hubbard) then
                      do im1=1,2*il+1
                         do im2=1,2*il+1
                            gtr_off(il,ix+im1,ix+im2,is,ip,ie) =  gtr_off(il,ix+im1,ix+im2,is,ip,ie) * &
                                  exp(2*conis*xphase(1,il,ip))/(2*il+1)
                            if(worker) gtr_offloc(il,ix+im1,ix+im2,is,ip,isize) = gtr_off(il,ix+im1,ix+im2,is,ip,ie)
                         enddo
                      enddo
                 endif

              enddo ! il
          enddo ! ip

        enddo  !ii
        enddo  !ie double energy loop

     enddo  ! is=1,2
	 if (worker) par_type = 2


       if (numprocs .gt. 1) then   ! if parallel calculation
            call seconds(wall_commst)

            length = (lx + 1) * (nphx + 1) * 2
            maxlen = (lx+1) * (nphx + 1) * 2
            ixl = length * isize
            iyl = maxlen * ne

            if (worker) then
               !-- Send pointers for gtr buffer to master
               !-- Send buffer
               call par_send_int_scalar(ixl,1,0,this_process)
               if (ixl .ne. 0) then  
                  call par_send_cmplx(gtrloc,ixl,0,this_process)
               endif

            else  ! master

               ! Set up map of gtr to processor
               do is = 1,numprocs-1
                  istart = is*iblock + 1
                  do ii=istart,ne,numprocs*iblock
                     do ie = ii, MIN(ii+iblock-1,ne)
                        mapMPI(ie) = is
                     enddo
                  enddo
               enddo
               do i = 1,numprocs-1
                  !-- Receive pointers for gtr buffer from i
                  !-- Receive buffer from i
                  call par_recv_int_scalar(ixl,1,i,i)
                  if (ixl .ne. 0) then
                     call par_recv_cmplx(gtrloc,ixl,i,i)
                     indx = 1
                     do j = 1,ne
                        if (mapMPI(j) .eq. i) then
                           do ip=0,nph
                              do il=0,lmaxph(ip)
                                 do is=1,2
                                    gtr(il,is,ip,j) = gtrloc(il,is,ip,indx)
                                 enddo
                              enddo
                           enddo
                           indx = indx + 1
                        endif
                     enddo
                  endif
               enddo
            endif  ! worker/master?
            call seconds(wall_commend)
            wall_comm = wall_comm + wall_commend - wall_commst
       endif

     if (numprocs .gt. 1) then
        call seconds(wall_commst)

        length_off = (lx+1)*((lx+1)**2)*((lx+1)**2)*2*(nphx+1)  
        maxlen_off = ((lx+1)**2)*((lx+1)**2)*2*(nphx+1)*(lx+1) 
        ixl_off =  length_off * isize 
        iyl_off =  maxlen_off * ne 

        if (worker) then
            !-- Send pointers for gtr buffer to master
            !-- Send buffer
            call par_send_int_scalar(ixl_off,1,0,this_process)
            if (ixl_off.ne.0) then
               call par_send_cmplx(gtr_offloc,ixl_off,0,this_process)
            endif
        else  ! master
            ! Set up map of gtr to processor
            do is = 1,numprocs-1
                istart = is*iblock + 1
                do ii=istart,ne,numprocs*iblock
                   do ie = ii, MIN(ii+iblock-1,ne)
                      mapMPI(ie) = is
                   enddo
                enddo
            enddo
            do i = 1,numprocs-1
                !-- Receive pointers for gtr buffer from i
                !-- Receive buffer from i
                call par_recv_int_scalar(ixl_off,1,i,i)
                if (ixl_off.ne.0) then
                   call par_recv_cmplx(gtr_offloc,ixl_off,i,i)
                   indx = 1
                   do j = 1,ne
                       if (mapMPI(j) .eq. i) then
                          do ip=0,nph
                             do il=0,lmaxph(ip)
                                if(il.eq.l_hubbard) then 
                                   do is=1,2
                                      do im1=1,(l_hubbard+1)**2
                                         do im2=1,(l_hubbard+1)**2
                                            gtr_off(il,im1,im2,is,ip,j) = gtr_offloc(il,im1,im2,is,ip,indx)
                                         enddo
                                      enddo
                                   enddo
                                endif
                             enddo
                          enddo
                          indx = indx + 1
                       endif
                   enddo
                endif
            enddo
        endif  ! worker or master?
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
     endif

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
        ! Set up map of gtr to processor
        do is = 1,numprocs-1
            istart = is*iblock + 1
            do ii=istart,ne,numprocs*iblock
               do ie = ii, MIN(ii+iblock-1,ne)
                  mapMPI(ie) = is
               enddo
            enddo
        enddo
        do i = 1,numprocs-1
            !-- Receive pointers for gtr buffer from i
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
  deallocate(gtrloc,gtr_mloc,gtr_offloc)
  deallocate(gg)
  deallocate(lcalc,xphase,mapMPI,ck)

  return
  end


