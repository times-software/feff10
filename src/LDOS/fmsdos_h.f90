  subroutine fmsdos_h(ifms, rclust, lfms,iph0,idwopt,tk,thetad,sigma2, &
        lmaxph, nat, iphat, ratdbl, inclus,                           &
        ne, ne1, ne3, nph, em, eref,eref_sp, iz, ph_sp,aph_sp,        &
        minv, rdirec, toler1, toler2,Vnlm,i_opt,xmu)

  use DimsMod, only: natx, nex, nphx=>nphu, lx, nspx=>nspu, ltot
  use controls,only : ispace,sprkkrpot,sprkkrklist,fullpot
  use constants
  use par
  use errorfile
  use hubbard_inp
  implicit none

  real*8 wall_commend, wall_commst
  integer, parameter :: iblock=1

  ! Input
  integer, intent(in) :: ne,ne1,ne3,ifms,lfms
  integer, intent(in) :: iph0,idwopt,nat,nph
  integer, intent(out) :: inclus
  real,    intent(in) :: rdirec,toler1,toler2,rclust
  integer, intent(in) :: iphat(natx)
  real*8,  intent(in) :: ratdbl(3,natx)
  complex*16, intent(in) :: em(nex),eref(nex)
  integer, intent(in) :: lmaxph(0:nphx)
  integer, intent(in) :: iz(0:nphx)
  real*8,  intent(in) :: tk,thetad,sigma2
  real*8, intent(in) :: Vnlm(0:lx,(lx+1)**2,2,0:nphx)

  complex*16, intent(inout) :: ph_sp(nex,ltot+1,2,0:nphx)
  integer, intent(inout) :: minv
  complex*16, intent(inout) :: aph_sp(nex,lx+1,(lx+1)**2,2,0:nphx)
  complex*16, intent(in) ::  eref_sp(nex,2)
  integer, intent(in) :: i_opt
  real*8, intent(in) :: xmu

  real :: rat(3,natx)
  real :: rpart,aipart
  integer :: ihole,iph,lmaxphpass(0:nphx) !KJ the last one to avoid sloppy interfaces

  character*30  fname
  character*512 slog
  complex, parameter :: conis=(0,1)
  
  ! Added to satisfy implicit none
  integer :: iat,idwopx
  integer :: ipot,imj,nk,ik,istart,isize,iverb
  integer :: nsp,ispin,j,maxl,ip,ios,isp
  integer :: is,i,indx,length,maxlen,ixl,iyl
  integer :: il,ix,im,ii,ie,ipp,ill
  real*8  :: wall_prep,wall_yprep
  integer :: i_for_next_report
  integer, parameter :: i_report_granularity=20

  !     fms staff
  integer, allocatable :: map(:)
  logical, allocatable :: lcalc(:)

  complex, allocatable, dimension(:,:,:) :: gg, gg_m
  complex, allocatable, dimension(:,:,:,:) :: gg_sp, gg_m_sp, gtr, gtrloc
  complex, allocatable, dimension(:,:,:,:,:) :: gtr_m, gtr_mloc
  complex, allocatable, dimension(:,:,:,:,:,:) :: gtr_off, gtr_offloc
  complex, allocatable :: xphase(:,:,:), xphase_m(:,:,:,:), ck(:)

  !J.K. Adding transformation matrices, TFrm and TFrmInv
  !For now, text with p states (3X3 matix). Later these will be passed as input.
  logical   UseTFrm(0:lx,0:nphx)
  complex, allocatable, dimension(:,:,:,:) :: TFrm, TFrmInv
  complex, allocatable, dimension(:,:,:,:,:) :: Trans, InvTrans
  complex*16 dck
  complex bmat(-1:1,-2:2,-4:4,-4:4)
  integer :: lkap(-1:1), lkapq(-2:2)
  integer im1,im2,imm,length_off,maxlen_off,ixl_off,iyl_off,length_m,maxlen_m,ixl_m,iyl_m

  save
 ! Allocate local variables
  allocate( lcalc(0:lx),  &
          gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx),           &
          xphase(nspx, -lx:lx, 0:nphx), ck(nspx),               &
          map(nex)     )
  allocate( gg_sp(nspx*(lx+1)**2, nspx*(lx+1)**2, 2, 0:nphx), &
          gg_m(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx), &
          gg_m_sp(nspx*(lx+1)**2, nspx*(lx+1)**2, 2, 0:nphx),&
          xphase_m(nspx,-lx:lx,(lx+1)**2,0:nphx))

  allocate(   TFrm(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx), TFrmInv(2*l_hubbard+1,2*l_hubbard+1,0:lx,0:nphx), &
  Trans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx),InvTrans(2*l_hubbard+1,2*l_hubbard+1,2,0:lx,0:nphx) )

  allocate(gtr(0:lx, 2, 0:nphx, nex), gtrloc(0:lx, 2, 0:nphx, nex))
  allocate(gtr_m(0:lx,nspx*(lx+1)**2,2,0:nphx,nex), &
           gtr_off(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2,0:nphx,nex), &
           gtr_offloc(0:lx,nspx*(lx+1)**2,nspx*(lx+1)**2,2, 0:nphx,nex), &
           gtr_mloc(0:lx,nspx*(lx+1)**2,2,0:nphx,nex) )

  UseTFrm(:,:) = .False.
  UseTFrm(l_hubbard,1) = .true.

  gtr = cmplx(0,0)
  gtrloc = cmplx(0,0)
  inclus=0

  gtr_m = cmplx(0,0)
  gtr_mloc=cmplx(0,0)
  gtr_off = cmplx(0,0)
  gtr_offloc = cmplx(0,0)

  if (rclust.le.0.0 .and. ispace.ne.0) goto 900

  write(slog,'(a,i5,a)') 'Using ', ne,' energy points in ldos_hubbard'
  if(master .and. iph0.eq.0) call wlog(slog)  !write this message only on first pass

  do iat=1,nat
    do j=1,3
       rat(j,iat) = real (ratdbl(j,iat))
    enddo
  enddo

  idwopx = -1

  call seconds(wall_prep)
  if (ispace.ne.0)  call yprep(iph0, nat, inclus, nph, iphat, rclust, rat )
  call seconds(wall_yprep)
  wall_yprep = wall_yprep - wall_prep

  if (inclus.gt.1 .or. ispace.eq.0) then

     ! call fms for a cluster around central atom
     if(ispace.ne.0) then
        write (slog,'("FMS for a cluster of ",i3," atoms around atom type ",i2)') inclus, iph0
     else
        write(slog,'("FMS for atom iph = ",i2)') iph0
     endif
     if(master) call wlog(slog)

     if(i_opt.eq.2) then
        open(file='transf.dat',unit=63,status='unknown',iostat=ios)
        open(file='Invtransf.dat',unit=64,status='unknown',iostat=ios)
        do is=1,2
           do ip=0,nphx
              do il=0,lx
                 do im1=1,2*l_hubbard+1
                    read(63,*)(Trans(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
                    read(64,*)(InvTrans(im1,im2,is,il,ip),im2=1,(2*l_hubbard+1))
                 enddo
              enddo
           enddo
        enddo
        close(63)
        close(64)
     endif

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
                 rpart  = dble( ph_sp(ie, abs(ill)+1, is, ipp))
                 aipart = dimag(ph_sp(ie, abs(ill)+1, is, ipp))
                 xphase(1, ill, ipp) = cmplx(rpart, aipart)

                 if(i_opt.eq.2) then
                    do imm = (abs(ill))**2+1,(abs(ill)+1)**2
                       rpart  = dble( aph_sp(ie, abs(ill)+1,imm,is,ipp))
                       aipart = dimag(aph_sp(ie, abs(ill)+1,imm,is,ipp))
                       xphase_m(1,ill, imm, ipp) = cmplx(rpart, aipart)
                    enddo
                 endif
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

          nsp = 1
          ispin = 0
          lcalc(0:lx) = .true.
          minv = 0

          if(i_opt.eq.1) then
             call fms(lfms, nsp, ispin, inclus, nph, ck, lmaxphpass, xphase,ie,iverb, minv, rdirec, toler1, toler2, lcalc, gg)
          elseif(i_opt.eq.2) then
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
             call fms_h(lfms, nsp,ispin,inclus,nph,ck,lmaxphpass,xphase_m,ie,iverb,minv,rdirec,toler1,toler2,lcalc,gg_m,  &
                     TFrm, TFrmInv, UseTFrm)
          endif

          do ip=0,nph
              do il=0,lmaxph(ip)
                 ix = il**2
                 do im1=1,2*il+1
                    do im2=1,2*il+1
                        if(i_opt.eq.1)    gg_sp(ix+im1,ix+im2,is,ip)= gg(ix+im1,ix+im2,ip)
                        if(i_opt.eq.2)    gg_m_sp(ix+im1,ix+im2,is,ip)=gg_m(ix+im1,ix+im2,ip)
                    enddo
                 enddo
              enddo
          enddo

          if (worker) isize = isize + 1

          do ip=0,nph
              do il=0,lmaxph(ip)
                 ix = il**2
                 do im=1,2*il+1
                    gtr(il,is,ip,ie)=gtr(il,is,ip,ie)+gg_sp(ix+im,ix+im,is,ip)
                    if(i_opt.eq.1) then
                       gtr_m(il,ix+im,is,ip,ie)=0.0*gtr_m(il,ix+im,is,ip,ie)+ gg_sp(ix+im,ix+im,is,ip)
                    elseif(i_opt.eq.2)then
                       gtr_m(il,ix+im,is,ip,ie)=0.0*gtr_m(il,ix+im,is,ip,ie)+  gg_m_sp(ix+im,ix+im,is,ip)
                    endif 
                 enddo

                 ! T.A. filling out the offdiagonal terms of the gg_matrix
                 if(i_opt.eq.1) then
                     if (il.eq.l_hubbard) then
                        ix = il**2
                        do im1=1,2*il+1
                           do im2=1,2*il+1
                                gtr_off(il,ix+im1,ix+im2,is,ip,ie)= gg_sp(ix+im1,ix+im2,is,ip)
                           enddo
                        enddo
                     endif
                 endif

                 gtr(il,is,ip,ie)= gtr(il,is,ip,ie)*  exp(2*conis*xphase(1, il,ip))/(2*il+1)
                 if (worker) gtrloc(il,is,ip,isize) = gtr(il,is,ip,ie)
                 do im=1,2*il+1
                    if(i_opt.eq.1) then
                       gtr_m(il,ix+im,is,ip,ie)= gtr_m(il,ix+im,is,ip,ie)*  exp(2*conis*xphase(1,il,ip))/(2*il+1)
                       if (worker) gtr_mloc(il,ix+im,is,ip,isize) =  gtr_m(il,ix+im,is,ip,ie)
                    elseif(i_opt.eq.2) then
                       gtr_m(il,ix+im,is,ip,ie)= gtr_m(il,ix+im,is,ip,ie)*   exp(2*conis*xphase_m(1,il,ix+im,ip))/(2*il+1)
                       if (worker) gtr_mloc(il,ix+im,is,ip,isize) =   gtr_m(il,ix+im,is,ip,ie)
                    endif
                 enddo

                 if(i_opt.eq.1) then
                     if(il.eq.l_hubbard) then
                          do im1=1,2*il+1
                             do im2=1,2*il+1
                                gtr_off(il,ix+im1,ix+im2,is,ip,ie) =   gtr_off(il,ix+im1,ix+im2,is,ip,ie)* &
                                      exp(2*conis*xphase(1,il,ip))/(2*il+1)
                                if(worker) gtr_offloc(il,ix+im1,ix+im2,is,ip,isize)=  gtr_off(il,ix+im1,ix+im2,is,ip,ie)
                             enddo
                          enddo
                     endif
                 endif

              enddo ! il
          enddo ! ip

!        call par_barrier

        enddo  !ii
        enddo  !ie double energy loop
      !KJ call par_barrier

     enddo  ! is=1,2
	 if (worker) par_type = 2

  endif  ! main "if" : do something or not?

  if(i_opt.eq.1) then
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
                        map(ie) = is
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
                        if (map(j) .eq. i) then
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
  endif

  if(i_opt.eq.1) then
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
                      map(ie) = is
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
                       if (map(j) .eq. i) then
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
                  map(ie) = is
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
                  if (map(j) .eq. i) then
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


!--   write gtr.bin
 900  continue
  if(i_opt.eq.1) then
     if (master) then
        write(fname,'("gtr", i2.2, ".bin")')  iph0
        open (unit=3, file=fname, status='unknown',   &  
              access='sequential', form='unformatted', iostat=ios)
        write(3) ne, ne1, ne3, nph, ifms
        write(3)((((gtr(il,is,ip,ie),il=0,lx),ip=0,nph),ie=1,ne),is=1,2)
        close (unit=3)

        write(fname,'("gtr_off", i2.2, ".bin")')  iph0
        open (unit=31, file=fname, status='unknown',  &  
             access='sequential', form='unformatted', iostat=ios)
        write(31) ne, ne1, ne3, nph, ifms
        write(31)((((((gtr_off(il,im1,im2,is,ip,ie),im1=1,(l_hubbard+1)**2) &
          ,im2=1,(l_hubbard+1)**2),ip=0,nph),ie=1,ne),is=1,2),il=0,lx)
        close (unit=31)

        if (iph0.eq.1 .and. lx.ge.3) then  !KJ I think this is just a debugging file
        open (unit=31, file='gtr_off.dat', form='formatted',access='append')
        write(31,*) 'NEW WRITE ******************'
        write(31,'(18(e12.5,1x,e12.5,3x))')((gtr_off(3,im1,im2,1,1,5),im1=10,16),im2=10,16)
        close (unit=31)
        endif
     endif
  endif
  if (master) then
     write(fname,'("gtr_m", i2.2, ".bin")')  iph0
     open (unit=13, file=fname, status='unknown',access='sequential', form='unformatted', iostat=ios)
     rewind(13)
     write(13) ne, ne1, ne3, nph, ifms
     write(13)(((((gtr_m(il,im,is,ip,ie), im=(il**2)+1,(il+1)**2),il=0,lx), ip=0,nph),ie=1,ne),is=1,2)
     close (unit=13)
  endif
  call par_barrier

  ! Deallocate local variables
  deallocate(gtr,gtrloc,gtr_m,gtr_mloc,gtr_off,gtr_offloc)
  deallocate(gg,gg_m,gg_m_sp,gg_sp)
  deallocate(lcalc,xphase,xphase_m,map,ck)
  deallocate(TFrm, TFrmInv, Trans, InvTrans)


  return
  end

