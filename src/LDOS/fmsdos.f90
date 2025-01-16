!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmsdos.f90,v $:
! $Revision: 1.14 $
! $Author: jorissen $
! $Date: 2012/05/15 22:57:34 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fmsdos(ifms, rclust, lfms,iph0,idwopt,tk,thetad,sigma2,   &
     &            lmaxph, nat, iphat, ratdbl, inclus,                &
     &            ne, ne1, ne3, nph, em, eref, iz, ph,               &
     &            minv, rdirec, toler1, toler2 )
  !     uses Bruce Ravel subroutine to do FMS for LDOS
  !     written by alexei ankudinov 06.2000 from fmstot.f

  use DimsMod, only: nphx=>nphu, nex, ltot, nspx=>nspu, lx, natx
  use controls,only : ispace,sprkkrpot,sprkkrklist,fullpot
  use kklist !KJ for debugging only
  use constants
  use par
  use errorfile
  implicit none

  real*8 wall_commend, wall_commst
  integer, parameter :: iblock=1

  ! Args
  !tk,thetad,sigma2,

  ! Input
  integer, intent(in) :: ne,ne1,ne3,ifms,lfms
  integer, intent(in) :: iph0,idwopt,nat,nph
  integer, intent(out) :: inclus
  real,    intent(in) :: rdirec,toler1,toler2,rclust
  integer, intent(in), dimension(natx) :: iphat
  real*8,  intent(in), dimension(3,natx) :: ratdbl
  complex*16, intent(in) :: em(nex),eref(nex)
  integer, intent(in) :: lmaxph(0:nphx)
  integer, intent(in) :: iz(0:nphx)
  real*8,  intent(in) :: tk,thetad,sigma2


  complex*16, intent(inout) :: ph(nex, ltot+1, 0:nphx)
  integer, intent(inout) :: minv

  real :: rat(3,natx)
  real :: rpart,aipart
  integer :: ihole,iph,lmaxphpass(0:nphx) !KJ the last one to avoid sloppy interfaces

  character*30  fname
  character*512 slog
    complex, parameter :: conis=(0,1)
  
  ! Added to satisfy implicit none
  integer :: iat,idwopx
  integer :: ipot,imj,nk,ik,istart,isize,iverb
  integer :: nsp,ispin,j,ip,ios,isp
  integer :: is,i,indx,length,maxlen,ixl,iyl
  integer :: il,ix,im,ii,ie,ipp,ill
  real*8  :: wall_prep,wall_yprep
  integer :: i_for_next_report
  integer, parameter :: i_report_granularity=20

  !     fms staff
  integer, allocatable :: map(:)
  logical, allocatable :: lcalc(:)

  complex, allocatable, dimension(:,:,:) :: gg, gtr, gtrloc
  complex, allocatable :: xphase(:,:,:), ck(:)

  complex*16 dck
  complex bmat(-1:1,-2:2,-4:4,-4:4)
  integer :: lkap(-1:1), lkapq(-2:2)
  
  logical setkgrid !KJ debugging variable

  !KJ following variables for reading SPRKKR phase shifts :
  real*8, allocatable :: pstab(:,:,:,:)
  complex*16, allocatable :: etab(:),ptab(:)

  integer sprne,sprnk,sprnmj
  character adsun
  integer nktab
  real*8 maxi(3),mini(3) !debug
  !KJ
  integer nlgtr  !KJ for m-dos

  save

  ! Allocate local variables
  allocate( lcalc(0:lx),  &
       &   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx),           &
       &   xphase(nspx, -lx:lx, 0:nphx), ck(nspx),               &
       &   pstab(nex,2*lx+1,3,0:nphx),                           &
       &   map(nex),etab(nex),ptab(nex)      )

  !KJ allocate gtr dynamically ; for mt-potential, sum over m
  ! for full potential, save m-dependent information.
  ! CC, for printing out lm-projected dos
  fullpot = .true.
  if(fullpot) then
     nlgtr=(lx+1)**2
  else
     nlgtr=lx
  endif

  allocate(gtr(0:nlgtr,0:nphx,nex),gtrloc(0:nlgtr,0:nphx,nex))

  gtr = cmplx(0,0)
  gtrloc = cmplx(0,0)
  lmaxphpass=lmaxph 
  inclus=0 !KJ 4-2012 cosmetics

  if (.not.(rclust.le.0.0 .and. ispace.ne.0)) then  !KJ 8/06 added ispace

     !     ne = 1
     write(slog,82) 'Using ', ne,' energy points'
82   format(a,i5,a)
     if(master .and. iph0.eq.0) call wlog(slog)  !write this message only on first pass

     do iat=1,nat
        do j=1,3
           rat(j,iat) = real (ratdbl(j,iat))
        enddo
     enddo

     !     transform to single precision
     idwopx = -1

     call seconds(wall_prep)
     if (ispace.ne.0) then 
        call yprep(iph0, nat, inclus, nph, iphat, rclust, rat)
     endif

     call seconds(wall_yprep)
     wall_yprep = wall_yprep - wall_prep

!!!!KJ What follows is a big block of code for reading k-points and phase shifts from SPRKKR
!!!!   Commenting out for now as it's likely outdated and may never have worked well.
!     ! for debugging :
!     if((sprkkrklist.eq.1).and.ispace.eq.0) then
!        mini=dble(0);maxi=dble(0)
!        ! read k-mesh from file :
!        OPEN(92,file='sprkkrkmesh.txt',form='formatted')
!        READ(92,*) NKTAB
!        if (nktab.gt.nkp) stop 'nkp too small'
!        nkp=nktab
!        DO I=1,NKTAB
!           read(92,*) BK(:,I) !,WEIGHT(I)
!           weight(i)=dble(1)/dble(nktab)
!           do j=1,3;if(bk(j,i).lt.mini(j)) mini(j)=bk(j,i)
!              if(bk(j,i).gt.maxi(j)) maxi(j)=bk(j,i);
!           enddo
!        ENDDO
!        CLOSE(92)
!        write(*,*) maxi;write(*,*) mini  !;stop
!     endif
!
!     if (sprkkrpot.eq.1.and.iph0.eq.0) then 
!        !KJ 8-06 read phase shifts from  SPRKKR code  
!        !KJ iph0 : make sure it happens only once
!        open(77,file='phasefeff.dat',form='formatted')
!        do iph=1,nph
!           do ie=1,ne
!              write(77,2078) em(ie),cdsqrt(2*(em(ie)-eref(ie))), ph(ie,:,iph)
!           enddo
!        enddo
!        close(77)
!2078    format(100e14.5)
!        if(.not.fullpot) then
!           open(78,file='phaseforfeff.dat',form='formatted', status='old',err=1475)
!           do iph=1,nph
!              read(78,*) adsun,ipot,sprne,sprnk,sprnmj
!              nk=sprnk
!              if(sprne.ne.ne) write(*,*) 'sprne and ne differ',sprne,ne
!              if(ipot.ne.iph) write(*,*) 'ipot and iph differ',ipot,iph
!              do ie=1,sprne
!                 read(78,2078) etab(ie),ptab(ie), ((pstab(ie,ik,imj,iph),ik=1,nk),imj=1,3)
!              enddo
!           enddo
!           close(78)
!           ! now copy them into the ph array :
!           call wlog('copying of phaseshifts is probably wrong.')
!           do ik=1,nk,2
!              ph(1:ne,1+(ik-1)/2,0:nph)=dcmplx(pstab(1:ne,ik,1,0:nph),dble(0))
!           enddo
!           ph(:,:,0)=ph(:,:,1)
!           call wlog('Just copying GS to core hole phase shifts.')
!
!           open(77,file='phasefeff2.dat',form='formatted')
!           do iph=1,nph
!              do ie=1,ne
!                 write(77,2078) em(ie),eref(ie), cdsqrt(2*(em(ie)-eref(ie))),ph(ie,:,iph)
!              enddo
!           enddo
!           close(77)
!        else
!           open(78,file='TFORFEFF1.DAT',status='old',err=1475)
!           close(78)
!        endif
!     endif
!
!     ! !KJ end my changes


     if (inclus.gt.1 .or. ispace.eq.0) then  !KJ 8/06 added "or ispace"

        !c      call fms for a cluster around central atom
        if(ispace.ne.0) then
          write (slog,35) inclus, iph0
35        format ('FMS for a cluster of ',i3,' atoms around atom type ',i2)
        else
          write(slog,36) iph0
36        format('FMS for atom iph = ',i2)
        endif
        if(master)call wlog (slog)
 !       call wlog (' Please, wait (updates every 20 points) ...')

        istart = this_process*iblock + 1
        isize = 0 !relic - not used
        i_for_next_report=1
        
        II_LOOP: do ii=istart,ne,numprocs*iblock
           IE_LOOP: do ie = ii, min(ii+iblock-1,ne)  ! This loop is inconsequential as iblock=1.  Hence ie=ii
              if (worker) par_type = 3
              dck=sqrt(2*(em(ie)-eref(ie)))
              rpart  = real( dble(dck))
              aipart = real(dimag(dck))
              ck(1) = cmplx(rpart, aipart)
              do ipp = 0,nph
                 do ill = -lmaxph(ipp), lmaxph(ipp)
                    rpart  = dble( ph(ie, abs(ill)+1, ipp))
                    aipart = dimag(ph(ie, abs(ill)+1, ipp))
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


              nsp = 1
              ispin = 0
              do ill = 0,lx
                 lcalc(ill) = .true.
              enddo
              minv = 0

              if (ispace.eq.0) then ! KJ
                 call fmskspace(ispin,ck,xphase,ie,em(ie)-eref(ie),gg,iverb,dble(0))  ! KJ
              else           ! KJ
                 call fms(lfms,nsp,ispin,inclus,nph,ck,lmaxphpass,xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc,gg) !KJ 11-2011 lmaxph->lmaxphpass
              endif      !  KJ

              if (worker) then
                 isize = isize + 1
              endif
              IP_LOOP: do ip=0,nph
                 IL_LOOP: do il=0,lmaxph(ip)
                    ix = il**2
                    do im=1,2*il+1
                       if(fullpot) then
                          gtr(ix+im,ip,ie)=gg(ix+im,ix+im,ip)*  exp(2*conis*xphase(1, il,ip))/(2*il+1)
                          if (worker) gtrloc(ix+im,ip,isize) = gtr(ix+im,ip,ie)
                       else
                          gtr(il,ip,ie)=gtr(il,ip,ie)+ gg(ix+im,ix+im,ip)
                       endif
                    enddo
                    if(.not.fullpot) gtr(il,ip,ie)= gtr(il,ip,ie)* exp(2*conis*xphase(1, il,ip))/(2*il+1)
                    if (worker.and.(.not.fullpot)) then 
                       gtrloc(il,ip,isize) = gtr(il,ip,ie)
                    endif
                 enddo IL_LOOP
              enddo IP_LOOP
           enddo IE_LOOP
        enddo II_LOOP

        if (worker) par_type = 2

     endif

     if (numprocs .gt. 1) then
        call seconds(wall_commst)
        length = (nlgtr + 1) * (nphx + 1)  !KJ
        maxlen = (nlgtr + 1) * (nphx + 1)  !KJ
        !        length = (lx + 1) * (nphx + 1)
        !        maxlen = (lx + 1) * (nphx + 1)
        ixl = length * isize
        iyl = maxlen * ne
        if (worker) then
           !-- Send pointers for gtr buffer to master
           call par_send_int_scalar (ixl,1,0,this_process)   !KJ added _scalar
           !-- Send buffer
           if (ixl .ne. 0)  call par_send_cmplx(gtrloc,ixl,0,this_process)
        else
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
              call par_recv_int_scalar(ixl,1,i,i)  !KJ added _scalar
              !-- Receive buffer from i
              if (ixl .ne. 0) then
                 call par_recv_cmplx(gtrloc,ixl,i,i)
                 indx = 1
                 do j = 1,ne
                    if (map(j) .eq. i) then
                       do ip=0,nph
                          do il=0,nlgtr !KJ lmaxph(ip)
                             gtr(il,ip,j) = gtrloc(il,ip,indx) 
                          enddo
                       enddo
                       indx = indx + 1
                    endif
                 enddo
              endif
           enddo
        endif
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
     endif
  endif

  !--   write gtr.bin

  if (master) then
     write(fname,920)  iph0
920  format('gtr', i2.2, '.bin')
     open (unit=3, file=fname, status='unknown', access='sequential', form='unformatted', iostat=ios)
     write(3) ne, ne1, ne3, nph, ifms
     write(3) (((gtr(il,ip,ie), il=0,nlgtr), ip=0,nph), ie=1,ne) !KJ
     !        write(3) (((gtr(il,ip,ie), il=0,lx), ip=0,nph), ie=1,ne)
     close (unit=3)
  endif
  call par_barrier

  ! Deallocate local variables
  deallocate(gtr,gtrloc)
  deallocate(lcalc,gg,xphase,pstab,map,ck,etab,ptab)
  fullpot = .false.
  return

!1475 continue
!  call wlog ('SPR-KKR must calculate phase shifts first.')
!  open(76,file='energygridldos.dat',form='formatted')
!  write(76,*) ne
!  do ie=1,ne
!     write(76,2078) ((em(ie)-eref(ie)),isp=1,nspx),em(ie), eref(ie)
!  enddo
!  close(76)
!  call wlog ('Energy grid written to energygridldos.dat .')
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

  ! Deallocate local variables
  !deallocate(gtr,gtrloc)
  !deallocate(lcalc,gg,xphase,pstab,map,ck,etab,ptab)
  !call WipeErrorfileAtFinish !I consider this a regular stop
  !stop

end subroutine fmsdos
