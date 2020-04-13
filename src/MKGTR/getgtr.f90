!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getgtr.f90,v $:
! $Revision: 1.19 $
! $Author: hebhop $
! $Date: 2013/01/09 21:32:49 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getgtr
  !     Calculates FMS contribution to absorption
  !     uses Bruce Ravel subroutine to do FMS in self-consistency loop
  !     notice that it can do FMS with polarization dependence and always include l-->l-1 channel.
  !     written by alexei ankudinov 06.1997
  use DimsMod, only: nphx=>nphu, ltot, nspx=>nspu, nex, lx
  use IOMod
  use ErrorMod
  use constants
  use par
  use fms_inp,only: minv,rclust=>rfms2,lmaxph,ipr3 !KJ this and following from fmstot
  use global_inp,only: ipol,ispin,le2,angks,ptz,do_nrixs
  use atoms_inp,only: nat,iphat,ratdbl=>rat
  use nrixs_inp,only: jmax,kfinmax,jinit
  use eels_inp,only: ipmin,ipmax,ipstep,eels
  implicit none


  integer, parameter :: npadx=8
  integer :: knd(8), lnd(8)
  logical ltrace
  integer ne, ne1, ne3,  nph, ihole, iz(0:nphx)
  integer  nsp, ie
  integer ip,nip !KJ 1-06 added this variable - just local index

  ! Added to satisfy implicit none
  integer :: is1,is2,iph,k1,k2,m1,m2 ! Loop indecies
  integer :: ind, ix1,ix2,ik0,lmaxp1
  integer :: kinit,linit,L1,L2,ms1,ms2
  integer :: ios
  integer nphtmp
  real*8 :: aa,bb,rnrmav,xmu,edge

  ! Work space
  complex*16, allocatable :: ph(:,:,:,:)
  integer, allocatable    :: lmax(:,:)
 
  ! Complex energy grid emg is decomposed into em and eref to have
  !  the same structure in phase.bin
  complex*16, allocatable  :: em(:), eref(:,:)
  character*6, allocatable  :: potlbl(:)

  ! FMS staff
  complex, allocatable, dimension(:,:,:,:) :: gg
  complex, allocatable, dimension(:,:) :: gtr
  !KJ commented out  1-06    dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
  complex*16, allocatable :: bmat(:,:,:,:,:,:,:) !KJ added last index
  complex*16, allocatable :: bmat0(:,:,:,:,:,:)  !KJ new variable 1-06      
  complex*16, allocatable :: rkk(:,:,:)
  complex*16, allocatable :: dum(:)
  complex bmatsmall(2*lx+1,2*lx+1)


  ! Allocate local variables
  allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(lmax(nex, 0:nphx))
  allocate(em(nex), eref(nex, nspx))
  allocate(potlbl(0:nphx))
  allocate(gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx, nex))
  allocate(gtr(ipmin:ipmax,nex))
  allocate(bmat(-lx:lx,0:1,8, -lx:lx,0:1,8,ipmin:ipmax))
  allocate(rkk(nex,8,nspx))
  allocate(bmat0(-lx:lx,0:1,8, -lx:lx,0:1,8))
  allocate(dum(nex))


  nip = 0

  gtr = cmplx(0)
  bmat=dcmplx(0)
  ph = dcmplx(0)
  lmax = 0
  em = dcmplx(0)
  eref = dcmplx(0)
  potlbl = ''
  gg = cmplx(0)
  bmat0=dcmplx(0)
  rkk = dcmplx(0)
  dum = dcmplx(0)


  nsp = 2 
  if (abs(ispin).eq.1 ) nsp = nspx
  !KJ new treatment:
  nsp=nspx

  !     need less data than rphbin.f provides, also dimensions of ph array are different.
  call rdxsph (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
  call setkap (ihole, kinit, linit)  

  !KJ 1-06  I added the do-loop around the call to bcoef for ELNES calcul.
  ltrace = .FALSE. ! Josh - bug fix for feff9.
  do ip=ipmin,ipmax,ipstep
     if (eels.eq.1) call iniptz(ptz,ip,2)  !KJ Only change ptz for ELNES
     call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, knd, lnd, bmat0)
     bmat(:,:,:,:,:,:,ip)=bmat0(:,:,:,:,:,:)
  enddo

  ! Read gg
  DO ie = 1, ne
     IF(.FALSE.) THEN
        nphtmp = nph
     ELSE
        nphtmp = 0
     END IF
     DO iph = 0, nphtmp
        CALL Read2D('gg.bin', gg(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,iph,ie), L1, L2)
        ! Check that bounds are correct.
        IF(L1.ne.nspx*(lmaxph(iph)+1)**2.or.L2.ne.nspx*(lmaxph(iph)+1)**2) CALL Error('Error when reading gg.bin')
     END DO
  END DO
  CALL CloseFl('gg.bin')
  do ie=1,ne
  CALL Write2D('gg.dat', gg(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,0,ie))
  enddo
  if(ipr3.ge.4) then  !KJ debugging output
     do ms1=0,1 ; do k1=1,8 ; do ms2=0,1 ; do k2=1,8
         bmatsmall(1:2*lx+1,1:2*lx+1)=cmplx(bmat(-lx:lx,ms2,k2,-lx:lx,ms1,k1,ipmin))
	     call Write2D('bmat.dat',bmatsmall(1:2*lx+1,1:2*lx+1))
     enddo ; enddo ; enddo ; enddo
     do ie=1,nspx 
        call Write2D('rkk.dat',cmplx(rkk(1:ne,1:8,ie)))
     enddo
  endif



  ! Form gtr
  do ie = 1, ne
     do k1 = 1,8
        do is1 = 1,nsp
           do k2 = 1,8
              do is2 = 1,nsp
                 ix1 = nsp * ( lnd(k1)**2 +  lnd(k1) )
                 ix2 = nsp * ( lnd(k2)**2 +  lnd(k2) )
                 ms1 = is1 - 1
                 ms2 = is2 - 1
                 if (lnd(k2).ge.0 .and. lnd(k1).ge.0) then
                    do m1=-lnd(k1), lnd(k1)
                       do m2=-lnd(k2), lnd(k2)
                          do ip=ipmin,ipmax,ipstep
                             gtr(ip,ie) = gtr(ip,ie) +                       &
                                  & gg(ix1+nsp*m1+is1,ix2+nsp*m2+is2,0,ie) * &
                                    bmat(m2,ms2,k2, m1,ms1,k1,ip) * &
                                    rkk(ie,k1,is1)*rkk(ie,k2,is2)
!                             write(99,*) gtr(ip,ie)
                          enddo  ! ip
                       enddo  ! m2
                    enddo  ! m1
                 endif
              enddo  ! is2
           enddo  ! k2
        enddo  ! is1
     enddo  ! k1
  enddo   ! ie



  !     write fms.bin
  open (unit=1, file='fms.bin', status='unknown', iostat=ios)
  !        write title line
  write(1,105) rclust*bohr
105 format('FMS rfms=', f7.4)
  write(1, 110) ne, ne1, ne3,  nph, npadx, nip  !KJ added nip 1-06
110 format(6(1x,i7))  !KJ changed 5 to 6
  do ip=ipmin,ipmax      !KJ I added the do loop
     do ie = 1, ne
        aa = dble ( real ( gtr(ip,ie) ) ) !KJ I added ip index
        bb = dble (aimag ( gtr(ip,ie) ) ) !KJ ditto
        dum(ie) = dcmplx (aa,bb)
!        write(98,*) dum(ie)        
     end do
     call wrpadx(1, npadx, dum, ne)
  enddo                  !KJ end of my loop
!  call wrpadc(1, npadx, gtr, (ipmax-ipmin+1)*ne)
  close (unit=1)
  open(99,file='gtr.dat',form='formatted',status='unknown')
  do ie=1,ne
     write(99,'(10f13.6)') em(ie),gtr(ipmin,ie)
  enddo
  close(99)

  ! Deallocate local variables
  deallocate(ph,lmax,em,eref,potlbl,bmat)
  deallocate(gg,gtr)
  deallocate(rkk,bmat0,dum)

  return

end subroutine getgtr
